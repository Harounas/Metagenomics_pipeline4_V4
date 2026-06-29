#!/usr/bin/env python3
"""
viral_classification_workflow.py

Workflow starting from assembled per-sample contigs:
  1. Merge per-sample contigs into one FASTA
  2. Kraken2 taxonomic classification on merged contigs
  3. DIAMOND blastx on merged contigs -> identify viral candidates
  4. Union viral contig IDs (Kraken2 + DIAMOND)
  5. Extract viral contigs -> viral_contigs_merged.fasta
  6. CD-HIT-EST clustering
  7. DIAMOND blastx on clustered representatives
  8. Annotate with virus names -> diamond_results_contig_with_sampleid.tsv
  9. Merge cluster info + DIAMOND annotations -> final TSV

Output TSV columns:
  query_id, clstr, clstr_size, length, clstr_rep, clstr_iden, clstr_cov,
  Sample_ID, contigs_len, virus, evalue, bitscore, pident, qcov
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

VIRAL_KEYWORDS = ("virus", "virinae", "viridae", "viricota", "viricetes", "phage")
KRAKEN_VIRAL_RANKS = {"F", "F1", "F2", "G", "G1", "G2", "S", "S1", "S2"}


# ---------------------------------------------------------------------------
# Step 1 – Merge per-sample contigs
# ---------------------------------------------------------------------------

def merge_contigs(output_dir: str, merged_fasta: str, min_length: int = 200) -> Path:
    """
    Scan output_dir for per-sample contigs and write a single merged FASTA.

    Looks for two layouts:
      - {output_dir}/{sample}_contigs.fasta        (flat)
      - {output_dir}/{sample}/contigs.fasta        (per-sample subdirectory)

    Renames every sequence to  {sample_id}|{original_contig_id}.
    Returns the Path of the merged FASTA.
    """
    base_dir = Path(output_dir)
    records = []
    seen_samples: set = set()

    for contigs_file in sorted(base_dir.glob("*_contigs.fasta")):
        sample_id = contigs_file.stem.replace("_contigs", "")
        seen_samples.add(sample_id)
        for rec in SeqIO.parse(contigs_file, "fasta"):
            if len(rec.seq) >= min_length:
                rec.id = f"{sample_id}|{rec.id}"
                rec.description = ""
                records.append(rec)

    for contigs_file in sorted(base_dir.glob("*/contigs.fasta")):
        sample_id = contigs_file.parent.name
        if sample_id in seen_samples:
            continue
        seen_samples.add(sample_id)
        for rec in SeqIO.parse(contigs_file, "fasta"):
            if len(rec.seq) >= min_length:
                rec.id = f"{sample_id}|{rec.id}"
                rec.description = ""
                records.append(rec)

    out_path = Path(merged_fasta)
    SeqIO.write(records, out_path, "fasta")
    print(f"Merged {len(records)} contigs (>={min_length} bp) "
          f"from {len(seen_samples)} samples -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Step 2 – Kraken2 on merged contigs
# ---------------------------------------------------------------------------

def run_kraken2_on_contigs(merged_fasta: str, kraken_db: str,
                            output_dir: str, threads: int = 8):
    """
    Run Kraken2 in single-end mode on the merged contigs FASTA.
    Returns (report_path, output_path).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    report = out_dir / "merged_contigs_kraken_report.txt"
    kout   = out_dir / "merged_contigs_kraken_output.txt"

    cmd = [
        "kraken2", "--db", kraken_db,
        "--report", str(report),
        "--output", str(kout),
        "--threads", str(threads),
        "--use-names",
        str(merged_fasta),
    ]
    print("Running Kraken2:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return report, kout


# ---------------------------------------------------------------------------
# Step 3 – Extract viral IDs from Kraken2 results
# ---------------------------------------------------------------------------

def extract_kraken_viral_ids(kraken_output: str, kraken_report: str) -> set:
    """
    Return the set of contig IDs classified as viral by Kraken2.

    Builds a viral taxid set from the report (rows whose rank is in
    KRAKEN_VIRAL_RANKS and whose name contains a VIRAL_KEYWORD), then
    collects contig IDs from the per-read output that matched those taxids.
    """
    viral_taxids: set = set()
    with open(kraken_report, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 6:
                continue
            rank  = row[3].strip()
            taxid = row[4].strip()
            name  = row[5].strip().lower()
            if rank in KRAKEN_VIRAL_RANKS and any(k in name for k in VIRAL_KEYWORDS):
                viral_taxids.add(taxid)

    viral_ids: set = set()
    with open(kraken_output, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 3:
                continue
            contig_id = row[1].strip()
            info      = row[2]
            for tid in viral_taxids:
                if f"taxid {tid}" in info or f"(taxid {tid})" in info:
                    viral_ids.add(contig_id)
                    break

    print(f"Kraken2 viral contig IDs: {len(viral_ids)}")
    return viral_ids


# ---------------------------------------------------------------------------
# Step 4 – Initial DIAMOND pass on merged contigs
# ---------------------------------------------------------------------------

def run_diamond_initial(merged_fasta: str, diamond_db: str,
                         output_dir: str, threads: int = 8) -> Path:
    """
    Run DIAMOND blastx on merged contigs.

    Output includes 'stitle' (subject title) so viral hits can be identified
    by keyword without rescanning the NR FASTA at this stage.
    Uses --sensitive and top-1 hit per query for speed.
    """
    out_file = Path(output_dir) / "diamond_initial.m8"
    cmd = [
        "diamond", "blastx",
        "--query",  str(merged_fasta),
        "--db",     diamond_db,
        "--out",    str(out_file),
        "--threads", str(threads),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "stitle",
        "--sensitive",
        "--max-target-seqs", "1",
        "--evalue", "1e-5",
    ]
    print("Running initial DIAMOND:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_file


# ---------------------------------------------------------------------------
# Step 5 – Extract viral IDs from DIAMOND initial results
# ---------------------------------------------------------------------------

def extract_diamond_viral_ids(diamond_m8: str,
                               filter_by_keywords: bool = False) -> set:
    """
    Return contig IDs from the initial DIAMOND output.

    When filter_by_keywords=True (mixed NR database), only IDs whose subject
    title contains a viral keyword are returned.  When False (default, for
    viral-focused databases such as nr_genomad.dmnd), all hits are treated as
    viral candidates — false positives are removed later by the annotation step.
    """
    col_names = [
        "query_id", "subject_id", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "stitle",
    ]
    df = pd.read_csv(diamond_m8, sep="\t", header=None, names=col_names)

    if filter_by_keywords:
        mask = df["stitle"].str.lower().str.contains(
            "|".join(VIRAL_KEYWORDS), na=False)
        viral_ids = set(df.loc[mask, "query_id"].unique())
    else:
        viral_ids = set(df["query_id"].unique())

    print(f"DIAMOND viral contig IDs: {len(viral_ids)}")
    return viral_ids


# ---------------------------------------------------------------------------
# Step 6 – Write merged viral contigs
# ---------------------------------------------------------------------------

def write_viral_contigs(merged_fasta: str, viral_ids: set,
                         output_fasta: str) -> Path:
    """Extract contigs matching viral_ids from merged_fasta and write to output_fasta."""
    records = [r for r in SeqIO.parse(merged_fasta, "fasta") if r.id in viral_ids]
    out_path = Path(output_fasta)
    SeqIO.write(records, out_path, "fasta")
    print(f"Wrote {len(records)} viral contigs -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Full workflow orchestrator
# ---------------------------------------------------------------------------

def run_full_workflow(
    output_dir:           str,
    kraken_db:            str,
    diamond_db:           str,
    nr_path:              str,
    threads:              int  = 32,
    min_length:           int  = 200,
    skip_existing:        bool = False,
    filter_by_keywords:   bool = False,
) -> str:
    """
    End-to-end viral classification workflow starting from per-sample contigs.

    Parameters
    ----------
    output_dir          Directory containing *_contigs.fasta files (phase 1 output).
    kraken_db           Kraken2 database path.
    diamond_db          DIAMOND database (.dmnd) path.
    nr_path             NR protein FASTA path used for virus-name annotation.
    threads             CPU threads for all tools.
    min_length          Minimum contig length (bp) to include.
    skip_existing       Skip steps whose output files already exist.
    filter_by_keywords  Filter initial DIAMOND hits by viral keywords in title.
                        Set True when diamond_db is the full NR; leave False for
                        viral-focused databases (e.g. nr_genomad.dmnd).

    Returns
    -------
    Path to the final TSV file.
    """
    from Metagenomics_pipeline4_V2.extract_contigs_diamond import (
        cluster_contigs,
        run_diamond,
        process_virus_contigs,
    )
    from Metagenomics_pipeline4_V2.process_clustered_contigs import (
        process_clustered_contigs,
    )

    out            = Path(output_dir)
    merged_fasta   = out / "merged_contigs.fasta"
    kraken_dir     = out / "kraken_output"
    viral_fasta    = out / "viral_contigs_merged.fasta"
    clustered_dir  = out / "clustered_output"
    clust_fasta    = clustered_dir / "clustered_contigs.fasta"
    clstr_file     = clustered_dir / "clustered_contigs.fasta.clstr"
    diamond_result = out / "results_clustered.m8"
    diamond_tsv    = out / "diamond_results_contig_with_sampleid.tsv"

    # ── 1. Merge per-sample contigs ────────────────────────────────────────
    if not skip_existing or not merged_fasta.exists():
        merge_contigs(output_dir, str(merged_fasta), min_length)
    else:
        print(f"[skip] Using existing merged contigs: {merged_fasta}")

    # ── 2. Kraken2 on merged contigs ──────────────────────────────────────
    kraken_report = kraken_dir / "merged_contigs_kraken_report.txt"
    kraken_output = kraken_dir / "merged_contigs_kraken_output.txt"
    if not skip_existing or not kraken_output.exists():
        kraken_report, kraken_output = run_kraken2_on_contigs(
            str(merged_fasta), kraken_db, str(kraken_dir), threads)
    else:
        print(f"[skip] Using existing Kraken2 output: {kraken_output}")

    # ── 3. Extract Kraken2 viral IDs ──────────────────────────────────────
    kraken_viral = extract_kraken_viral_ids(str(kraken_output), str(kraken_report))

    # ── 4. DIAMOND initial pass ───────────────────────────────────────────
    initial_m8 = out / "diamond_initial.m8"
    if not skip_existing or not initial_m8.exists():
        initial_m8 = run_diamond_initial(
            str(merged_fasta), diamond_db, str(out), threads)
    else:
        print(f"[skip] Using existing initial DIAMOND: {initial_m8}")

    # ── 5. Extract DIAMOND viral IDs ──────────────────────────────────────
    diamond_viral = extract_diamond_viral_ids(
        str(initial_m8), filter_by_keywords=filter_by_keywords)

    # ── 6. Write merged viral contigs ─────────────────────────────────────
    all_viral = kraken_viral | diamond_viral
    print(f"Total unique viral contigs (Kraken2 + DIAMOND): {len(all_viral)}")
    if not skip_existing or not viral_fasta.exists():
        write_viral_contigs(str(merged_fasta), all_viral, str(viral_fasta))

    if not viral_fasta.exists() or viral_fasta.stat().st_size == 0:
        print("No viral contigs found — cannot continue.")
        return ""

    # ── 7. CD-HIT-EST clustering ──────────────────────────────────────────
    if not skip_existing or not clust_fasta.exists():
        cluster_contigs(
            viral_fasta, str(clustered_dir),
            final_output="clustered_contigs.fasta",
            threads=threads,
        )
    else:
        print(f"[skip] Using existing clustered contigs: {clust_fasta}")

    # ── 8. DIAMOND on clustered contigs ───────────────────────────────────
    if not skip_existing or not diamond_result.exists():
        run_diamond(diamond_db, str(clust_fasta), str(diamond_result), threads)
    else:
        print(f"[skip] Using existing clustered DIAMOND: {diamond_result}")

    # ── 9. Annotate DIAMOND results with virus names ───────────────────────
    process_virus_contigs(nr_path, str(diamond_result), output_dir)

    # ── 10. Build final TSV ────────────────────────────────────────────────
    final_tsv = process_clustered_contigs(
        str(clstr_file), str(diamond_tsv), output_dir)

    print(f"\nWorkflow complete! Final TSV: {final_tsv}")
    return final_tsv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Viral classification workflow: merged contigs -> Kraken2 "
                    "+ DIAMOND -> CD-HIT -> annotated TSV")

    parser.add_argument("--output_dir",  required=True,
                        help="Directory with per-sample *_contigs.fasta "
                             "(phase 1 output); also used for all workflow outputs")
    parser.add_argument("--kraken_db",   required=True,
                        help="Path to Kraken2 database")
    parser.add_argument("--diamond_db",  required=True,
                        help="Path to DIAMOND database (.dmnd)")
    parser.add_argument("--nr_path",     required=True,
                        help="Path to NR protein FASTA for virus-name annotation")
    parser.add_argument("--threads",     type=int, default=32,
                        help="CPU threads (default: 32)")
    parser.add_argument("--min_length",  type=int, default=200,
                        help="Minimum contig length in bp (default: 200)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip steps whose output files already exist")
    parser.add_argument("--filter_by_keywords", action="store_true",
                        help="Filter initial DIAMOND hits by viral keywords "
                             "(use when diamond_db is full NR, not a viral-only DB)")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    run_full_workflow(
        output_dir          = args.output_dir,
        kraken_db           = args.kraken_db,
        diamond_db          = args.diamond_db,
        nr_path             = args.nr_path,
        threads             = args.threads,
        min_length          = args.min_length,
        skip_existing       = args.skip_existing,
        filter_by_keywords  = args.filter_by_keywords,
    )


if __name__ == "__main__":
    main()
