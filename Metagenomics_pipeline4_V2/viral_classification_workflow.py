#!/usr/bin/env python3
"""
viral_classification_workflow.py

Workflow starting from assembled per-sample contigs:
  1.  Merge per-sample contigs into one FASTA
  2.  Kraken2 taxonomic classification on merged contigs
  3.  DIAMOND blastx (initial pass) on merged contigs -> viral candidates
  4.  Union viral contig IDs (Kraken2 + DIAMOND)
  5.  Extract viral contigs -> viral_contigs_merged.fasta
  6.  CD-HIT-EST clustering
  7.  DIAMOND blastx on clustered contigs (with stitle for virus names)
  8.  Build diamond_results_contig_with_sampleid.tsv from stitle
  9.  Merge cluster info + DIAMOND annotations -> final TSV

Output TSV columns:
  query_id, clstr, clstr_size, length, clstr_rep, clstr_iden, clstr_cov,
  Sample_ID, contigs_len, virus, evalue, bitscore, pident, qcov
"""

import argparse
import csv
import os
import re
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
    """Run Kraken2 in single-end mode on the merged contigs FASTA."""
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
    """Return contig IDs classified as viral by Kraken2."""
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
# Step 4 – Initial DIAMOND pass (viral candidate detection)
# ---------------------------------------------------------------------------

def run_diamond_initial(merged_fasta: str, diamond_db: str,
                         output_dir: str, threads: int = 8) -> Path:
    """
    Run DIAMOND blastx (standard outfmt 6, top-1 hit) on merged contigs.
    All hits are treated as viral candidates; non-viral contigs are
    filtered out later in the annotation step.
    """
    out_file = Path(output_dir) / "diamond_initial.m8"
    cmd = [
        "diamond", "blastx",
        "--query",           str(merged_fasta),
        "--db",              diamond_db,
        "--out",             str(out_file),
        "--threads",         str(threads),
        "--outfmt",          "6",
        "--sensitive",
        "--max-target-seqs", "1",
        "--evalue",          "1e-5",
    ]
    print("Running initial DIAMOND:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_file


# ---------------------------------------------------------------------------
# Step 5 – Extract viral IDs from initial DIAMOND results
# ---------------------------------------------------------------------------

def extract_diamond_viral_ids(diamond_m8: str) -> set:
    """Return all query contig IDs that got any DIAMOND hit."""
    col_names = [
        "query_id", "subject_id", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    df = pd.read_csv(diamond_m8, sep="\t", header=None, names=col_names)
    viral_ids = set(df["query_id"].unique())
    print(f"DIAMOND viral contig IDs: {len(viral_ids)}")
    return viral_ids


# ---------------------------------------------------------------------------
# Step 6 – Write merged viral contigs
# ---------------------------------------------------------------------------

def write_viral_contigs(merged_fasta: str, viral_ids: set,
                         output_fasta: str) -> Path:
    """Extract contigs matching viral_ids and write to output_fasta."""
    records = [r for r in SeqIO.parse(merged_fasta, "fasta") if r.id in viral_ids]
    out_path = Path(output_fasta)
    SeqIO.write(records, out_path, "fasta")
    print(f"Wrote {len(records)} viral contigs -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Step 7 – Final DIAMOND on clustered contigs (with stitle)
# ---------------------------------------------------------------------------

def run_diamond_with_stitle(query_fasta: str, diamond_db: str,
                             output_file: str, threads: int = 8) -> Path:
    """
    Run DIAMOND blastx with stitle included in the output.
    Passes all outfmt fields as a single string (required by DIAMOND).
    """
    out_path = Path(output_file)
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    cmd = [
        "diamond", "blastx",
        "--query",   str(query_fasta),
        "--db",      diamond_db,
        "--out",     str(out_path),
        "--threads", str(threads),
        "--outfmt",  outfmt,
        "--sensitive",
        "--evalue",  "1e-5",
    ]
    print("Running final DIAMOND (with stitle):", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_path


# ---------------------------------------------------------------------------
# Step 8 – Build diamond_results_contig_with_sampleid.tsv from stitle
# ---------------------------------------------------------------------------

def build_diamond_tsv(diamond_m8_with_stitle: str,
                       output_dir: str) -> Path:
    """
    Parse DIAMOND output (with stitle) and produce the annotated TSV
    expected by process_clustered_contigs.

    Extracts virus name from stitle using bracket pattern [Virus name]
    or the full stitle when no brackets are present.
    Selects best hit per query by bitscore.
    """
    col_names = [
        "query_id", "subject_id", "pident", "aln_len", "mismatches", "gaps",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
    ]
    df = pd.read_csv(diamond_m8_with_stitle, sep="\t", header=None,
                     names=col_names)

    # Extract virus name from stitle: prefer text inside [...] containing a keyword
    def parse_virus(stitle: str) -> str:
        matches = re.findall(r'\[([^\]]+)\]', str(stitle))
        for m in matches:
            if any(k in m.lower() for k in VIRAL_KEYWORDS):
                return m
        # Fallback: return the full stitle up to first bracket or truncated
        return str(stitle).split("[")[0].strip() or str(stitle)

    df["virus"] = df["stitle"].apply(parse_virus)

    # Keep best hit per query (highest bitscore)
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    best = df.loc[df.groupby("query_id")["bitscore"].idxmax()].copy()

    # Extract Sample_ID and contig length from query_id
    best["Sample_ID"] = best["query_id"].str.split("|").str[0]
    tail = best["query_id"].str.split("|").str[-1]
    best["contigs_len"] = pd.to_numeric(
        tail.str.extract(r'length_(\d+)', expand=False), errors="coerce")
    fallback = pd.to_numeric(tail.str.split("_").str[3], errors="coerce")
    best["contigs_len"] = best["contigs_len"].fillna(fallback)

    # Compute query coverage
    best["qstart"]      = pd.to_numeric(best["qstart"],  errors="coerce")
    best["qend"]        = pd.to_numeric(best["qend"],    errors="coerce")
    best["aln_len"]     = pd.to_numeric(best["aln_len"], errors="coerce")
    best["contigs_len"] = pd.to_numeric(best["contigs_len"], errors="coerce")
    best["qcov"] = ((best["qend"] - best["qstart"] + 1) / best["contigs_len"]) * 100

    out_path = Path(output_dir) / "diamond_results_contig_with_sampleid.tsv"
    best.to_csv(out_path, sep="\t", index=False)
    print(f"Diamond TSV written -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Full workflow orchestrator
# ---------------------------------------------------------------------------

def run_full_workflow(
    output_dir:    str,
    kraken_db:     str,
    diamond_db:    str,
    threads:       int  = 32,
    min_length:    int  = 200,
    skip_existing: bool = False,
) -> str:
    """
    End-to-end viral classification workflow starting from per-sample contigs.
    No NR FASTA required — virus names are extracted directly from DIAMOND stitle.
    """
    from Metagenomics_pipeline4_V2.extract_contigs_diamond import cluster_contigs
    from Metagenomics_pipeline4_V2.process_clustered_contigs import process_clustered_contigs

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
        print(f"[skip] merged_contigs.fasta exists")

    # ── 2. Kraken2 on merged contigs ──────────────────────────────────────
    kraken_report = kraken_dir / "merged_contigs_kraken_report.txt"
    kraken_output = kraken_dir / "merged_contigs_kraken_output.txt"
    if not skip_existing or not kraken_output.exists():
        kraken_report, kraken_output = run_kraken2_on_contigs(
            str(merged_fasta), kraken_db, str(kraken_dir), threads)
    else:
        print(f"[skip] Kraken2 output exists")

    # ── 3. Extract Kraken2 viral IDs ──────────────────────────────────────
    kraken_viral = extract_kraken_viral_ids(str(kraken_output), str(kraken_report))

    # ── 4. DIAMOND initial pass ───────────────────────────────────────────
    initial_m8 = out / "diamond_initial.m8"
    if not skip_existing or not initial_m8.exists():
        initial_m8 = run_diamond_initial(
            str(merged_fasta), diamond_db, str(out), threads)
    else:
        print(f"[skip] diamond_initial.m8 exists")

    # ── 5. Extract DIAMOND viral IDs ──────────────────────────────────────
    diamond_viral = extract_diamond_viral_ids(str(initial_m8))

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
        print(f"[skip] clustered_contigs.fasta exists")

    # ── 8. Final DIAMOND with stitle ──────────────────────────────────────
    if not skip_existing or not diamond_result.exists():
        run_diamond_with_stitle(
            str(clust_fasta), diamond_db, str(diamond_result), threads)
    else:
        print(f"[skip] results_clustered.m8 exists")

    # ── 9. Build annotated TSV from stitle ────────────────────────────────
    if not skip_existing or not diamond_tsv.exists():
        build_diamond_tsv(str(diamond_result), output_dir)
    else:
        print(f"[skip] diamond_results_contig_with_sampleid.tsv exists")

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
                        help="Directory with per-sample *_contigs.fasta files "
                             "(phase 1 output); also used for all workflow outputs")
    parser.add_argument("--kraken_db",   required=True,
                        help="Path to Kraken2 database")
    parser.add_argument("--diamond_db",  required=True,
                        help="Path to DIAMOND database (.dmnd)")
    parser.add_argument("--threads",     type=int, default=32,
                        help="CPU threads (default: 32)")
    parser.add_argument("--min_length",  type=int, default=200,
                        help="Minimum contig length in bp (default: 200)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip steps whose output files already exist")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    run_full_workflow(
        output_dir    = args.output_dir,
        kraken_db     = args.kraken_db,
        diamond_db    = args.diamond_db,
        threads       = args.threads,
        min_length    = args.min_length,
        skip_existing = args.skip_existing,
    )


if __name__ == "__main__":
    main()
