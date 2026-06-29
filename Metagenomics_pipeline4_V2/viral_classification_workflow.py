#!/usr/bin/env python3
"""
viral_classification_workflow.py

Workflow starting from assembled per-sample contigs:
  1.  Merge per-sample contigs into one FASTA
  2.  Kraken2 taxonomic classification on merged contigs -> viral contig IDs
  3.  geNomad end-to-end on merged contigs -> viral contigs FASTA
  4.  Merge viral contigs (Kraken2 IDs + geNomad FASTA, deduplicated)
  5.  CD-HIT-EST clustering
  6.  DIAMOND blastx on clustered contigs (with stitle for virus names)
  7.  Build diamond_results_contig_with_sampleid.tsv from stitle
  8.  Merge cluster info + DIAMOND annotations -> final TSV

Output TSV columns:
  query_id, clstr, clstr_size, length, clstr_rep, clstr_iden, clstr_cov,
  Sample_ID, contigs_len, virus, evalue, bitscore, pident, qcov
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

VIRAL_KEYWORDS   = ("virus", "virinae", "viridae", "viricota", "viricetes", "phage")
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

    Renames every sequence to {sample_id}|{original_contig_id}.
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
# Step 3 – geNomad on merged contigs
# ---------------------------------------------------------------------------

def run_genomad_on_contigs(merged_fasta: str, genomad_db: str,
                            output_dir: str, threads: int = 8,
                            min_score: float = 0.5,
                            splits: int = 8,
                            genomad_min_length: int = 1000) -> Path:
    """
    Run geNomad end-to-end on merged contigs.

    Pre-filters contigs to genomad_min_length before running to reduce
    MMseqs2 memory. geNomad recommends >=1000 bp (2500 bp for best accuracy).
    Returns path to the viral contigs FASTA produced by geNomad.
    """
    genomad_bin = shutil.which("genomad")
    if genomad_bin is None:
        raise FileNotFoundError(
            "genomad executable not found in PATH.\n"
            "Make sure the genomad conda environment is activated:\n"
            "  eval \"$(mamba shell hook --shell bash)\" && mamba activate genomad")

    out_dir = Path(output_dir) / "genomad_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Pre-filter to genomad_min_length to reduce protein count for MMseqs2
    filtered_fasta = out_dir / "contigs_for_genomad.fasta"
    if not filtered_fasta.exists():
        records = [r for r in SeqIO.parse(merged_fasta, "fasta")
                   if len(r.seq) >= genomad_min_length]
        SeqIO.write(records, filtered_fasta, "fasta")
        print(f"Filtered {len(records)} contigs (>={genomad_min_length} bp) "
              f"for geNomad -> {filtered_fasta}")
    else:
        print(f"[skip] contigs_for_genomad.fasta exists")

    fasta_stem  = filtered_fasta.stem           # "contigs_for_genomad"
    summary_dir = out_dir / f"{fasta_stem}_summary"
    virus_fasta = summary_dir / f"{fasta_stem}_virus.fna"

    cmd = [
        genomad_bin, "end-to-end",
        str(filtered_fasta),
        str(out_dir),
        genomad_db,
        "--min-score", str(min_score),
        "--threads",   str(threads),
        "--splits",    str(splits),
        "--restart",
    ]
    print("Running geNomad:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    if not virus_fasta.exists():
        raise FileNotFoundError(
            f"geNomad viral FASTA not found: {virus_fasta}\n"
            f"Contents of {out_dir}: {list(out_dir.iterdir())}")

    return virus_fasta


def extract_genomad_viral_ids(genomad_virus_fasta: str) -> set:
    """Return set of contig IDs from the geNomad viral FASTA."""
    ids = {rec.id for rec in SeqIO.parse(genomad_virus_fasta, "fasta")}
    print(f"geNomad viral contig IDs: {len(ids)}")
    return ids


# ---------------------------------------------------------------------------
# Step 4 – Merge Kraken2 + geNomad viral contigs (deduplicated)
# ---------------------------------------------------------------------------

def write_merged_viral_contigs(merged_fasta: str,
                                viral_ids: set,
                                output_fasta: str) -> Path:
    """
    Extract contigs from merged_fasta whose IDs are in viral_ids and
    write them to output_fasta. Deduplication is by contig ID.
    """
    seen: set = set()
    records = []
    for rec in SeqIO.parse(merged_fasta, "fasta"):
        if rec.id in viral_ids and rec.id not in seen:
            seen.add(rec.id)
            records.append(rec)

    out_path = Path(output_fasta)
    SeqIO.write(records, out_path, "fasta")
    print(f"Wrote {len(records)} merged viral contigs -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Step 5 – CD-HIT-EST clustering  (reuses existing function)
# ---------------------------------------------------------------------------
# cluster_contigs() from extract_contigs_diamond.py is called in the orchestrator


# ---------------------------------------------------------------------------
# Step 6 – Final DIAMOND on clustered contigs (with stitle)
# ---------------------------------------------------------------------------

def run_diamond_with_stitle(query_fasta: str, diamond_db: str,
                             output_file: str, threads: int = 8) -> Path:
    """
    Run DIAMOND blastx with stitle in the output so virus names can be
    extracted without scanning an NR FASTA file.
    outfmt fields are passed as separate list items (required by subprocess).
    """
    out_path = Path(output_file)
    cmd = [
        "diamond", "blastx",
        "--query",   str(query_fasta),
        "--db",      diamond_db,
        "--out",     str(out_path),
        "--threads", str(threads),
        "--outfmt",  "6",
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
        "--sensitive",
        "--evalue",  "1e-5",
    ]
    print("Running final DIAMOND:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_path


# ---------------------------------------------------------------------------
# Step 7 – Build annotated TSV from DIAMOND stitle
# ---------------------------------------------------------------------------

def build_diamond_tsv(diamond_m8_with_stitle: str, output_dir: str) -> Path:
    """
    Parse DIAMOND output (with stitle) and produce the annotated TSV
    expected by process_clustered_contigs.

    - Extracts virus name from stitle (prefers bracketed [Virus name])
    - Selects best hit per query by bitscore
    - Computes query coverage from qstart/qend/contig length
    """
    col_names = [
        "query_id", "subject_id", "pident", "aln_len", "mismatches", "gaps",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
    ]
    df = pd.read_csv(diamond_m8_with_stitle, sep="\t", header=None,
                     names=col_names)

    def parse_virus(stitle: str) -> str:
        matches = re.findall(r'\[([^\]]+)\]', str(stitle))
        for m in matches:
            if any(k in m.lower() for k in VIRAL_KEYWORDS):
                return m
        return str(stitle).split("[")[0].strip() or str(stitle)

    df["virus"]    = df["stitle"].apply(parse_virus)
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    best = df.loc[df.groupby("query_id")["bitscore"].idxmax()].copy()

    # Sample_ID and contig length from query_id  ({sample}|NODE_x_length_y_cov_z)
    best["Sample_ID"] = best["query_id"].str.split("|").str[0]
    tail = best["query_id"].str.split("|").str[-1]
    best["contigs_len"] = pd.to_numeric(
        tail.str.extract(r'length_(\d+)', expand=False), errors="coerce")
    fallback = pd.to_numeric(tail.str.split("_").str[3], errors="coerce")
    best["contigs_len"] = best["contigs_len"].fillna(fallback)

    # Query coverage
    for col in ("qstart", "qend", "aln_len", "contigs_len"):
        best[col] = pd.to_numeric(best[col], errors="coerce")
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
    genomad_db:    str,
    threads:       int   = 32,
    min_length:    int   = 200,
    min_score:     float = 0.5,
    splits:              int   = 8,
    genomad_min_length:  int   = 1000,
    skip_existing:       bool  = False,
) -> str:
    """
    End-to-end workflow:
      Kraken2 + geNomad -> merge viral contigs -> CD-HIT -> DIAMOND -> TSV
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

    kraken_viral = extract_kraken_viral_ids(str(kraken_output), str(kraken_report))

    # ── 3. geNomad on merged contigs ──────────────────────────────────────
    genomad_virus_fasta = out / "genomad_out" / "merged_contigs_summary" / "merged_contigs_virus.fna"
    if not skip_existing or not genomad_virus_fasta.exists():
        genomad_virus_fasta = run_genomad_on_contigs(
            str(merged_fasta), genomad_db, str(out), threads,
            min_score, splits, genomad_min_length)
    else:
        print(f"[skip] geNomad viral FASTA exists")

    genomad_viral = extract_genomad_viral_ids(str(genomad_virus_fasta))

    # ── 4. Merge Kraken2 + geNomad viral contigs ──────────────────────────
    all_viral = kraken_viral | genomad_viral
    print(f"Total unique viral contigs (Kraken2 + geNomad): {len(all_viral)}")
    if not skip_existing or not viral_fasta.exists():
        write_merged_viral_contigs(str(merged_fasta), all_viral, str(viral_fasta))

    if not viral_fasta.exists() or viral_fasta.stat().st_size == 0:
        print("No viral contigs found — cannot continue.")
        return ""

    # ── 5. CD-HIT-EST clustering ──────────────────────────────────────────
    if not skip_existing or not clust_fasta.exists():
        cluster_contigs(
            viral_fasta, str(clustered_dir),
            final_output="clustered_contigs.fasta",
            threads=threads,
        )
    else:
        print(f"[skip] clustered_contigs.fasta exists")

    # ── 6. Final DIAMOND on clustered contigs ─────────────────────────────
    diamond_done = diamond_result.exists() and diamond_result.stat().st_size > 0
    if not skip_existing or not diamond_done:
        run_diamond_with_stitle(
            str(clust_fasta), diamond_db, str(diamond_result), threads)
    else:
        print(f"[skip] results_clustered.m8 exists ({diamond_result.stat().st_size} bytes)")

    # ── 7. Build annotated TSV from stitle ────────────────────────────────
    diamond_tsv_done = diamond_tsv.exists() and diamond_tsv.stat().st_size > 0
    if not skip_existing or not diamond_tsv_done:
        build_diamond_tsv(str(diamond_result), output_dir)
    else:
        print(f"[skip] diamond_results_contig_with_sampleid.tsv exists")

    # ── 8. Build final TSV (cluster info + DIAMOND annotations) ───────────
    final_tsv = process_clustered_contigs(
        str(clstr_file), str(diamond_tsv), output_dir)

    print(f"\nWorkflow complete! Final TSV: {final_tsv}")
    return final_tsv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Viral classification: Kraken2 + geNomad -> CD-HIT -> DIAMOND -> TSV")

    parser.add_argument("--output_dir",  required=True,
                        help="Directory with per-sample *_contigs.fasta files")
    parser.add_argument("--kraken_db",   required=True,
                        help="Kraken2 database path")
    parser.add_argument("--diamond_db",  required=True,
                        help="DIAMOND database (.dmnd)")
    parser.add_argument("--genomad_db",  required=True,
                        help="geNomad database path")
    parser.add_argument("--threads",     type=int,   default=32)
    parser.add_argument("--min_length",  type=int,   default=200,
                        help="Minimum contig length in bp (default: 200)")
    parser.add_argument("--min_score",   type=float, default=0.5,
                        help="geNomad minimum virus score (default: 0.5)")
    parser.add_argument("--splits",             type=int, default=8,
                        help="geNomad MMseqs2 splits to reduce memory usage (default: 8)")
    parser.add_argument("--genomad_min_length", type=int, default=1000,
                        help="Min contig length (bp) fed to geNomad (default: 1000)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip steps whose output files already exist")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    run_full_workflow(
        output_dir    = args.output_dir,
        kraken_db     = args.kraken_db,
        diamond_db    = args.diamond_db,
        genomad_db    = args.genomad_db,
        threads       = args.threads,
        min_length    = args.min_length,
        min_score     = args.min_score,
        splits              = args.splits,
        genomad_min_length  = args.genomad_min_length,
        skip_existing       = args.skip_existing,
    )


if __name__ == "__main__":
    main()
