import os
import shutil
import subprocess
import tempfile
import logging

def run_spades(forward, reverse, base_name, output_dir, threads=8):
    """
    Runs MetaSPAdes for de novo assembly.
    SPAdes writes configs and temp files to a local tmpdir to avoid NFS I/O
    errors, then copies contigs.fasta to the final NFS output directory.
    """
    sample_outdir = os.path.join(output_dir, base_name)
    contigs_output = os.path.join(sample_outdir, "contigs.fasta")

    if os.path.exists(contigs_output):
        logging.info(f"[SKIP] SPAdes assembly already exists for {base_name}.")
        return contigs_output

    os.makedirs(sample_outdir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"spades_{base_name}_") as tmp_spades:
        cmd = [
            "metaspades.py",
            "-1", forward,
            "-2", reverse,
            "-o", tmp_spades,
            "-t", str(threads)
        ]

        logging.info(f"[RUN] Running MetaSPAdes for {base_name}: {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SPAdes failed for {base_name} (exit {e.returncode})") from e

        tmp_contigs = os.path.join(tmp_spades, "contigs.fasta")
        if not os.path.exists(tmp_contigs):
            raise RuntimeError(f"SPAdes finished but no contigs.fasta found for {base_name}")

        shutil.copy2(tmp_contigs, contigs_output)
        logging.info(f"[DONE] Assembly complete for {base_name}: {contigs_output}")

    return contigs_output
