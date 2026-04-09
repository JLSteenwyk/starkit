#!/usr/bin/env python3
"""
Download annotated GenBank files for benchmark genomes from NCBI.

Reads accessions from scripts/benchmark_accessions.txt and downloads
GBFF files via the NCBI Datasets API.

Usage:
    python scripts/download_benchmark_genomes.py
"""

import os
import subprocess
import sys
import zipfile
import tempfile
import shutil

ACCESSION_FILE = "scripts/benchmark_accessions.txt"
OUTPUT_DIR = "benchmarking/genomes"


def download_genome(accession, output_dir):
    """Download a single genome's GBFF file from NCBI."""
    output_path = os.path.join(output_dir, f"{accession}.gbk")

    if os.path.exists(output_path) and os.path.getsize(output_path) > 1000:
        return output_path, "already exists"

    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{accession}/download?include_annotation_type=GENOME_GBFF"
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = os.path.join(tmpdir, "dataset.zip")

        # Download
        result = subprocess.run(
            ["curl", "-L", "-s", "-o", zip_path, url],
            capture_output=True, timeout=300,
        )
        if result.returncode != 0:
            return None, f"curl failed: {result.stderr.decode()[:100]}"

        if not os.path.exists(zip_path) or os.path.getsize(zip_path) < 100:
            return None, "download too small"

        # Extract GBFF
        try:
            with zipfile.ZipFile(zip_path) as zf:
                gbff_files = [n for n in zf.namelist() if n.endswith(".gbff")]
                if not gbff_files:
                    return None, "no .gbff in zip"

                with zf.open(gbff_files[0]) as src, open(output_path, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        except zipfile.BadZipFile:
            return None, "bad zip file"

    size_mb = os.path.getsize(output_path) / 1024 / 1024
    return output_path, f"downloaded ({size_mb:.1f} MB)"


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(ACCESSION_FILE) as f:
        entries = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t", 1)
            acc = parts[0]
            org = parts[1] if len(parts) > 1 else ""
            entries.append((acc, org))

    print(f"Downloading {len(entries)} genomes to {OUTPUT_DIR}/")
    print()

    success = 0
    failed = []

    for i, (acc, org) in enumerate(entries, 1):
        path, status = download_genome(acc, OUTPUT_DIR)
        if path:
            success += 1
            print(f"  [{i:>3}/{len(entries)}] {acc} {org[:35]:<35} {status}")
        else:
            failed.append((acc, org, status))
            print(f"  [{i:>3}/{len(entries)}] {acc} {org[:35]:<35} FAILED: {status}")

    print(f"\nDone: {success} downloaded, {len(failed)} failed")
    if failed:
        print("\nFailed genomes:")
        for acc, org, reason in failed:
            print(f"  {acc} {org}: {reason}")


if __name__ == "__main__":
    main()
