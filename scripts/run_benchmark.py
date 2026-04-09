#!/usr/bin/env python3
"""
Run StarKIT on all benchmark genomes and collect results.

Usage:
    python scripts/run_benchmark.py
"""

import csv
import os
import sys
import time

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starkit.starkit import run
from starkit.write import write_tsv, write_fasta
from starkit.report import generate_report

GENOME_DIR = "benchmarking/genomes"
RESULTS_DIR = "benchmarking/results"
SUMMARY_FILE = "benchmarking/benchmark_summary.tsv"


def run_starkit_on_genome(gbk_path, output_prefix):
    """Run StarKIT on a single genome, return the StarKITRun object."""
    try:
        starkit_run = run(
            input_file=gbk_path,
            output_prefix=output_prefix,
            gff_file=None,
            evalue=1e-10,
            min_size=0,
            max_size=0,
            evidence="all",
            use_log=False,
            quiet=True,
            no_homology=False,
        )
        write_tsv(starkit_run, output_prefix)
        write_fasta(starkit_run, output_prefix)
        generate_report(starkit_run, output_prefix)
        return starkit_run
    except Exception as e:
        print(f"    ERROR: {e}")
        return None


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load accession list
    accessions = {}
    with open("scripts/benchmark_accessions.txt") as f:
        for line in f:
            parts = line.strip().split("\t", 1)
            accessions[parts[0]] = parts[1] if len(parts) > 1 else ""

    # Find downloaded genomes
    gbk_files = sorted(
        f for f in os.listdir(GENOME_DIR) if f.endswith(".gbk")
    )

    print(f"Running StarKIT on {len(gbk_files)} genomes...")
    print()

    summary_rows = []

    for i, gbk_file in enumerate(gbk_files, 1):
        acc = gbk_file.replace(".gbk", "")
        org = accessions.get(acc, "")
        gbk_path = os.path.join(GENOME_DIR, gbk_file)
        output_prefix = os.path.join(RESULTS_DIR, acc)

        # Skip if already processed
        tsv_path = f"{output_prefix}.tsv"
        if os.path.exists(tsv_path) and os.path.getsize(tsv_path) > 0:
            # Read existing results
            with open(tsv_path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
            n_ships = len(rows)
            n_high = sum(1 for r in rows if r.get("evidence_level") == "high")
            n_medium = sum(1 for r in rows if r.get("evidence_level") == "medium")
            n_novel = sum(1 for r in rows if r.get("novelty") == "new")
            print(f"  [{i:>3}/{len(gbk_files)}] {acc} {org[:40]:<40} CACHED: {n_ships} ships")
            summary_rows.append({
                "accession": acc,
                "organism": org,
                "total_starships": n_ships,
                "high_confidence": n_high,
                "medium_confidence": n_medium,
                "novel": n_novel,
                "runtime": "cached",
                "status": "ok",
            })
            continue

        start = time.time()
        print(f"  [{i:>3}/{len(gbk_files)}] {acc} {org[:40]:<40} ", end="", flush=True)

        result = run_starkit_on_genome(gbk_path, output_prefix)
        elapsed = time.time() - start

        if result:
            n_ships = len(result.starships)
            n_high = sum(1 for s in result.starships if s.evidence_level.value == "high")
            n_medium = sum(1 for s in result.starships if s.evidence_level.value == "medium")
            n_novel = sum(1 for s in result.starships if s.is_novel)
            print(f"{n_ships} ships ({n_high}H/{n_medium}M) in {elapsed:.1f}s")
            summary_rows.append({
                "accession": acc,
                "organism": org,
                "total_starships": n_ships,
                "high_confidence": n_high,
                "medium_confidence": n_medium,
                "novel": n_novel,
                "runtime": f"{elapsed:.1f}",
                "status": "ok",
            })
        else:
            print(f"FAILED in {elapsed:.1f}s")
            summary_rows.append({
                "accession": acc,
                "organism": org,
                "total_starships": 0,
                "high_confidence": 0,
                "medium_confidence": 0,
                "novel": 0,
                "runtime": f"{elapsed:.1f}",
                "status": "failed",
            })

    # Write summary
    with open(SUMMARY_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "accession", "organism", "total_starships", "high_confidence",
            "medium_confidence", "novel", "runtime", "status",
        ])
        writer.writeheader()
        writer.writerows(summary_rows)

    # Print summary stats
    total_genomes = len(summary_rows)
    ok = [r for r in summary_rows if r["status"] == "ok"]
    total_ships = sum(r["total_starships"] for r in ok)
    total_high = sum(r["high_confidence"] for r in ok)
    total_medium = sum(r["medium_confidence"] for r in ok)
    total_novel = sum(r["novel"] for r in ok)

    print()
    print(f"{'='*60}")
    print(f"BENCHMARK SUMMARY")
    print(f"{'='*60}")
    print(f"Genomes processed: {len(ok)}/{total_genomes}")
    print(f"Total Starships:   {total_ships}")
    print(f"  High confidence: {total_high}")
    print(f"  Medium:          {total_medium}")
    print(f"  Novel:           {total_novel}")
    print(f"Summary written to {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
