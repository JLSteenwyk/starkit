#!/usr/bin/env python3
"""
Export Starship nucleotide sequences from the Starbase SQLite database as a
FASTA reference file for use in homology searches.

Usage:
    python scripts/export_starship_fasta.py

Outputs:
    starkit/data/starships_ref.fasta
"""

import os
import sqlite3
import sys
from collections import Counter

DB_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "starbase_v0.1.0-pre.sqlite",
)
OUTPUT_FASTA = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    os.pardir,
    "starkit",
    "data",
    "starships_ref.fasta",
)
OUTPUT_FASTA = os.path.normpath(OUTPUT_FASTA)


def export_starship_reference_fasta(conn, output_path):
    """
    Extract all Starship nucleotide sequences from the ships table,
    join with family information, and write to a FASTA file.

    FASTA headers: >{ship_id}|{family_name}|{sequence_length}bp

    Returns a dict of statistics.
    """
    cursor = conn.cursor()

    # For ships appearing in multiple joined_ships rows (different families),
    # we pick one family per ship using MIN(fn.familyName) to be deterministic.
    # Ships with no joined_ships entry get NULL family via LEFT JOIN.
    cursor.execute("""
        SELECT
            s.id,
            s.sequence,
            s.sequence_length,
            MIN(fn.familyName) AS family_name
        FROM ships s
        LEFT JOIN joined_ships js ON s.id = js.ship_id
        LEFT JOIN family_names fn ON js.ship_family_id = fn.id
        WHERE s.sequence IS NOT NULL
          AND s.sequence != ''
        GROUP BY s.id
        ORDER BY s.id
    """)

    rows = cursor.fetchall()

    family_counts = Counter()
    size_bins = Counter()  # for size distribution
    total = 0

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as fasta_out:
        for ship_id, sequence, sequence_length, family_name in rows:
            # Normalize family name
            if family_name is None or family_name == "None" or family_name.strip() == "":
                family_name = "unclassified"

            # Use recorded length if available, otherwise compute it
            seq_clean = sequence.strip()
            if sequence_length is None:
                sequence_length = len(seq_clean)

            # Write FASTA entry
            header = f">{ship_id}|{family_name}|{sequence_length}bp"
            fasta_out.write(f"{header}\n")

            # Write sequence in 80-character lines
            for i in range(0, len(seq_clean), 80):
                fasta_out.write(seq_clean[i : i + 80] + "\n")

            family_counts[family_name] += 1
            total += 1

            # Bin sizes
            if sequence_length < 10000:
                size_bins["<10kb"] += 1
            elif sequence_length < 50000:
                size_bins["10-50kb"] += 1
            elif sequence_length < 100000:
                size_bins["50-100kb"] += 1
            elif sequence_length < 200000:
                size_bins["100-200kb"] += 1
            elif sequence_length < 500000:
                size_bins["500kb-1Mb"] += 1
            else:
                size_bins[">=500kb"] += 1

    stats = {
        "total": total,
        "family_counts": dict(family_counts.most_common()),
        "size_distribution": dict(
            sorted(size_bins.items(), key=lambda x: ["<10kb", "10-50kb", "50-100kb", "100-200kb", "500kb-1Mb", ">=500kb"].index(x[0]) if x[0] in ["<10kb", "10-50kb", "50-100kb", "100-200kb", "500kb-1Mb", ">=500kb"] else 99)
        ),
        "output_path": output_path,
    }
    return stats


def print_stats(stats):
    """Print summary statistics."""
    print(f"  Total sequences: {stats['total']}")
    print(f"  Output: {stats['output_path']}")
    print()

    print("  Per-family counts:")
    for family, count in stats["family_counts"].items():
        print(f"    {family}: {count}")
    print()

    print("  Size distribution:")
    for size_bin, count in stats["size_distribution"].items():
        print(f"    {size_bin}: {count}")


def main():
    db_path = sys.argv[1] if len(sys.argv) > 1 else DB_PATH
    output_path = sys.argv[2] if len(sys.argv) > 2 else OUTPUT_FASTA

    print("Starship Reference FASTA Exporter")
    print(f"  Database: {db_path}")
    print(f"  Output:   {output_path}")
    print()

    if not os.path.exists(db_path):
        print(f"Error: Database not found at {db_path}")
        sys.exit(1)

    conn = sqlite3.connect(db_path)

    print("Exporting Starship nucleotide sequences...")
    stats = export_starship_reference_fasta(conn, output_path)
    print_stats(stats)

    conn.close()
    print()
    print("Done!")


if __name__ == "__main__":
    main()
