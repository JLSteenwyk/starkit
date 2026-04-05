#!/usr/bin/env python3
"""
Build HMM profiles for StarKIT from the Starbase SQLite database.

Downloads captain protein sequences grouped by family, aligns them with
MUSCLE (or MAFFT), and builds HMM profiles with pyhmmer.

Usage:
    python scripts/build_hmm_profiles.py scripts/starbase_v0.1.0-pre.sqlite

Outputs:
    starkit/data/hmm/captains.hmm       — all-captain HMM for detection
    starkit/data/families/<family>.hmm   — per-family HMMs for classification
    starkit/data/boundaries/tir_pwms.json — TIR position-weight matrices
"""

import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, OrderedDict, defaultdict

import pyhmmer


DB_PATH = sys.argv[1] if len(sys.argv) > 1 else "scripts/starbase_v0.1.0-pre.sqlite"
OUTPUT_HMM_DIR = os.path.join("starkit", "data", "hmm")
OUTPUT_FAMILY_DIR = os.path.join("starkit", "data", "families")
OUTPUT_BOUNDARY_DIR = os.path.join("starkit", "data", "boundaries")
OUTPUT_STARSHIP_FASTA = os.path.join("starkit", "data", "starships_ref.fasta")
MIN_SEQUENCES_FOR_HMM = 3  # minimum sequences to build a family HMM


def get_connection(db_path):
    if not os.path.exists(db_path):
        print(f"Error: Database not found at {db_path}")
        sys.exit(1)
    return sqlite3.connect(db_path)


def extract_captains_by_family(conn):
    """Extract captain protein sequences grouped by family name."""
    cursor = conn.cursor()

    # Get captains with family assignments
    cursor.execute("""
        SELECT c.captainID, c.sequence, fn.familyName
        FROM captains c
        JOIN joined_ships js ON c.ship_id = js.ship_id
        JOIN family_names fn ON js.ship_family_id = fn.id
        WHERE c.sequence IS NOT NULL
          AND c.sequence != ''
          AND fn.familyName IS NOT NULL
          AND fn.familyName != 'None'
    """)

    families = defaultdict(list)
    all_captains = []

    for captain_id, sequence, family_name in cursor.fetchall():
        # Clean sequence: remove gaps, stops, whitespace
        seq = sequence.strip().replace("-", "").replace("*", "").replace(" ", "")
        if len(seq) < 100:  # skip very short sequences
            continue

        families[family_name].append((captain_id, seq))
        all_captains.append((captain_id, seq))

    return families, all_captains


def extract_captains_unassigned(conn):
    """Extract captain sequences without family assignments."""
    cursor = conn.cursor()

    cursor.execute("""
        SELECT c.captainID, c.sequence
        FROM captains c
        LEFT JOIN joined_ships js ON c.ship_id = js.ship_id
        WHERE c.sequence IS NOT NULL
          AND c.sequence != ''
          AND (js.ship_family_id IS NULL
               OR js.ship_family_id IN (
                   SELECT id FROM family_names WHERE familyName IS NULL OR familyName = 'None'
               ))
    """)

    captains = []
    for captain_id, sequence in cursor.fetchall():
        seq = sequence.strip().replace("-", "").replace("*", "").replace(" ", "")
        if len(seq) >= 100:
            captains.append((captain_id, seq))

    return captains


def write_fasta(sequences, output_path):
    """Write sequences as FASTA file. sequences = list of (id, seq) tuples."""
    with open(output_path, "w") as f:
        for seq_id, seq in sequences:
            f.write(f">{seq_id}\n{seq}\n")


def align_sequences(fasta_path, aligned_path):
    """Align sequences using MUSCLE or MAFFT."""
    # Try MUSCLE first
    try:
        subprocess.run(
            ["muscle", "-align", fasta_path, "-output", aligned_path],
            capture_output=True, check=True, timeout=300,
        )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass

    # Try MUSCLE v3 syntax
    try:
        subprocess.run(
            ["muscle", "-in", fasta_path, "-out", aligned_path],
            capture_output=True, check=True, timeout=300,
        )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass

    # Fallback to MAFFT
    try:
        with open(aligned_path, "w") as out:
            subprocess.run(
                ["mafft", "--auto", fasta_path],
                stdout=out, stderr=subprocess.DEVNULL, check=True, timeout=300,
            )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass

    print(f"  Warning: Neither MUSCLE nor MAFFT available. Skipping alignment.")
    return False


def build_hmm_from_alignment(aligned_path, hmm_name):
    """Build an HMM from a multiple sequence alignment using pyhmmer."""
    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.MSAFile(aligned_path, digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()

    msa.name = hmm_name.encode()

    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)

    return hmm


def build_hmms_for_families(families, output_dir):
    """Build per-family HMM profiles."""
    os.makedirs(output_dir, exist_ok=True)

    hmms = []
    for family_name, sequences in sorted(families.items()):
        if len(sequences) < MIN_SEQUENCES_FOR_HMM:
            print(f"  Skipping {family_name}: only {len(sequences)} sequences "
                  f"(need {MIN_SEQUENCES_FOR_HMM})")
            continue

        print(f"  Building HMM for {family_name} ({len(sequences)} sequences)...")

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, f"{family_name}.fasta")
            aligned_path = os.path.join(tmpdir, f"{family_name}_aligned.fasta")

            write_fasta(sequences, fasta_path)

            if not align_sequences(fasta_path, aligned_path):
                continue

            try:
                hmm = build_hmm_from_alignment(aligned_path, family_name)
                hmms.append(hmm)
                print(f"    -> Built HMM: {hmm.M} match states")
            except Exception as e:
                print(f"    -> Error building HMM: {e}")
                continue

    # Write all family HMMs to a single file
    if hmms:
        output_path = os.path.join(output_dir, "families.hmm")
        with open(output_path, "wb") as f:
            for hmm in hmms:
                hmm.write(f)
        print(f"  Wrote {len(hmms)} family HMMs to {output_path}")

    return hmms


def build_detection_hmm(all_captains, output_dir):
    """Build a single HMM from all captain sequences for detection."""
    os.makedirs(output_dir, exist_ok=True)

    print(f"  Aligning {len(all_captains)} captain sequences...")

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "all_captains.fasta")
        aligned_path = os.path.join(tmpdir, "all_captains_aligned.fasta")

        write_fasta(all_captains, fasta_path)

        if not align_sequences(fasta_path, aligned_path):
            return None

        try:
            hmm = build_hmm_from_alignment(aligned_path, "captain_all")
            output_path = os.path.join(output_dir, "captains.hmm")
            with open(output_path, "wb") as f:
                hmm.write(f)
            print(f"    -> Built detection HMM: {hmm.M} match states")
            print(f"    -> Wrote to {output_path}")
            return hmm
        except Exception as e:
            print(f"    -> Error: {e}")
            return None


def build_tir_pwms(conn, output_dir):
    """Build TIR position-weight matrices from Starbase boundary data."""
    os.makedirs(output_dir, exist_ok=True)

    cursor = conn.cursor()

    # Get all non-trivial TIR sequences (length > 5, not just dots)
    cursor.execute("""
        SELECT upTIR, downTIR
        FROM starship_features
        WHERE upTIR IS NOT NULL AND LENGTH(upTIR) > 5
          AND upTIR NOT LIKE '%.%'
    """)

    up_tirs = []
    down_tirs = []
    for up, down in cursor.fetchall():
        up = up.strip().upper()
        down = down.strip().upper() if down else ""
        if len(up) >= 8 and all(c in "ACGTN" for c in up):
            up_tirs.append(up)
        if len(down) >= 8 and all(c in "ACGTN" for c in down):
            down_tirs.append(down)

    print(f"  Found {len(up_tirs)} upstream TIRs, {len(down_tirs)} downstream TIRs")

    pwms = []

    # Build PWMs from TIR sequences grouped by length
    for label, tir_seqs in [("upstream", up_tirs), ("downstream", down_tirs)]:
        # Group by length
        by_length = defaultdict(list)
        for seq in tir_seqs:
            by_length[len(seq)].append(seq)

        # Build PWM for each common length group (>= 10 sequences)
        for length, seqs in sorted(by_length.items()):
            if len(seqs) < 10:
                continue

            # Build position frequency matrix
            matrix = []
            for pos in range(length):
                counts = Counter(seq[pos] for seq in seqs if pos < len(seq))
                total = sum(counts.values())
                freq = [
                    counts.get("A", 0) / total,
                    counts.get("C", 0) / total,
                    counts.get("G", 0) / total,
                    counts.get("T", 0) / total,
                ]
                matrix.append(freq)

            # Calculate min_score as 60% of max possible score
            max_score = sum(max(row) for row in matrix)
            min_score = max_score * 0.6

            pwms.append({
                "name": f"tir_{label}_{length}bp",
                "length": length,
                "count": len(seqs),
                "matrix": matrix,
                "min_score": round(min_score, 3),
            })

    output_path = os.path.join(output_dir, "tir_pwms.json")
    with open(output_path, "w") as f:
        json.dump(pwms, f, indent=2)

    print(f"  Built {len(pwms)} TIR PWMs, wrote to {output_path}")

    # Also save DR (direct repeat / TSD) data for reference
    cursor.execute("""
        SELECT upDR
        FROM starship_features
        WHERE upDR IS NOT NULL AND LENGTH(upDR) >= 3
          AND upDR NOT LIKE '%.%'
    """)
    drs = [row[0].strip().upper() for row in cursor.fetchall()
           if all(c in "ACGTN" for c in row[0].strip().upper())]
    print(f"  Found {len(drs)} direct repeats (TSDs) for reference")

    return pwms


def export_starship_reference_fasta(conn, output_path):
    """
    Extract all Starship nucleotide sequences from the ships table,
    join with family information, and write a FASTA reference file
    suitable for homology searches.

    FASTA headers: >{ship_id}|{family_name}|{sequence_length}bp

    Ships appearing in multiple joined_ships rows (different families) are
    deduplicated by choosing one family deterministically (MIN).  Ships
    without a family assignment are labelled 'unclassified'.

    Returns a dict with summary statistics.
    """
    cursor = conn.cursor()

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
    size_bins = Counter()
    total = 0

    SIZE_BIN_ORDER = ["<10kb", "10-50kb", "50-100kb", "100-200kb", "200-500kb", ">=500kb"]

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with open(output_path, "w") as fasta_out:
        for ship_id, sequence, sequence_length, family_name in rows:
            if family_name is None or family_name == "None" or family_name.strip() == "":
                family_name = "unclassified"

            seq_clean = sequence.strip()
            if sequence_length is None:
                sequence_length = len(seq_clean)

            header = f">{ship_id}|{family_name}|{sequence_length}bp"
            fasta_out.write(f"{header}\n")
            for i in range(0, len(seq_clean), 80):
                fasta_out.write(seq_clean[i : i + 80] + "\n")

            family_counts[family_name] += 1
            total += 1

            if sequence_length < 10000:
                size_bins["<10kb"] += 1
            elif sequence_length < 50000:
                size_bins["10-50kb"] += 1
            elif sequence_length < 100000:
                size_bins["50-100kb"] += 1
            elif sequence_length < 200000:
                size_bins["100-200kb"] += 1
            elif sequence_length < 500000:
                size_bins["200-500kb"] += 1
            else:
                size_bins[">=500kb"] += 1

    stats = {
        "total": total,
        "family_counts": OrderedDict(family_counts.most_common()),
        "size_distribution": OrderedDict(
            (b, size_bins.get(b, 0)) for b in SIZE_BIN_ORDER if size_bins.get(b, 0) > 0
        ),
        "output_path": output_path,
    }

    print(f"  Total sequences: {stats['total']}")
    print(f"  Output: {stats['output_path']}")
    print(f"  Per-family counts:")
    for family, count in stats["family_counts"].items():
        print(f"    {family}: {count}")
    print(f"  Size distribution:")
    for size_bin, count in stats["size_distribution"].items():
        print(f"    {size_bin}: {count}")

    return stats


def main():
    print(f"StarKIT HMM Profile Builder")
    print(f"Database: {DB_PATH}")
    print()

    conn = get_connection(DB_PATH)

    # Step 1: Extract captain sequences
    print("Step 1: Extracting captain sequences from Starbase...")
    families, all_captains = extract_captains_by_family(conn)
    unassigned = extract_captains_unassigned(conn)

    print(f"  Found {len(all_captains)} captains in {len(families)} families")
    print(f"  Found {len(unassigned)} unassigned captains")
    for family, seqs in sorted(families.items(), key=lambda x: -len(x[1])):
        print(f"    {family}: {len(seqs)}")
    print()

    # Include unassigned in all_captains for detection HMM
    all_for_detection = all_captains + unassigned

    # Step 2: Build detection HMM (all captains)
    print("Step 2: Building detection HMM (all captains)...")
    build_detection_hmm(all_for_detection, OUTPUT_HMM_DIR)
    print()

    # Step 3: Build per-family HMMs
    print("Step 3: Building per-family HMMs...")
    build_hmms_for_families(families, OUTPUT_FAMILY_DIR)
    print()

    # Step 4: Build TIR PWMs
    print("Step 4: Building TIR position-weight matrices...")
    build_tir_pwms(conn, OUTPUT_BOUNDARY_DIR)
    print()

    # Step 5: Export FASTA files for reference
    fasta_dir = os.path.join("scripts", "fasta_exports")
    os.makedirs(fasta_dir, exist_ok=True)

    write_fasta(all_for_detection, os.path.join(fasta_dir, "all_captains.fasta"))
    for family, seqs in families.items():
        write_fasta(seqs, os.path.join(fasta_dir, f"captains_{family}.fasta"))
    print(f"Step 5: Exported FASTA files to {fasta_dir}/")
    print()

    # Step 6: Export Starship nucleotide reference FASTA
    print("Step 6: Exporting Starship nucleotide reference FASTA...")
    export_starship_reference_fasta(conn, OUTPUT_STARSHIP_FASTA)
    print()

    conn.close()
    print()
    print("Done!")


if __name__ == "__main__":
    main()
