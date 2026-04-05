"""
Starship boundary detection module.

Uses a tiered approach:
  Tier 1: Homology-based boundaries (from reference Starship alignments)
  Tier 2: Family-specific DR motif scanning
  Tier 3: Captain position + family size prior (estimated)
"""

import json
import logging
import os
from typing import Optional

from .helpers import reverse_complement
from .models import CaptainHit, TIR
from .settings import (
    DEFAULT_UPSTREAM_SCAN,
    MIN_TIR_LENGTH,
    MIN_TIR_IDENTITY,
    TIR_SCAN_WINDOW,
)

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
#  Reference data loading
# --------------------------------------------------------------------------- #

def load_tir_pwms(pwm_file):
    """Load TIR position-weight matrices from a JSON file."""
    if not os.path.exists(pwm_file):
        return []
    try:
        with open(pwm_file) as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return []


def load_family_reference(ref_file):
    """Load family-specific DR motifs and size distributions."""
    if not os.path.exists(ref_file):
        return {}
    try:
        with open(ref_file) as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return {}


# --------------------------------------------------------------------------- #
#  DR / TSD motif scanning
# --------------------------------------------------------------------------- #

def find_motif_positions(sequence, motif):
    """Find all positions of a motif in a sequence (case-insensitive).

    Returns list of 0-based start positions.
    """
    positions = []
    seq_upper = sequence.upper()
    motif_upper = motif.upper()
    start = 0
    while True:
        idx = seq_upper.find(motif_upper, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def find_dr_pair(sequence, captain_start_in_seq, motif, min_size, max_size):
    """Search for a pair of DR motifs flanking a Starship.

    Scans upstream of the captain for the 5' DR, then searches downstream
    for a matching 3' DR that gives a valid element size.

    Returns (five_prime_pos, three_prime_pos) in sequence-local coords,
    or (None, None) if no valid pair is found.
    """
    positions = find_motif_positions(sequence, motif)
    if not positions:
        return None, None

    # 5' DR: closest motif upstream of captain (within 10kb)
    upstream_positions = [
        p for p in positions
        if p < captain_start_in_seq
        and captain_start_in_seq - p <= DEFAULT_UPSTREAM_SCAN
    ]
    if not upstream_positions:
        return None, None

    five_prime = max(upstream_positions)  # closest upstream

    # 3' DR: find a matching motif downstream that gives valid size
    downstream_positions = [
        p for p in positions
        if p > five_prime + min_size and p <= five_prime + max_size
    ]
    if not downstream_positions:
        return five_prime, None

    three_prime = downstream_positions[0]
    return five_prime, three_prime


# --------------------------------------------------------------------------- #
#  TIR detection
# --------------------------------------------------------------------------- #

def find_tirs_pwm(sequence, pwm_data):
    """Score sequence against TIR position-weight matrices."""
    if not pwm_data:
        return []

    seq_upper = sequence.upper()
    base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    results = []

    for pwm in pwm_data:
        matrix = pwm["matrix"]
        pwm_len = len(matrix)
        min_score = pwm.get("min_score", 0)

        for i in range(len(seq_upper) - pwm_len + 1):
            score = 0
            for j, row in enumerate(matrix):
                base = seq_upper[i + j]
                if base in base_idx:
                    score += row[base_idx[base]]
            if score >= min_score:
                results.append((i, i + pwm_len, score))

    return results


def find_tirs_denovo(upstream_seq, downstream_seq,
                     min_length=MIN_TIR_LENGTH,
                     min_identity=MIN_TIR_IDENTITY,
                     scan_window=TIR_SCAN_WINDOW):
    """De novo inverted repeat detection via k-mer seeding + extension."""
    up = upstream_seq[:scan_window].upper()
    down_raw = downstream_seq[-scan_window:] if len(downstream_seq) >= scan_window else downstream_seq
    down = down_raw.upper()
    rc_down = reverse_complement(down)

    k = min_length
    if len(up) < k or len(rc_down) < k:
        return None, None

    kmer_index = {}
    for i in range(len(rc_down) - k + 1):
        kmer = rc_down[i:i + k]
        if "N" in kmer:
            continue
        kmer_index.setdefault(kmer, []).append(i)

    best_match = None
    best_length = 0

    for i in range(len(up) - k + 1):
        kmer = up[i:i + k]
        if "N" in kmer or kmer not in kmer_index:
            continue
        for j in kmer_index[kmer]:
            length = k
            while (i + length < len(up) and j + length < len(rc_down)
                   and up[i + length] == rc_down[j + length]
                   and up[i + length] != "N"):
                length += 1
            if length > best_length:
                best_length = length
                best_match = (i, j, length)

    if best_match is None or best_length < min_length:
        return None, None

    i, j, length = best_match
    tir_seq = up[i:i + length]

    down_offset = len(downstream_seq) - scan_window + j if len(downstream_seq) >= scan_window else j

    tir_left = (i, i + length, tir_seq)
    tir_right = (down_offset, down_offset + length, reverse_complement(tir_seq))

    return tir_left, tir_right


# --------------------------------------------------------------------------- #
#  Main boundary detection (tiered)
# --------------------------------------------------------------------------- #

def define_boundaries(captain_hit, record, min_size, max_size,
                      pwm_data=None, family_ref=None, homology_hit=None):
    """Define Starship boundaries using a tiered approach.

    Tier 1: Use homology alignment coordinates if available.
    Tier 2: Scan for family-specific DR motifs.
    Tier 3: Captain position + family size prior (estimated).
    """
    contig_len = len(record.seq)
    seq = str(record.seq)
    captain_start = captain_hit.start

    # ------------------------------------------------------------------- #
    # Tier 1: Homology-based boundaries
    # ------------------------------------------------------------------- #
    if homology_hit is not None:
        start = homology_hit.start
        end = homology_hit.end
        truncated = (start <= 10000) or (end >= contig_len - 1000)

        tir_left, tir_right, tsd = _find_structural_features(
            seq, start, end, captain_start, pwm_data, family_ref, captain_hit,
        )

        logger.info(
            f"Tier 1 (homology): {record.id}:{start}-{end} "
            f"({end - start:,} bp)"
        )
        return dict(
            start=start, end=end,
            tir_left=tir_left, tir_right=tir_right, tsd=tsd,
            truncated=truncated, boundary_method="homology",
        )

    # ------------------------------------------------------------------- #
    # Tier 2: Family-specific DR scanning
    # ------------------------------------------------------------------- #
    if family_ref:
        family_name = captain_hit.hmm_name
        fam_data = family_ref.get(family_name, {})
        dr_motif = fam_data.get("dr_motif")

        if dr_motif:
            five_prime, three_prime = find_dr_pair(
                seq, captain_start, dr_motif, min_size, max_size,
            )

            if five_prime is not None and three_prime is not None:
                start = five_prime
                end = three_prime + len(dr_motif)
                truncated = (start <= 10000) or (end >= contig_len - 1000)
                tsd = seq[five_prime:five_prime + len(dr_motif)].upper()

                tir_left_result = None
                tir_right_result = None
                boundary_up = seq[start:start + TIR_SCAN_WINDOW]
                boundary_down = seq[max(0, end - TIR_SCAN_WINDOW):end]

                tl, tr = find_tirs_denovo(boundary_up, boundary_down)
                if tl is not None:
                    tir_left_result = TIR(
                        start=start + tl[0], end=start + tl[1], sequence=tl[2],
                    )
                if tr is not None:
                    down_base = max(0, end - TIR_SCAN_WINDOW)
                    tir_right_result = TIR(
                        start=down_base + tr[0], end=down_base + tr[1],
                        sequence=tr[2],
                    )

                logger.info(
                    f"Tier 2 (DR '{dr_motif}'): {record.id}:{start}-{end} "
                    f"({end - start:,} bp)"
                )
                return dict(
                    start=start, end=end,
                    tir_left=tir_left_result, tir_right=tir_right_result,
                    tsd=tsd, truncated=truncated,
                    boundary_method="dr_motif",
                )

    # ------------------------------------------------------------------- #
    # Tier 3: Captain position + family size prior (estimated)
    # ------------------------------------------------------------------- #
    median_size = 65000  # default across all families
    if family_ref:
        family_name = captain_hit.hmm_name
        fam_data = family_ref.get(family_name, {})
        median_size = fam_data.get("size_median", 65000)

    start = max(0, captain_start - 2000)
    end = min(contig_len, start + median_size)

    if end - start < min_size:
        end = min(contig_len, start + min_size)
    if end - start > max_size:
        end = start + max_size

    truncated = (start <= 10000) or (end >= contig_len - 1000)

    logger.info(
        f"Tier 3 (estimated, median={median_size:,}): "
        f"{record.id}:{start}-{end} ({end - start:,} bp)"
    )
    return dict(
        start=start, end=end,
        tir_left=None, tir_right=None, tsd=None,
        truncated=truncated, boundary_method="estimated",
    )


def _find_structural_features(seq, start, end, captain_start,
                              pwm_data, family_ref, captain_hit):
    """Try to find TIRs and TSDs within known boundaries."""
    tir_left = None
    tir_right = None
    tsd = None

    if family_ref and captain_hit:
        family_name = captain_hit.hmm_name
        fam_data = family_ref.get(family_name, {})
        dr_motif = fam_data.get("dr_motif")
        if dr_motif:
            flank_5 = seq[max(0, start - 20):start + 50].upper()
            flank_3 = seq[max(0, end - 50):end + 20].upper()
            if dr_motif in flank_5 and dr_motif in flank_3:
                tsd = dr_motif

    boundary_up = seq[start:start + TIR_SCAN_WINDOW]
    boundary_down = seq[max(0, end - TIR_SCAN_WINDOW):end]

    tl, tr = find_tirs_denovo(boundary_up, boundary_down)
    if tl is not None:
        tir_left = TIR(start=start + tl[0], end=start + tl[1], sequence=tl[2])
    if tr is not None:
        down_base = max(0, end - TIR_SCAN_WINDOW)
        tir_right = TIR(start=down_base + tr[0], end=down_base + tr[1], sequence=tr[2])

    return tir_left, tir_right, tsd


# Backwards compatibility alias
def find_tsd_motifs(sequence, motif="TTAC"):
    """Scan a DNA sequence for all occurrences of a motif."""
    return find_motif_positions(sequence, motif)
