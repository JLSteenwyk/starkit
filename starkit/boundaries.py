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
import re
from typing import Optional

from .helpers import reverse_complement
from .models import CaptainHit, TIR
from .settings import (
    DEFAULT_UPSTREAM_SCAN,
    MIN_TIR_LENGTH,
    MIN_TIR_IDENTITY,
    TIR_SCAN_WINDOW,
    TIR_SEED_K,
)

logger = logging.getLogger(__name__)

_MYB_PATTERNS = re.compile(r"\bmyb\b|SANT|Myb_DNA-bind|PF13837|PF00249|IPR017930", re.IGNORECASE)


def _is_myb_feature(feature):
    for key in ("product", "note", "db_xref", "inference"):
        for val in feature.qualifiers.get(key, []):
            if _MYB_PATTERNS.search(val):
                return True
    return False

def find_myb_boundary(record, captain_start, captain_end, captain_strand, max_size):
    best_feature = None
    best_distance = -1
    for feature in record.features:
        if feature.type != "CDS" or not _is_myb_feature(feature):
            continue
        feat_start = int(feature.location.start)
        feat_end = int(feature.location.end)
        if captain_strand >= 0:
            if feat_start <= captain_start:
                continue
            distance = feat_end - captain_start
        else:
            if feat_end >= captain_end:
                continue
            distance = captain_end - feat_start
        if distance > max_size or distance <= 0:
            continue
        if distance > best_distance:
            best_distance = distance
            best_feature = feature
    if best_feature is None:
        return None, None
    feat_start = int(best_feature.location.start)
    feat_end = int(best_feature.location.end)
    boundary_position = feat_end if captain_strand >= 0 else feat_start
    return boundary_position, best_feature


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
        if p > five_prime + (min_size or 0)
        and (not max_size or p <= five_prime + max_size)
    ]
    if not downstream_positions:
        return five_prime, None

    three_prime = downstream_positions[0]
    return five_prime, three_prime


def find_dr_from_library(sequence, captain_start_in_seq, dr_library,
                         family_size_median, size_tolerance=0.6):
    """Search for known DR sequences from a family-specific library.

    For each known DR, finds occurrences upstream of the captain (5' end),
    then searches downstream for the same DR or close variants to define
    the 3' boundary. Candidates are scored by DR length and proximity to
    the family median size.

    Parameters
    ----------
    sequence : str
        Full contig sequence.
    captain_start_in_seq : int
        Captain gene start position (0-based).
    dr_library : set of str
        Known DR sequences for this family from Starbase.
    family_size_median : int
        Median element size for this family.
    size_tolerance : float
        Fraction of median size to allow as deviation (default 0.6 = ±60%).

    Returns
    -------
    tuple of (five_prime_pos, three_prime_pos, five_dr_seq, three_dr_seq)
        or (None, None, None, None) if no pair found.
    """
    if not dr_library:
        return None, None, None, None

    seq_upper = sequence.upper()
    min_size = max(5000, int(family_size_median * (1 - size_tolerance)))
    max_size = int(family_size_median * (1 + size_tolerance))

    up_start = max(0, captain_start_in_seq - DEFAULT_UPSTREAM_SCAN)
    upstream = seq_upper[up_start:captain_start_in_seq]

    best = None
    best_score = 0

    for dr_seq in dr_library:
        if len(dr_seq) < 4:
            continue

        # Find occurrences of this DR upstream of captain
        pos = 0
        while True:
            idx = upstream.find(dr_seq, pos)
            if idx == -1:
                break
            five_prime_pos = up_start + idx

            # Build a set of sequences to search for at the 3' end:
            # the same DR plus variants (±1bp from start/end)
            search_variants = {dr_seq}
            if len(dr_seq) >= 5:
                search_variants.add(dr_seq[1:])      # trim first base
                search_variants.add(dr_seq[:-1])      # trim last base
            # Also add other library DRs that share a core with this one
            for other_dr in dr_library:
                if len(other_dr) >= 4 and (
                    dr_seq.startswith(other_dr) or other_dr.startswith(dr_seq)
                    or dr_seq.endswith(other_dr) or other_dr.endswith(dr_seq)
                ):
                    search_variants.add(other_dr)

            search_start = five_prime_pos + min_size
            search_end = min(len(seq_upper), five_prime_pos + max_size)
            if search_start >= search_end:
                pos = idx + 1
                continue

            downstream = seq_upper[search_start:search_end]

            for variant in search_variants:
                dpos = 0
                while True:
                    didx = downstream.find(variant, dpos)
                    if didx == -1:
                        break
                    three_prime_pos = search_start + didx + len(variant)
                    element_size = three_prime_pos - five_prime_pos

                    size_ratio = element_size / family_size_median
                    match_bonus = 1.0 if variant == dr_seq else 0.8
                    score = len(dr_seq) * match_bonus * (1.0 - abs(1.0 - size_ratio) * 0.3)

                    if score > best_score:
                        best_score = score
                        best = (five_prime_pos, three_prime_pos, dr_seq, variant)

                    dpos = didx + 1

            pos = idx + 1

    if best is None:
        return None, None, None, None

    return best


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
    """De novo inverted repeat detection with mismatch-tolerant extension.

    Uses short k-mer seeds (6bp) to find candidate matches between the
    upstream boundary and the reverse-complement of the downstream boundary,
    then extends with mismatch tolerance (up to 30% mismatches). This catches
    imperfect TIRs that have diverged over time.

    Parameters
    ----------
    upstream_seq : str
        Sequence at the 5' boundary region.
    downstream_seq : str
        Sequence at the 3' boundary region.
    min_length : int
        Minimum total TIR length to report (default 8).
    min_identity : float
        Minimum fraction of matching bases in the TIR (default 0.70).
    scan_window : int
        How many bp to scan at each boundary (default 1000).
    """
    up = upstream_seq[:scan_window].upper()
    down_raw = downstream_seq[-scan_window:] if len(downstream_seq) >= scan_window else downstream_seq
    down = down_raw.upper()
    rc_down = reverse_complement(down)

    k = TIR_SEED_K
    if len(up) < k or len(rc_down) < k:
        return None, None

    # Build k-mer index from reverse-complemented downstream
    kmer_index = {}
    for i in range(len(rc_down) - k + 1):
        kmer = rc_down[i:i + k]
        if "N" in kmer:
            continue
        kmer_index.setdefault(kmer, []).append(i)

    best_match = None
    best_score = 0  # score = matches - mismatches

    for i in range(len(up) - k + 1):
        kmer = up[i:i + k]
        if "N" in kmer or kmer not in kmer_index:
            continue
        for j in kmer_index[kmer]:
            # Extend with mismatch tolerance
            matches = k
            mismatches = 0
            length = k
            max_mismatches_in_row = 0
            consecutive_mismatches = 0

            while i + length < len(up) and j + length < len(rc_down):
                a = up[i + length]
                b = rc_down[j + length]
                if a == "N" or b == "N":
                    break
                if a == b:
                    matches += 1
                    consecutive_mismatches = 0
                else:
                    mismatches += 1
                    consecutive_mismatches += 1
                    # Stop if 3+ mismatches in a row (left the TIR region)
                    if consecutive_mismatches >= 3:
                        # Back up to remove trailing mismatches
                        length -= (consecutive_mismatches - 1)
                        mismatches -= (consecutive_mismatches - 1)
                        break
                length += 1

                # Stop if identity drops below threshold
                if length > min_length and matches / length < min_identity:
                    length -= 1
                    if up[i + length] != rc_down[j + length]:
                        mismatches -= 1
                    else:
                        matches -= 1
                    break

            total_len = length
            if total_len < min_length:
                continue
            identity = matches / total_len if total_len > 0 else 0
            if identity < min_identity:
                continue

            # Score: prefer longer TIRs with higher identity
            score = matches + (identity * total_len * 0.5)
            if score > best_score:
                best_score = score
                best_match = (i, j, total_len, matches, mismatches)

    if best_match is None:
        return None, None

    i, j, length, matches, mismatches = best_match
    tir_seq_up = up[i:i + length]

    # Map j (index in rc_down) back to downstream_seq coordinates.
    # rc_down is the reverse complement of down (last scan_window bp).
    # Position j in rc_down corresponds to position (len(down)-1-j) in down,
    # counted from the start of the extracted window.
    down_window_start = max(0, len(downstream_seq) - scan_window)
    down_len = len(down)
    # The right TIR in the original sequence spans:
    #   start = down_window_start + (down_len - j - length)
    #   end   = down_window_start + (down_len - j)
    right_start = down_window_start + (down_len - j - length)
    right_end = down_window_start + (down_len - j)
    right_seq = downstream_seq[right_start:right_end].upper()

    tir_left = (i, i + length, tir_seq_up)
    tir_right = (right_start, right_end, right_seq)

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
    # Tier 2: Family-specific DR library scanning
    # ------------------------------------------------------------------- #
    if family_ref:
        family_name = captain_hit.hmm_name
        fam_data = family_ref.get(family_name, {})
        dr_library = fam_data.get("dr_library")
        family_size_median = fam_data.get("size_median", 65000)

        if dr_library:
            result = find_dr_from_library(
                seq, captain_start, set(dr_library),
                family_size_median,
            )
            five_pos, three_pos, five_dr, three_dr = result

            if five_pos is not None and three_pos is not None:
                start = five_pos
                end = three_pos
                truncated = (start <= 10000) or (end >= contig_len - 1000)
                tsd = five_dr

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
                    f"Tier 2 (DR library '{five_dr}'/'{three_dr}'): "
                    f"{record.id}:{start}-{end} ({end - start:,} bp)"
                )
                return dict(
                    start=start, end=end,
                    tir_left=tir_left_result, tir_right=tir_right_result,
                    tsd=tsd, truncated=truncated,
                    boundary_method="dr_motif",
                )

    # ------------------------------------------------------------------- #
    # Tier 2.5: MYB/SANT transcription factor at 3' boundary
    # ------------------------------------------------------------------- #
    search_max = max_size if max_size else 200000
    myb_pos, myb_feat = find_myb_boundary(
        record, captain_hit.start, captain_hit.end,
        captain_hit.strand, search_max,
    )
    if myb_pos is not None:
        if captain_hit.strand >= 0:
            start = max(0, captain_hit.start - 2000)
            end = min(contig_len, myb_pos)
        else:
            start = max(0, myb_pos)
            end = min(contig_len, captain_hit.end + 2000)
        element_size = end - start
        effective_min = min_size if min_size else 5000
        if element_size >= effective_min:
            truncated = (start <= 10000) or (end >= contig_len - 1000)
            logger.info(
                f"Tier 2.5 (MYB TF): {record.id}:{start}-{end} ({element_size:,} bp)"
            )
            return dict(
                start=start, end=end,
                tir_left=None, tir_right=None, tsd=None,
                truncated=truncated, boundary_method="myb_tf",
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

    if min_size and end - start < min_size:
        end = min(contig_len, start + min_size)
    if max_size and end - start > max_size:
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
        dr_library = fam_data.get("dr_library", [])
        # Check if any known DR is AT the boundary (within 5bp of start/end)
        for dr_seq in dr_library:
            if len(dr_seq) < 4:
                continue
            near_5 = seq[max(0, start - 5):start + len(dr_seq) + 5].upper()
            near_3 = seq[max(0, end - len(dr_seq) - 5):end + 5].upper()
            if dr_seq in near_5 and dr_seq in near_3:
                tsd = dr_seq
                break

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
