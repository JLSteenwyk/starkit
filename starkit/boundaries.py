"""
Boundary detection for Starship transposable elements.

Identifies the 5' and 3' boundaries of Starship elements using:
- Target site duplication (TSD) motifs (consensus: TTAC)
- Terminal inverted repeats (TIRs) via PWM scoring or de novo detection
- Captain gene proximity constraints
"""

import json
import logging
import os
import re

from starkit.models import CaptainHit, TIR
from starkit.settings import (
    DEFAULT_UPSTREAM_SCAN,
    DEFAULT_TSD_MOTIF,
    MIN_TIR_LENGTH,
    MIN_TIR_IDENTITY,
    TIR_SCAN_WINDOW,
)

try:
    from starkit.helpers import reverse_complement
except ImportError:

    def reverse_complement(seq):
        """Return the reverse complement of a DNA sequence."""
        complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(complement)[::-1]


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# TSD motif scanning
# ---------------------------------------------------------------------------


def find_tsd_motifs(sequence, motif=DEFAULT_TSD_MOTIF):
    """Scan a DNA sequence for all occurrences of a TSD motif.

    Parameters
    ----------
    sequence : str
        DNA sequence to scan (case-insensitive).
    motif : str, optional
        Short motif to search for (default ``"TTAC"``).

    Returns
    -------
    list of int
        0-based positions where *motif* starts in *sequence*.
    """
    sequence_upper = sequence.upper()
    motif_upper = motif.upper()
    positions = []
    start = 0
    while True:
        idx = sequence_upper.find(motif_upper, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


# ---------------------------------------------------------------------------
# PWM-based TIR detection
# ---------------------------------------------------------------------------


def find_tirs_pwm(sequence, pwm_data):
    """Score a sequence against TIR position-weight matrices.

    Parameters
    ----------
    sequence : str
        DNA sequence to scan.
    pwm_data : list of dict or None
        Each dict must contain:
        - ``"name"``  : str – identifier for the PWM
        - ``"matrix"``: list of list – one inner list per position,
          columns ordered [A, C, G, T]
        - ``"min_score"``: float – minimum score threshold

    Returns
    -------
    list of tuple
        ``(start, end, score)`` for every position whose score meets or
        exceeds the PWM threshold.  *start* is 0-based, *end* is
        exclusive.
    """
    if not pwm_data:
        return []

    base_index = {"A": 0, "C": 1, "G": 2, "T": 3}
    seq_upper = sequence.upper()
    hits = []

    for pwm in pwm_data:
        matrix = pwm["matrix"]
        pwm_len = len(matrix)
        min_score = pwm["min_score"]

        if pwm_len > len(seq_upper):
            continue

        for i in range(len(seq_upper) - pwm_len + 1):
            score = 0.0
            for j, row in enumerate(matrix):
                base = seq_upper[i + j]
                col = base_index.get(base)
                if col is None:
                    # ambiguous base – treat as zero contribution
                    continue
                score += row[col]

            if score >= min_score:
                hits.append((i, i + pwm_len, score))

    return hits


# ---------------------------------------------------------------------------
# De novo TIR detection
# ---------------------------------------------------------------------------


def find_tirs_denovo(
    upstream_seq,
    downstream_seq,
    min_length=MIN_TIR_LENGTH,
    min_identity=MIN_TIR_IDENTITY,
    scan_window=TIR_SCAN_WINDOW,
):
    """De novo TIR detection using k-mer sliding window.

    Starship TIRs are inverted repeats: the downstream TIR is the
    reverse complement of the upstream TIR.  We slide a window of
    *min_length* bp across the first *scan_window* bp of
    *upstream_seq* and compare against the reverse complement of the
    last *scan_window* bp of *downstream_seq*.

    Parameters
    ----------
    upstream_seq : str
        Sequence around the candidate 5' boundary.
    downstream_seq : str
        Sequence around the candidate 3' boundary.
    min_length : int
        Minimum TIR length (k-mer size for initial seeding).
    min_identity : float
        Minimum fraction of identical bases between matched TIRs.
    scan_window : int
        Number of bases at each boundary to scan.

    Returns
    -------
    tuple
        ``(tir_left, tir_right)`` where each is
        ``(start_offset, end_offset, sequence)`` or ``(None, None)``
        if no qualifying pair is found.
    """
    up = upstream_seq[:scan_window].upper()
    down = downstream_seq[-scan_window:].upper() if len(downstream_seq) >= scan_window else downstream_seq.upper()
    down_rc = reverse_complement(down)

    if len(up) < min_length or len(down_rc) < min_length:
        return (None, None)

    # Build a set of k-mers from the downstream reverse complement for
    # fast lookup, mapping kmer -> list of positions in down_rc.
    kmer_index = {}
    for j in range(len(down_rc) - min_length + 1):
        kmer = down_rc[j : j + min_length]
        if "N" in kmer:
            continue
        kmer_index.setdefault(kmer, []).append(j)

    best_left = None
    best_right = None
    best_identity = 0.0
    best_match_len = 0

    for i in range(len(up) - min_length + 1):
        seed = up[i : i + min_length]
        if "N" in seed:
            continue
        if seed not in kmer_index:
            continue

        for j in kmer_index[seed]:
            # Extend the match rightward as far as bases are
            # unambiguous and identical.
            match_len = min_length
            while (
                i + match_len < len(up)
                and j + match_len < len(down_rc)
            ):
                base_up = up[i + match_len]
                base_dn = down_rc[j + match_len]
                if base_up == base_dn and base_up != "N":
                    match_len += 1
                else:
                    break

            # Compute identity over the extended region
            matches = sum(
                1
                for k in range(match_len)
                if up[i + k] == down_rc[j + k]
            )
            identity = matches / match_len if match_len > 0 else 0.0

            if identity >= min_identity and match_len > best_match_len:
                best_match_len = match_len
                best_identity = identity

                tir_seq_left = up[i : i + match_len]
                tir_seq_right = reverse_complement(down_rc[j : j + match_len])

                best_left = (i, i + match_len, tir_seq_left)
                # Map j back to the original downstream coordinate.
                # down_rc[j] corresponds to downstream position
                # (len(down) - 1 - j) counting from the start of the
                # scanned window, but since down is the *last*
                # scan_window of downstream_seq we compute the offset
                # relative to the start of downstream_seq.
                down_window_start = max(0, len(downstream_seq) - scan_window)
                # In the original downstream orientation the TIR runs
                # from the position that maps to the *end* of the
                # reverse-complement match back to its start.
                right_end_in_down = len(down) - j
                right_start_in_down = len(down) - (j + match_len)
                best_right = (
                    down_window_start + right_start_in_down,
                    down_window_start + right_end_in_down,
                    tir_seq_right,
                )

    if best_left is not None:
        return (best_left, best_right)

    return (None, None)


# ---------------------------------------------------------------------------
# PWM loading
# ---------------------------------------------------------------------------


def load_tir_pwms(pwm_file):
    """Load TIR position-weight matrices from a JSON file.

    The JSON file should contain a list of objects each with keys
    ``"name"``, ``"matrix"``, and ``"min_score"``.

    Parameters
    ----------
    pwm_file : str
        Path to the JSON file.

    Returns
    -------
    list of dict
        Parsed PWM data, or an empty list if the file does not exist
        or cannot be parsed.
    """
    if not os.path.isfile(pwm_file):
        logger.debug("PWM file not found: %s", pwm_file)
        return []
    try:
        with open(pwm_file, "r") as fh:
            data = json.load(fh)
        if not isinstance(data, list):
            logger.warning("PWM file does not contain a list: %s", pwm_file)
            return []
        return data
    except (json.JSONDecodeError, OSError) as exc:
        logger.warning("Failed to load PWM file %s: %s", pwm_file, exc)
        return []


# ---------------------------------------------------------------------------
# Main boundary definition
# ---------------------------------------------------------------------------


def define_boundaries(captain_hit, record, min_size, max_size, pwm_data=None):
    """Predict the 5' and 3' boundaries of a Starship element.

    Parameters
    ----------
    captain_hit : CaptainHit
        Detected captain gene.
    record : Bio.SeqRecord.SeqRecord
        SeqRecord for the contig that contains *captain_hit*.
    min_size : int
        Minimum allowed element size (bp).
    max_size : int
        Maximum allowed element size (bp).
    pwm_data : list of dict or None
        Pre-loaded TIR PWM data (see :func:`load_tir_pwms`).

    Returns
    -------
    dict
        Keys: ``start``, ``end``, ``tir_left``, ``tir_right``, ``tsd``,
        ``truncated``.
    """
    contig_seq = str(record.seq).upper()
    contig_len = len(contig_seq)

    captain_start = captain_hit.start
    captain_end = captain_hit.end
    strand = captain_hit.strand  # +1 or -1

    # ----------------------------------------------------------------
    # Determine scan direction: captain is always at the 5' end of the
    # element, so the element extends *downstream* of the captain on
    # the + strand and *upstream* on the - strand.
    # ----------------------------------------------------------------
    if strand >= 0:
        # Captain is on + strand -> 5' boundary is upstream of captain_start
        scan_5p_start = max(0, captain_start - DEFAULT_UPSTREAM_SCAN)
        scan_5p_end = captain_start
        scan_3p_start = captain_end
        scan_3p_end = min(contig_len, captain_end + max_size)
    else:
        # Captain is on - strand -> 5' boundary is downstream of captain_end
        scan_5p_start = captain_end
        scan_5p_end = min(contig_len, captain_end + DEFAULT_UPSTREAM_SCAN)
        scan_3p_start = max(0, captain_start - max_size)
        scan_3p_end = captain_start

    # ----------------------------------------------------------------
    # Step (a): Find candidate 5' boundary via TSD motif
    # ----------------------------------------------------------------
    region_5p = contig_seq[scan_5p_start:scan_5p_end]
    tsd_5p_positions = find_tsd_motifs(region_5p)

    candidate_5p = None
    tsd_seq = None

    if strand >= 0:
        # We want the closest TTAC *upstream* of the captain, i.e. the
        # rightmost hit in the 5' scan region.
        if tsd_5p_positions:
            best_offset = tsd_5p_positions[-1]  # closest to captain
            candidate_5p = scan_5p_start + best_offset
    else:
        # On - strand the 5' boundary is downstream of captain_end; we
        # want the leftmost TTAC in that region.
        if tsd_5p_positions:
            best_offset = tsd_5p_positions[0]
            candidate_5p = scan_5p_start + best_offset

    # ----------------------------------------------------------------
    # Step (b): Find candidate 3' boundary via matching TSD
    # ----------------------------------------------------------------
    region_3p = contig_seq[scan_3p_start:scan_3p_end]
    tsd_3p_positions = find_tsd_motifs(region_3p)

    candidate_3p = None

    if candidate_5p is not None and tsd_3p_positions:
        for offset_3p in tsd_3p_positions:
            abs_3p = scan_3p_start + offset_3p
            if strand >= 0:
                element_size = abs_3p - candidate_5p
            else:
                element_size = candidate_5p - abs_3p

            if min_size <= element_size <= max_size:
                candidate_3p = abs_3p
                tsd_seq = DEFAULT_TSD_MOTIF
                # Take the first valid match (closest to captain on 3' side)
                break

    # Fallback: if no matching TSD found, use size-based heuristic
    if candidate_5p is None:
        candidate_5p = scan_5p_start  # contig edge or scan limit

    if candidate_3p is None:
        # Place 3' boundary at median expected distance, clipped to contig
        if strand >= 0:
            candidate_3p = min(contig_len, candidate_5p + min_size)
        else:
            candidate_3p = max(0, candidate_5p - min_size)
        tsd_seq = None

    # Ensure start < end regardless of strand
    element_start = min(candidate_5p, candidate_3p)
    element_end = max(candidate_5p, candidate_3p)

    # ----------------------------------------------------------------
    # Step (c): TIR detection
    # ----------------------------------------------------------------
    tir_left = None
    tir_right = None

    # Extract boundary regions for TIR scanning
    boundary_5p_start = max(0, element_start)
    boundary_5p_end = min(contig_len, element_start + TIR_SCAN_WINDOW)
    boundary_3p_start = max(0, element_end - TIR_SCAN_WINDOW)
    boundary_3p_end = min(contig_len, element_end)

    upstream_region = contig_seq[boundary_5p_start:boundary_5p_end]
    downstream_region = contig_seq[boundary_3p_start:boundary_3p_end]

    # Try PWM-based detection first
    pwm_hits_left = find_tirs_pwm(upstream_region, pwm_data)
    pwm_hits_right = find_tirs_pwm(downstream_region, pwm_data)

    if pwm_hits_left and pwm_hits_right:
        # Take highest-scoring hit from each side
        best_left = max(pwm_hits_left, key=lambda h: h[2])
        best_right = max(pwm_hits_right, key=lambda h: h[2])

        tir_left = TIR(
            start=boundary_5p_start + best_left[0],
            end=boundary_5p_start + best_left[1],
            sequence=upstream_region[best_left[0] : best_left[1]],
        )
        tir_right = TIR(
            start=boundary_3p_start + best_right[0],
            end=boundary_3p_start + best_right[1],
            sequence=downstream_region[best_right[0] : best_right[1]],
        )
    else:
        # Fallback to de novo detection
        denovo_result = find_tirs_denovo(
            upstream_region,
            downstream_region,
            min_length=MIN_TIR_LENGTH,
            min_identity=MIN_TIR_IDENTITY,
            scan_window=TIR_SCAN_WINDOW,
        )
        left_hit, right_hit = denovo_result

        if left_hit is not None and right_hit is not None:
            tir_left = TIR(
                start=boundary_5p_start + left_hit[0],
                end=boundary_5p_start + left_hit[1],
                sequence=left_hit[2],
            )
            tir_right = TIR(
                start=boundary_3p_start + right_hit[0],
                end=boundary_3p_start + right_hit[1],
                sequence=right_hit[2],
            )

    # ----------------------------------------------------------------
    # Step (d): Validation – truncation flags
    # ----------------------------------------------------------------
    truncated = False
    # Captain within 10kb of contig start
    if captain_start < DEFAULT_UPSTREAM_SCAN:
        truncated = True
    # Predicted end within 1kb of contig end
    if element_end > contig_len - 1000:
        truncated = True
    # Predicted start at contig start
    if element_start == 0:
        truncated = True

    # ----------------------------------------------------------------
    # Step (e): Return result
    # ----------------------------------------------------------------
    return {
        "start": element_start,
        "end": element_end,
        "tir_left": tir_left,
        "tir_right": tir_right,
        "tsd": tsd_seq,
        "truncated": truncated,
    }
