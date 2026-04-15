"""
Confidence scoring module.

Assigns a composite confidence score and an evidence level to each
Starship prediction based on four components: captain gene detection,
homology to known Starships, boundary quality, and structural features.
"""

import logging
import math
from typing import List

from starkit.models import EvidenceLevel, StarshipResult

logger = logging.getLogger(__name__)

# Component weights
W_CAPTAIN = 0.3
W_HOMOLOGY = 0.3
W_BOUNDARY = 0.2
W_STRUCTURAL = 0.2


def compute_evidence_level(result: StarshipResult) -> EvidenceLevel:
    """Determine the evidence level for a Starship prediction.

    HIGH:   captain + homology OR captain + structural features OR
            captain + MYB-at-opposite-edge (strong canonical marker)
    MEDIUM: captain alone, strong homology alone, or partial structural
    LOW:    weak homology only, no captain, no structural features,
            captain orientation violation (demoted), Tier 3 boundaries
    """
    # Tier 3 (estimated) boundaries are unreliable — cap at LOW
    if result.boundary_method == "estimated":
        return EvidenceLevel.LOW

    has_captain = result.captain.evalue < 1.0  # real captain, not placeholder
    has_homology = result.homology_identity >= 0.80 and result.homology_coverage >= 0.50
    has_both_tirs = result.tir_left is not None and result.tir_right is not None
    has_any_tir = result.tir_left is not None or result.tir_right is not None
    has_tsd = result.tsd is not None
    has_structural = has_any_tir or has_tsd
    has_myb_edge = getattr(result, "myb_at_opposite_edge", False)
    has_duf3723 = getattr(result, "has_duf3723_marker", False)
    orientation_bad = getattr(result, "captain_orientation_flag", False)

    # A canonical MYB at the opposite edge is a strong structural marker on
    # its own and can rescue an element to HIGH even without TIRs/TSD.
    strong_marker = has_myb_edge or has_duf3723

    # Orientation violation demotes the prediction by one level.
    if orientation_bad:
        if has_captain and has_homology and strong_marker:
            return EvidenceLevel.MEDIUM
        if has_captain:
            return EvidenceLevel.MEDIUM
        return EvidenceLevel.LOW

    if has_captain and has_homology:
        return EvidenceLevel.HIGH
    if has_captain and has_structural:
        return EvidenceLevel.HIGH
    if has_captain and strong_marker:
        return EvidenceLevel.HIGH
    if has_captain:
        return EvidenceLevel.MEDIUM
    if has_homology and result.homology_coverage >= 0.80:
        return EvidenceLevel.MEDIUM
    if has_homology:
        return EvidenceLevel.MEDIUM
    return EvidenceLevel.LOW


def compute_confidence_score(result: StarshipResult) -> float:
    """Compute a composite confidence score in [0, 1].

    Four weighted components plus a marker-gene bonus:

    - Captain (0.3): -log10(evalue) / 50, capped at 1.0.
      A strong HMM hit (e-value 1e-100) scores 1.0.
      No captain (evalue=1.0) scores 0.0.

    - Homology (0.3): identity * coverage.
      100% identity at 100% coverage = 1.0.
      No homology match = 0.0.

    - Boundary (0.2): quality of boundary determination.
      homology-based = 1.0, DR-motif = 0.7, estimated = 0.3.

    - Structural (0.2): TIR and TSD presence.
      Both TIRs + TSD = 1.0, TIRs only = 0.6, TSD only = 0.4,
      any single TIR = 0.3, none = 0.0.

    Marker bonuses (applied after the weighted sum):
      + MYB/SANT at opposite edge with inverted orientation: +0.10
      + MYB/SANT present (but not at the canonical edge): +0.03
      + DUF3723 present: +0.05

    Penalties:
      - Captain orientation flag set: -0.15
    """
    # Captain component
    if result.captain.evalue < 1.0:
        captain_score = min(1.0, -math.log10(max(result.captain.evalue, 1e-300)) / 50)
    else:
        captain_score = 0.0

    # Homology component
    homology_score = result.homology_identity * result.homology_coverage

    # Boundary component
    boundary_map = {"homology": 1.0, "dr_motif": 0.7, "myb_tf": 0.55, "estimated": 0.3}
    boundary_score = boundary_map.get(result.boundary_method, 0.3)

    # Structural component
    has_both_tirs = result.tir_left is not None and result.tir_right is not None
    has_any_tir = result.tir_left is not None or result.tir_right is not None
    has_tsd = result.tsd is not None

    if has_both_tirs and has_tsd:
        structural_score = 1.0
    elif has_both_tirs:
        structural_score = 0.6
    elif has_tsd:
        structural_score = 0.4
    elif has_any_tir:
        structural_score = 0.3
    else:
        structural_score = 0.0

    score = (
        W_CAPTAIN * captain_score
        + W_HOMOLOGY * homology_score
        + W_BOUNDARY * boundary_score
        + W_STRUCTURAL * structural_score
    )

    # Marker bonuses for Starship-specific auxiliary genes
    if getattr(result, "myb_at_opposite_edge", False):
        score += 0.10
    elif getattr(result, "has_myb_marker", False):
        score += 0.03
    if getattr(result, "has_duf3723_marker", False):
        score += 0.05

    # Penalty for orientation violations
    if getattr(result, "captain_orientation_flag", False):
        score -= 0.15

    score = max(0.0, min(1.0, score))
    return round(score, 4)


def score_starships(
    starship_results: List[StarshipResult],
) -> List[StarshipResult]:
    """Score confidence and evidence level for each Starship prediction.

    Updates each StarshipResult in place.
    """
    for result in starship_results:
        result.evidence_level = compute_evidence_level(result)
        result.confidence_score = compute_confidence_score(result)

    level_counts = {}
    for result in starship_results:
        level = result.evidence_level.value
        level_counts[level] = level_counts.get(level, 0) + 1

    logger.info(
        f"Scored {len(starship_results)} Starship(s): "
        + ", ".join(f"{k}={v}" for k, v in sorted(level_counts.items()))
    )

    return starship_results
