"""
Confidence scoring module.

Assigns a composite confidence score and an evidence level to each
Starship prediction based on the structural features detected (captain
gene, TIRs, TSD) and element size.
"""

import logging
import math
from typing import List

from starkit.models import EvidenceLevel, StarshipResult
from starkit.settings import (
    CONFIDENCE_CAPTAIN_WEIGHT,
    CONFIDENCE_TIR_WEIGHT,
    CONFIDENCE_TSD_WEIGHT,
    CONFIDENCE_SIZE_WEIGHT,
    MEDIAN_STARSHIP_SIZE,
    STARSHIP_SIZE_STD,
)

logger = logging.getLogger(__name__)


def compute_evidence_level(result: StarshipResult) -> EvidenceLevel:
    """Determine the evidence level for a Starship prediction.

    Parameters
    ----------
    result : StarshipResult
        A Starship prediction with structural features already detected.

    Returns
    -------
    EvidenceLevel
        - ``HIGH``   -- captain + both TIRs + TSD
        - ``MEDIUM`` -- captain + partial evidence (one TIR, or TSD only,
          or both TIRs without TSD)
        - ``LOW``    -- captain only (no TIRs, no TSD)
    """
    has_both_tirs = result.tir_left is not None and result.tir_right is not None
    has_any_tir = result.tir_left is not None or result.tir_right is not None
    has_tsd = result.tsd is not None

    if has_both_tirs and has_tsd:
        return EvidenceLevel.HIGH
    elif has_any_tir or has_tsd:
        return EvidenceLevel.MEDIUM
    else:
        return EvidenceLevel.LOW


def compute_confidence_score(result: StarshipResult) -> float:
    """Compute a composite confidence score in [0, 1] for a Starship.

    The score is a weighted sum of four components:

    - **Captain** (weight 0.4): ``-log10(evalue) / 50``, capped at 1.0.
    - **TIR** (weight 0.3): 1.0 if both TIRs found, 0.5 if one, 0.0 if none.
    - **TSD** (weight 0.2): 1.0 if a target-site duplication is present.
    - **Size** (weight 0.1): Gaussian centred on the median Starship size.

    Parameters
    ----------
    result : StarshipResult
        A Starship prediction with structural features already detected.

    Returns
    -------
    float
        Confidence score rounded to four decimal places.
    """
    # Captain score: -log10(evalue) / 50, capped at 1.0
    captain_score = min(
        1.0,
        -math.log10(max(result.captain.evalue, 1e-300)) / 50,
    )

    # TIR score
    if result.tir_left is not None and result.tir_right is not None:
        tir_score = 1.0
    elif result.tir_left is not None or result.tir_right is not None:
        tir_score = 0.5
    else:
        tir_score = 0.0

    # TSD score
    tsd_score = 1.0 if result.tsd is not None else 0.0

    # Size score: Gaussian centred on median Starship size
    size_diff = abs(result.size - MEDIAN_STARSHIP_SIZE)
    size_score = math.exp(-(size_diff ** 2) / (2 * STARSHIP_SIZE_STD ** 2))

    score = (
        CONFIDENCE_CAPTAIN_WEIGHT * captain_score
        + CONFIDENCE_TIR_WEIGHT * tir_score
        + CONFIDENCE_TSD_WEIGHT * tsd_score
        + CONFIDENCE_SIZE_WEIGHT * size_score
    )

    return round(score, 4)


def score_starships(
    starship_results: List[StarshipResult],
) -> List[StarshipResult]:
    """Score confidence and evidence level for each Starship prediction.

    This is the main entry point for confidence scoring. It updates each
    :class:`StarshipResult` in place with its computed ``evidence_level``
    and ``confidence_score``.

    Parameters
    ----------
    starship_results : list of StarshipResult
        Starship predictions to score.

    Returns
    -------
    list of StarshipResult
        The input list, updated in place.
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
