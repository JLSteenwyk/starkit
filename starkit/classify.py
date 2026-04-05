"""
Starship family classification module.

Classifies Starship elements into families by scoring their captain proteins
against per-family HMM profiles using pyhmmer.
"""

import os
import logging
from typing import List, Optional, Tuple

import pyhmmer
import pyhmmer.easel
import pyhmmer.plan7

from starkit.exceptions import StarKITException
from starkit.models import CaptainHit, StarshipResult
from starkit.settings import FAMILY_HMM_DIR

logger = logging.getLogger(__name__)


def load_family_hmms(
    family_hmm_dir: str,
) -> List[pyhmmer.plan7.HMM]:
    """Load per-family HMM profiles from a directory.

    Parameters
    ----------
    family_hmm_dir : str
        Path to directory containing ``.hmm`` family profile files.

    Returns
    -------
    list of pyhmmer.plan7.HMM
        Loaded HMM profile objects, one per family.

    Raises
    ------
    StarKITException
        If no ``.hmm`` files are found in *family_hmm_dir*.
    """
    hmms: List[pyhmmer.plan7.HMM] = []

    if not os.path.isdir(family_hmm_dir):
        raise StarKITException(
            f"Family HMM directory not found: {family_hmm_dir}"
        )

    hmm_files = sorted(
        f for f in os.listdir(family_hmm_dir) if f.endswith(".hmm")
    )

    for filename in hmm_files:
        filepath = os.path.join(family_hmm_dir, filename)
        with pyhmmer.plan7.HMMFile(filepath) as hmm_file:
            for hmm in hmm_file:
                hmms.append(hmm)

    if not hmms:
        raise StarKITException(
            f"No family HMM profiles found in {family_hmm_dir}"
        )

    logger.info(f"Loaded {len(hmms)} family HMM profile(s)")
    return hmms


def classify_captain(
    captain_hit: CaptainHit,
    records: list,
    family_hmms: List[pyhmmer.plan7.HMM],
) -> Tuple[str, float]:
    """Score a captain protein against each family HMM.

    Extracts the captain's protein sequence from the genome records and runs
    pyhmmer.hmmsearch with the family HMMs against that single protein.

    Parameters
    ----------
    captain_hit : CaptainHit
        The captain gene hit to classify.
    records : list of Bio.SeqRecord.SeqRecord
        Annotated genome records containing CDS features with translations.
    family_hmms : list of pyhmmer.plan7.HMM
        Family HMM profiles to score against.

    Returns
    -------
    tuple of (str, float)
        ``(family_name, score)`` for the best-scoring family, or
        ``("unclassified", 0.0)`` if no family HMM hits the captain.
    """
    # Find the captain protein sequence from the records
    protein_seq: Optional[str] = None
    protein_id = captain_hit.protein_id

    for record in records:
        if record.id != captain_hit.contig_id:
            continue
        for feature in record.features:
            if feature.type != "CDS":
                continue
            # Match by location to the captain feature
            if (
                int(feature.location.start) == captain_hit.start
                and int(feature.location.end) == captain_hit.end
            ):
                translations = feature.qualifiers.get("translation", [])
                if translations and translations[0]:
                    protein_seq = translations[0]
                    break
        if protein_seq is not None:
            break

    if protein_seq is None:
        logger.warning(
            f"Could not extract protein sequence for captain {protein_id}"
        )
        return ("unclassified", 0.0)

    # Build a single-sequence database for hmmsearch
    alphabet = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(
        name=protein_id.encode(),
        sequence=protein_seq,
    ).digitize(alphabet)

    # Score against all family HMMs
    best_family = "unclassified"
    best_score = 0.0

    for top_hits in pyhmmer.hmmsearch(family_hmms, [seq]):
        family_name = top_hits.query_name.decode()
        for hit in top_hits:
            if hit.included and hit.score > best_score:
                best_score = hit.score
                best_family = family_name

    if best_family != "unclassified":
        logger.info(
            f"Captain {protein_id} classified as {best_family} "
            f"(score={best_score:.1f})"
        )
    else:
        logger.info(f"Captain {protein_id} could not be classified")

    return (best_family, best_score)


def classify_starships(
    starship_results: List[StarshipResult],
    records: list,
    family_hmm_dir: str = FAMILY_HMM_DIR,
) -> List[StarshipResult]:
    """Classify each Starship element into a family.

    This is the main entry point for family classification. For each
    StarshipResult, the captain protein is scored against per-family HMM
    profiles and the best-matching family is assigned.

    Parameters
    ----------
    starship_results : list of StarshipResult
        Starship predictions to classify.
    records : list of Bio.SeqRecord.SeqRecord
        Annotated genome records containing CDS features with translations.
    family_hmm_dir : str, optional
        Path to directory containing family HMM profiles.
        Defaults to :data:`starkit.settings.FAMILY_HMM_DIR`.

    Returns
    -------
    list of StarshipResult
        The input list, with ``captain_family`` and ``family_score``
        updated in place.
    """
    family_hmms = load_family_hmms(family_hmm_dir)

    for result in starship_results:
        family_name, score = classify_captain(
            result.captain, records, family_hmms
        )
        result.captain_family = family_name
        result.family_score = score

    classified = sum(
        1 for r in starship_results if r.captain_family != "unclassified"
    )
    logger.info(
        f"Classified {classified}/{len(starship_results)} Starship(s) "
        f"into families"
    )

    return starship_results
