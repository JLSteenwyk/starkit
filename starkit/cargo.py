"""
Cargo gene extraction and annotation module.

Cargo genes are the functional payload carried within Starship transposable
elements. This module identifies CDS features within Starship boundaries
and tags known auxiliary gene types using Pfam HMM domain searches.
"""

import logging
import os
from typing import Dict, List, Optional

import pyhmmer
import pyhmmer.easel
import pyhmmer.plan7

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from starkit.models import CargoGene
from starkit.settings import DATA_DIR

logger = logging.getLogger(__name__)

AUXILIARY_HMM_DIR = os.path.join(DATA_DIR, "auxiliary_hmms")

# Map HMM profile names to cargo tags
_HMM_TAG_MAP = {
    "DUF3435": "tyr",
    "NB-ARC": "nlr",
    "Ferric_reduct": "fre",
    "Vip3A_N": "d37",          # DUF3723 family
    "DUF3723": "d37",
    "Aminotran_1_2": "plp",
    "Myb_DNA-bind_4": "myb",
    "Myb_DNA-binding": "myb",
}

# Default e-value threshold for auxiliary gene detection
_AUX_EVALUE = 1e-5


def load_auxiliary_hmms() -> List[pyhmmer.plan7.HMM]:
    """Load auxiliary gene HMM profiles from the data directory."""
    if not os.path.isdir(AUXILIARY_HMM_DIR):
        return []

    hmms = []
    for fname in sorted(os.listdir(AUXILIARY_HMM_DIR)):
        if not fname.endswith(".hmm"):
            continue
        path = os.path.join(AUXILIARY_HMM_DIR, fname)
        try:
            with pyhmmer.plan7.HMMFile(path) as hf:
                for hmm in hf:
                    hmms.append(hmm)
        except Exception as e:
            logger.warning(f"Could not load auxiliary HMM {fname}: {e}")

    if hmms:
        logger.info(f"Loaded {len(hmms)} auxiliary gene HMM profile(s)")
    return hmms


def tag_proteins_by_hmm(
    proteins: List[tuple],
    hmm_profiles: List[pyhmmer.plan7.HMM],
    evalue_threshold: float = _AUX_EVALUE,
) -> Dict[str, str]:
    """Search proteins against auxiliary HMMs and return a tag map.

    Parameters
    ----------
    proteins : list of (protein_id, protein_sequence_str)
    hmm_profiles : list of pyhmmer HMM objects
    evalue_threshold : float

    Returns
    -------
    dict mapping protein_id -> tag string
    """
    if not proteins or not hmm_profiles:
        return {}

    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    for pid, seq in proteins:
        try:
            sequences.append(
                pyhmmer.easel.TextSequence(
                    name=pid.encode() if isinstance(pid, str) else pid,
                    sequence=seq,
                ).digitize(alphabet)
            )
        except Exception:
            continue

    if not sequences:
        return {}

    tags = {}  # protein_id -> (tag, evalue) — keep best hit per protein

    for top_hits in pyhmmer.hmmsearch(hmm_profiles, sequences):
        query_name = top_hits.query.name
        hmm_name = query_name.decode() if isinstance(query_name, bytes) else str(query_name)
        tag = _HMM_TAG_MAP.get(hmm_name, "")
        if not tag:
            continue

        for hit in top_hits:
            if hit.included and hit.evalue <= evalue_threshold:
                hit_name = hit.name
                pid = hit_name.decode() if isinstance(hit_name, bytes) else str(hit_name)
                # Keep the best (lowest e-value) tag per protein
                if pid not in tags or hit.evalue < tags[pid][1]:
                    tags[pid] = (tag, hit.evalue)

    return {pid: tag for pid, (tag, _) in tags.items()}


def extract_cargo(
    record: SeqRecord,
    start: int,
    end: int,
    captain_feature: Optional[SeqFeature],
    hmm_tags: Optional[Dict[str, str]] = None,
) -> List[CargoGene]:
    """Collect all CDS features from *record* that fall within [start, end].

    The captain gene is excluded by comparing feature locations.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        Annotated genome record (contig) containing CDS features.
    start : int
        Start coordinate of the Starship element boundary.
    end : int
        End coordinate of the Starship element boundary.
    captain_feature : Bio.SeqFeature.SeqFeature or None
        The captain gene feature to exclude from the cargo list.
    hmm_tags : dict or None
        Pre-computed protein_id -> tag mapping from HMM search.

    Returns
    -------
    list of CargoGene
        Cargo genes found within the Starship boundaries.
    """
    cargo_genes: List[CargoGene] = []

    if captain_feature is not None:
        captain_start = int(captain_feature.location.start)
        captain_end = int(captain_feature.location.end)
    else:
        captain_start = None
        captain_end = None

    for feature in record.features:
        if feature.type != "CDS":
            continue

        feat_start = int(feature.location.start)
        feat_end = int(feature.location.end)

        # Check that the feature falls within the Starship boundaries
        if feat_start < start or feat_end > end:
            continue

        # Exclude the captain gene by matching its location
        if captain_start is not None and feat_start == captain_start and feat_end == captain_end:
            continue

        gene_id = feature.qualifiers.get(
            "protein_id",
            [feature.qualifiers.get("ID", ["unknown"])[0]],
        )[0]
        product = feature.qualifiers.get("product", ["hypothetical protein"])[0]
        strand = feature.location.strand or 1

        # Get tag from HMM search results
        tag = ""
        if hmm_tags:
            tag = hmm_tags.get(gene_id, "")

        cargo_genes.append(
            CargoGene(
                gene_id=gene_id,
                product=product,
                start=feat_start,
                end=feat_end,
                strand=strand,
                tag=tag,
            )
        )

    logger.info(
        f"Extracted {len(cargo_genes)} cargo gene(s) from "
        f"{record.id}:{start}-{end}"
    )
    return cargo_genes
