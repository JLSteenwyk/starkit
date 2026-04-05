"""
Cargo gene extraction module.

Cargo genes are the functional payload carried within Starship transposable
elements. This module identifies all CDS features that fall within a
Starship's boundaries, excluding the captain gene itself.
"""

import logging
from typing import List

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from starkit.models import CargoGene

logger = logging.getLogger(__name__)


def extract_cargo(
    record: SeqRecord,
    start: int,
    end: int,
    captain_feature: SeqFeature,
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
    captain_feature : Bio.SeqFeature.SeqFeature
        The captain gene feature to exclude from the cargo list.

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

        cargo_genes.append(
            CargoGene(
                gene_id=gene_id,
                product=product,
                start=feat_start,
                end=feat_end,
                strand=strand,
            )
        )

    logger.info(
        f"Extracted {len(cargo_genes)} cargo gene(s) from "
        f"{record.id}:{start}-{end}"
    )
    return cargo_genes
