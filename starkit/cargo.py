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


def sixframe_marker_search(
    record,
    region_start: int,
    region_end: int,
    hmm_profiles: List[pyhmmer.plan7.HMM],
    evalue_threshold: float = 1e-5,
) -> List[CargoGene]:
    """Search a region for MYB/SANT and DUF3723 marker genes via six-frame
    translation. Used to find markers missed by the gene predictor.

    Only the marker-relevant HMMs are used (MYB/SANT, DUF3723, Vip3A_N).
    Returns CargoGene objects for any hits found.
    """
    if not hmm_profiles:
        return []

    # Only use marker-relevant HMMs
    marker_tags = {"myb", "d37"}
    marker_hmms = []
    for hmm in hmm_profiles:
        hname = hmm.name
        hname_str = hname.decode() if isinstance(hname, bytes) else str(hname)
        tag = _HMM_TAG_MAP.get(hname_str, "")
        if tag in marker_tags:
            marker_hmms.append(hmm)

    if not marker_hmms:
        return []

    alphabet = pyhmmer.easel.Alphabet.amino()
    seq = record.seq[region_start:region_end]
    rc_seq = seq.reverse_complement()

    orf_sequences = []
    orf_metadata = {}  # orf_id -> (gs, ge, strand)
    min_aa_len = 50  # allow smaller ORFs for marker genes

    for frame in range(3):
        # Forward strand
        subseq = seq[frame:]
        subseq = subseq[:len(subseq) - len(subseq) % 3]
        if len(subseq) >= min_aa_len * 3:
            protein = str(subseq.translate(table=1))
            aa_pos = 0
            for orf_idx, orf in enumerate(protein.split("*")):
                if len(orf) >= min_aa_len:
                    gs = region_start + frame + aa_pos * 3
                    ge = gs + len(orf) * 3
                    orf_id = f"mk_{record.id}_{frame}f_{gs}"
                    orf_sequences.append(
                        pyhmmer.easel.TextSequence(
                            name=orf_id.encode(), sequence=orf,
                        ).digitize(alphabet)
                    )
                    orf_metadata[orf_id] = (gs, ge, 1)
                aa_pos += len(orf) + 1

        # Reverse strand
        subseq_rc = rc_seq[frame:]
        subseq_rc = subseq_rc[:len(subseq_rc) - len(subseq_rc) % 3]
        if len(subseq_rc) >= min_aa_len * 3:
            protein_rc = str(subseq_rc.translate(table=1))
            aa_pos = 0
            for orf_idx, orf in enumerate(protein_rc.split("*")):
                if len(orf) >= min_aa_len:
                    rc_start = frame + aa_pos * 3
                    rc_end = rc_start + len(orf) * 3
                    ge = region_end - rc_start
                    gs = region_end - rc_end
                    orf_id = f"mk_{record.id}_{frame}r_{gs}"
                    orf_sequences.append(
                        pyhmmer.easel.TextSequence(
                            name=orf_id.encode(), sequence=orf,
                        ).digitize(alphabet)
                    )
                    orf_metadata[orf_id] = (gs, ge, -1)
                aa_pos += len(orf) + 1

    if not orf_sequences:
        return []

    # Search with marker HMMs
    marker_genes = []
    seen_regions = []  # (start, end) tuples to dedupe overlapping ORFs

    hits = []
    for top_hits in pyhmmer.hmmsearch(marker_hmms, orf_sequences):
        query_name = top_hits.query.name
        hmm_name = query_name.decode() if isinstance(query_name, bytes) else str(query_name)
        tag = _HMM_TAG_MAP.get(hmm_name, "")
        if not tag:
            continue
        for hit in top_hits:
            if hit.included and hit.evalue <= evalue_threshold:
                hit_name = hit.name
                orf_id = hit_name.decode() if isinstance(hit_name, bytes) else str(hit_name)
                meta = orf_metadata.get(orf_id)
                if meta is None:
                    continue
                gs, ge, strand = meta
                hits.append({
                    "orf_id": orf_id,
                    "start": gs,
                    "end": ge,
                    "strand": strand,
                    "evalue": hit.evalue,
                    "tag": tag,
                })

    # Deduplicate overlapping hits (keep best e-value per region)
    hits.sort(key=lambda h: h["evalue"])
    for h in hits:
        overlapping = False
        for s, e in seen_regions:
            if h["start"] < e and h["end"] > s:
                overlapping = True
                break
        if overlapping:
            continue
        seen_regions.append((h["start"], h["end"]))
        marker_genes.append(
            CargoGene(
                gene_id=f"sixframe_{h['orf_id']}",
                product=f"six-frame predicted {h['tag'].upper()}",
                start=h["start"],
                end=h["end"],
                strand=h["strand"],
                tag=h["tag"],
            )
        )

    return marker_genes


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
