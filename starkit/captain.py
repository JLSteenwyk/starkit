"""
Captain gene detection module.

Captain genes are tyrosine recombinases that define the 5' boundary of
Starship transposable elements. This module detects them by searching
protein sequences against curated HMM profiles using pyhmmer.
"""

import os
from typing import List, Tuple

import pyhmmer
import pyhmmer.easel
import pyhmmer.plan7

from starkit.exceptions import StarKITException
from starkit.helpers import get_protein_sequences
from starkit.models import CaptainHit
from starkit.settings import CAPTAIN_HMM_DIR


def load_hmm_profiles(hmm_dir: str) -> List[pyhmmer.plan7.HMM]:
    """Load all HMM profiles from the given directory.

    Parameters
    ----------
    hmm_dir : str
        Path to directory containing ``.hmm`` profile files.

    Returns
    -------
    list of pyhmmer.plan7.HMM
        Loaded HMM profile objects.

    Raises
    ------
    StarKITException
        If no ``.hmm`` files are found in *hmm_dir*.
    """
    hmms: List[pyhmmer.plan7.HMM] = []

    hmm_files = sorted(
        f for f in os.listdir(hmm_dir) if f.endswith(".hmm")
    )

    for filename in hmm_files:
        filepath = os.path.join(hmm_dir, filename)
        with pyhmmer.plan7.HMMFile(filepath) as hmm_file:
            for hmm in hmm_file:
                hmms.append(hmm)

    if not hmms:
        raise StarKITException("No captain HMM profiles found")

    return hmms


def search_captains(
    proteins: List[Tuple],
    hmm_profiles: List[pyhmmer.plan7.HMM],
    evalue_threshold: float,
) -> List[CaptainHit]:
    """Search protein sequences against captain HMM profiles.

    Parameters
    ----------
    proteins : list of tuple
        Each tuple is ``(protein_id, protein_sequence_str, contig_id, feature)``
        as returned by :func:`starkit.helpers.get_protein_sequences`.
    hmm_profiles : list of pyhmmer.plan7.HMM
        HMM profiles loaded by :func:`load_hmm_profiles`.
    evalue_threshold : float
        Maximum E-value for a hit to be retained.

    Returns
    -------
    list of CaptainHit
        Detected captain gene hits passing the threshold.
    """
    if not proteins:
        return []

    # Build a lookup from protein_id -> (contig_id, feature) for mapping
    # hits back to their genomic context.
    protein_lookup = {
        pid: (contig_id, feature)
        for pid, _, contig_id, feature in proteins
    }

    # Convert protein sequences to digitised pyhmmer sequences.
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = [
        pyhmmer.easel.TextSequence(
            name=pid.encode(),
            sequence=seq,
        ).digitize(alphabet)
        for pid, seq, _, _ in proteins
    ]

    captain_hits: List[CaptainHit] = []

    for top_hits in pyhmmer.hmmsearch(hmm_profiles, sequences):
        hmm_name = top_hits.query_name.decode()
        for hit in top_hits:
            if hit.included and hit.evalue <= evalue_threshold:
                protein_id = hit.name.decode()
                contig_id, feature = protein_lookup[protein_id]

                location = feature.location
                start = int(location.start)
                end = int(location.end)
                strand = int(location.strand)

                captain_hits.append(
                    CaptainHit(
                        feature=feature,
                        contig_id=contig_id,
                        start=start,
                        end=end,
                        strand=strand,
                        evalue=hit.evalue,
                        score=hit.score,
                        hmm_name=hmm_name,
                        protein_id=protein_id,
                    )
                )

    return captain_hits


def detect_captains(
    records,
    hmm_dir: str,
    evalue_threshold: float,
) -> List[CaptainHit]:
    """Detect captain genes in annotated genome records.

    This is the main entry point for captain detection. It extracts protein
    sequences from the genome records, loads captain HMM profiles, and
    searches the proteins against the profiles.

    Parameters
    ----------
    records : list of Bio.SeqRecord.SeqRecord
        Annotated genome records (contigs) with CDS features.
    hmm_dir : str
        Path to directory containing captain HMM profile files.
    evalue_threshold : float
        Maximum E-value for a hit to be retained.

    Returns
    -------
    list of CaptainHit
        Detected captain gene hits.
    """
    proteins = get_protein_sequences(records)
    hmm_profiles = load_hmm_profiles(hmm_dir)
    return search_captains(proteins, hmm_profiles, evalue_threshold)
