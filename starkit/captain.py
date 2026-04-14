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
        query_name = top_hits.query.name
        hmm_name = query_name.decode() if isinstance(query_name, bytes) else str(query_name)
        for hit in top_hits:
            if hit.included and hit.evalue <= evalue_threshold:
                hit_name = hit.name
                protein_id = hit_name.decode() if isinstance(hit_name, bytes) else str(hit_name)
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


def sixframe_captain_search(
    record,
    region_start: int,
    region_end: int,
    hmm_profiles: List[pyhmmer.plan7.HMM],
    evalue_threshold: float,
) -> List[CaptainHit]:
    """Search for captain genes via six-frame translation of a genomic region.

    Used as a fallback for homology-only predictions where no annotated
    captain was found. Translates the region in all 6 reading frames and
    searches the resulting ORFs against captain HMM profiles.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        The contig containing the region.
    region_start, region_end : int
        Coordinates of the region to search.
    hmm_profiles : list of pyhmmer.plan7.HMM
        Captain HMM profiles.
    evalue_threshold : float
        Maximum E-value threshold.

    Returns
    -------
    list of CaptainHit
        Captain hits found in the six-frame translations.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    seq = record.seq[region_start:region_end]
    rc_seq = seq.reverse_complement()

    min_orf_len = 300  # minimum ORF length in nucleotides (100 aa)

    alphabet = pyhmmer.easel.Alphabet.amino()
    orf_sequences = []
    orf_metadata = []  # (orf_id, contig_start, contig_end, strand)

    for frame in range(3):
        # Forward strand
        subseq = seq[frame:]
        # Trim to multiple of 3
        subseq = subseq[:len(subseq) - len(subseq) % 3]
        if len(subseq) >= min_orf_len:
            protein = str(subseq.translate(table=1))
            # Split on stop codons, keep ORFs >= 100 aa
            pos = 0
            for orf_idx, orf in enumerate(protein.split("*")):
                orf_nt_start = frame + pos * 3
                if len(orf) >= 100:
                    orf_id = f"6f_{record.id}_{region_start}_{frame}f_{orf_idx}"
                    genomic_start = region_start + orf_nt_start
                    genomic_end = genomic_start + len(orf) * 3
                    orf_sequences.append(
                        pyhmmer.easel.TextSequence(
                            name=orf_id.encode(),
                            sequence=orf,
                        ).digitize(alphabet)
                    )
                    orf_metadata.append((orf_id, genomic_start, genomic_end, 1))
                pos += len(orf) + 1  # +1 for the stop codon

        # Reverse strand
        subseq_rc = rc_seq[frame:]
        subseq_rc = subseq_rc[:len(subseq_rc) - len(subseq_rc) % 3]
        if len(subseq_rc) >= min_orf_len:
            protein_rc = str(subseq_rc.translate(table=1))
            pos = 0
            for orf_idx, orf in enumerate(protein_rc.split("*")):
                orf_nt_start = frame + pos * 3
                if len(orf) >= 100:
                    orf_id = f"6f_{record.id}_{region_start}_{frame}r_{orf_idx}"
                    # Map reverse strand coordinates back to forward strand
                    rc_start = orf_nt_start
                    rc_end = rc_start + len(orf) * 3
                    genomic_end = region_end - rc_start
                    genomic_start = region_end - rc_end
                    orf_sequences.append(
                        pyhmmer.easel.TextSequence(
                            name=orf_id.encode(),
                            sequence=orf,
                        ).digitize(alphabet)
                    )
                    orf_metadata.append((orf_id, genomic_start, genomic_end, -1))
                pos += len(orf) + 1

    if not orf_sequences:
        return []

    # Build lookup
    meta_lookup = {oid: (gs, ge, strand) for oid, gs, ge, strand in orf_metadata}

    hits = []
    for top_hits in pyhmmer.hmmsearch(hmm_profiles, orf_sequences):
        query_name = top_hits.query.name
        hmm_name = query_name.decode() if isinstance(query_name, bytes) else str(query_name)
        for hit in top_hits:
            if hit.included and hit.evalue <= evalue_threshold:
                hit_name = hit.name
                orf_id = hit_name.decode() if isinstance(hit_name, bytes) else str(hit_name)
                if orf_id not in meta_lookup:
                    continue
                genomic_start, genomic_end, strand = meta_lookup[orf_id]

                # Create a synthetic feature for this ORF
                feature = SeqFeature(
                    location=SimpleLocation(genomic_start, genomic_end, strand=strand),
                    type="CDS",
                    qualifiers={"note": ["six-frame translation captain"]},
                )

                hits.append(
                    CaptainHit(
                        feature=feature,
                        contig_id=record.id,
                        start=genomic_start,
                        end=genomic_end,
                        strand=strand,
                        evalue=hit.evalue,
                        score=hit.score,
                        hmm_name=hmm_name,
                        protein_id=f"sixframe_{orf_id}",
                    )
                )

    # Keep only the best hit per region (avoid reporting multiple overlapping ORFs)
    if hits:
        hits.sort(key=lambda h: h.evalue)
        return [hits[0]]

    return []


def _sixframe_orfs_from_seq(
    record_id: str,
    seq_start: int,
    full_seq,
    alphabet,
    min_aa_len: int = 100,
):
    """Generate ORFs from six-frame translation of a sequence.

    Yields (orf_id, digitised_sequence, genomic_start, genomic_end, strand).
    """
    length = len(full_seq)
    rc_full = full_seq.reverse_complement()

    for frame in range(3):
        # Forward strand
        subseq = full_seq[frame:]
        subseq = subseq[:len(subseq) - len(subseq) % 3]
        if len(subseq) >= min_aa_len * 3:
            protein = str(subseq.translate(table=1))
            aa_pos = 0
            for orf_idx, orf in enumerate(protein.split("*")):
                if len(orf) >= min_aa_len:
                    genomic_start = seq_start + frame + aa_pos * 3
                    genomic_end = genomic_start + len(orf) * 3
                    orf_id = f"6f_{record_id}_{frame}f_{genomic_start}"
                    digitised = pyhmmer.easel.TextSequence(
                        name=orf_id.encode(),
                        sequence=orf,
                    ).digitize(alphabet)
                    yield orf_id, digitised, genomic_start, genomic_end, 1
                aa_pos += len(orf) + 1  # +1 for the stop codon

        # Reverse strand
        subseq_rc = rc_full[frame:]
        subseq_rc = subseq_rc[:len(subseq_rc) - len(subseq_rc) % 3]
        if len(subseq_rc) >= min_aa_len * 3:
            protein_rc = str(subseq_rc.translate(table=1))
            aa_pos = 0
            for orf_idx, orf in enumerate(protein_rc.split("*")):
                if len(orf) >= min_aa_len:
                    rc_start = frame + aa_pos * 3
                    rc_end = rc_start + len(orf) * 3
                    genomic_end = seq_start + length - rc_start
                    genomic_start = seq_start + length - rc_end
                    orf_id = f"6f_{record_id}_{frame}r_{genomic_start}"
                    digitised = pyhmmer.easel.TextSequence(
                        name=orf_id.encode(),
                        sequence=orf,
                    ).digitize(alphabet)
                    yield orf_id, digitised, genomic_start, genomic_end, -1
                aa_pos += len(orf) + 1


def full_genome_captain_search(
    records,
    hmm_profiles: List[pyhmmer.plan7.HMM],
    evalue_threshold: float = 1e-15,
    existing_captains: List[CaptainHit] = None,
) -> List[CaptainHit]:
    """Search for captain genes across the full genome via six-frame translation.

    This supplements annotation-based captain detection by finding captains
    that were missed by the gene predictor. Uses a stricter e-value threshold
    than annotated captain search to control false positives.

    Parameters
    ----------
    records : list of Bio.SeqRecord.SeqRecord
        Genome contigs.
    hmm_profiles : list of pyhmmer.plan7.HMM
        Captain HMM profiles.
    evalue_threshold : float
        E-value threshold (default 1e-15, stricter than annotation search).
    existing_captains : list of CaptainHit
        Already-detected captains (from annotation). Six-frame hits
        overlapping these are discarded (annotation takes precedence).

    Returns
    -------
    list of CaptainHit
        Novel captain hits from six-frame translation, not overlapping
        any existing annotated captain.
    """
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    if existing_captains is None:
        existing_captains = []

    # Index existing captains by contig for fast overlap checks
    existing_by_contig = {}
    for cap in existing_captains:
        existing_by_contig.setdefault(cap.contig_id, []).append((cap.start, cap.end))

    alphabet = pyhmmer.easel.Alphabet.amino()
    all_orf_sequences = []
    all_orf_metadata = {}  # orf_id -> (contig_id, start, end, strand)

    for record in records:
        seq = record.seq
        for orf_id, digitised, gs, ge, strand in _sixframe_orfs_from_seq(
            record.id, 0, seq, alphabet, min_aa_len=100,
        ):
            all_orf_sequences.append(digitised)
            all_orf_metadata[orf_id] = (record.id, gs, ge, strand)

    if not all_orf_sequences:
        return []

    # Run the HMM search across all ORFs at once
    all_hits = []
    for top_hits in pyhmmer.hmmsearch(hmm_profiles, all_orf_sequences):
        query_name = top_hits.query.name
        hmm_name = query_name.decode() if isinstance(query_name, bytes) else str(query_name)
        for hit in top_hits:
            if hit.included and hit.evalue <= evalue_threshold:
                hit_name = hit.name
                orf_id = hit_name.decode() if isinstance(hit_name, bytes) else str(hit_name)
                meta = all_orf_metadata.get(orf_id)
                if meta is None:
                    continue
                contig_id, gs, ge, strand = meta

                # Skip if this ORF overlaps an existing annotated captain
                overlaps_existing = False
                for ex_s, ex_e in existing_by_contig.get(contig_id, []):
                    if gs < ex_e and ge > ex_s:
                        overlaps_existing = True
                        break
                if overlaps_existing:
                    continue

                all_hits.append({
                    "orf_id": orf_id,
                    "contig_id": contig_id,
                    "start": gs,
                    "end": ge,
                    "strand": strand,
                    "evalue": hit.evalue,
                    "score": hit.score,
                    "hmm_name": hmm_name,
                })

    # Deduplicate overlapping six-frame hits (keep best e-value per region)
    all_hits.sort(key=lambda h: h["evalue"])
    kept = []
    used_regions = {}  # contig_id -> list of (start, end)
    for h in all_hits:
        regions = used_regions.setdefault(h["contig_id"], [])
        overlap = False
        for s, e in regions:
            if h["start"] < e and h["end"] > s:
                overlap = True
                break
        if not overlap:
            regions.append((h["start"], h["end"]))
            kept.append(h)

    # Convert to CaptainHit objects
    captain_hits = []
    for h in kept:
        feature = SeqFeature(
            location=SimpleLocation(h["start"], h["end"], strand=h["strand"]),
            type="CDS",
            qualifiers={"note": ["six-frame translation captain"]},
        )
        captain_hits.append(
            CaptainHit(
                feature=feature,
                contig_id=h["contig_id"],
                start=h["start"],
                end=h["end"],
                strand=h["strand"],
                evalue=h["evalue"],
                score=h["score"],
                hmm_name=h["hmm_name"],
                protein_id=f"sixframe_{h['orf_id']}",
            )
        )

    return captain_hits
