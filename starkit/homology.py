"""
Nucleotide homology search for Starship detection.

Uses mappy (Python minimap2 bindings) to align known Starship sequences
from Starbase against the input genome, detecting elements that may have
degraded or missing captain genes.
"""

import logging
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import mappy

logger = logging.getLogger(__name__)


@dataclass
class HomologyHit:
    """A region in the target genome matching a known Starship."""
    contig_id: str
    start: int
    end: int
    query_name: str        # reference Starship ID
    query_family: str      # family from reference header
    query_length: int      # length of the reference Starship
    aligned_length: int    # total aligned bases in target
    identity: float        # percent identity (0-1)
    coverage: float        # fraction of query covered (0-1)
    strand: int            # 1 for forward, -1 for reverse


def load_reference_starships(ref_fasta_path: str) -> str:
    """Verify the reference FASTA exists and return its path."""
    if not os.path.exists(ref_fasta_path):
        logger.warning(f"Reference Starship FASTA not found: {ref_fasta_path}")
        return None
    return ref_fasta_path


def parse_ref_header(header: str) -> Tuple[str, str, int]:
    """Parse reference FASTA header: >ship_id|family|lengthbp
    Returns (ship_id, family, length)."""
    parts = header.split("|")
    ship_id = parts[0] if len(parts) > 0 else header
    family = parts[1] if len(parts) > 1 else "unclassified"
    length = 0
    if len(parts) > 2:
        try:
            length = int(parts[2].replace("bp", ""))
        except ValueError:
            length = 0
    return ship_id, family, length


def search_homology(
    genome_fasta_path: str,
    ref_fasta_path: str,
    min_identity: float = 0.80,
    min_coverage: float = 0.50,
) -> List[HomologyHit]:
    """
    Search for Starship homologs in the genome using minimap2.

    Args:
        genome_fasta_path: Path to the input genome FASTA (we'll create a temp FASTA from SeqRecords)
        ref_fasta_path: Path to reference Starship FASTA
        min_identity: Minimum alignment identity (0-1)
        min_coverage: Minimum fraction of reference Starship covered

    Returns:
        List of HomologyHit objects passing thresholds
    """
    # Index the genome for nucleotide-to-nucleotide alignment
    # Use asm20 preset which is good for cross-species/divergent alignments
    aligner = mappy.Aligner(genome_fasta_path, preset="asm20", best_n=5)
    if not aligner:
        logger.warning("Failed to build minimap2 index for genome")
        return []

    hits = []

    # Read reference sequences and align each against the genome
    for name, seq, qual in mappy.fastx_read(ref_fasta_path):
        ship_id, family, query_length = parse_ref_header(name)
        if query_length == 0:
            query_length = len(seq)

        # Collect all alignments for this reference Starship
        alignments = list(aligner.map(seq))
        if not alignments:
            continue

        # For each alignment, compute identity and coverage
        for aln in alignments:
            identity = aln.mlen / aln.blen if aln.blen > 0 else 0
            # Coverage = how much of the query (reference Starship) was aligned
            coverage = (aln.q_en - aln.q_st) / query_length if query_length > 0 else 0

            if identity >= min_identity and coverage >= min_coverage:
                strand = 1 if aln.strand == 1 else -1
                hits.append(HomologyHit(
                    contig_id=aln.ctg,
                    start=aln.r_st,
                    end=aln.r_en,
                    query_name=ship_id,
                    query_family=family,
                    query_length=query_length,
                    aligned_length=aln.r_en - aln.r_st,
                    identity=round(identity, 4),
                    coverage=round(coverage, 4),
                    strand=strand,
                ))

    logger.info(f"Found {len(hits)} homology hit(s) passing thresholds")
    return hits


def merge_overlapping_hits(hits: List[HomologyHit], merge_distance: int = 5000) -> List[HomologyHit]:
    """
    Merge overlapping or nearby homology hits on the same contig into
    non-redundant regions. Keep the best-scoring hit's metadata.

    Hits within merge_distance bp of each other on the same contig are merged.
    """
    if not hits:
        return []

    # Group by contig
    by_contig = {}
    for hit in hits:
        by_contig.setdefault(hit.contig_id, []).append(hit)

    merged = []
    for contig_id, contig_hits in by_contig.items():
        # Sort by start position
        contig_hits.sort(key=lambda h: h.start)

        current = contig_hits[0]
        for hit in contig_hits[1:]:
            # If overlapping or within merge_distance, merge
            if hit.start <= current.end + merge_distance:
                # Extend the region
                new_end = max(current.end, hit.end)
                # Keep metadata from the hit with better coverage
                if hit.coverage > current.coverage:
                    best = hit
                else:
                    best = current
                current = HomologyHit(
                    contig_id=contig_id,
                    start=min(current.start, hit.start),
                    end=new_end,
                    query_name=best.query_name,
                    query_family=best.query_family,
                    query_length=best.query_length,
                    aligned_length=new_end - min(current.start, hit.start),
                    identity=best.identity,
                    coverage=best.coverage,
                    strand=best.strand,
                )
            else:
                merged.append(current)
                current = hit
        merged.append(current)

    logger.info(f"Merged into {len(merged)} non-redundant region(s)")
    return merged


def filter_novel_hits(
    homology_hits: List[HomologyHit],
    captain_starships: list,
    overlap_threshold: float = 0.5,
) -> List[HomologyHit]:
    """
    Filter out homology hits that overlap with already-detected captain-based
    Starships. Only return truly novel hits.

    Args:
        homology_hits: Merged homology hits
        captain_starships: List of StarshipResult from captain-based detection
        overlap_threshold: Minimum reciprocal overlap fraction to consider redundant

    Returns:
        Homology hits that don't overlap with captain-based detections
    """
    novel = []

    for hit in homology_hits:
        is_redundant = False
        for starship in captain_starships:
            if hit.contig_id != starship.contig_id:
                continue

            # Calculate overlap
            overlap_start = max(hit.start, starship.start)
            overlap_end = min(hit.end, starship.end)
            overlap_len = max(0, overlap_end - overlap_start)

            hit_len = hit.end - hit.start
            starship_len = starship.end - starship.start

            # Reciprocal overlap check
            if hit_len > 0 and starship_len > 0:
                overlap_frac_hit = overlap_len / hit_len
                overlap_frac_starship = overlap_len / starship_len
                if overlap_frac_hit >= overlap_threshold or overlap_frac_starship >= overlap_threshold:
                    is_redundant = True
                    break

        if not is_redundant:
            novel.append(hit)

    logger.info(f"{len(novel)} novel homology-only region(s) after filtering")
    return novel


def detect_by_homology(
    genome_fasta_path: str,
    ref_fasta_path: str,
    captain_starships: list = None,
    min_identity: float = 0.80,
    min_coverage: float = 0.50,
) -> List[HomologyHit]:
    """
    Main entry point for homology-based Starship detection.

    Returns novel HomologyHit regions not already found by captain detection.
    """
    if captain_starships is None:
        captain_starships = []

    ref_path = load_reference_starships(ref_fasta_path)
    if ref_path is None:
        return []

    hits = search_homology(genome_fasta_path, ref_path, min_identity, min_coverage)
    merged = merge_overlapping_hits(hits)
    novel = filter_novel_hits(merged, captain_starships)

    return novel
