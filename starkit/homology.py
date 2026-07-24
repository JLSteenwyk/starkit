"""
Nucleotide homology search for Starship detection.

Uses mappy (Python minimap2 bindings) to align known Starship reference
sequences against the input genome, detecting elements that may have degraded
or missing captain genes.
"""

import logging
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import mappy

from .references import parse_ref_header

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
    query_start: int = 0   # start position in reference (query) coords
    query_end: int = 0     # end position in reference (query) coords
    fragments_merged: int = 1  # how many raw hits were merged into this


def load_reference_starships(ref_fasta_path: str) -> str:
    """Verify the reference FASTA exists and return its path."""
    if not os.path.exists(ref_fasta_path):
        logger.warning(f"Reference Starship FASTA not found: {ref_fasta_path}")
        return None
    return ref_fasta_path


def search_homology(
    genome_fasta_path: str,
    ref_fasta_path: str,
    min_identity: float = 0.80,
    min_coverage: float = 0.50,
    min_fragment_coverage: float = 0.15,
) -> Tuple[List[HomologyHit], List[HomologyHit]]:
    """
    Search for Starship homologs in the genome using minimap2.

    Returns two lists:
      - hits: full hits passing both identity AND coverage thresholds
      - fragments: partial hits passing identity but below full coverage
        (used by the fragment merging pass to reconstruct large
        rearranged elements from their conserved edges)
    """
    aligner = mappy.Aligner(genome_fasta_path, preset="asm20", best_n=5)
    if not aligner:
        logger.warning("Failed to build minimap2 index for genome")
        return [], []

    hits = []
    fragments = []

    for name, seq, qual in mappy.fastx_read(ref_fasta_path):
        ship_id, family, query_length = parse_ref_header(name)
        if query_length == 0:
            query_length = len(seq)

        alignments = list(aligner.map(seq))
        if not alignments:
            continue

        for aln in alignments:
            identity = aln.mlen / aln.blen if aln.blen > 0 else 0
            coverage = (aln.q_en - aln.q_st) / query_length if query_length > 0 else 0

            if identity < min_identity:
                continue

            strand = 1 if aln.strand == 1 else -1
            hit = HomologyHit(
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
                query_start=aln.q_st,
                query_end=aln.q_en,
            )

            if coverage >= min_coverage:
                hits.append(hit)
            elif coverage >= min_fragment_coverage:
                fragments.append(hit)

    logger.info(
        f"Found {len(hits)} full homology hit(s) and {len(fragments)} "
        f"fragment hit(s)"
    )
    return hits, fragments


def merge_fragment_pairs(
    fragments: List[HomologyHit],
    max_element_size: int = 700000,
    min_combined_coverage: float = 0.60,
) -> List[HomologyHit]:
    """
    Merge fragment hits that cover different parts of the SAME reference
    Starship on the SAME contig into a single reconstructed element.

    This addresses cases where a large Starship has rearranged internals
    but conserved edges: minimap2 reports two separate partial alignments
    (one covering reference bases 1-60kb, another covering 300-400kb),
    and we want to merge them into one element spanning both.

    The two fragments must:
      - map to the same reference Starship
      - cover non-overlapping parts of the reference (different q-coords)
      - combined cover at least min_combined_coverage of the reference
      - be on the same contig within max_element_size bp of each other
      - be in the same strand orientation
    """
    if not fragments:
        return []

    # Group fragments by (contig, reference)
    groups: dict = {}
    for f in fragments:
        key = (f.contig_id, f.query_name)
        groups.setdefault(key, []).append(f)

    merged = []
    for key, frags in groups.items():
        if len(frags) < 2:
            continue

        contig_id, ref_name = key
        frags.sort(key=lambda f: f.start)

        # For each pair, check if they can be merged
        used = set()
        for i in range(len(frags)):
            if i in used:
                continue
            a = frags[i]

            # Greedy: absorb every compatible partner into a
            cluster = [a]
            cluster_indices = {i}

            for j in range(i + 1, len(frags)):
                if j in used or j in cluster_indices:
                    continue
                b = frags[j]

                # Must be on same strand
                if a.strand != b.strand:
                    continue

                # Genomic distance check — don't merge hits too far apart
                genomic_span = max(x.end for x in cluster) - min(x.start for x in cluster + [b])
                if genomic_span > max_element_size:
                    continue

                # Check that b's query coords don't heavily overlap any member
                # of the cluster (we want non-redundant fragments of the same
                # reference)
                b_overlaps_cluster = False
                for c in cluster:
                    q_overlap_start = max(c.query_start, b.query_start)
                    q_overlap_end = min(c.query_end, b.query_end)
                    q_overlap = max(0, q_overlap_end - q_overlap_start)
                    c_q_len = c.query_end - c.query_start
                    b_q_len = b.query_end - b.query_start
                    if c_q_len > 0 and q_overlap / c_q_len > 0.5:
                        b_overlaps_cluster = True
                        break
                    if b_q_len > 0 and q_overlap / b_q_len > 0.5:
                        b_overlaps_cluster = True
                        break

                if b_overlaps_cluster:
                    continue

                cluster.append(b)
                cluster_indices.add(j)

            if len(cluster) < 2:
                continue

            # Compute combined coverage
            ref_len = a.query_length
            if ref_len <= 0:
                continue

            # Union of query ranges
            q_ranges = sorted((c.query_start, c.query_end) for c in cluster)
            union = []
            for qs, qe in q_ranges:
                if union and qs <= union[-1][1]:
                    union[-1] = (union[-1][0], max(union[-1][1], qe))
                else:
                    union.append((qs, qe))
            combined_covered = sum(qe - qs for qs, qe in union)
            combined_coverage = combined_covered / ref_len

            if combined_coverage < min_combined_coverage:
                continue

            # Build merged hit
            merged_start = min(c.start for c in cluster)
            merged_end = max(c.end for c in cluster)
            best_identity = max(c.identity for c in cluster)

            merged.append(HomologyHit(
                contig_id=contig_id,
                start=merged_start,
                end=merged_end,
                query_name=ref_name,
                query_family=cluster[0].query_family,
                query_length=ref_len,
                aligned_length=merged_end - merged_start,
                identity=round(best_identity, 4),
                coverage=round(combined_coverage, 4),
                strand=cluster[0].strand,
                query_start=min(c.query_start for c in cluster),
                query_end=max(c.query_end for c in cluster),
                fragments_merged=len(cluster),
            ))

            for idx in cluster_indices:
                used.add(idx)

    if merged:
        logger.info(
            f"Reconstructed {len(merged)} element(s) from fragment pairs"
        )
    return merged


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
                    query_start=best.query_start,
                    query_end=best.query_end,
                    fragments_merged=max(
                        current.fragments_merged, hit.fragments_merged
                    ),
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

    hits, fragments = search_homology(
        genome_fasta_path, ref_path, min_identity, min_coverage,
    )
    # Reconstruct large elements from fragment pairs
    reconstructed = merge_fragment_pairs(fragments)
    all_hits = hits + reconstructed
    merged = merge_overlapping_hits(all_hits)
    novel = filter_novel_hits(merged, captain_starships)

    return novel
