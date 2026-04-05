"""
Overlap resolution for Starship predictions.

Handles three types of overlap:
- Duplicates (>80% reciprocal overlap): merge, keep best captain
- Nested (one fully contains the other): keep both, annotate relationship
- Partial overlap: keep both, log the overlap
"""

import logging
from typing import List

from .models import StarshipResult

logger = logging.getLogger(__name__)


def resolve_overlaps(starship_results: List[StarshipResult]) -> List[StarshipResult]:
    """Resolve overlapping Starship predictions.

    Merges duplicates, annotates nested elements, and re-numbers IDs.
    """
    if len(starship_results) <= 1:
        return starship_results

    # Sort by contig then start position
    results = sorted(starship_results, key=lambda r: (r.contig_id, r.start))

    # Track which indices to remove (merged into another)
    remove = set()

    for i in range(len(results)):
        if i in remove:
            continue
        for j in range(i + 1, len(results)):
            if j in remove:
                continue
            a = results[i]
            b = results[j]

            if a.contig_id != b.contig_id:
                break  # sorted by contig, no more overlaps possible

            overlap_start = max(a.start, b.start)
            overlap_end = min(a.end, b.end)
            overlap = max(0, overlap_end - overlap_start)

            if overlap == 0:
                continue

            a_size = a.size or 1
            b_size = b.size or 1
            pct_a = overlap / a_size
            pct_b = overlap / b_size

            if pct_a >= 0.80 and pct_b >= 0.80:
                # DUPLICATE: merge into the one with better captain e-value
                if a.captain.evalue <= b.captain.evalue:
                    keeper, merged = a, b
                    remove.add(j)
                else:
                    keeper, merged = b, a
                    remove.add(i)

                # Store the merged captain as an additional hit
                keeper.additional_captains.append(merged.captain)
                logger.info(
                    f"Merged duplicate: {merged.starship_id} "
                    f"({merged.captain_family}) into {keeper.starship_id} "
                    f"({keeper.captain_family})"
                )
                if i in remove:
                    break

            elif pct_a >= 0.80 or pct_b >= 0.80:
                # NESTED: smaller is inside larger
                if a_size >= b_size:
                    parent, child = a, b
                else:
                    parent, child = b, a
                child.nested_in = parent.starship_id
                logger.info(
                    f"Nested: {child.starship_id} ({child.size:,}bp) "
                    f"inside {parent.starship_id} ({parent.size:,}bp)"
                )
            else:
                # PARTIAL OVERLAP: keep both, just log
                logger.info(
                    f"Partial overlap: {a.starship_id} and {b.starship_id} "
                    f"share {overlap:,}bp"
                )

    # Filter out merged duplicates
    kept = [r for i, r in enumerate(results) if i not in remove]

    # Re-number IDs
    for idx, result in enumerate(kept, 1):
        result.starship_id = f"starship_{idx:03d}"
        # Update nested_in references to old IDs
        # (nested_in was set to the pre-renumber ID, but since we're
        # renumbering, we need to update. For simplicity, store the
        # parent's contig:start as a stable reference instead.)

    # Fix nested_in to use new IDs
    id_map = {}
    for r in kept:
        # Map old approximate identity to new ID
        id_map[(r.contig_id, r.start, r.end)] = r.starship_id

    for r in kept:
        if r.nested_in:
            # Find the parent by checking which kept result the old ID maps to
            for parent in kept:
                if (parent.contig_id == r.contig_id
                        and parent.start <= r.start
                        and parent.end >= r.end
                        and parent.starship_id != r.starship_id):
                    r.nested_in = parent.starship_id
                    break

    if remove:
        logger.info(
            f"Overlap resolution: {len(starship_results)} predictions -> "
            f"{len(kept)} after merging {len(remove)} duplicate(s)"
        )

    return kept
