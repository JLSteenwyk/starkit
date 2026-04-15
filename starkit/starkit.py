#!/usr/bin/env python

import logging
import os
import sys
import tempfile
import time

from Bio import SeqIO

from .args_processing import process_args
from .boundaries import define_boundaries, load_tir_pwms, load_family_reference
from .captain import detect_captains, load_hmm_profiles, full_genome_captain_search
from .cargo import extract_cargo, load_auxiliary_hmms, tag_proteins_by_hmm, sixframe_marker_search
from .classify import classify_starships
from .confidence import score_starships
from .files import load_genome
from .helpers import compute_genome_stats
from .homology import detect_by_homology, HomologyHit
from .logger import logger, log_file_logger
from .models import CaptainHit, EvidenceLevel, StarshipResult, StarKITRun
from .parser import create_parser
from .report import generate_report
from .settings import (
    CAPTAIN_HMM_DIR, BOUNDARY_DATA_DIR, FAMILY_HMM_DIR, STARSHIP_REF_FASTA,
)
from .version import __version__ as current_version
from .dedup import resolve_overlaps
from .write import write_user_args, write_output_stats, write_tsv, write_fasta, write_bed


def _write_temp_fasta(records):
    """Write SeqRecords to a temporary FASTA file for mappy. Returns path."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False,
    )
    SeqIO.write(records, tmp, "fasta")
    tmp.close()
    return tmp.name


def _homology_hit_to_starship(hit, record, idx, cargo_genes):
    """Convert a HomologyHit into a StarshipResult with a synthetic captain."""
    # Create a placeholder CaptainHit for homology-only predictions
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    placeholder_feature = SeqFeature(
        location=SimpleLocation(hit.start, min(hit.start + 1000, hit.end), strand=hit.strand),
        type="CDS",
        qualifiers={"note": ["homology-based prediction, no captain gene detected"]},
    )
    placeholder_captain = CaptainHit(
        feature=placeholder_feature,
        contig_id=hit.contig_id,
        start=hit.start,
        end=min(hit.start + 1000, hit.end),
        strand=hit.strand,
        evalue=1.0,  # no HMM hit
        score=0.0,
        hmm_name="homology",
        protein_id=f"homology_{hit.query_name}",
    )

    return StarshipResult(
        starship_id=f"starship_{idx:03d}",
        contig_id=hit.contig_id,
        region=record,
        start=hit.start,
        end=hit.end,
        captain=placeholder_captain,
        captain_family=hit.query_family if hit.query_family != "unclassified" else "unclassified",
        family_score=0.0,
        tir_left=None,
        tir_right=None,
        tsd=None,
        cargo_genes=cargo_genes,
        confidence_score=0.0,
        evidence_level=EvidenceLevel.LOW,
        truncated=False,
        boundary_method="homology",
        homology_identity=hit.identity,
        homology_coverage=hit.coverage,
        homology_reference_id=hit.query_name,
        homology_reference_family=hit.query_family,
    )


def _find_best_homology_for_captain(captain, homology_hits):
    """Find the best homology hit that overlaps a captain gene location."""
    best = None
    best_coverage = 0
    for hit in homology_hits:
        if hit.contig_id != captain.contig_id:
            continue
        # Captain must be inside the homology region
        if hit.start <= captain.start and hit.end >= captain.end:
            if hit.coverage > best_coverage:
                best = hit
                best_coverage = hit.coverage
    return best


def _find_nearby_homology_end(captain, homology_hits, fragment_hits=None, max_distance=700000):
    """Find the 3' end position of any homology hit/fragment nearby the
    captain (on the same contig, within max_distance). Returns the farthest
    valid 3' position as a boundary hint for DR scanning, or None.
    """
    candidates = []
    all_hits = list(homology_hits or [])
    if fragment_hits:
        all_hits.extend(fragment_hits)

    for hit in all_hits:
        if hit.contig_id != captain.contig_id:
            continue
        # Look for hits that extend downstream of the captain
        if captain.strand >= 0:
            if hit.end >= captain.start and hit.end - captain.start <= max_distance:
                candidates.append(hit.end)
        else:
            if hit.start <= captain.end and captain.end - hit.start <= max_distance:
                candidates.append(hit.start)

    if not candidates:
        return None
    # Return the farthest boundary — real large elements get a hint toward
    # their actual 3' position.
    if captain.strand >= 0:
        return max(candidates)
    return min(candidates)


def _homology_hit_used_by_captain(hit, captain_starships):
    """Check if a homology hit region is already covered by a captain prediction."""
    for sr in captain_starships:
        if hit.contig_id != sr.contig_id:
            continue
        overlap_start = max(hit.start, sr.start)
        overlap_end = min(hit.end, sr.end)
        overlap = max(0, overlap_end - overlap_start)
        hit_len = hit.end - hit.start
        sr_len = sr.end - sr.start
        if hit_len > 0 and sr_len > 0:
            # Consider used if substantial overlap from either perspective
            if overlap / hit_len >= 0.5 or overlap / sr_len >= 0.5:
                return True
    return False


def run(
    input_file,
    output_prefix,
    gff_file,
    evalue,
    min_size,
    max_size,
    evidence,
    use_log,
    quiet,
    no_homology=False,
    library=None,
    run_mode="relaxed",
):
    """Core pipeline: parse input -> detect captains -> homology search ->
    tiered boundary detection -> extract cargo -> classify -> score."""

    # Step 1: Parse & normalize input
    records = load_genome(input_file, gff_file)
    genome_stats = compute_genome_stats(records)

    # Step 2: Detect captain genes from annotated proteins
    captain_hits = detect_captains(records, CAPTAIN_HMM_DIR, evalue)

    # Step 2a: Full-genome six-frame captain search to catch unannotated captains.
    # Gene predictors sometimes miss captain genes — search the raw sequence
    # in all 6 reading frames. Uses stricter e-value to control false positives.
    try:
        captain_hmm_profiles = load_hmm_profiles(CAPTAIN_HMM_DIR)
        sixframe_novel = full_genome_captain_search(
            records,
            captain_hmm_profiles,
            evalue_threshold=1e-15,
            existing_captains=captain_hits,
        )
        if sixframe_novel:
            logger.info(
                f"Full-genome six-frame search found {len(sixframe_novel)} "
                f"additional captain(s) not present in annotation"
            )
            captain_hits.extend(sixframe_novel)
    except Exception as e:
        logger.warning(f"Full-genome captain search failed: {e}")

    # Step 2b: Tag all proteins with auxiliary gene HMMs (once for the whole genome)
    from .helpers import get_protein_sequences
    all_proteins = get_protein_sequences(records)
    aux_hmms = load_auxiliary_hmms()
    hmm_tags = tag_proteins_by_hmm(
        [(pid, seq) for pid, seq, _, _ in all_proteins],
        aux_hmms,
    ) if aux_hmms else {}

    # Step 2c: Classify captains early (needed for family-specific DR scanning)
    #          Store results to apply to StarshipResults later without re-running.
    captain_classifications = {}  # protein_id -> (family_name, family_score)
    family_hmm_dir = FAMILY_HMM_DIR
    family_hmms = None
    if os.path.isdir(family_hmm_dir) and os.listdir(family_hmm_dir):
        from .classify import classify_captain, load_family_hmms
        family_hmms = load_family_hmms(family_hmm_dir)
        for captain in captain_hits:
            family_name, family_score = classify_captain(
                captain, records, family_hmms,
            )
            captain.hmm_name = family_name
            captain_classifications[captain.protein_id] = (family_name, family_score)

    # Load reference data for boundary detection
    pwm_file = os.path.join(BOUNDARY_DATA_DIR, "tir_pwms.json")
    pwm_data = load_tir_pwms(pwm_file)
    family_ref_file = os.path.join(FAMILY_HMM_DIR, "family_reference.json")
    family_ref = load_family_reference(family_ref_file)

    records_by_id = {rec.id: rec for rec in records}

    # Step 3: Run homology search (needed BEFORE boundary detection for Tier 1)
    homology_hits = []
    fragment_hits = []
    genome_fasta = None
    combined_ref = None
    try:
        # Build the reference FASTA (built-in + optional custom library)
        ref_sources = []
        if os.path.exists(STARSHIP_REF_FASTA):
            ref_sources.append(STARSHIP_REF_FASTA)
        if library and os.path.exists(library):
            ref_sources.append(library)
            logger.info(f"Using custom library: {library}")

        if not no_homology and ref_sources:
            genome_fasta = _write_temp_fasta(records)

            # Concatenate reference sources if multiple
            if len(ref_sources) == 1:
                ref_fasta = ref_sources[0]
            else:
                combined_ref = tempfile.NamedTemporaryFile(
                    mode="w", suffix=".fasta", delete=False,
                )
                for src in ref_sources:
                    with open(src) as f:
                        combined_ref.write(f.read())
                combined_ref.close()
                ref_fasta = combined_ref.name

            from .homology import (
                search_homology,
                merge_overlapping_hits,
                merge_fragment_pairs,
            )
            full_hits, raw_fragment_hits = search_homology(
                genome_fasta, ref_fasta,
                min_identity=0.80, min_coverage=0.50,
            )
            fragment_hits = raw_fragment_hits  # keep for boundary hinting
            # Reconstruct large elements from fragment pairs (addresses
            # cases where the element middle has rearranged but the
            # edges match a known reference)
            reconstructed = merge_fragment_pairs(fragment_hits)
            if reconstructed:
                logger.info(
                    f"Reconstructed {len(reconstructed)} element(s) from "
                    f"conserved edges of rearranged Starships"
                )
            all_hits = full_hits + reconstructed
            homology_hits = merge_overlapping_hits(all_hits, merge_distance=1000)
    finally:
        if genome_fasta and os.path.exists(genome_fasta):
            os.unlink(genome_fasta)
        if combined_ref and os.path.exists(combined_ref.name):
            os.unlink(combined_ref.name)

    # Step 4: Tiered boundary detection + cargo extraction for each captain
    starship_results = []
    for idx, captain in enumerate(captain_hits, 1):
        record = records_by_id.get(captain.contig_id)
        if record is None:
            logger.warning(
                f"Contig '{captain.contig_id}' not found for captain "
                f"'{captain.protein_id}'; skipping"
            )
            continue

        # Find best overlapping homology hit for Tier 1
        best_homology = _find_best_homology_for_captain(captain, homology_hits)

        # Find a nearby homology/fragment 3' position as a boundary hint
        # for Tier 2 DR scanning (helps pick the right DR pair when many
        # copies of the same motif exist)
        homology_end_hint = _find_nearby_homology_end(
            captain, homology_hits, fragment_hits,
        )

        # Tiered boundary detection
        boundary = define_boundaries(
            captain, record, min_size, max_size,
            pwm_data=pwm_data,
            family_ref=family_ref,
            homology_hit=best_homology,
            homology_end_hint=homology_end_hint,
        )

        cargo_genes = extract_cargo(
            record, boundary["start"], boundary["end"], captain.feature,
            hmm_tags=hmm_tags,
        )

        result = StarshipResult(
            starship_id=f"starship_{idx:03d}",
            contig_id=captain.contig_id,
            region=record,
            start=boundary["start"],
            end=boundary["end"],
            captain=captain,
            tir_left=boundary["tir_left"],
            tir_right=boundary["tir_right"],
            tsd=boundary["tsd"],
            cargo_genes=cargo_genes,
            truncated=boundary["truncated"],
            boundary_method=boundary.get("boundary_method", "estimated"),
            homology_identity=best_homology.identity if best_homology else 0.0,
            homology_coverage=best_homology.coverage if best_homology else 0.0,
            homology_reference_id=best_homology.query_name if best_homology else "",
            homology_reference_family=best_homology.query_family if best_homology else "",
        )

        # Flag captain orientation issues.
        # Canonical: captain at one END of the element, pointing INWARD toward
        # the rest of the element. For + strand captains, they belong at the
        # 5' end (near start), pointing right into the cargo. For - strand
        # captains, they belong at the 3' end (near end), pointing left.
        #
        # We flag if EITHER:
        #   - the captain is not within the first/last 20% of the element
        #     (position check), OR
        #   - the captain points outward instead of inward (direction check).
        element_size = result.end - result.start
        if element_size > 0:
            # Position check: how far from each edge is the captain?
            dist_from_5p = captain.start - result.start   # 0 = flush at 5' edge
            dist_from_3p = result.end - captain.end       # 0 = flush at 3' edge
            edge_threshold = element_size * 0.20
            near_5p = dist_from_5p <= edge_threshold
            near_3p = dist_from_3p <= edge_threshold

            if captain.strand >= 0:
                # + strand captain: must be near 5' end and pointing right (+)
                if not near_5p:
                    result.captain_orientation_flag = True
            else:
                # - strand captain: must be near 3' end and pointing left (-)
                if not near_3p:
                    result.captain_orientation_flag = True

            # Additionally: check merged additional captains. If any of the
            # merged captains has an orientation inconsistent with the primary,
            # flag the element as having overlapping/conflicting captains.
            for extra in (result.additional_captains or []):
                ex_5p = extra.start - result.start <= edge_threshold
                ex_3p = result.end - extra.end <= edge_threshold
                if extra.strand >= 0 and not ex_5p:
                    result.captain_orientation_flag = True
                    break
                if extra.strand < 0 and not ex_3p:
                    result.captain_orientation_flag = True
                    break
                # Also flag if extra captain points in opposite direction
                # from primary — two captains facing away from each other
                # indicates a split/fragmented element, not a canonical one.
                if extra.strand != captain.strand:
                    result.captain_orientation_flag = True
                    break

        # Flag captain truncation (protein < 400 aa equivalent)
        captain_nt_len = captain.end - captain.start
        captain_aa_len = captain_nt_len / 3
        if captain_aa_len < 400:
            result.captain_truncated_flag = True

        starship_results.append(result)

    # Step 5: Add novel homology-only Starships (no captain detected)
    # For each, attempt six-frame translation to find unannotated captains
    from .captain import sixframe_captain_search
    try:
        sixframe_hmms = load_hmm_profiles(CAPTAIN_HMM_DIR)
    except Exception:
        sixframe_hmms = None

    next_idx = len(starship_results) + 1
    homology_added = 0
    sixframe_found = 0
    for hit in homology_hits:
        hit_size = hit.end - hit.start
        if min_size and hit_size < min_size:
            continue
        if max_size and hit_size > max_size:
            continue
        if _homology_hit_used_by_captain(hit, starship_results):
            continue

        record = records_by_id.get(hit.contig_id)
        if record is None:
            continue

        # Try six-frame translation to find unannotated captain
        sixframe_captain = None
        if sixframe_hmms:
            sixframe_hits = sixframe_captain_search(
                record, hit.start, hit.end, sixframe_hmms, evalue,
            )
            if sixframe_hits:
                sixframe_captain = sixframe_hits[0]
                sixframe_found += 1
                logger.info(
                    f"Six-frame captain found: {sixframe_captain.protein_id} "
                    f"(e={sixframe_captain.evalue:.2e}) in {hit.contig_id}:"
                    f"{hit.start:,}-{hit.end:,}"
                )

        if sixframe_captain:
            # Upgrade to captain-based prediction
            cargo_genes = extract_cargo(
                record, hit.start, hit.end, sixframe_captain.feature,
                hmm_tags=hmm_tags,
            )
            result = StarshipResult(
                starship_id=f"starship_{next_idx:03d}",
                contig_id=hit.contig_id,
                region=record,
                start=hit.start,
                end=hit.end,
                captain=sixframe_captain,
                captain_family=hit.query_family if hit.query_family != "unclassified" else "unclassified",
                family_score=0.0,
                cargo_genes=cargo_genes,
                boundary_method="homology",
                homology_identity=hit.identity,
                homology_coverage=hit.coverage,
                homology_reference_id=hit.query_name,
                homology_reference_family=hit.query_family,
            )
            # Classify the six-frame captain
            if family_hmms is not None:
                fname, fscore = classify_captain(sixframe_captain, records, family_hmms)
                sixframe_captain.hmm_name = fname
                captain_classifications[sixframe_captain.protein_id] = (fname, fscore)
        else:
            cargo_genes = extract_cargo(record, hit.start, hit.end, None, hmm_tags=hmm_tags)
            result = _homology_hit_to_starship(hit, record, next_idx, cargo_genes)

        starship_results.append(result)
        next_idx += 1
        homology_added += 1

    if homology_added:
        logger.info(
            f"Added {homology_added} Starship(s) from homology search"
            f" ({sixframe_found} with six-frame captains)"
        )

    if not starship_results:
        logger.info("No Starships predicted.")
        return StarKITRun(
            input_file=input_file,
            genome_stats=genome_stats,
            starships=[],
            parameters=dict(
                output_prefix=output_prefix, evalue=evalue,
                min_size=min_size, max_size=max_size, evidence=evidence,
            ),
            version=current_version,
        )

    # Step 6: Apply family classifications from Step 2b (no re-computation)
    for result in starship_results:
        pid = result.captain.protein_id
        if pid in captain_classifications:
            result.captain_family, result.family_score = captain_classifications[pid]

    # Step 7: Resolve overlapping predictions
    starship_results = resolve_overlaps(starship_results)

    # Step 7a: Annotate Starship-associated marker genes (MYB/SANT, DUF3723).
    # These are the two most Starship-specific auxiliary genes and serve as
    # confidence markers. A MYB/SANT gene at the opposite edge from the
    # captain with inverted orientation strongly supports element boundaries.
    #
    # Supplement annotated cargo with de novo six-frame marker search — some
    # gene predictors miss MYB/SANT or DUF3723 genes, so we scan the element
    # region directly for these markers.
    for r in starship_results:
        existing_myb_duf = [c for c in r.cargo_genes if c.tag in ("myb", "d37")]
        # Only run six-frame search if markers are missing from annotation,
        # to keep runtime low.
        if not existing_myb_duf and aux_hmms:
            denovo_markers = sixframe_marker_search(
                r.region, r.start, r.end, aux_hmms,
            )
            # Don't duplicate annotated genes at the same position
            existing_positions = {
                (c.start, c.end) for c in r.cargo_genes
            }
            for marker in denovo_markers:
                overlaps = any(
                    marker.start < end and marker.end > start
                    for start, end in existing_positions
                )
                if not overlaps:
                    r.cargo_genes.append(marker)
                    logger.info(
                        f"De novo {marker.tag} marker in {r.starship_id}: "
                        f"{r.contig_id}:{marker.start:,}-{marker.end:,}"
                    )

        myb_genes = [c for c in r.cargo_genes if c.tag == "myb"]
        duf_genes = [c for c in r.cargo_genes if c.tag == "d37"]
        r.has_myb_marker = bool(myb_genes)
        r.has_duf3723_marker = bool(duf_genes)

        # Check MYB/SANT position + orientation relative to the captain.
        # Canonical: captain at one edge, MYB/SANT at the OPPOSITE edge,
        # pointing in the OPPOSITE direction.
        if myb_genes:
            element_size = r.end - r.start
            edge_threshold = element_size * 0.25 if element_size > 0 else 0
            captain_near_5p = (r.captain.start - r.start) <= edge_threshold
            for myb in myb_genes:
                myb_near_5p = (myb.start - r.start) <= edge_threshold
                myb_near_3p = (r.end - myb.end) <= edge_threshold
                opposite_edge = (
                    (captain_near_5p and myb_near_3p)
                    or (not captain_near_5p and myb_near_5p)
                )
                opposite_strand = (
                    (r.captain.strand >= 0 and myb.strand < 0)
                    or (r.captain.strand < 0 and myb.strand >= 0)
                )
                if opposite_edge and opposite_strand:
                    r.myb_at_opposite_edge = True
                    break

    # Step 7b: Re-check orientation flags now that additional_captains are
    # populated by resolve_overlaps. If any merged captain points in a
    # different direction than the primary, flag the element.
    for r in starship_results:
        if r.captain_orientation_flag:
            continue  # already flagged
        primary_strand = r.captain.strand
        element_size = r.end - r.start
        edge_threshold = element_size * 0.20 if element_size > 0 else 0
        for extra in (r.additional_captains or []):
            ex_5p = extra.start - r.start <= edge_threshold
            ex_3p = r.end - extra.end <= edge_threshold
            # Extra captain not at a canonical edge, or facing opposite way
            if extra.strand != primary_strand:
                r.captain_orientation_flag = True
                break
            if extra.strand >= 0 and not ex_5p:
                r.captain_orientation_flag = True
                break
            if extra.strand < 0 and not ex_3p:
                r.captain_orientation_flag = True
                break

    # Step 8: Score confidence
    score_starships(starship_results)

    # Filter by evidence level
    if evidence != "all":
        level_map = {
            "high": [EvidenceLevel.HIGH],
            "medium": [EvidenceLevel.HIGH, EvidenceLevel.MEDIUM],
            "low": [EvidenceLevel.HIGH, EvidenceLevel.MEDIUM, EvidenceLevel.LOW],
        }
        allowed = level_map.get(evidence, [e for e in EvidenceLevel])
        starship_results = [r for r in starship_results if r.evidence_level in allowed]

    # Strict mode: remove predictions that violate canonical structure
    # Captain must be at the element edge with correct orientation
    if run_mode == "strict":
        before = len(starship_results)
        strict_results = []
        for r in starship_results:
            # Homology-only predictions without a real captain are excluded
            if r.captain.evalue >= 1.0:
                logger.info(
                    f"Strict: removing {r.starship_id} (no captain gene)"
                )
                continue
            # Captain must be at element edge with canonical orientation
            if r.captain_orientation_flag:
                logger.info(
                    f"Strict: removing {r.starship_id} (captain orientation)"
                )
                continue
            strict_results.append(r)
        starship_results = strict_results
        removed = before - len(starship_results)
        if removed:
            logger.info(f"Strict mode: removed {removed} prediction(s)")

    return StarKITRun(
        input_file=input_file,
        genome_stats=genome_stats,
        starships=starship_results,
        parameters=dict(
            output_prefix=output_prefix, evalue=evalue,
            min_size=min_size, max_size=max_size, evidence=evidence,
            mode=run_mode,
        ),
        version=current_version,
    )


def execute(
    input_file,
    output_prefix,
    gff_file,
    evalue,
    min_size,
    max_size,
    evidence,
    use_log,
    quiet,
    **kwargs,
):
    """High-level orchestration: logging setup, run pipeline, write output."""
    if use_log:
        log_file_logger.setLevel(logging.DEBUG)
        log_file_logger.propagate = False
        fh = logging.FileHandler(f"{output_prefix}.log", mode="w")
        fh.setLevel(logging.DEBUG)
        log_file_logger.addHandler(fh)

    if quiet:
        logger.disabled = True

    start_time = time.time()

    # Display run parameters
    write_user_args(
        input_file, output_prefix, evalue, min_size, max_size,
        evidence, gff_file, use_log, quiet,
    )

    # Run the pipeline
    no_homology = kwargs.get("no_homology", False)
    library = kwargs.get("library", None)
    run_mode = kwargs.get("run_mode", "relaxed")
    starkit_run = run(
        input_file, output_prefix, gff_file,
        evalue, min_size, max_size, evidence,
        use_log, quiet, no_homology=no_homology,
        library=library, run_mode=run_mode,
    )

    # Write output
    elapsed = time.time() - start_time
    write_tsv(starkit_run, output_prefix)
    write_fasta(starkit_run, output_prefix)
    write_bed(starkit_run, output_prefix)
    generate_report(starkit_run, output_prefix, runtime=elapsed)

    # Display stats
    write_output_stats(starkit_run, start_time)


def main(argv=None):
    """Parse arguments and execute the pipeline."""
    parser = create_parser()
    args = parser.parse_args(argv)
    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])
