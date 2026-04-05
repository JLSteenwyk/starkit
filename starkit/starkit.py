#!/usr/bin/env python

import logging
import os
import sys
import tempfile
import time

from Bio import SeqIO

from .args_processing import process_args
from .boundaries import define_boundaries, load_tir_pwms, load_family_reference
from .captain import detect_captains
from .cargo import extract_cargo
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
from .write import write_user_args, write_output_stats, write_tsv


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
):
    """Core pipeline: parse input -> detect captains -> homology search ->
    tiered boundary detection -> extract cargo -> classify -> score."""

    # Step 1: Parse & normalize input
    records = load_genome(input_file, gff_file)
    genome_stats = compute_genome_stats(records)

    # Step 2: Detect captain genes
    captain_hits = detect_captains(records, CAPTAIN_HMM_DIR, evalue)

    # Step 2b: Classify captains early (needed for family-specific DR scanning)
    family_hmm_dir = FAMILY_HMM_DIR
    if captain_hits and os.path.isdir(family_hmm_dir) and os.listdir(family_hmm_dir):
        from .classify import classify_captain, load_family_hmms
        family_hmms = load_family_hmms(family_hmm_dir)
        for captain in captain_hits:
            family_name, family_score = classify_captain(
                captain, records, family_hmms,
            )
            captain.hmm_name = family_name

    # Load reference data for boundary detection
    pwm_file = os.path.join(BOUNDARY_DATA_DIR, "tir_pwms.json")
    pwm_data = load_tir_pwms(pwm_file)
    family_ref_file = os.path.join(FAMILY_HMM_DIR, "family_reference.json")
    family_ref = load_family_reference(family_ref_file)

    records_by_id = {rec.id: rec for rec in records}

    # Step 3: Run homology search (needed BEFORE boundary detection for Tier 1)
    homology_hits = []
    genome_fasta = None
    try:
        if not no_homology and os.path.exists(STARSHIP_REF_FASTA):
            genome_fasta = _write_temp_fasta(records)
            from .homology import search_homology, merge_overlapping_hits
            raw_hits = search_homology(
                genome_fasta, STARSHIP_REF_FASTA,
                min_identity=0.80, min_coverage=0.50,
            )
            # DON'T merge aggressively — keep individual high-coverage hits
            # Only merge hits from the SAME reference Starship
            homology_hits = merge_overlapping_hits(raw_hits, merge_distance=1000)
    finally:
        if genome_fasta and os.path.exists(genome_fasta):
            os.unlink(genome_fasta)

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

        # Tiered boundary detection
        boundary = define_boundaries(
            captain, record, min_size, max_size,
            pwm_data=pwm_data,
            family_ref=family_ref,
            homology_hit=best_homology,
        )

        cargo_genes = extract_cargo(
            record, boundary["start"], boundary["end"], captain.feature,
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
        )
        starship_results.append(result)

    # Step 5: Add novel homology-only Starships (no captain detected)
    next_idx = len(starship_results) + 1
    homology_added = 0
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

        cargo_genes = extract_cargo(record, hit.start, hit.end, None)
        result = _homology_hit_to_starship(hit, record, next_idx, cargo_genes)
        starship_results.append(result)
        next_idx += 1
        homology_added += 1

    if homology_added:
        logger.info(f"Added {homology_added} Starship(s) from homology search")

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

    # Step 6: Classify families
    family_hmm_dir = FAMILY_HMM_DIR
    if os.path.isdir(family_hmm_dir) and os.listdir(family_hmm_dir):
        captain_based = [r for r in starship_results if r.captain.hmm_name != "homology"]
        if captain_based:
            classify_starships(captain_based, records, family_hmm_dir)

    # Step 7: Score confidence
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

    return StarKITRun(
        input_file=input_file,
        genome_stats=genome_stats,
        starships=starship_results,
        parameters=dict(
            output_prefix=output_prefix, evalue=evalue,
            min_size=min_size, max_size=max_size, evidence=evidence,
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
    starkit_run = run(
        input_file, output_prefix, gff_file,
        evalue, min_size, max_size, evidence,
        use_log, quiet, no_homology=no_homology,
    )

    # Write output
    write_tsv(starkit_run, output_prefix)
    generate_report(starkit_run, output_prefix)

    # Display stats
    write_output_stats(starkit_run, start_time)


def main(argv=None):
    """Parse arguments and execute the pipeline."""
    parser = create_parser()
    args = parser.parse_args(argv)
    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])
