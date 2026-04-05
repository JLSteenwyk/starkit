#!/usr/bin/env python

import logging
import sys
import time

from .args_processing import process_args
from .boundaries import define_boundaries, load_tir_pwms
from .captain import detect_captains
from .cargo import extract_cargo
from .classify import classify_starships
from .confidence import score_starships
from .files import load_genome
from .helpers import compute_genome_stats
from .logger import logger, log_file_logger
from .models import StarshipResult, StarKITRun
from .parser import create_parser
from .report import generate_report
from .settings import CAPTAIN_HMM_DIR, BOUNDARY_DATA_DIR, FAMILY_HMM_DIR
from .version import __version__ as current_version
from .write import write_user_args, write_output_stats, write_tsv

import os


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
):
    """Core pipeline: parse input → detect captains → find boundaries →
    extract cargo → classify → score confidence → return results."""

    # Step 1: Parse & normalize input
    records = load_genome(input_file, gff_file)
    genome_stats = compute_genome_stats(records)

    # Step 2: Detect captain genes
    captain_hits = detect_captains(records, CAPTAIN_HMM_DIR, evalue)

    if not captain_hits:
        logger.info("No captain genes detected. No Starships predicted.")
        return StarKITRun(
            input_file=input_file,
            genome_stats=genome_stats,
            starships=[],
            parameters=dict(
                output_prefix=output_prefix,
                evalue=evalue,
                min_size=min_size,
                max_size=max_size,
                evidence=evidence,
            ),
            version=current_version,
        )

    # Load TIR PWMs for boundary detection
    pwm_file = os.path.join(BOUNDARY_DATA_DIR, "tir_pwms.json")
    pwm_data = load_tir_pwms(pwm_file)

    # Build a dict of records by id for fast lookup
    records_by_id = {rec.id: rec for rec in records}

    # Step 3-4: Define boundaries and extract cargo for each captain
    starship_results = []
    for idx, captain in enumerate(captain_hits, 1):
        record = records_by_id.get(captain.contig_id)
        if record is None:
            logger.warning(
                f"Contig '{captain.contig_id}' not found for captain "
                f"'{captain.protein_id}'; skipping"
            )
            continue

        # Step 3: Define boundaries
        boundary = define_boundaries(captain, record, min_size, max_size, pwm_data)

        # Step 4: Extract cargo genes
        cargo_genes = extract_cargo(
            record, boundary["start"], boundary["end"], captain.feature
        )

        starship_id = f"starship_{idx:03d}"

        result = StarshipResult(
            starship_id=starship_id,
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
        )
        starship_results.append(result)

    # Step 5: Classify families
    family_hmm_dir = FAMILY_HMM_DIR
    if os.path.isdir(family_hmm_dir) and os.listdir(family_hmm_dir):
        classify_starships(starship_results, records, family_hmm_dir)

    # Step 6: Score confidence
    score_starships(starship_results)

    # Filter by evidence level
    if evidence != "all":
        from .models import EvidenceLevel

        level_map = {
            "high": [EvidenceLevel.HIGH],
            "medium": [EvidenceLevel.HIGH, EvidenceLevel.MEDIUM],
            "low": [EvidenceLevel.HIGH, EvidenceLevel.MEDIUM, EvidenceLevel.LOW],
        }
        allowed = level_map.get(evidence, [e for e in EvidenceLevel])
        starship_results = [r for r in starship_results if r.evidence_level in allowed]

    starkit_run = StarKITRun(
        input_file=input_file,
        genome_stats=genome_stats,
        starships=starship_results,
        parameters=dict(
            output_prefix=output_prefix,
            evalue=evalue,
            min_size=min_size,
            max_size=max_size,
            evidence=evidence,
        ),
        version=current_version,
    )

    return starkit_run


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
    starkit_run = run(
        input_file, output_prefix, gff_file,
        evalue, min_size, max_size, evidence,
        use_log, quiet,
    )

    # Step 7: Write output
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
