#!/usr/bin/env python

import logging
import os
import sys
import tempfile
import time

from Bio import SeqIO

from .args_processing import process_args
from .boundaries import define_boundaries, load_tir_pwms
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
    )


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
    """Core pipeline: parse input -> detect captains -> homology search ->
    find boundaries -> extract cargo -> classify -> score -> return results."""

    # Step 1: Parse & normalize input
    records = load_genome(input_file, gff_file)
    genome_stats = compute_genome_stats(records)

    # Step 2: Detect captain genes
    captain_hits = detect_captains(records, CAPTAIN_HMM_DIR, evalue)

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

    # Step 5: Homology-based detection
    genome_fasta = None
    try:
        if os.path.exists(STARSHIP_REF_FASTA):
            genome_fasta = _write_temp_fasta(records)
            novel_hits = detect_by_homology(
                genome_fasta, STARSHIP_REF_FASTA,
                captain_starships=starship_results,
                min_identity=0.80,
                min_coverage=0.50,
            )

            # Convert novel homology hits to StarshipResults
            next_idx = len(starship_results) + 1
            for hit in novel_hits:
                # Filter by size
                hit_size = hit.end - hit.start
                if hit_size < min_size or hit_size > max_size:
                    continue

                record = records_by_id.get(hit.contig_id)
                if record is None:
                    continue

                cargo_genes = extract_cargo(record, hit.start, hit.end, None)

                result = _homology_hit_to_starship(
                    hit, record, next_idx, cargo_genes,
                )
                starship_results.append(result)
                next_idx += 1

            if novel_hits:
                logger.info(
                    f"Added {next_idx - len(captain_hits) - 1} Starship(s) "
                    f"from homology search"
                )
        else:
            logger.info(
                "Reference Starship FASTA not found; skipping homology search"
            )
    finally:
        if genome_fasta and os.path.exists(genome_fasta):
            os.unlink(genome_fasta)

    if not starship_results:
        logger.info("No Starships predicted.")
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

    # Step 6: Classify families (for captain-based hits that aren't classified yet)
    family_hmm_dir = FAMILY_HMM_DIR
    if os.path.isdir(family_hmm_dir) and os.listdir(family_hmm_dir):
        # Only classify captain-based results (homology ones already have families)
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
