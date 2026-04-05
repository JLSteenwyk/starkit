import logging
import time

from .models import EvidenceLevel, StarKITRun

logger = logging.getLogger(__name__)


def write_user_args(
    input_file,
    output_prefix,
    evalue,
    min_size,
    max_size,
    evidence,
    gff_file,
    use_log,
    quiet,
):
    """Print the run parameters to stdout using the logger."""
    logger.info(f"\n  Input file: {input_file}")
    logger.info(f"  Output prefix: {output_prefix}")
    logger.info(f"  GFF3 file: {gff_file}")
    logger.info(f"  E-value threshold: {evalue}")
    min_str = f"{min_size}" if min_size else "none"
    max_str = f"{max_size}" if max_size else "none"
    logger.info(f"  Size range: {min_str} - {max_str} bp")
    logger.info(f"  Evidence filter: {evidence}")


def write_output_stats(starkit_run: StarKITRun, start_time: float):
    """Print summary statistics about the completed run."""
    starships = starkit_run.starships
    total = len(starships)

    high = sum(1 for s in starships if s.evidence_level == EvidenceLevel.HIGH)
    medium = sum(1 for s in starships if s.evidence_level == EvidenceLevel.MEDIUM)
    low = sum(1 for s in starships if s.evidence_level == EvidenceLevel.LOW)

    output_prefix = starkit_run.parameters.get("output_prefix", starkit_run.input_file + ".starkit")
    elapsed = time.time() - start_time

    logger.info(f"\n  Starships found: {total}")
    logger.info(f"    High confidence: {high}")
    logger.info(f"    Medium confidence: {medium}")
    logger.info(f"    Low confidence: {low}")
    logger.info("")
    logger.info("  Output files:")
    logger.info(f"    {output_prefix}.tsv")
    logger.info(f"    {output_prefix}.html")
    logger.info(f"\n  Execution time: {elapsed:.2f}s")


def write_tsv(starkit_run: StarKITRun, output_prefix: str):
    """Write Starship results to a TSV file at {output_prefix}.tsv."""
    output_file = f"{output_prefix}.tsv"

    header = "\t".join([
        "starship_id",
        "contig",
        "start",
        "end",
        "size",
        "captain_gene",
        "captain_evalue",
        "family",
        "family_score",
        "evidence_level",
        "confidence_score",
        "tir_left",
        "tir_right",
        "tsd_sequence",
        "cargo_gene_count",
        "truncated",
    ])

    with open(output_file, "w") as f:
        f.write(header + "\n")

        for result in starkit_run.starships:
            tir_left = (
                f"{result.tir_left.start}-{result.tir_left.end}"
                if result.tir_left
                else "NA"
            )
            tir_right = (
                f"{result.tir_right.start}-{result.tir_right.end}"
                if result.tir_right
                else "NA"
            )
            tsd_sequence = result.tsd or "NA"

            row = "\t".join([
                result.starship_id,
                result.contig_id,
                str(result.start),
                str(result.end),
                str(result.size),
                result.captain.protein_id,
                f"{result.captain.evalue:.2e}",
                result.captain_family,
                f"{result.family_score:.1f}",
                result.evidence_level.value,
                f"{result.confidence_score:.4f}",
                tir_left,
                tir_right,
                tsd_sequence,
                str(len(result.cargo_genes)),
                str(result.truncated),
            ])
            f.write(row + "\n")

    logger.info(f"  Results written to {output_file}")
