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
    logger.info(f"    {output_prefix}.fasta")
    logger.info(f"    {output_prefix}.bed")
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
        "boundary_method",
        "tir_left",
        "tir_right",
        "dr_sequence",
        "cargo_gene_count",
        "truncated",
        "captain_orientation_flag",
        "captain_truncated_flag",
        "novelty",
        "nested_in",
        "classification_status",
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

            novelty = "new" if result.is_novel else "known"

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
                result.boundary_method,
                tir_left,
                tir_right,
                tsd_sequence,
                str(len(result.cargo_genes)),
                str(result.truncated),
                str(result.captain_orientation_flag),
                str(result.captain_truncated_flag),
                novelty,
                result.nested_in or "NA",
                result.classification_status,
            ])
            f.write(row + "\n")

    logger.info(f"  Results written to {output_file}")


def write_fasta(starkit_run: StarKITRun, output_prefix: str):
    """Write predicted Starship nucleotide sequences to a FASTA file.

    Each Starship is extracted from its contig and written as a separate
    entry. Header format: >starship_id contig:start-end family=X size=Xbp
    """
    output_file = f"{output_prefix}.fasta"

    with open(output_file, "w") as f:
        for result in starkit_run.starships:
            seq = str(result.region.seq[result.start:result.end])
            header = (
                f">{result.starship_id} "
                f"{result.contig_id}:{result.start}-{result.end} "
                f"family={result.captain_family} "
                f"size={result.size}bp "
                f"evidence={result.evidence_level.value}"
            )
            f.write(header + "\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    logger.info(f"  Sequences written to {output_file}")


def write_bed(starkit_run: StarKITRun, output_prefix: str):
    """Write Starship features to a BED-like file at {output_prefix}.bed.

    Each Starship produces rows for:
      - captain gene (tag="cap")
      - cargo genes (tag=".")
      - TIR left/right if present (tag="flank", featureID has |TIR|up / |TIR|down)
      - DR left/right if present (tag="flank", featureID has |DR|up / |DR|down)
    """
    output_file = f"{output_prefix}.bed"

    header = "\t".join([
        "contigID", "begin", "end", "featureID",
        "tag", "strand", "starshipID", "annotation",
    ])

    def _strand_str(strand_val):
        if strand_val == 1:
            return "+"
        elif strand_val == -1:
            return "-"
        return "."

    with open(output_file, "w") as f:
        f.write(header + "\n")

        for result in starkit_run.starships:
            sid = result.starship_id
            contig = result.contig_id
            cap = result.captain
            cap_strand = _strand_str(cap.strand)

            # Collect all rows for this Starship, then sort by start position
            rows = []

            # DR left (upstream) — sits just outside the TIR on the upstream side
            if result.tsd is not None and result.tir_left is not None:
                dr_seq = result.tsd.lower()
                dr_start = result.tir_left.start - len(dr_seq)
                dr_end = result.tir_left.start
                rows.append((dr_start, dr_end, f"{cap.protein_id}|DR|up",
                             "flank", cap_strand, dr_seq))

            # TIR left
            if result.tir_left is not None:
                tir = result.tir_left
                rows.append((tir.start, tir.end, f"{cap.protein_id}|TIR|up",
                             "flank", cap_strand, tir.sequence))

            # Captain gene
            rows.append((cap.start, cap.end, cap.protein_id,
                         "cap", cap_strand, "."))

            # Cargo genes
            for cargo in result.cargo_genes:
                cargo_strand = _strand_str(cargo.strand)
                annotation = cargo.product if cargo.product != "hypothetical protein" else "."
                rows.append((cargo.start, cargo.end, cargo.gene_id,
                             ".", cargo_strand, annotation))

            # TIR right
            if result.tir_right is not None:
                tir = result.tir_right
                rows.append((tir.start, tir.end, f"{cap.protein_id}|TIR|down",
                             "flank", cap_strand, tir.sequence))

            # DR right (downstream) — sits just outside the TIR on the downstream side
            if result.tsd is not None and result.tir_right is not None:
                dr_seq = result.tsd.lower()
                dr_start = result.tir_right.end
                dr_end = result.tir_right.end + len(dr_seq)
                rows.append((dr_start, dr_end, f"{cap.protein_id}|DR|down",
                             "flank", cap_strand, dr_seq))

            # Sort all rows by start coordinate and write
            rows.sort(key=lambda r: r[0])
            for start, end, feat_id, tag, strand, annot in rows:
                line = "\t".join([
                    contig, str(start), str(end), feat_id,
                    tag, strand, sid, annot,
                ])
                f.write(line + "\n")

    logger.info(f"  BED features written to {output_file}")
