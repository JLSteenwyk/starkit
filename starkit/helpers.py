import logging

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

COMPLEMENT_MAP = str.maketrans("ATCGatcg", "TAGCtagc")


def get_protein_sequences(
    records: list[SeqRecord],
) -> list[tuple[str, str, str, SeqFeature]]:
    """Extract all CDS features with protein translations from SeqRecords.

    Returns a list of tuples:
        (protein_id, protein_sequence_str, contig_id, feature)
    """
    proteins = []
    unnamed_count = 0

    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue

            translations = feature.qualifiers.get("translation", [])
            if not translations or not translations[0]:
                continue

            protein_seq = translations[0]

            # Determine a protein identifier
            protein_id = None
            for key in ("protein_id", "locus_tag", "ID", "Name", "gene"):
                vals = feature.qualifiers.get(key, [])
                if vals and vals[0]:
                    protein_id = vals[0]
                    break

            if protein_id is None:
                unnamed_count += 1
                protein_id = f"{record.id}_CDS_{unnamed_count}"

            proteins.append((protein_id, protein_seq, record.id, feature))

    logger.info(f"Extracted {len(proteins)} protein sequence(s)")
    return proteins


def reverse_complement(seq_str: str) -> str:
    """Return the reverse complement of a DNA sequence string."""
    return seq_str.translate(COMPLEMENT_MAP)[::-1]


def compute_genome_stats(records: list[SeqRecord]) -> dict:
    """Compute basic genome statistics from a list of SeqRecords.

    Returns a dict with keys:
        contig_count  - number of contigs/records
        total_length  - sum of all contig lengths in bp
        gc_content    - GC fraction (0.0-1.0) across the entire genome
    """
    contig_count = len(records)
    total_length = 0
    gc_count = 0

    for record in records:
        seq_str = str(record.seq).upper()
        total_length += len(seq_str)
        gc_count += seq_str.count("G") + seq_str.count("C")

    gc_content = gc_count / total_length if total_length > 0 else 0.0

    contig_lengths = {record.id: len(record.seq) for record in records}

    return {
        "contig_count": contig_count,
        "total_length": total_length,
        "gc_content": gc_content,
        "contig_lengths": contig_lengths,
    }
