import logging
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .exceptions import InvalidInputFileFormat, MissingGFFError, NoProteinsError

logger = logging.getLogger(__name__)

GENBANK_EXTENSIONS = {".gbk", ".gb", ".genbank"}
FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}


def detect_input_format(input_file: str) -> str:
    """Determine whether input file is GenBank or FASTA based on extension.

    Returns 'genbank' or 'fasta'.
    Raises InvalidInputFileFormat if the extension is not recognized.
    """
    _, ext = os.path.splitext(input_file.lower())
    if ext in GENBANK_EXTENSIONS:
        return "genbank"
    elif ext in FASTA_EXTENSIONS:
        return "fasta"
    else:
        raise InvalidInputFileFormat(
            f"Unrecognized file extension '{ext}'. "
            f"Accepted GenBank extensions: {sorted(GENBANK_EXTENSIONS)}. "
            f"Accepted FASTA extensions: {sorted(FASTA_EXTENSIONS)}."
        )


def parse_genbank(input_file: str) -> list[SeqRecord]:
    """Read a GenBank file and return a list of SeqRecord objects."""
    records = list(SeqIO.parse(input_file, "genbank"))
    if not records:
        raise InvalidInputFileFormat(
            f"No records found in GenBank file: {input_file}"
        )
    logger.info(f"Parsed {len(records)} contig(s) from GenBank file")
    return records


def _parse_gff_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3 column 9 (semicolon-delimited key=value pairs) into a dict."""
    attributes = {}
    for item in attr_string.strip().rstrip(";").split(";"):
        item = item.strip()
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key.strip()] = value.strip()
    return attributes


def parse_fasta_gff(fasta_file: str, gff_file: str) -> list[SeqRecord]:
    """Read a FASTA file and GFF3 annotation, returning SeqRecords with CDS features.

    CDS features are added to the corresponding SeqRecord based on contig/seqid.
    Protein translations are extracted from the nucleotide sequence if not
    otherwise available.
    """
    # Read FASTA sequences into a dict keyed by record id
    records_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    if not records_dict:
        raise InvalidInputFileFormat(
            f"No sequences found in FASTA file: {fasta_file}"
        )

    # Parse GFF3 and add CDS features to the appropriate SeqRecord
    with open(gff_file) as gff:
        for line in gff:
            line = line.strip()
            # Skip comments and blank lines
            if not line or line.startswith("#"):
                continue
            # GFF3 directive to stop parsing
            if line.startswith("##FASTA"):
                break

            cols = line.split("\t")
            if len(cols) != 9:
                continue

            seqid = cols[0]
            feature_type = cols[2]

            # We only care about CDS features
            if feature_type != "CDS":
                continue

            # Skip features on contigs not in our FASTA
            if seqid not in records_dict:
                logger.warning(
                    f"GFF3 references seqid '{seqid}' not found in FASTA; skipping"
                )
                continue

            try:
                start = int(cols[3]) - 1  # GFF3 is 1-based, BioPython is 0-based
                end = int(cols[4])
            except ValueError:
                continue

            strand_str = cols[6]
            if strand_str == "+":
                strand = 1
            elif strand_str == "-":
                strand = -1
            else:
                strand = 0

            # Parse attributes from column 9
            attrs = _parse_gff_attributes(cols[8])

            qualifiers = {}
            if "ID" in attrs:
                qualifiers["ID"] = [attrs["ID"]]
            if "Name" in attrs:
                qualifiers["Name"] = [attrs["Name"]]
            if "protein_id" in attrs:
                qualifiers["protein_id"] = [attrs["protein_id"]]
            if "product" in attrs:
                qualifiers["product"] = [attrs["product"]]

            # Build the location and translate CDS to protein
            location = SimpleLocation(start, end, strand=strand)
            record = records_dict[seqid]
            try:
                nuc_seq = location.extract(record.seq)
                protein_seq = nuc_seq.translate(table=1, to_stop=True)
                qualifiers["translation"] = [str(protein_seq)]
            except Exception as e:
                logger.warning(
                    f"Could not translate CDS at {seqid}:{start+1}-{end}: {e}"
                )

            # Create the SeqFeature after qualifiers are fully populated
            feature = SeqFeature(location=location, type="CDS", qualifiers=qualifiers)
            record.features.append(feature)

    records = list(records_dict.values())
    logger.info(
        f"Parsed {len(records)} contig(s) from FASTA + GFF3"
    )
    return records


def _has_protein_sequences(records: list[SeqRecord]) -> bool:
    """Check that at least one CDS feature with a translation exists."""
    for record in records:
        for feature in record.features:
            if feature.type == "CDS":
                translations = feature.qualifiers.get("translation", [])
                if translations and translations[0]:
                    return True
    return False


def load_genome(input_file: str, gff_file: str = None) -> list[SeqRecord]:
    """Main entry point for loading genome data.

    Detects the input format, calls the appropriate parser, validates that
    protein sequences exist, and returns a list of SeqRecord objects.
    """
    fmt = detect_input_format(input_file)

    if fmt == "genbank":
        records = parse_genbank(input_file)
    elif fmt == "fasta":
        if gff_file is None:
            raise MissingGFFError(
                "A GFF3 annotation file (--gff) is required when the input is FASTA."
            )
        records = parse_fasta_gff(input_file, gff_file)
    else:
        raise InvalidInputFileFormat(f"Unsupported format: {fmt}")

    if not _has_protein_sequences(records):
        raise NoProteinsError(
            "No CDS features with protein translations found in the input. "
            "Ensure the genome annotation includes CDS features with translations."
        )

    return records
