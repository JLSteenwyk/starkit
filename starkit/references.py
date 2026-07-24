"""Utilities for interpreting Starship reference-library records."""

from typing import Tuple


_NAMED_FAMILIES = {
    "Arwing",
    "Enterprise",
    "Galactica",
    "Hephaestus",
    "Moya",
    "Mystery",
    "Phoenix",
    "Prometheus",
    "Serenity",
    "Tardis",
    "Voyager",
}


def parse_ref_header(header: str) -> Tuple[str, str, int]:
    """Return ``(reference_id, family, length)`` for a FASTA header.

    Starbase records use ``ship_id|family|lengthbp``. The curated
    Aspergillus library uses a pipe to separate the element ID from an
    orientation/haplotype label, so that complete header is retained as the
    reference ID rather than treating the label as a family. Vogan Starship
    records encode the family at the start of the record name.
    """
    parts = header.split("|")

    # Aspergillus representative-element headers look like
    # "aspacu8_s00007|-_navis0001hap0435". The second field is not a family.
    if len(parts) == 2 and parts[1].startswith(("+_", "-_")):
        return header, "unclassified", 0

    reference_id = parts[0] if parts else header
    family = parts[1] if len(parts) > 1 else "unclassified"
    length = 0
    if len(parts) > 2:
        try:
            length = int(parts[2].removesuffix("bp"))
        except ValueError:
            length = 0

    if family == "unclassified":
        for marker in ("-fam_", "_fam_"):
            if marker in reference_id:
                candidate = reference_id.split(marker, 1)[0]
                if candidate in _NAMED_FAMILIES:
                    family = candidate
                break
        else:
            candidate = reference_id.split("_", 1)[0]
            if candidate in _NAMED_FAMILIES:
                family = candidate

    return reference_id, family, length
