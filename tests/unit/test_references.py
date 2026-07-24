import os

from starkit.references import parse_ref_header
from starkit.settings import STARSHIP_REF_FASTAS


def test_builtin_reference_libraries_include_curated_fastas():
    expected_names = {
        "starships_ref.fasta",
        "aspergillus_representative_elements.fna",
        "Vogan_Regions_260610.fna",
        "Vogan_Starships_260610.fna",
    }

    assert {os.path.basename(path) for path in STARSHIP_REF_FASTAS} == expected_names
    assert all(os.path.isfile(path) for path in STARSHIP_REF_FASTAS[1:])


def test_parse_starbase_reference_header():
    assert parse_ref_header("42|Prometheus|32787bp") == (
        "42",
        "Prometheus",
        32787,
    )


def test_parse_aspergillus_reference_header():
    header = "aspacu8_s00007|-_navis0001hap0435"
    assert parse_ref_header(header) == (header, "unclassified", 0)


def test_parse_vogan_starship_family_headers():
    assert parse_ref_header("Prometheus-fam_Rheap_MPI-PUGE-AT-0058") == (
        "Prometheus-fam_Rheap_MPI-PUGE-AT-0058",
        "Prometheus",
        0,
    )
    assert parse_ref_header("Tardis_fam_Pencan_IBT18980") == (
        "Tardis_fam_Pencan_IBT18980",
        "Tardis",
        0,
    )
    assert parse_ref_header("Phoenix_Mpha") == (
        "Phoenix_Mpha",
        "Phoenix",
        0,
    )
