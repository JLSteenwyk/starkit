import pytest

from starkit.files import detect_input_format
from starkit.exceptions import InvalidInputFileFormat


def test_detect_genbank():
    assert detect_input_format("genome.gbk") == "genbank"
    assert detect_input_format("genome.gb") == "genbank"
    assert detect_input_format("genome.genbank") == "genbank"


def test_detect_fasta():
    assert detect_input_format("genome.fasta") == "fasta"
    assert detect_input_format("genome.fa") == "fasta"
    assert detect_input_format("genome.fna") == "fasta"


def test_detect_invalid():
    with pytest.raises(InvalidInputFileFormat):
        detect_input_format("genome.txt")


def test_detect_case_insensitive():
    assert detect_input_format("genome.GBK") == "genbank"
    assert detect_input_format("genome.FASTA") == "fasta"
