import sys
import pytest

from starkit.parser import create_parser


def test_parser_basic(monkeypatch):
    """Test that parser accepts a positional input argument."""
    monkeypatch.setattr(sys, "argv", ["starkit", "test.gbk"])
    parser = create_parser()
    args = parser.parse_args(["test.gbk"])
    assert args.input == "test.gbk"


def test_parser_with_gff(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["starkit", "test.fa", "--gff", "test.gff3"])
    parser = create_parser()
    args = parser.parse_args(["test.fa", "--gff", "test.gff3"])
    assert args.input == "test.fa"
    assert args.gff == "test.gff3"


def test_parser_with_options(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["starkit", "test.gbk", "-e", "1e-5"])
    parser = create_parser()
    args = parser.parse_args(["test.gbk", "-e", "1e-5", "--min-size", "20000"])
    assert args.evalue == 1e-5
    assert args.min_size == 20000


def test_parser_exits_with_no_args(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["starkit"])
    with pytest.raises(SystemExit):
        create_parser()
