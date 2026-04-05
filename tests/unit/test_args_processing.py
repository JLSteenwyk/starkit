import os
import tempfile
import pytest
from types import SimpleNamespace

from starkit.args_processing import process_args


def test_process_args_defaults():
    with tempfile.NamedTemporaryFile(suffix=".gbk", delete=False) as f:
        f.write(b"test")
        tmpfile = f.name

    try:
        args = SimpleNamespace(
            input=tmpfile,
            output=None,
            gff=None,
            evalue=None,
            min_size=None,
            max_size=None,
            evidence=None,
            log=False,
            quiet=False,
        )
        result = process_args(args)
        assert result["input_file"] == tmpfile
        assert result["output_prefix"] == f"{tmpfile}.starkit"
        assert result["evalue"] == 1e-10
        assert result["min_size"] == 15000
        assert result["max_size"] == 700000
        assert result["evidence"] == "all"
    finally:
        os.unlink(tmpfile)


def test_process_args_custom():
    with tempfile.NamedTemporaryFile(suffix=".gbk", delete=False) as f:
        f.write(b"test")
        tmpfile = f.name

    try:
        args = SimpleNamespace(
            input=tmpfile,
            output="custom_output",
            gff=None,
            evalue=1e-5,
            min_size=20000,
            max_size=500000,
            evidence="high",
            log=True,
            quiet=False,
        )
        result = process_args(args)
        assert result["output_prefix"] == "custom_output"
        assert result["evalue"] == 1e-5
        assert result["min_size"] == 20000
        assert result["evidence"] == "high"
    finally:
        os.unlink(tmpfile)
