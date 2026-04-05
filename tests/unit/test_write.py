import os
import tempfile

from starkit.write import write_tsv


def test_write_tsv(sample_starkit_run):
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test_output")
        write_tsv(sample_starkit_run, prefix)

        tsv_path = f"{prefix}.tsv"
        assert os.path.exists(tsv_path)

        with open(tsv_path) as f:
            lines = f.readlines()

        # Header + 1 data row
        assert len(lines) == 2

        header = lines[0].strip().split("\t")
        assert "starship_id" in header
        assert "contig" in header
        assert "cargo_gene_count" in header

        data = lines[1].strip().split("\t")
        assert data[0] == "starship_001"
        assert data[1] == "contig_1"
