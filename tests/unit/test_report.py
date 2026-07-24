import os
import tempfile

from starkit.report import generate_report, generate_svg_diagram


def test_generate_svg_diagram(sample_starship_result):
    svg = generate_svg_diagram(sample_starship_result)
    assert "<svg" in svg
    assert "</svg>" in svg
    # Should contain captain gene color
    assert "#c77c11" in svg or "c77c11" in svg


def test_generate_report(sample_starkit_run):
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test_output")
        path = generate_report(sample_starkit_run, prefix)

        assert os.path.exists(path)
        assert path.endswith(".html")

        with open(path) as f:
            html = f.read()

        assert "StarKIT" in html
        assert "starship_001" in html
        assert "<table" in html


def test_generate_report_for_non_starbase_reference(sample_starkit_run):
    result = sample_starkit_run.starships[0]
    result.homology_identity = 0.95
    result.homology_coverage = 0.90
    result.homology_reference_id = "Prometheus-fam_Rheap_MPI-PUGE-AT-0058"
    result.homology_reference_family = "Prometheus"

    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test_output")
        path = generate_report(sample_starkit_run, prefix)

        with open(path) as f:
            html = f.read()

        assert "matches a reference Starship" in html
        assert result.homology_reference_id in html
        assert "starbase.serve.scilifelab.se/ships/" not in html
