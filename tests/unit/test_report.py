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
