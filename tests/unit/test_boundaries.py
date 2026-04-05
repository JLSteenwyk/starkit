from starkit.boundaries import find_tsd_motifs, find_tirs_denovo, load_tir_pwms


def test_find_tsd_motifs_basic():
    seq = "AAAATTACGGGTTACAAA"
    positions = find_tsd_motifs(seq, "TTAC")
    assert 4 in positions
    assert 11 in positions
    assert len(positions) == 2


def test_find_tsd_motifs_case_insensitive():
    seq = "AAAAttacGGG"
    positions = find_tsd_motifs(seq, "TTAC")
    assert 4 in positions


def test_find_tsd_motifs_none():
    seq = "AAAAAAGGGGGCCCCC"
    positions = find_tsd_motifs(seq, "TTAC")
    assert len(positions) == 0


def test_find_tirs_denovo_match():
    # Upstream has "ATGCATGCATGC" starting at position 0
    # Downstream has the reverse complement at the end
    tir_seq = "ATGCATGCATGC"
    from starkit.helpers import reverse_complement
    rc = reverse_complement(tir_seq)

    upstream = tir_seq + "N" * 488  # pad to 500bp
    downstream = "N" * 488 + rc

    tir_l, tir_r = find_tirs_denovo(
        upstream, downstream,
        min_length=10, min_identity=0.70, scan_window=500,
    )
    assert tir_l is not None
    assert tir_r is not None
    assert len(tir_l[2]) >= 10


def test_find_tirs_denovo_no_match():
    upstream = "A" * 500
    downstream = "C" * 500
    tir_l, tir_r = find_tirs_denovo(upstream, downstream, min_length=10)
    assert tir_l is None
    assert tir_r is None


def test_load_tir_pwms_missing_file():
    result = load_tir_pwms("/nonexistent/path.json")
    assert result == []
