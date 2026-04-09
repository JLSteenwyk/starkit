from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from starkit.confidence import compute_evidence_level, compute_confidence_score, score_starships
from starkit.models import EvidenceLevel, StarshipResult, CaptainHit, TIR


def _make_result(tir_left=None, tir_right=None, tsd=None, evalue=1e-50,
                 size=100000, homology_identity=0.0, homology_coverage=0.0,
                 boundary_method="homology"):
    record = SeqRecord(Seq("A" * 200000), id="c1")
    feature = SeqFeature(location=SimpleLocation(0, 100, strand=1), type="CDS")
    captain = CaptainHit(
        feature=feature, contig_id="c1", start=0, end=100,
        strand=1, evalue=evalue, score=50.0, hmm_name="test", protein_id="p1",
    )
    return StarshipResult(
        starship_id="s1", contig_id="c1", region=record,
        start=0, end=size, captain=captain,
        tir_left=tir_left, tir_right=tir_right, tsd=tsd,
        boundary_method=boundary_method,
        homology_identity=homology_identity,
        homology_coverage=homology_coverage,
    )


def test_evidence_high_captain_plus_homology():
    result = _make_result(evalue=1e-50, homology_identity=0.95, homology_coverage=0.90)
    assert compute_evidence_level(result) == EvidenceLevel.HIGH


def test_evidence_high_captain_plus_structural():
    tir_l = TIR(start=0, end=50, sequence="ATGC" * 12)
    result = _make_result(tir_left=tir_l, tsd="AAGCC")
    assert compute_evidence_level(result) == EvidenceLevel.HIGH


def test_evidence_medium_captain_only():
    result = _make_result(evalue=1e-50)
    assert compute_evidence_level(result) == EvidenceLevel.MEDIUM


def test_evidence_medium_strong_homology():
    result = _make_result(evalue=1.0, homology_identity=0.95, homology_coverage=0.90)
    assert compute_evidence_level(result) == EvidenceLevel.MEDIUM


def test_evidence_low_nothing():
    result = _make_result(evalue=1.0)
    assert compute_evidence_level(result) == EvidenceLevel.LOW


def test_confidence_score_range():
    result = _make_result()
    score = compute_confidence_score(result)
    assert 0.0 <= score <= 1.0


def test_confidence_captain_plus_homology_is_high():
    result = _make_result(
        evalue=1e-100, homology_identity=1.0, homology_coverage=1.0,
        boundary_method="homology",
    )
    score = compute_confidence_score(result)
    assert score >= 0.8


def test_confidence_homology_only_is_reasonable():
    result = _make_result(
        evalue=1.0, homology_identity=0.95, homology_coverage=0.90,
        boundary_method="homology",
    )
    score = compute_confidence_score(result)
    # Should be ~0.3*0 + 0.3*0.855 + 0.2*1.0 + 0.2*0 = 0.4565
    assert score >= 0.4


def test_confidence_ordering():
    """Captain+homology > homology-only > captain-only > nothing."""
    both = _make_result(evalue=1e-100, homology_identity=1.0, homology_coverage=1.0,
                        boundary_method="homology")
    captain = _make_result(evalue=1e-50, boundary_method="estimated")
    homology = _make_result(evalue=1.0, homology_identity=0.95, homology_coverage=0.90,
                            boundary_method="homology")
    nothing = _make_result(evalue=1.0, boundary_method="estimated")

    assert compute_confidence_score(both) > compute_confidence_score(homology)
    assert compute_confidence_score(homology) > compute_confidence_score(captain)
    assert compute_confidence_score(captain) > compute_confidence_score(nothing)


def test_score_starships():
    result = _make_result(evalue=1e-50, homology_identity=0.9, homology_coverage=0.8)
    scored = score_starships([result])
    assert scored[0].evidence_level == EvidenceLevel.HIGH
    assert 0.0 <= scored[0].confidence_score <= 1.0


def test_tier3_estimated_always_low():
    """Estimated boundaries (Tier 3) should always get LOW evidence."""
    result = _make_result(
        evalue=1e-100, homology_identity=0.95, homology_coverage=0.90,
        boundary_method="estimated",
    )
    assert compute_evidence_level(result) == EvidenceLevel.LOW
