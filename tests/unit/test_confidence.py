import math
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from starkit.confidence import compute_evidence_level, compute_confidence_score, score_starships
from starkit.models import EvidenceLevel, StarshipResult, CaptainHit, TIR


def _make_result(tir_left=None, tir_right=None, tsd=None, evalue=1e-50, size=100000):
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
    )


def test_evidence_high():
    tir_l = TIR(start=0, end=50, sequence="ATGC" * 12)
    tir_r = TIR(start=99950, end=100000, sequence="GCAT" * 12)
    result = _make_result(tir_left=tir_l, tir_right=tir_r, tsd="TTACAAAAAAA")
    assert compute_evidence_level(result) == EvidenceLevel.HIGH


def test_evidence_medium_one_tir():
    tir_l = TIR(start=0, end=50, sequence="ATGC" * 12)
    result = _make_result(tir_left=tir_l, tir_right=None, tsd=None)
    assert compute_evidence_level(result) == EvidenceLevel.MEDIUM


def test_evidence_medium_tsd_only():
    result = _make_result(tsd="TTACAAAAAAA")
    assert compute_evidence_level(result) == EvidenceLevel.MEDIUM


def test_evidence_low():
    result = _make_result()
    assert compute_evidence_level(result) == EvidenceLevel.LOW


def test_confidence_score_range():
    result = _make_result()
    score = compute_confidence_score(result)
    assert 0.0 <= score <= 1.0


def test_confidence_high_is_higher():
    tir_l = TIR(start=0, end=50, sequence="ATGC" * 12)
    tir_r = TIR(start=99950, end=100000, sequence="GCAT" * 12)
    high = _make_result(tir_left=tir_l, tir_right=tir_r, tsd="TTAC", evalue=1e-100)
    low = _make_result(evalue=1e-11)
    assert compute_confidence_score(high) > compute_confidence_score(low)


def test_score_starships():
    result = _make_result()
    scored = score_starships([result])
    assert scored[0].evidence_level == EvidenceLevel.LOW
    assert 0.0 <= scored[0].confidence_score <= 1.0
