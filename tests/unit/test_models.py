from starkit.models import EvidenceLevel, StarshipResult, StarKITRun, CaptainHit, TIR, CargoGene


def test_evidence_level_values():
    assert EvidenceLevel.HIGH.value == "high"
    assert EvidenceLevel.MEDIUM.value == "medium"
    assert EvidenceLevel.LOW.value == "low"


def test_starship_result_size(sample_starship_result):
    assert sample_starship_result.size == 16000  # 20000 - 4000


def test_starship_result_defaults():
    """Test that defaults work for optional fields."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    record = SeqRecord(Seq("ATGC" * 100), id="test")
    feature = SeqFeature(location=SimpleLocation(0, 100, strand=1), type="CDS")
    captain = CaptainHit(
        feature=feature, contig_id="test", start=0, end=100,
        strand=1, evalue=1e-10, score=50.0, hmm_name="test", protein_id="p1",
    )
    result = StarshipResult(
        starship_id="s1", contig_id="test", region=record,
        start=0, end=1000, captain=captain,
    )
    assert result.captain_family == "unclassified"
    assert result.family_score == 0.0
    assert result.tir_left is None
    assert result.tir_right is None
    assert result.tsd is None
    assert result.cargo_genes == []
    assert result.confidence_score == 0.0
    assert result.evidence_level == EvidenceLevel.LOW
    assert result.truncated is False


def test_starkit_run(sample_starkit_run):
    assert sample_starkit_run.version == "0.1.0"
    assert len(sample_starkit_run.starships) == 1
    assert sample_starkit_run.genome_stats["contig_count"] == 1
