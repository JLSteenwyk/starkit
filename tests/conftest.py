import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from starkit.models import (
    CaptainHit,
    CargoGene,
    EvidenceLevel,
    StarshipResult,
    StarKITRun,
    TIR,
)


SAMPLES_DIR = os.path.join(os.path.dirname(__file__), "integration", "samples")


@pytest.fixture
def sample_genbank():
    """Path to the test GenBank file."""
    return os.path.join(SAMPLES_DIR, "test_genome.gbk")


@pytest.fixture
def simple_record():
    """A minimal SeqRecord with CDS features for unit testing."""
    # Create a 200kb contig with some sequence
    seq = Seq("ATGCTTAC" * 25000)  # 200kb
    record = SeqRecord(seq, id="contig_1", name="contig_1")

    # Add a captain gene CDS at position 5000-6500 (+ strand)
    captain_quals = {
        "protein_id": ["captain_001"],
        "product": ["tyrosine recombinase"],
        "translation": ["MAKTLRKHRYGSILD" + "A" * 100],
    }
    captain_feature = SeqFeature(
        location=SimpleLocation(5000, 6500, strand=1),
        type="CDS",
        qualifiers=captain_quals,
    )
    record.features.append(captain_feature)

    # Add cargo gene CDS features downstream
    for i, (start, end, product) in enumerate([
        (8000, 9500, "hypothetical protein"),
        (10000, 12000, "transporter"),
        (13000, 14500, "oxidoreductase"),
        (15000, 16000, "methyltransferase"),
    ]):
        quals = {
            "protein_id": [f"cargo_{i:03d}"],
            "product": [product],
            "translation": ["M" + "A" * 50],
        }
        feature = SeqFeature(
            location=SimpleLocation(start, end, strand=1 if i % 2 == 0 else -1),
            type="CDS",
            qualifiers=quals,
        )
        record.features.append(feature)

    return record


@pytest.fixture
def sample_captain_hit(simple_record):
    """A CaptainHit for the captain gene in simple_record."""
    captain_feature = simple_record.features[0]
    return CaptainHit(
        feature=captain_feature,
        contig_id="contig_1",
        start=5000,
        end=6500,
        strand=1,
        evalue=1e-50,
        score=200.0,
        hmm_name="captain_family_A",
        protein_id="captain_001",
    )


@pytest.fixture
def sample_starship_result(simple_record, sample_captain_hit):
    """A StarshipResult for testing."""
    return StarshipResult(
        starship_id="starship_001",
        contig_id="contig_1",
        region=simple_record,
        start=4000,
        end=20000,
        captain=sample_captain_hit,
        captain_family="Voyager",
        family_score=150.0,
        tir_left=TIR(start=4000, end=4050, sequence="ATGCATGCATGC" * 4),
        tir_right=TIR(start=19950, end=20000, sequence="GCATGCATGCAT" * 4),
        tsd="TTACAAGCCGTA",
        cargo_genes=[
            CargoGene("cargo_000", "hypothetical protein", 8000, 9500, 1),
            CargoGene("cargo_001", "transporter", 10000, 12000, -1),
            CargoGene("cargo_002", "oxidoreductase", 13000, 14500, 1),
            CargoGene("cargo_003", "methyltransferase", 15000, 16000, -1),
        ],
        confidence_score=0.85,
        evidence_level=EvidenceLevel.HIGH,
        truncated=False,
    )


@pytest.fixture
def sample_starkit_run(sample_starship_result):
    """A StarKITRun for testing output modules."""
    return StarKITRun(
        input_file="test_genome.gbk",
        genome_stats={
            "contig_count": 1,
            "total_length": 200000,
            "gc_content": 0.50,
        },
        starships=[sample_starship_result],
        parameters={
            "output_prefix": "test_output",
            "evalue": 1e-10,
            "min_size": 15000,
            "max_size": 700000,
            "evidence": "all",
        },
        version="0.1.0",
    )
