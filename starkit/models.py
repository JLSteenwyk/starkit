from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


class EvidenceLevel(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class CaptainHit:
    feature: SeqFeature
    contig_id: str
    start: int
    end: int
    strand: int
    evalue: float
    score: float
    hmm_name: str
    protein_id: str


@dataclass
class TIR:
    start: int
    end: int
    sequence: str


@dataclass
class CargoGene:
    gene_id: str
    product: str
    start: int
    end: int
    strand: int


@dataclass
class StarshipResult:
    starship_id: str
    contig_id: str
    region: SeqRecord
    start: int
    end: int
    captain: CaptainHit
    captain_family: str = "unclassified"
    family_score: float = 0.0
    tir_left: Optional[TIR] = None
    tir_right: Optional[TIR] = None
    tsd: Optional[str] = None
    cargo_genes: list = field(default_factory=list)
    confidence_score: float = 0.0
    evidence_level: EvidenceLevel = EvidenceLevel.LOW
    truncated: bool = False
    boundary_method: str = "estimated"  # "homology", "dr_motif", or "estimated"
    homology_identity: float = 0.0     # 0-1, from homology alignment
    homology_coverage: float = 0.0     # 0-1, from homology alignment
    nested_in: Optional[str] = None    # starship_id of parent if nested
    additional_captains: list = field(default_factory=list)  # extra CaptainHits merged

    @property
    def size(self) -> int:
        return self.end - self.start

    @property
    def is_novel(self) -> bool:
        return not (self.homology_identity >= 0.80 and self.homology_coverage >= 0.50)


@dataclass
class StarKITRun:
    input_file: str
    genome_stats: dict
    starships: list
    parameters: dict
    version: str
