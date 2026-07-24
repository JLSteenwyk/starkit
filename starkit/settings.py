import os

# Default thresholds
DEFAULT_EVALUE = 1e-10
DEFAULT_MIN_SIZE = 0
DEFAULT_MAX_SIZE = 0  # 0 = no limit
DEFAULT_UPSTREAM_SCAN = 10000
DEFAULT_TSD_MOTIF = "TTAC"
DEFAULT_TSD_CONSENSUS_REGEX = r"TTAC.{7}A"

# Family classification
WEAK_CLASSIFICATION_THRESHOLD = 200  # family HMM score below this = "weakly classified"

# De novo TIR detection
MIN_TIR_LENGTH = 8
MIN_TIR_IDENTITY = 0.70
TIR_SCAN_WINDOW = 1000
TIR_SEED_K = 6  # k-mer seed size for initial matching

# Data paths (relative to package)
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
CAPTAIN_HMM_DIR = os.path.join(DATA_DIR, "hmm")
BOUNDARY_DATA_DIR = os.path.join(DATA_DIR, "boundaries")
FAMILY_HMM_DIR = os.path.join(DATA_DIR, "families")
STARSHIP_REF_FASTA = os.path.join(DATA_DIR, "starships_ref.fasta")
ASPERGILLUS_REF_FASTA = os.path.join(
    DATA_DIR, "aspergillus_representative_elements.fna",
)
VOGAN_REGIONS_REF_FASTA = os.path.join(
    DATA_DIR, "Vogan_Regions_260610.fna",
)
VOGAN_STARSHIPS_REF_FASTA = os.path.join(
    DATA_DIR, "Vogan_Starships_260610.fna",
)

# Reference libraries searched by default. STARSHIP_REF_FASTA is generated
# from Starbase and may be absent from lightweight installations; the three
# curated .fna libraries are shipped with StarKIT.
STARSHIP_REF_FASTAS = (
    STARSHIP_REF_FASTA,
    ASPERGILLUS_REF_FASTA,
    VOGAN_REGIONS_REF_FASTA,
    VOGAN_STARSHIPS_REF_FASTA,
)
