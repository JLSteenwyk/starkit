import os

# Default thresholds
DEFAULT_EVALUE = 1e-10
DEFAULT_MIN_SIZE = 0
DEFAULT_MAX_SIZE = 0  # 0 = no limit
DEFAULT_UPSTREAM_SCAN = 10000
DEFAULT_TSD_MOTIF = "TTAC"
DEFAULT_TSD_CONSENSUS_REGEX = r"TTAC.{7}A"

# Evidence level thresholds
CONFIDENCE_CAPTAIN_WEIGHT = 0.4
CONFIDENCE_TIR_WEIGHT = 0.3
CONFIDENCE_TSD_WEIGHT = 0.2
CONFIDENCE_SIZE_WEIGHT = 0.1

# De novo TIR detection
MIN_TIR_LENGTH = 10
MIN_TIR_IDENTITY = 0.70
TIR_SCAN_WINDOW = 500

# Starship size distribution (for confidence scoring)
MEDIAN_STARSHIP_SIZE = 110000
STARSHIP_SIZE_STD = 80000

# Data paths (relative to package)
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
CAPTAIN_HMM_DIR = os.path.join(DATA_DIR, "hmm")
BOUNDARY_DATA_DIR = os.path.join(DATA_DIR, "boundaries")
FAMILY_HMM_DIR = os.path.join(DATA_DIR, "families")
STARSHIP_REF_FASTA = os.path.join(DATA_DIR, "starships_ref.fasta")
