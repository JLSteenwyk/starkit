# StarKIT Design Spec

**Date:** 2026-04-05
**Author:** Jacob L. Steenwyk
**Status:** Draft

## Overview

StarKIT is a command-line tool for predicting Starship giant mobile genetic elements in single fungal genomes. It follows the ClipKIT architectural pattern: one main command, internal modules for each pipeline step, dataclasses for results.

Starships are large DNA transposons (15–700+ kb) endemic to Pezizomycotina fungi. They are defined by a captain gene (tyrosine recombinase) at the 5' end, terminal inverted repeats (TIRs) at the boundaries, and target site duplications (TSDs) flanking the insertion. They carry diverse cargo genes including secondary metabolite clusters, drug resistance genes, and virulence factors.

## Input

Two input modes:

1. **GenBank file:** `starkit input.gbk`
2. **FASTA + GFF3:** `starkit input.fasta --gff annotations.gff3`

Both are normalized internally into BioPython SeqRecord objects with gene features. The `--gff` flag is required when the input is FASTA format; an error is raised if FASTA is provided without GFF3.

## Pipeline

### Step 1: Parse & Normalize Input (`files.py`)

- Detect input format (GenBank vs. FASTA)
- If FASTA + GFF3, merge into SeqRecord objects with features
- Validate that protein sequences are available (or can be extracted/translated from CDS features)
- Return a list of SeqRecord objects (one per contig/chromosome)

### Step 2: Detect Captain Genes (`captain.py`)

- Extract all protein sequences from annotated CDS features
- Search against bundled captain HMM profiles using pyhmmer
- HMM profiles are built from curated captain sequences in Starbase, grouped by family
- Filter hits by e-value threshold (default: 1e-10, configurable via `--evalue`)
- Captain genes encode tyrosine recombinases with the conserved RKHRY catalytic pentad
- Return a list of captain hit locations with associated metadata (e-value, HMM score, family match)

### Step 3: Define Starship Boundaries (`boundaries.py`)

This is the most complex step. For each captain gene hit:

**3a. Find the 5' boundary:**
- The captain gene is always at the 5' end of a Starship, typically within 2–7 kb of the boundary
- Scan upstream of the captain (up to 10 kb) for:
  - `TTAC` target site motif (Starship insertion site, 98% conserved)
  - TIR sequences using Starbase-derived position-weight matrices (PWMs)
- If at a contig edge, flag as potentially truncated

**3b. Find the 3' boundary:**
- Scan downstream of the captain up to `--max-size` (default: 700,000 bp)
- Search for:
  - Matching TIR (inverted complement of the 5' TIR)
  - Matching TSD/SDR: copy of the `TTAC(N7)A` motif flanking the 3' boundary
  - Gene density drop-off as a secondary signal
- If at a contig edge, flag as potentially truncated

**3c. De novo TIR detection (fallback):**
- When Starbase-derived PWMs fail to find TIRs (novel families, degraded repeats)
- Take flanking sequence on the left side of the search window
- Reverse-complement flanking sequence on the right side
- Align or use k-mer matching to find the best inverted repeat pair
- Apply minimum length and identity thresholds (to be calibrated from Starbase data)

**3d. Validate boundaries:**
- Check for matching TSD/SDR sequences flanking both boundaries
- Confirm predicted size is within plausible range (`--min-size` to `--max-size`)
- Handle edge cases:
  - Contig edges: flag as potentially truncated, report partial prediction
  - Multiple captains in proximity: resolve as nested or tandem insertions using TIR evidence
  - Degraded/absent TIRs: report captain-only prediction with low confidence

**Empirical data from literature informing this step:**
- TSDs are variable length across families: 4 bp (Defiant), 6 bp (Phoenix), 12 bp (Voyager)
- Consensus target site: `TTAC(N7)A` (from 132 insertion sites, PNAS 2023)
- TIR lengths/sequences are poorly characterized; Starbase data will be used to establish family-specific models
- Starship sizes: 27–393 kb observed (avg ~110 kb), 15 kb minimum cutoff

### Step 4: Extract Cargo Genes (`cargo.py`)

- Collect all annotated genes (CDS features) within each predicted Starship boundary
- Pull existing annotations directly from the input GenBank/GFF3 (no re-annotation)
- Exclude the captain gene from the cargo list
- Return cargo gene lists per Starship with: gene ID, product annotation, start, end, strand

### Step 5: Classify Starship Families (`classify.py`)

- Score each captain protein against per-family HMM profiles built from Starbase
- Assign the best-matching family based on highest HMM score above a threshold
- If no family scores above threshold, assign "unclassified"
- Return family name and score per Starship

### Step 6: Score Confidence (`confidence.py`)

Composite confidence score (0–1) and evidence level per prediction:

| Evidence Level | Criteria |
|----------------|----------|
| High | Captain + both TIRs + TSD detected |
| Medium | Captain + partial boundary evidence (one TIR, or TSD only) |
| Low | Captain gene only, no boundary features detected |

The composite score incorporates:
- Captain HMM e-value (stronger hit = higher score)
- TIR detection quality (both found, one found, none)
- TSD presence and match to `TTAC(N7)A` consensus
- Size plausibility relative to known Starship distribution

### Step 7: Write Output (`write.py` + `report.py`)

Two output files produced:

**TSV summary table** (`<prefix>.tsv`):
- One row per predicted Starship
- Columns: starship_id, contig, start, end, size, captain_gene, captain_evalue, family, family_score, evidence_level, confidence_score, tir_left, tir_right, tsd_sequence, cargo_gene_count, truncated

**HTML report** (`<prefix>.html`):
- Self-contained single file (all CSS/JS inline)
- Header: StarKIT banner, input file info, genome stats, run parameters, summary counts
- Summary table: sortable, one row per Starship
- Expandable detail sections per Starship:
  - Linear gene diagram (inline SVG): captain gene (amber), TIR markers (blue), cargo gene arrows (strand-colored), scale bar
  - Details table: captain HMM e-value, family score, TSD sequence, TIR sequences
  - Cargo gene table: gene ID, product, start, end, strand

## Architecture

Follows ClipKIT pattern: single entry point, internal modules, no exposed subcommands.

### Project Structure

```
StarKIT/
├── starkit/
│   ├── __init__.py          # imports api.starkit
│   ├── __main__.py          # entry point: calls main()
│   ├── starkit.py           # CLI orchestrator: main() → execute() → run()
│   ├── parser.py            # argparse with ASCII banner + grouped help
│   ├── args_processing.py   # argument validation & defaults
│   ├── version.py           # version string
│   ├── settings.py          # constants (default thresholds, data paths)
│   ├── exceptions.py        # StarKITException
│   ├── logger.py            # dual logging (stdout + file)
│   ├── helpers.py           # shared utilities
│   ├── files.py             # GenBank/FASTA+GFF3 parsing, format detection
│   ���── write.py             # output formatting (TSV)
│   ├── captain.py           # captain gene detection (pyhmmer)
│   ├── boundaries.py        # TIR/TSD scanning, boundary definition
│   ├── cargo.py             # cargo gene extraction
│   ├── classify.py          # Starship family classification
│   ├── confidence.py        # confidence scoring + evidence levels
│   ├── report.py            # HTML report generation
│   └── data/                # bundled reference data
│       ├── hmm/             # captain HMM profiles (from Starbase)
│       ├── boundaries/      # TIR PWMs (from Starbase)
│       └── families/        # per-family HMM profiles for classification
├── tests/
│   ├── unit/
│   └── integration/
│       └── samples/         # test GenBank files
├── setup.py
├── Makefile
└── requirements.txt
```

### Entry Point Flow

```
starkit input.gbk [options]
    → starkit.py:main(argv)
        → parser.py:create_parser()
        → args_processing.py:process_args()
        → starkit.py:execute()
            → files.py: parse input → SeqRecords
            → captain.py: detect captains
            → boundaries.py: define boundaries
            → cargo.py: extract cargo
            → classify.py: classify families
            → confidence.py: score predictions
            → write.py: write TSV
            → report.py: write HTML
        → starkit.py:run() returns StarKITRun
```

### Key Data Structures

```python
class EvidenceLevel(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"

@dataclass
class StarshipResult:
    region: SeqRecord          # source contig
    start: int                 # Starship start coordinate
    end: int                   # Starship end coordinate
    captain: Feature           # captain gene feature
    captain_evalue: float      # captain HMM e-value
    captain_family: str        # classified family name
    family_score: float        # HMM score for classification
    tir_left: tuple | None     # (start, end, sequence) or None
    tir_right: tuple | None    # (start, end, sequence) or None
    tsd: str | None            # TSD sequence or None
    cargo_genes: list          # list of gene Features within boundaries
    confidence_score: float    # composite 0–1
    evidence_level: EvidenceLevel
    truncated: bool            # True if at contig edge

@dataclass
class StarKITRun:
    input_file: str
    genome_stats: dict         # contig count, total length, etc.
    starships: list            # list of StarshipResult
    parameters: dict           # run parameters (thresholds, etc.)
    version: str
```

## CLI Interface

```
starkit <input> [options]

Required:
  input                    Genome file (GenBank or FASTA)

Optional:
  --gff, -g               GFF3 annotation file (required if input is FASTA)
  --output, -o             Output prefix (default: <input>.starkit)
  --evalue, -e             Captain HMM e-value threshold (default: 1e-10)
  --min-size               Minimum Starship size in bp (default: 15000)
  --max-size               Maximum Starship size in bp (default: 700000)
  --evidence, -ev          Minimum evidence level to report: high, medium, low, all (default: all)
  --log                    Write log file
  --quiet, -q              Suppress stdout messages
  --version, -v            Show version
```

## HTML Report Design

Self-contained single HTML file. Color scheme matches ClipKIT documentation (jlsteenwyk.com/ClipKIT):

| Element | Color |
|---------|-------|
| Captain genes | `#c77c11` (amber) |
| TIR markers | `#0b81d5` (primary blue) |
| Cargo genes (+ strand) | `#4B77BE` (secondary blue) |
| Cargo genes (- strand) | `#6BAADB` (lighter blue) |
| Headings/labels | `#444` |
| Table headers/panels | `#EEE` background |
| Confidence: high | `#0b81d5` |
| Confidence: medium | `#c77c11` |
| Confidence: low | `#999` |
| Font family | Oxygen, Garamond, Georgia, serif |

Structure:
- Header with banner, input info, genome stats, run parameters, summary counts
- Sortable summary table (one row per Starship)
- Expandable detail sections per Starship with:
  - Inline SVG linear gene diagram (captain, TIRs, cargo arrows, scale bar)
  - Details table (e-values, scores, TSD/TIR sequences)
  - Cargo gene table (ID, product, coordinates, strand)

## Dependencies

```
biopython>=1.82     # sequence/annotation parsing
pyhmmer>=0.10.0     # HMM search (no external HMMER needed)
numpy>=1.24.0       # numerical operations
jinja2>=3.1.0       # HTML report templating
```

## Reference Data

All bundled in `starkit/data/`, built from Starbase (https://github.com/FungAGE/starbase):

- **Captain HMM profiles** (`data/hmm/`): built from curated captain protein sequences in Starbase, used for detection (Step 2)
- **TIR position-weight matrices** (`data/boundaries/`): built from curated Starship boundary sequences in Starbase, used for boundary scanning (Step 3)
- **Family HMM profiles** (`data/families/`): per-family captain HMMs for classification (Step 5)

Data preparation is a separate preprocessing step (not part of the StarKIT tool). Scripts for building these profiles from Starbase will live in a `scripts/` directory and are not installed with the package.

## Benchmarking

Benchmarking is separate from the StarKIT tool. It will:
- Compare StarKIT predictions against curated Starships in Starbase
- Compute precision, recall, F1 at the Starship level
- Evaluate boundary accuracy
- Results stored in a directory added to `.gitignore`

## Out of Scope for v1

- Re-annotation of cargo genes (users provide annotated genomes)
- Phylogeny-based family classification (v1 uses HMM-score similarity)
- Multi-genome comparative analysis
- Functional categorization of cargo genes
- Comparison to other tools (e.g., starfish)
