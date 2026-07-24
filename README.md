# StarKIT

StarKIT is a beta command-line tool for identifying Starship transposable
elements in a single fungal genome. It combines captain-gene detection,
homology to known Starships, boundary inference, and cargo annotation to
produce TSV, FASTA, BED, and HTML results.

## Local installation for beta testers

StarKIT requires Git and Python 3.9 or newer. We recommend installing it in a
dedicated virtual environment.

```bash
git clone https://github.com/JLSteenwyk/starkit.git
cd starkit

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

On Windows PowerShell, activate the environment with:

```powershell
.\.venv\Scripts\Activate.ps1
```

Confirm that the installation works:

```bash
starkit --version
starkit --help
```

## Quick start

Run StarKIT on an annotated GenBank genome:

```bash
starkit genome.gbk -o genome_results
```

FASTA input requires a matching GFF3 annotation:

```bash
starkit genome.fasta --gff genome.gff3 -o genome_results
```

StarKIT is under active beta development. Please report the command used,
console output, and relevant input details when sharing a problem.
