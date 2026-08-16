# Bakta Distributed Annotation

Distributed bacterial genome annotation with Bakta and auditable protein prioritization for functional refinement.

## Overview

This repository contains two components external to Bakta, used in a didactic and reproducible bacterial genome annotation workflow:

1. **`anot00_bakta_01_fleet.py`** — organizes samples, creates a shared coordinator–worker queue, and runs Bakta in a distributed manner.
2. **`anot00_bakta_02_selecionados.py`** — reads Bakta outputs, checks consistency across files, and prioritizes proteins for subsequent functional refinement.

> Bakta performs the genome annotation itself. Task distribution, auditing, and priority classification are functions of this external pipeline and are not native Bakta features.

## Repository structure

```text
bakta-distributed-annotation/
├── README.md
├── .gitignore
├── scripts/
│   ├── anot00_bakta_01_fleet.py
│   └── anot00_bakta_02_selecionados.py
└── docs/
    └── tutorial.md
```

## Main requirements

- Linux
- Python 3
- Bakta
- a database compatible with the installed Bakta version
- Conda or Mamba is recommended for environment installation

The second script uses only the Python standard library and requires Python 3.9 or later.

## Quick start

First, check the command-line help for each script:

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

Basic workflow:

```text
NCBI genomes
      ↓
anot00_bakta_01_fleet.py init
      ↓
job queue
      ↓
anot00_bakta_01_fleet.py work
      ↓
Bakta
      ↓
output00_bakta
      ↓
anot00_bakta_02_selecionados.py
      ↓
proteins prioritized for downstream analyses
```

Installation and execution instructions are available in [`docs/tutorial.md`](docs/tutorial.md).

## Expected inputs

The first script treats each immediate subdirectory of the input root as one sample and recursively searches for:

- `*_genomic.fna`
- `*_genomic.gbff`

NCBI FASTA files containing only CDS or RNA sequences are explicitly ignored.

## Main outputs

The selection stage classifies eligible CDS into three priority levels:

- `high`
- `medium`
- `low`

Pseudogenes and sORFs are handled separately by default. The `<sample>.priority_high.faa` file was designed as the input contract for the next curated-search stage.

## Reproducibility and auditing

The workflow includes identity checks between JSON and FAA files, SHA-256 hashes, atomic file writing, and per-sample manifests to improve traceability and reproducibility.

## Bakta

Bakta is developed by Schwengers and collaborators. Please consult the official Bakta project and cite the original publication when using the software in scientific work.

## License

The license for this repository has not yet been defined. Until a `LICENSE` file is added, no reuse license should be assumed for these scripts.
