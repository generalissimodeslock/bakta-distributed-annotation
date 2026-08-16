# Bakta Distributed Annotation

[![Python checks](https://github.com/generalissimodeslock/bakta-distributed-annotation/actions/workflows/python-checks.yml/badge.svg)](https://github.com/generalissimodeslock/bakta-distributed-annotation/actions/workflows/python-checks.yml)
[![License: GPL-3.0-only](https://img.shields.io/badge/License-GPL--3.0--only-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)

Distributed bacterial genome annotation with **Bakta**, followed by an auditable protein-prioritization step for downstream functional refinement.

**Documentation:** English · [Português (Brasil)](docs/tutorial.pt-BR.md)

## Overview

This repository contains two components external to Bakta that were designed for a reproducible bacterial genome-annotation workflow:

1. **`anot00_bakta_01_fleet.py`** — discovers valid NCBI genome inputs, creates a shared coordinator–worker queue, runs Bakta across multiple workers, records provenance, and monitors execution.
2. **`anot00_bakta_02_selecionados.py`** — audits Bakta JSON/FAA outputs, checks sequence and identifier consistency, separates special protein classes, and assigns eligible CDS to `high`, `medium`, or `low` priority for subsequent functional refinement.

> **Important:** Bakta performs the genome annotation. Queue management, distributed execution, cross-file auditing, and protein-priority classification are implemented by this external workflow and are not native Bakta features.

## Repository structure

```text
bakta-distributed-annotation/
├── README.md
├── LICENSE
├── environment.yml
├── .gitignore
├── .github/
│   └── workflows/
│       └── python-checks.yml
├── docs/
│   ├── tutorial.md
│   └── tutorial.pt-BR.md
├── scripts/
│   ├── README.md
│   ├── anot00_bakta_01_fleet.py
│   └── anot00_bakta_02_selecionados.py
└── tests/
    └── test_smoke.py
```

## Requirements

- Linux
- Python 3.9 or later
- Bakta
- a Bakta database compatible with the installed Bakta version
- Conda or Mamba recommended for environment management

The selector script uses only the Python standard library. The reference implementation targets Bakta 1.12.x outputs.

## Reproducible environment

The easiest setup is to use the supplied environment file:

```bash
mamba env create -f environment.yml
mamba activate bakta_env
```

Conda can be used instead of Mamba:

```bash
conda env create -f environment.yml
conda activate bakta_env
```

The Bakta database is installed separately because it is large and versioned independently from the Conda environment.

## Quick start

Check the command-line interfaces first:

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

A minimal workflow is:

```text
NCBI genome folders
        ↓
anot00_bakta_01_fleet.py init
        ↓
shared job queue
        ↓
anot00_bakta_01_fleet.py work
        ↓
Bakta annotation
        ↓
output00_bakta/
        ↓
anot00_bakta_02_selecionados.py --dry-run
        ↓
anot00_bakta_02_selecionados.py
        ↓
output00_bakta_selecao/
        ↓
prioritized proteins for downstream analysis
```

Full installation and execution instructions are available in [`docs/tutorial.md`](docs/tutorial.md).

## Input contract

The fleet script treats each immediate subdirectory of the input root as one sample and recursively searches for one primary genomic FASTA and its corresponding GenBank flat file:

```text
*_genomic.fna
*_genomic.gbff
```

NCBI FASTA files containing only CDS or RNA sequences are explicitly excluded. When multiple assembly candidates are found, the script applies deterministic selection and consistency checks instead of silently choosing unrelated files.

## Distributed execution

Initialize the queue:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root /path/to/etapa00_genomas_ncbi \
  --output-root /path/to/output00_bakta \
  init
```

Start a worker:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /path/to/output00_bakta \
  work \
  --id worker-01 \
  --cpus 4 \
  --db /path/to/bakta/db \
  --skip-plot \
  --until-empty
```

Monitor the queue:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /path/to/output00_bakta \
  status --watch
```

All workers must see the same shared queue and input/output filesystem. Each worker should use its own functional Bakta installation and preferably a local identical copy of the Bakta database.

## Protein prioritization

Preview the samples without writing outputs:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /path/to/output00_bakta \
  --output-root /path/to/output00_bakta_selecao \
  --dry-run
```

Run the selection stage:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /path/to/output00_bakta \
  --output-root /path/to/output00_bakta_selecao
```

The default policy separates pseudogenes and sORFs from the conventional search set and classifies eligible CDS into three mutually exclusive priority levels:

- **high** — missing, hypothetical, or uncharacterized functional description;
- **medium** — generic, family-, domain-, DUF-, or UPF-level description;
- **low** — uncertain wording or informative annotations supported by weak identity/coverage metrics.

Flags such as short protein length or missing EC/COG/GO identifiers are audit markers; they are not automatically equivalent to an annotation error.

## Main outputs

For each sample, the selection stage can produce:

```text
<sample>.targets_refined.tsv
<sample>.selection_audit.tsv
<sample>.priority_high.faa
<sample>.priority_medium.faa
<sample>.priority_low.faa
<sample>.sorf_separate.faa
<sample>.pseudogene_separate.faa
<sample>.target_report.tsv
<sample>.target_manifest.json
```

`<sample>.priority_high.faa` is the intended contract for the next curated functional-search stage.

## Reproducibility and auditing

The workflow includes several safeguards intended to make processing traceable and reviewable:

- atomic job claiming in the shared queue;
- explicit `pending`, `running`, `done`, and `failed` states;
- Bakta version, database path, worker identity, CPU count, parameters, and timestamps recorded in job metadata;
- verification of required Bakta output files;
- JSON/FAA locus-set consistency checks;
- JSON/FAA amino-acid sequence identity checks;
- validation of the official `hypotheticals.faa` subset;
- mutually exclusive priority FASTA files;
- atomic output writing;
- SHA-256 fingerprints and per-sample manifests;
- safe reruns that skip outputs only when the saved manifest and fingerprints remain valid.

## Automated checks

GitHub Actions runs the repository checks against Python 3.9, 3.10, 3.11, 3.12, and 3.13. The workflow verifies that both scripts compile, that their command-line help can be opened, and that the smoke tests pass.

The smoke tests cover core behavior that does not require downloading the Bakta database, including genomic-FASTA recognition, sample-name sanitization, accession normalization, and representative `high`/`medium`/`low` selection rules.

Run them locally with:

```bash
python3 -m unittest discover -s tests -v
```

## Optional completion e-mail

Notification settings can be supplied through command-line options. For portable deployments, configure recipient/sender information explicitly rather than relying on environment-specific defaults. Passwords are read from an environment variable rather than being written into the command line or repository.

```bash
export BAKTA_GMAIL_APP_PASSWORD="application-password"

python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-01 \
  --email-to recipient@example.org \
  --smtp-user sender@example.org \
  --notify-email \
  --until-empty
```

Secrets should never be committed to the repository.

## Scope

This repository intentionally covers the **Bakta annotation stage and the immediate post-Bakta protein selection/audit stage**. Later functional-refinement stages, such as Swiss-Prot searches, belong to the broader annotation pipeline and are not included here.

## Bakta

Bakta is developed by Schwengers and collaborators. This repository does not redistribute the Bakta database or Bakta source code. When using Bakta in scientific work, consult the official project and cite the appropriate Bakta publication.

Official project: https://github.com/oschwengers/bakta

## License

The repository-specific scripts and documentation are released under the **GNU General Public License v3.0 only (`GPL-3.0-only`)**. See [`LICENSE`](LICENSE).
