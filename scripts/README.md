# Scripts

This directory contains the two external workflow components described in the repository documentation.

## `anot00_bakta_01_fleet.py`

Distributed Bakta orchestration. It discovers NCBI genome/GBFF pairs, initializes a shared queue, atomically claims jobs, runs Bakta, records provenance and progress, validates required output files, and moves jobs to `done` or `failed`.

Main commands:

```bash
python3 anot00_bakta_01_fleet.py init
python3 anot00_bakta_01_fleet.py work --id worker-01 --cpus 4 --until-empty
python3 anot00_bakta_01_fleet.py status --watch
```

## `anot00_bakta_02_selecionados.py`

Post-Bakta auditing and protein prioritization. It checks JSON/FAA identifier and sequence consistency, validates the official hypothetical-protein subset, separates sORFs and pseudogenes by default, and writes mutually exclusive `high`, `medium`, and `low` priority FASTA files plus audit tables and a manifest.

Main commands:

```bash
python3 anot00_bakta_02_selecionados.py --dry-run
python3 anot00_bakta_02_selecionados.py
```

Use `--help` on either script for the complete command-line interface.
