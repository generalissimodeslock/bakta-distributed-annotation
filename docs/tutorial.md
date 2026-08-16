# Basic installation and execution tutorial

This tutorial covers the essential commands required to install Bakta and run the two scripts in this repository.

## 1. Create the software environment

The repository includes an `environment.yml` file. With Mamba:

```bash
mamba env create -f environment.yml
mamba activate bakta_env
```

With Conda:

```bash
conda env create -f environment.yml
conda activate bakta_env
```

Verify the installation:

```bash
bakta --version
bakta --help
```

## 2. Install a Bakta database

List compatible database versions:

```bash
bakta_db list
```

Create a database directory:

```bash
mkdir -p ~/bakta_db
```

Download the full database:

```bash
bakta_db download --output ~/bakta_db --type full
```

Or the reduced database:

```bash
bakta_db download --output ~/bakta_db --type light
```

The fleet script resolves the database location in this order:

```text
--db
↓
BAKTA_DB environment variable
↓
script default path
```

For example:

```bash
export BAKTA_DB=/path/to/bakta/db
```

## 3. Inspect the command-line interfaces

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

## 4. Organize the input samples

The fleet script expects one immediate subdirectory per sample. Within each sample directory it recursively searches for:

```text
*_genomic.fna
*_genomic.gbff
```

NCBI FASTA files containing only CDS or RNA sequences are excluded.

Example:

```text
etapa00_genomas_ncbi/
├── sample_01/
│   ├── assembly_genomic.fna
│   └── assembly_genomic.gbff
├── sample_02/
└── sample_03/
```

The scripts contain default paths from the environment in which the workflow was developed. On another system, pass explicit `--input-root`, `--output-root`, and `--db` paths.

## 5. Initialize the shared queue

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root /path/to/etapa00_genomas_ncbi \
  --output-root /path/to/output00_bakta \
  init
```

The queue is created under `output00_bakta/bakta_queue/` with `pending`, `running`, `done`, `failed`, and `logs` directories.

## 6. Monitor queue status

One-time status:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /path/to/output00_bakta \
  status
```

Continuous monitoring:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /path/to/output00_bakta \
  status --watch
```

## 7. Test a single job

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /path/to/output00_bakta \
  work \
  --id worker-01 \
  --cpus 4 \
  --db /path/to/bakta/db \
  --once
```

`--once` claims at most one job and then exits, which is useful for validating a new worker.

## 8. Process the queue until completion

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

Start additional workers with distinct identifiers such as `worker-02` and `worker-03`.

All workers must see the same shared input/output filesystem and queue. Each worker must have a functional Bakta installation and access to a compatible database. For read-heavy database workloads, local identical database copies are preferable to a single network-hosted database when practical.

## 9. Important Bakta-stage options

```text
--cpus N
--db PATH
--complete
--compliant
--keep-contig-headers
--skip-plot
--gram +|-|?
--translation-table 11|4|25
--tmp-dir PATH
--once
--until-empty
```

Use `--complete` only when **all sequences in the FASTA are complete replicons**.

## 10. Optional completion e-mail

The public documentation recommends providing notification settings explicitly:

```bash
export BAKTA_EMAIL_TO="recipient@example.org"
export BAKTA_SMTP_USER="sender@example.org"
export BAKTA_GMAIL_APP_PASSWORD="application-password"
```

Then run:

```bash
python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-01 \
  --notify-email \
  --until-empty
```

Do not commit passwords, `.env` files, or other secrets.

## 11. Preview the post-Bakta selection stage

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /path/to/output00_bakta \
  --output-root /path/to/output00_bakta_selecao \
  --dry-run
```

`--dry-run` validates discoverability of samples and required inputs without producing new selection outputs.

## 12. Run the post-Bakta selection stage

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /path/to/output00_bakta \
  --output-root /path/to/output00_bakta_selecao
```

## 13. Process a single sample

Batch mode by sample name:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /path/to/output00_bakta \
  --output-root /path/to/output00_bakta_selecao \
  --sample sample_01
```

Direct mode:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-dir /path/to/output00_bakta/sample_01 \
  --outdir /path/to/output00_bakta_selecao/sample_01
```

## 14. Default selection policy

```text
short protein:            < 90 aa  → descriptive audit flag
weak identity:            < 90%    → low priority when applicable
weak query coverage:      < 80%    → low priority when applicable
weak subject coverage:    < 80%    → low priority when applicable
```

Pseudogenes and sORFs are separated by default. Missing, hypothetical, or uncharacterized products are assigned high priority; generic family/domain/DUF/UPF descriptions are assigned medium priority; uncertain wording or weak inference metrics may lead to low priority.

Audit flags such as short protein length or missing EC/COG/GO identifiers do not, by themselves, mean that an annotation is wrong.

## 15. Main selection outputs

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

## 16. Minimal workflow summary

```text
NCBI genome folders
        ↓
queue initialization
        ↓
distributed Bakta execution
        ↓
Bakta output directories
        ↓
selection dry-run
        ↓
selection/audit
        ↓
priority FASTA files + audit tables + manifest
```

## Conceptual note

Bakta performs genome annotation. The distributed task queue, file-integrity auditing, and downstream protein-priority classification are implemented by the scripts in this repository and are not native Bakta functions.
