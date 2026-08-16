# Installation and usage tutorial

This tutorial describes the complete workflow for installing and using the two scripts in the **Bakta Distributed Annotation** repository. It is intended to let a new user reproduce the workflow either on a single computer or across multiple workers connected to a shared queue.

The repository contains two main stages:

1. `anot00_bakta_01_fleet.py` — discovers samples, creates the task queue, runs Bakta, and monitors workers.
2. `anot00_bakta_02_selecionados.py` — audits Bakta outputs and prioritizes proteins for downstream functional refinement.

> Bakta performs genome annotation. The distributed queue, cross-file auditing, and protein-priority classification are external functions implemented by the scripts in this repository.

---

## 1. Requirements

Before starting, the computer should have:

- Linux;
- Git;
- Conda or Mamba;
- internet access to install the environment and download the Bakta database;
- enough disk space for the Bakta database and analysis outputs.

The environment supplied by this repository uses:

```text
Python >= 3.9
Bakta 1.12.x
```

The Bakta database is **not included in the repository** because it is large and independently versioned.

---

## 2. Clone the repository

In a terminal:

```bash
git clone https://github.com/generalissimodeslock/bakta-distributed-annotation.git
cd bakta-distributed-annotation
```

From this point onward, the commands in this tutorial assume that the terminal is inside the repository directory.

The main structure should look like:

```text
bakta-distributed-annotation/
├── README.md
├── LICENSE
├── environment.yml
├── docs/
├── scripts/
│   ├── anot00_bakta_01_fleet.py
│   └── anot00_bakta_02_selecionados.py
└── tests/
    └── test_smoke.py
```

---

## 3. Create the software environment

### With Mamba

```bash
mamba env create -f environment.yml
mamba activate bakta_env
```

### With Conda

```bash
conda env create -f environment.yml
conda activate bakta_env
```

Confirm that the environment is active:

```bash
python --version
bakta --version
```

You can also inspect Bakta help:

```bash
bakta --help
```

---

## 4. Verify the scripts

Check script versions:

```bash
python3 scripts/anot00_bakta_01_fleet.py --version
python3 scripts/anot00_bakta_02_selecionados.py --version
```

Inspect all command-line options:

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

These commands do not run analyses; they only confirm that Python can load the scripts.

---

## 5. Install the Bakta database

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

Or download the reduced database:

```bash
bakta_db download --output ~/bakta_db --type light
```

After the download, identify the actual database directory created by `bakta_db` and use that path in the commands below.

`anot00_bakta_01_fleet.py` resolves the database in this order:

```text
1. --db
2. BAKTA_DB environment variable
3. default path stored in the script
```

To expose the path in the current shell:

```bash
export BAKTA_DB=/actual/path/to/bakta/database
```

On another computer, it is recommended to pass `--db` explicitly or define `BAKTA_DB` on that worker.

---

## 6. Prepare the samples

The fleet script treats each immediate subdirectory of the input root as one sample.

Example:

```text
etapa00_genomas_ncbi/
├── Acinetobacter_schindleri_HZE30-1/
│   ├── GCF_xxxxx_genomic.fna
│   └── GCF_xxxxx_genomic.gbff
├── Delftia_tsuruhatensis_CM13/
│   ├── GCF_yyyyy_genomic.fna
│   └── GCF_yyyyy_genomic.gbff
└── another_sample/
    ├── assembly_genomic.fna
    └── assembly_genomic.gbff
```

Inside each sample directory, the script recursively searches for:

```text
*_genomic.fna
*_genomic.gbff
```

NCBI FASTA files containing only CDS or RNA sequences are excluded, such as:

```text
*_cds_from_genomic.fna
*_rna_from_genomic.fna
```

The `GBFF` file is used to recover metadata, whereas the `FNA` contains the assembly passed to Bakta.

---

## 7. Define working directories

The scripts preserve default paths from the environment in which the pipeline was originally developed. On another system, explicitly providing your own paths is safer.

This tutorial uses:

```text
/path/to/etapa00_genomas_ncbi
/path/to/output00_bakta
/path/to/output00_bakta_selecao
```

Replace `/path/to/` with locations that exist on your system.

You may also define shell convenience variables:

```bash
INPUT_ROOT=/path/to/etapa00_genomas_ncbi
BAKTA_OUTPUT=/path/to/output00_bakta
SELECTION_OUTPUT=/path/to/output00_bakta_selecao
```

These variables are only shell shortcuts; the scripts do not depend on them.

---

# Part I — Running Bakta

## 8. Create the task queue

The queue must be initialized before starting workers.

Using explicit paths:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root /path/to/etapa00_genomas_ncbi \
  --output-root /path/to/output00_bakta \
  init
```

Or, if you defined the convenience variables:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  init
```

The queue is created inside:

```text
output00_bakta/
└── bakta_queue/
    ├── pending/
    ├── running/
    ├── done/
    ├── failed/
    └── logs/
```

Each valid sample generates a small JSON work record, or `job`.

### Rebuilding the queue

The command supports:

```bash
--overwrite
```

For example:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  init --overwrite
```

This clears the queue states before recreating them. **Do not use it routinely**; use it only when you intentionally want to rebuild the queue.

---

## 9. Check queue status

One-time status:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status
```

Continuous monitoring:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status --watch
```

By default, `--watch` refreshes every five seconds.

Use another interval with:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status --watch --interval 10
```

The panel reports how many jobs are in:

```text
pending
running
done
failed
```

It also displays active workers, the sample currently being processed, and an estimated progress value parsed from Bakta messages.

---

## 10. Test a single job

Before releasing the entire queue, run one sample only:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --once
```

The option:

```text
--once
```

makes the worker claim at most one job and exit afterward.

This is the recommended way to validate a new installation before processing all samples.

---

## 11. Process the full queue on one computer

Distributed execution is **not required**. A single computer can process the entire queue:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --until-empty
```

The option:

```text
--until-empty
```

keeps the worker active until no jobs remain in `pending` or `running`.

---

## 12. Run multiple computers

For distributed execution, all workers must be able to access:

- the same input root;
- the same `output00_bakta` root;
- the same `bakta_queue` directory.

Each computer runs the same command with a different worker identifier and, if needed, a different local database path.

### Worker 01

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --until-empty
```

### Worker 02

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-02 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --until-empty
```

Worker identifiers must be distinct. They may be names such as `worker-01`, `worker-02`, or IP addresses.

The Bakta database does not have to live at the same physical path on every machine: each worker can point `--db` to its own compatible local copy.

For read-heavy databases, identical local copies are generally preferable to a single network-hosted database.

---

## 13. Main Bakta-stage options

The most useful `work` options are:

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
--no-dashboard
--dashboard-interval SECONDS
```

### `--complete`

Use only when **all sequences in the FASTA are complete replicons**.

Do not use it indiscriminately for fragmented WGS assemblies.

### `--skip-plot`

Disables Bakta PNG/SVG plot generation. This is useful when annotation files are the main goal and unnecessary graphical outputs should be avoided.

### `--gram`

Accepted values:

```text
+
-
?
```

The default is `?`.

### `--translation-table`

Accepted values:

```text
11
4
25
```

The default is `11`.

### `--tmp-dir`

Allows Bakta temporary files to be directed to a fast local disk:

```bash
--tmp-dir /path/to/local/tmp
```

---

## 14. Worker dashboard

During execution, each worker has its own terminal dashboard.

Disable it with:

```bash
--no-dashboard
```

Change its refresh interval with:

```bash
--dashboard-interval 2
```

The value is expressed in seconds.

---

## 15. Optional completion e-mail

The public repository contains no personal e-mail address or password embedded in the code.

To enable notifications, define the settings explicitly.

Example using the default environment-variable names:

```bash
export BAKTA_EMAIL_TO="recipient@example.org"
export BAKTA_SMTP_USER="sender@example.org"
export BAKTA_GMAIL_APP_PASSWORD="application-password"
```

Then choose **one** worker to be responsible for notification:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --notify-email \
  --until-empty
```

The recipient and SMTP user can also be provided directly with:

```text
--email-to
--smtp-user
--smtp-host
--smtp-port
--smtp-password-env
```

The password should remain in an environment variable and **must never be committed to GitHub**.

The repository `.gitignore` already excludes `.env` files and secret directories, but this does not replace careful handling by the user.

---

## 16. Verify annotation completion

When a job finishes, it is moved to:

```text
done/
```

or:

```text
failed/
```

Logs are stored under:

```text
output00_bakta/bakta_queue/logs/
```

For a final check:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status
```

A job is only considered successful when Bakta exits appropriately and the required output files are present.

---

# Part II — Protein audit and selection

## 17. Inputs used by the second stage

`anot00_bakta_02_selecionados.py` primarily uses:

```text
<sample>.json
<sample>.faa
```

When available, it also uses:

```text
<sample>.inference.tsv
<sample>.hypotheticals.faa
<sample>.hypotheticals.tsv
```

The JSON file provides annotation metadata and the FAA file is treated as the contractual source of protein sequences.

Before classifying proteins, the script checks identifier and sequence consistency across the relevant files.

---

## 18. Run `--dry-run` first

Before creating selection outputs:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --dry-run
```

`--dry-run`:

- discovers samples;
- identifies the primary JSON and FAA files;
- reports optional auxiliary files;
- shows the output directories that would be used;
- does not create final selection outputs.

This should be the first command used in the second stage.

---

## 19. Run selection for all samples

After a successful `--dry-run`:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT"
```

For each sample, the program reports counts assigned to:

```text
high
medium
low
sorf
pseudogene
```

At the end, it also reports how many samples completed, were skipped, or failed.

---

## 20. Process one or several samples

### One sample

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --sample sample_01
```

### Multiple samples

Repeat `--sample`:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --sample sample_01 \
  --sample sample_02
```

### Sample list file

A file containing one sample name per line can also be used:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --samples-file samples.txt
```

Lines beginning with `#` may be used as comments in the list.

---

## 21. Direct mode for a single Bakta directory

To analyze one specific Bakta output directory directly:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-dir /path/to/output00_bakta/sample_01 \
  --outdir /path/to/output00_bakta_selecao/sample_01
```

If the Bakta file prefix cannot be inferred automatically, provide:

```bash
--prefix PREFIX
```

---

## 22. Default classification policy

Default values are:

```text
short protein:            < 90 amino acids
weak identity:            < 90%
weak query coverage:      < 80%
weak subject coverage:    < 80%
```

These thresholds can be changed with:

```text
--short-aa-threshold
--weak-identity
--weak-query-cov
--weak-subject-cov
```

The final three options receive values between `0` and `1`.

For example, `0.90` represents 90%.

### `high` priority

Generally includes proteins with missing, hypothetical, or uncharacterized descriptions.

### `medium` priority

Generally includes generic family-, domain-, DUF-, or UPF-level descriptions.

### `low` priority

Generally includes uncertain wording or informative annotations supported by identity or coverage below the configured thresholds.

### Not selected

Proteins with informative products and evidence considered adequate remain in the audit tables but are not written to priority FASTA files.

---

## 23. Audit flags are not automatically errors

The script may record flags such as:

```text
short_protein
gene_missing
ec_missing
cog_missing
go_missing
inference_missing
edge_feature
truncated_...
```

These flags describe properties that may be useful during review.

They **do not mean by themselves that the annotation is wrong**, and they are not necessarily selection criteria.

---

## 24. sORFs and pseudogenes

By default:

- sORFs are preserved separately;
- pseudogenes are preserved separately;
- neither group is automatically mixed into conventional protein targets.

This behavior can be changed with:

```text
--include-sorf-targets
--include-pseudogenes
```

Use these options deliberately because they change the conventional classification set.

---

## 25. Options that require prior review

The selector is designed to stop when certain integrity inconsistencies are detected.

The following options exist to continue **only after the underlying issue has been consciously reviewed**:

```text
--allow-id-mismatch
--allow-sequence-mismatch
--allow-hypothetical-mismatch
```

They mean, respectively:

- accept differences between JSON and FAA locus sets;
- accept sequence differences between JSON and FAA;
- accept inconsistencies involving the official `hypotheticals.faa` subset.

Do not use these options merely to suppress an error message.

Determine the cause of the inconsistency first.

---

## 26. Using `--force`

Both scripts expose a `--force` option in different contexts.

### In the fleet script

```text
work --force
```

passes `--force` to Bakta and allows an existing output directory to be overwritten.

### In the selector

```text
--force
```

allows known selection-stage output files to be replaced.

Under normal conditions, the selector supports safe reruns: valid outputs can be recognized from the manifest and input fingerprints and are skipped instead of being regenerated.

---

## 27. Stop the batch on the first error

By default, the selector can continue to other samples after one sample fails.

To stop immediately:

```bash
--fail-fast
```

This is especially useful during validation or development.

---

## 28. Selection outputs

For each sample, the second stage can produce:

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

### `targets_refined.tsv`

Contains selected proteins, original Bakta information, inference metrics, flags, and the priority rationale.

### `selection_audit.tsv`

Records the decision taken for audited proteins, including proteins that were not selected or were separated.

### Priority FASTA files

```text
priority_high.faa
priority_medium.faa
priority_low.faa
```

contain sequences intended for downstream searches according to the assigned priority.

### `target_report.tsv`

Summarizes general counts and frequencies of important flags.

### `target_manifest.json`

Records parameters, versions, input files, SHA-256 hashes, counts, warnings, and output fingerprints.

It is the main traceability record for the selection stage.

---

## 29. Run repository tests

The basic tests can be executed locally without running a complete Bakta annotation.

First check Python syntax:

```bash
python -m py_compile \
  scripts/anot00_bakta_01_fleet.py \
  scripts/anot00_bakta_02_selecionados.py
```

Then run the tests:

```bash
python -m unittest discover -s tests -v
```

The tests verify, among other points:

- recognition of the primary genomic FASTA;
- rejection of CDS- and RNA-only FASTA files;
- sample-name sanitization;
- accession normalization;
- representative `high`, `medium`, and `low` classifications;
- non-selection of an informative annotation with strong evidence.

GitHub Actions automatically runs these checks under Python 3.9, 3.10, 3.11, 3.12, and 3.13 on repository updates.

---

## 30. Recommended minimal sequence

For users who only want the normal workflow, the sequence is:

```bash
# 1. Enter the repository
cd bakta-distributed-annotation

# 2. Activate the environment
mamba activate bakta_env

# 3. Initialize the queue
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  init

# 4. Test one job
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work --id worker-01 --cpus 4 --db "$BAKTA_DB" --skip-plot --once

# 5. Process the remaining queue
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work --id worker-01 --cpus 4 --db "$BAKTA_DB" --skip-plot --until-empty

# 6. Preview selection without writing results
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --dry-run

# 7. Run audit and selection
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT"
```

---

## 31. Conceptual workflow

```text
NCBI genome folders
        ↓
queue initialization
        ↓
pending → running → done/failed
        ↓
Bakta execution by one or more workers
        ↓
output00_bakta
        ↓
selection dry-run
        ↓
JSON × FAA × auxiliary-file auditing
        ↓
high / medium / low classification
        ↓
FASTAs + tables + report + manifest
        ↓
downstream functional refinement
```

---

## 32. Repository scope

This repository ends with Bakta annotation and the immediate post-Bakta selection/audit stage.

The file:

```text
<sample>.priority_high.faa
```

was designed as the input contract for the next curated stage in the broader pipeline, but later stages are outside the scope of this repository.

For a high-level description of the project, see the repository [`README.md`](../README.md).
