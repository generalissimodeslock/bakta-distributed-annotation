#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
"""Etapa 01 do pipeline Bakta: anota genomas bacterianos obtidos do NCBI.

Cada subdiretório imediato de ``etapa00_genomas_ncbi`` é tratado como uma
amostra. Dentro dele, o script procura recursivamente um único
``*_genomic.fna`` e seu ``*_genomic.gbff`` correspondente. FASTA de CDS e RNA
do NCBI são ignorados explicitamente.

A execução é distribuída por uma fila coordinator–worker compartilhada. Cada
resultado é gravado diretamente em ``output00_bakta/<nome_da_pasta>/`` e usa
o mesmo nome seguro como prefixo. O caminho do banco obrigatório pode ser
informado por ``--db``, pela variável ``BAKTA_DB`` ou pelo caminho padrão.
"""

import argparse
import datetime as _dt
import glob
import hashlib
import json
import os
import queue as _queue
import random
import re
import shlex
import shutil
import smtplib
import socket
import ssl
import subprocess
import sys
import threading
import time
from email.message import EmailMessage
from typing import Optional, Tuple

PIPELINE_NAME = "anot00_bakta"
PIPELINE_STAGE = "01_fleet"
SCRIPT_NAME = "anot00_bakta_01_fleet.py"
SCRIPT_VERSION = "1.1.4"
INPUT_SOURCE = "NCBI"

DEFAULT_INPUT_ROOT = "/home/bioinfo/anotacao_coordinator_worker/etapa00_genomas_ncbi"
DEFAULT_OUTPUT_ROOT = "/home/bioinfo/anotacao_coordinator_worker/output00_bakta"
DEFAULT_DB_PATH = "/home/bioinfo/coordinator_worker_databank/bakta/db"
QUEUE_FOLDER_NAME = "bakta_queue"

DEFAULT_EMAIL_TO = os.environ.get("BAKTA_EMAIL_TO", "")
DEFAULT_SMTP_HOST = "smtp.gmail.com"
DEFAULT_SMTP_PORT = 587
DEFAULT_SMTP_USER = os.environ.get("BAKTA_SMTP_USER", "")
DEFAULT_SMTP_PASSWORD_ENV = "BAKTA_GMAIL_APP_PASSWORD"
EMAIL_SENT_MARKER = "notification_email.sent.json"
EMAIL_LOCK_MARKER = "notification_email.lock"
EMAIL_FAILURE_MARKER = "notification_email.failed.json"
EMAIL_RETRY_SECONDS = 300

ALLOWED_RE = re.compile(r'[^A-Za-z0-9._-]')
ALNUM_RE = re.compile(r'[^A-Za-z0-9]')
IPV4_RE = re.compile(r"^\d{1,3}(\.\d{1,3}){3}$")
ASSEMBLY_ACCESSION_RE = re.compile(r"(?<![A-Za-z0-9])(GC[AF]_\d+\.\d+)(?!\d)", re.I)

EXCLUDED_GENOMIC_FASTA_SUFFIXES = (
    "_cds_from_genomic.fna",
    "_rna_from_genomic.fna",
)
EXCLUDED_GENOMIC_FASTA_BASENAMES = {
    "cds_from_genomic.fna",
    "rna_from_genomic.fna",
}
IGNORED_SAMPLE_FOLDER_NAMES = {
    "genomas_referencia_ncbi",
    QUEUE_FOLDER_NAME,
}

def now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")

def sanitize(s: str) -> str:
    return ALLOWED_RE.sub('', s)

def only_alnum(s: str) -> str:
    return ALNUM_RE.sub('', s)

def get_local_ip() -> str:
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("1.1.1.1", 80))
            return s.getsockname()[0]
    except Exception:
        try:
            return socket.gethostbyname(socket.gethostname())
        except Exception:
            return "0.0.0.0"

def resolve_db_path(cli_db: Optional[str]) -> str:
    """Resolve o banco: --db > BAKTA_DB > caminho padrão do pipeline."""
    raw = cli_db or os.environ.get("BAKTA_DB") or DEFAULT_DB_PATH
    return os.path.abspath(os.path.expanduser(raw))

def validate_runtime(db_path: str) -> Optional[str]:
    """Retorna uma mensagem de erro se Bakta ou o banco não estiverem prontos."""
    if shutil.which("bakta") is None:
        return "executável 'bakta' não encontrado no PATH"
    if not os.path.isdir(db_path):
        return f"diretório do banco Bakta não encontrado: {db_path}"
    try:
        if not any(os.scandir(db_path)):
            return f"diretório do banco Bakta está vazio: {db_path}"
    except OSError as exc:
        return f"não foi possível acessar o banco Bakta {db_path}: {exc}"
    return None

def shorten_strain(s: str) -> str:
    if "-" in s and len(s) >= 20:
        s = s.split("-", 1)[0]
    if len(s) > 30:
        s = s[:30]
    return s

def normalize_accession(acc: str) -> str:
    if not acc:
        return ""
    acc = acc.strip()
    acc = acc.split(".", 1)[0]
    if acc.startswith("NZ_"):
        acc = acc[3:]
    return acc

def extract_assembly_accession(*values: str) -> str:
    """Obtém GCF_/GCA_ com versão a partir do nome da pasta ou do arquivo."""
    for value in values:
        if not value:
            continue
        match = ASSEMBLY_ACCESSION_RE.search(value)
        if match:
            return match.group(1).upper()
    return ""

def assembly_accession_version(path: str) -> Optional[Tuple[str, int]]:
    """Retorna (accession-base, versão) para GCF_/GCA_ presente no caminho."""
    accession = extract_assembly_accession(os.path.basename(path))
    if not accession or "." not in accession:
        return None
    base, version_text = accession.rsplit(".", 1)
    try:
        return base, int(version_text)
    except ValueError:
        return None

def safe_sample_name(folder_name: str) -> str:
    """Converte o nome da pasta em um prefixo aceito pelo Bakta."""
    sample = sanitize(folder_name.replace(" ", "_"))
    sample = sample.strip("._-")
    if not sample:
        raise ValueError(f"Nome de pasta inválido para amostra: {folder_name!r}")
    return sample[:160]

def is_primary_genomic_fna(path: str) -> bool:
    """Aceita FASTA genômico e exclui explicitamente FASTA de CDS/RNA."""
    name = os.path.basename(path).lower()
    if name in EXCLUDED_GENOMIC_FASTA_BASENAMES:
        return False
    return name.endswith("_genomic.fna") and not name.endswith(
        EXCLUDED_GENOMIC_FASTA_SUFFIXES
    )

def file_sha256(path: str, chunk_size: int = 1024 * 1024) -> str:
    """Calcula SHA-256 para confirmar se duas cópias FASTA são idênticas."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()

def choose_duplicate_copy(paths: list, sample_dir: str) -> str:
    """Escolhe uma entre cópias byte a byte idênticas; rejeita divergências."""
    hashes = {file_sha256(path) for path in paths}
    if len(hashes) != 1:
        names = ", ".join(os.path.relpath(path, sample_dir) for path in paths)
        raise ValueError(
            "mais de um arquivo para a mesma versão, mas com conteúdos "
            f"diferentes: {names}"
        )
    return min(
        paths,
        key=lambda path: (
            len(os.path.relpath(path, sample_dir).split(os.sep)),
            len(os.path.relpath(path, sample_dir)),
            os.path.relpath(path, sample_dir),
        ),
    )

def read_first_match(path: str, patterns) -> str:
    with open(path, "rt", errors="replace") as fh:
        for line in fh:
            for rx, group in patterns:
                m = rx.search(line)
                if m:
                    return m.group(group).strip()
    return ""

def parse_gbff(gbff_path: str) -> Tuple[str, str, str, str]:
    """
    Retorna: (genus, species, accession_sem_versao, strain_ou_isolate)
    Usa apenas o PRIMEIRO ACCESSION encontrado no arquivo.
    """
    acc = read_first_match(gbff_path, [(re.compile(r"^ACCESSION\s+(\S+)"), 1)])
    acc = acc.split(".", 1)[0] if acc else ""

    organism = read_first_match(
        gbff_path,
        [
            (re.compile(r'/organism="([^"]+)"'), 1),
            (re.compile(r"^  ORGANISM\s+(.+)$"), 1),
        ],
    )

    genus = ""
    species = ""
    if organism:
        parts = organism.split()
        if len(parts) >= 1:
            genus = parts[0]
        if len(parts) >= 2:
            species = parts[1]

    strain = read_first_match(
        gbff_path,
        [
            (re.compile(r'/strain="([^"]+)"'), 1),
            (re.compile(r'/isolate="([^"]+)"'), 1),
        ],
    )

    return genus, species, acc, strain

def build_prefix(genus: str, species: str, acc: str, fallback_base: str) -> str:
    """
    prefix:
      gênero completo + "_" + espécie (minúscula) + "_" + accession
    com a regra de remover NZ_ do accession.

    Ex.:
      Acidimicrobiales + bacterium + CP064997 -> Acidimicrobiales_bacterium_CP064997
    """
    genus = sanitize(genus)
    species = sanitize(species.lower())
    acc = normalize_accession(acc)
    acc = sanitize(acc)

    if genus and species and acc:
        prefix = f"{genus}_{species}_{acc}"
    elif genus and species:
        prefix = f"{genus}_{species}"
    elif acc:
        prefix = f"{sanitize(fallback_base)}_{acc}"
    else:
        prefix = sanitize(fallback_base)

    if len(prefix) > 160:
        prefix = prefix[:160]
    return prefix

def build_locustag(genus: str, species: str, acc: str, base: str) -> str:
    """
    locustag (até 12 caracteres, alfanumérico, começa com letra):
    1 letra do gênero + 5 letras da espécie + 6 primeiros chars do accession (sem NZ_)
    Se faltar algo, completa com hash do base.
    """
    g = (genus[:1] or "X")
    g = only_alnum(g).upper()
    if not g:
        g = "X"

    sp = only_alnum((species or "").lower())
    sp5 = (sp[:5] if sp else "sp")
    sp5 = only_alnum(sp5).upper()

    accn = only_alnum(normalize_accession(acc)).upper()
    acc6 = accn[:6] if accn else ""

    core = (g + sp5 + acc6)
    core = only_alnum(core).upper()

    if not core or not core[0].isalpha():
        core = "X" + core

    if len(core) < 3:
        h = hashlib.md5(base.encode("utf-8")).hexdigest().upper()
        core = (core + h)[:3]

    if len(core) < 12:
        h = hashlib.md5(base.encode("utf-8")).hexdigest().upper()
        need = 12 - len(core)
        core = core + only_alnum(h)[:need]

    if len(core) > 12:
        core = core[:12]

    return core

def find_sample_inputs(sample_dir: str) -> Tuple[str, str, str]:
    """Localiza o genoma e prefere GCF_/GCA_ oficial e a versão mais nova."""
    candidates = sorted(
        path
        for path in glob.glob(
            os.path.join(sample_dir, "**", "*_genomic.fna"), recursive=True
        )
        if os.path.isfile(path) and is_primary_genomic_fna(path)
    )
    if not candidates:
        raise ValueError("nenhum arquivo *_genomic.fna principal encontrado")

    selection_note = ""
    parsed = [(path, assembly_accession_version(path)) for path in candidates]
    recognized = [(path, value) for path, value in parsed if value is not None]

    if recognized:
        accession_bases = {value[0] for _path, value in recognized}
        if len(accession_bases) != 1:
            names = ", ".join(
                os.path.relpath(path, sample_dir) for path, _value in recognized
            )
            raise ValueError(
                "mais de uma montagem oficial GCF_/GCA_ diferente encontrada: "
                + names
            )

        highest_version = max(value[1] for _path, value in recognized)
        newest = [
            path
            for path, value in recognized
            if value[1] == highest_version
        ]
        fna = newest[0] if len(newest) == 1 else choose_duplicate_copy(
            newest, sample_dir
        )

        ignored = [
            os.path.relpath(path, sample_dir)
            for path in candidates
            if path != fna
        ]
        if ignored:
            selection_note = (
                f"selecionado {os.path.relpath(fna, sample_dir)}; "
                "cópia(s), alias ou versão(ões) anterior(es) ignorada(s): "
                + ", ".join(ignored)
            )
    elif len(candidates) == 1:
        fna = candidates[0]
    else:
        fna = choose_duplicate_copy(candidates, sample_dir)
        ignored = [
            os.path.relpath(path, sample_dir)
            for path in candidates
            if path != fna
        ]
        selection_note = (
            f"selecionado {os.path.relpath(fna, sample_dir)}; "
            "cópia(s) idêntica(s) ignorada(s): " + ", ".join(ignored)
        )

    gbff = os.path.splitext(fna)[0] + ".gbff"
    if not os.path.isfile(gbff):
        raise ValueError(
            "GBFF correspondente não encontrado: "
            + os.path.relpath(gbff, sample_dir)
        )
    return os.path.abspath(fna), os.path.abspath(gbff), selection_note

def init_queue(input_root: str, output_root: str, overwrite: bool) -> None:
    input_root = os.path.abspath(os.path.expanduser(input_root))
    output_root = os.path.abspath(os.path.expanduser(output_root))
    queue_dir = os.path.join(output_root, QUEUE_FOLDER_NAME)
    pending = os.path.join(queue_dir, "pending")
    running = os.path.join(queue_dir, "running")
    done = os.path.join(queue_dir, "done")
    failed = os.path.join(queue_dir, "failed")
    logs = os.path.join(queue_dir, "logs")

    for d in (pending, running, done, failed, logs):
        os.makedirs(d, exist_ok=True)

    if overwrite:
        for folder in (pending, running, done, failed):
            for f in glob.glob(os.path.join(folder, "*.json")):
                os.remove(f)

    if not os.path.isdir(input_root):
        print(f"[init] Diretório de entrada não encontrado: {input_root}", file=sys.stderr)
        return

    sample_dirs = sorted(
        path
        for path in glob.glob(os.path.join(input_root, "*"))
        if os.path.isdir(path) and not os.path.basename(path).startswith(".")
    )
    if not sample_dirs:
        print(f"[init] Nenhuma pasta de amostra encontrada em {input_root}", file=sys.stderr)
        return

    created = 0
    skipped = 0
    invalid = 0
    ignored_folders = 0
    for sample_dir in sample_dirs:
        folder_name = os.path.basename(sample_dir)
        if folder_name in IGNORED_SAMPLE_FOLDER_NAMES:
            ignored_folders += 1
            print(
                f"[init] PULADA {folder_name}: pasta agregadora/reservada, "
                "não é uma amostra individual"
            )
            continue
        try:
            sample = safe_sample_name(folder_name)
            if sample == QUEUE_FOLDER_NAME:
                raise ValueError(f"nome reservado: {QUEUE_FOLDER_NAME}")
            fna, gbff, selection_note = find_sample_inputs(sample_dir)
        except ValueError as exc:
            invalid += 1
            print(f"[init] IGNORADA {folder_name}: {exc}", file=sys.stderr)
            continue

        if selection_note:
            print(f"[init] {folder_name}: {selection_note}")

        job_name = f"{sample}.json"
        known_jobs = [
            os.path.join(folder, job_name)
            for folder in (pending, running, done, failed)
        ]
        if any(os.path.exists(path) for path in known_jobs):
            skipped += 1
            continue

        assembly_accession = extract_assembly_accession(
            folder_name, os.path.basename(fna), os.path.basename(gbff)
        )
        job_path = os.path.join(pending, job_name)
        job = {
            "pipeline": PIPELINE_NAME,
            "stage": PIPELINE_STAGE,
            "script": SCRIPT_NAME,
            "script_version": SCRIPT_VERSION,
            "input_source": INPUT_SOURCE,
            "sample": sample,
            "source_folder_name": folder_name,
            "sample_dir": os.path.abspath(sample_dir),
            "fna": fna,
            "gbff": gbff,
            "input_selection_note": selection_note,
            "assembly_accession": assembly_accession,
            "output_dir": os.path.join(output_root, folder_name),
            "created": now_iso(),
        }
        with open(job_path, "wt") as fh:
            json.dump(job, fh, indent=2, sort_keys=True)
        created += 1

    if created:
        for marker_name in (EMAIL_SENT_MARKER, EMAIL_FAILURE_MARKER):
            marker_path = os.path.join(queue_dir, marker_name)
            try:
                os.remove(marker_path)
            except FileNotFoundError:
                pass

    print(
        f"[init] pastas={len(sample_dirs)} criados={created} "
        f"existentes={skipped} invalidos={invalid} "
        f"pastas_puladas={ignored_folders} pending={pending}"
    )

def claim_one_job(pending_dir: str, running_dir: str) -> Optional[str]:
    jobs = glob.glob(os.path.join(pending_dir, "*.json"))
    if not jobs:
        return None
    random.shuffle(jobs)
    for job in jobs:
        dst = os.path.join(running_dir, os.path.basename(job))
        try:
            os.replace(job, dst)
            return dst
        except (FileNotFoundError, OSError):
            continue
    return None

def queue_counts(queue_dir: str) -> dict:
    """Conta os estados atuais da fila compartilhada."""
    return {
        state: len(glob.glob(os.path.join(queue_dir, state, "*.json")))
        for state in ("pending", "running", "done", "failed")
    }

def queue_execution_summary(queue_dir: str) -> dict:
    """Resume início da fila e quantidade de jobs processados por IP."""
    earliest_start = ""
    per_ip = {}

    for state in ("pending", "running", "done", "failed"):
        state_dir = os.path.join(queue_dir, state)
        for path in glob.glob(os.path.join(state_dir, "*.json")):
            job = read_json_safe(path)
            started = str(job.get("start") or "").strip()
            if started and (not earliest_start or started < earliest_start):
                earliest_start = started

            if state not in ("done", "failed"):
                continue

            owner_ip = str(job.get("owner_ip") or "desconhecido").strip()
            if not owner_ip:
                owner_ip = "desconhecido"
            row = per_ip.setdefault(
                owner_ip,
                {"total": 0, "done": 0, "failed": 0},
            )
            row["total"] += 1
            row[state] += 1

    ordered_per_ip = {
        ip: per_ip[ip]
        for ip in sorted(
            per_ip,
            key=lambda value: (-per_ip[value]["total"], value),
        )
    }
    return {
        "started_at": earliest_start,
        "per_ip": ordered_per_ip,
    }


def format_ip_work_summary(per_ip: dict) -> str:
    """Formata o resumo por IP para o corpo textual do e-mail."""
    if not per_ip:
        return "  - nenhum IP registrado\n"

    lines = []
    for ip, values in per_ip.items():
        lines.append(
            f"  - {ip}: {values.get('total', 0)} job(s) "
            f"({values.get('done', 0)} concluído(s), "
            f"{values.get('failed', 0)} falha(s))"
        )
    return "\n".join(lines) + "\n"


def write_json_atomic(path: str, payload: dict) -> None:
    """Grava um JSON por substituição atômica no mesmo sistema de arquivos."""
    temporary = f"{path}.tmp.{os.getpid()}"
    try:
        with open(temporary, "wt", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


PROGRESS_STAGE_RULES = (
    # Regras alinhadas às mensagens reais emitidas pelo Bakta 1.12.x.
    (re.compile(r"annotation successfully finished", re.I), 100.0, "Finalizado"),
    (re.compile(r"export annotation results|human readable tsv|gff3|genbank|embl|machine readable json|genome and annotation summary", re.I), 98.0, "Gravando resultados"),
    (re.compile(r"genome statistics|annotation summary", re.I), 97.0, "Gerando resumo"),
    (re.compile(r"improve annotations", re.I), 96.0, "Aprimorando anotações"),
    (re.compile(r"select features and create locus tags", re.I), 95.0, "Selecionando features e locus tags"),
    (re.compile(r"apply feature overlap filters", re.I), 94.0, "Resolvendo sobreposições"),
    (re.compile(r"detect oriCs?/oriVs?|detect oriTs?", re.I), 92.0, "Detectando origens de replicação"),
    (re.compile(r"detect gaps", re.I), 90.0, "Detectando gaps"),
    (re.compile(r"detect & annotate sORF|sORFs?", re.I), 88.0, "Anotando sORFs"),
    (re.compile(r"revise special cases", re.I), 84.0, "Revisando casos especiais"),
    (re.compile(r"analyze hypothetical proteins|Pfam hits|proteins statistics", re.I), 81.0, "Analisando proteínas hipotéticas"),
    (re.compile(r"detect pseudogenes|pseudogene candidates|verified", re.I), 76.0, "Detectando pseudogenes"),
    (re.compile(r"combine annotations and mark hypotheticals", re.I), 70.0, "Combinando anotações"),
    (re.compile(r"conduct expert systems|amrfinder|protein sequences", re.I), 64.0, "Aplicando sistemas especialistas"),
    (re.compile(r"lookup annotations", re.I), 58.0, "Consultando anotações"),
    (re.compile(r"found PSCs?|found PSCCs?", re.I), 52.0, "Buscando clusters proteicos"),
    (re.compile(r"detected IPSs?", re.I), 45.0, "Identificando sequências proteicas"),
    (re.compile(r"predict & annotate CDSs?|predicted:", re.I), 38.0, "Predizendo e anotando CDS"),
    (re.compile(r"predict CRISPR arrays?", re.I), 32.0, "Detectando CRISPR"),
    (re.compile(r"predict ncRNA regions?", re.I), 28.0, "Anotando regiões de ncRNA"),
    (re.compile(r"predict ncRNAs?", re.I), 24.0, "Anotando ncRNAs"),
    (re.compile(r"predict rRNAs?", re.I), 20.0, "Anotando rRNAs"),
    (re.compile(r"predict tmRNAs?", re.I), 16.0, "Anotando tmRNA"),
    (re.compile(r"predict tRNAs?", re.I), 12.0, "Anotando tRNAs"),
    (re.compile(r"start annotation", re.I), 8.0, "Iniciando anotação"),
    (re.compile(r"parse genome sequences|imported:|contigs:|chromosomes:|plasmids:", re.I), 5.0, "Lendo sequências"),
    (re.compile(r"Bakta v|options and arguments|db:|database", re.I), 3.0, "Validando configurações"),
)


def format_duration(seconds: float) -> str:
    """Formata segundos como HH:MM:SS."""
    seconds = max(0, int(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def progress_bar(percent: float, width: int = 30) -> str:
    """Cria uma barra textual de progresso."""
    percent = max(0.0, min(100.0, float(percent)))
    filled = int(round(width * percent / 100.0))
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def terminal_clear() -> None:
    """Limpa o terminal somente quando a saída é interativa."""
    if sys.stdout.isatty():
        print("\033[2J\033[H", end="")


def compact_text(value: str, limit: int) -> str:
    """Encurta texto longo preservando início e fim."""
    value = " ".join(str(value or "").split())
    if len(value) <= limit:
        return value
    if limit < 10:
        return value[:limit]
    left = (limit - 3) // 2
    right = limit - 3 - left
    return value[:left] + "..." + value[-right:]


def read_json_safe(path: str) -> dict:
    """Lê JSON sem interromper o monitor se ele estiver sendo atualizado."""
    try:
        with open(path, "rt", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError, json.JSONDecodeError):
        return {}


def infer_bakta_progress(
    line: str,
    current_percent: float,
    current_stage: str,
) -> Tuple[float, str]:
    """Infere uma etapa aproximada a partir das mensagens emitidas pelo Bakta."""
    clean = " ".join((line or "").split())
    if not clean:
        return current_percent, current_stage
    best_percent = current_percent
    best_stage = current_stage
    for rx, percent, stage in PROGRESS_STAGE_RULES:
        if rx.search(clean) and percent > best_percent:
            best_percent = percent
            best_stage = stage
    return best_percent, best_stage


def queue_progress(queue_dir: str) -> dict:
    """Calcula progresso exato por jobs e estimado incluindo jobs em execução."""
    counts = queue_counts(queue_dir)
    total = sum(counts.values())
    completed = counts["done"] + counts["failed"]
    exact_percent = (100.0 * completed / total) if total else 0.0

    running_progress = []
    running_dir = os.path.join(queue_dir, "running")
    for path in glob.glob(os.path.join(running_dir, "*.json")):
        job = read_json_safe(path)
        progress = job.get("progress") or {}
        try:
            percent = float(progress.get("estimated_percent", 0.0))
        except (TypeError, ValueError):
            percent = 0.0
        running_progress.append(max(0.0, min(100.0, percent)))

    estimated_units = completed + sum(p / 100.0 for p in running_progress)
    estimated_percent = (100.0 * estimated_units / total) if total else 0.0

    return {
        "counts": counts,
        "total": total,
        "completed": completed,
        "exact_percent": exact_percent,
        "estimated_percent": estimated_percent,
    }


def render_worker_dashboard(
    job: dict,
    queue_dir: str,
    elapsed_seconds: float,
    stage: str,
    estimated_percent: float,
    last_message: str,
    cpus: int,
    db_path: str,
    log_path: str,
    process_pid: Optional[int],
) -> None:
    """Mostra o painel do job que este worker está executando."""
    width = max(80, min(140, shutil.get_terminal_size((110, 30)).columns))
    q = queue_progress(queue_dir)
    counts = q["counts"]
    terminal_clear()
    print("=" * width)
    print(" BAKTA FLEET — MONITOR DO WORKER")
    print("=" * width)
    print(f"Amostra        : {job.get('sample', '?')}")
    print(f"Arquivo        : {job.get('fna', '?')}")
    print(f"Worker         : {job.get('owner_id', '?')}")
    print(f"Hostname       : {job.get('owner_hostname', '?')}")
    print(f"IP             : {job.get('owner_ip', '?')}")
    print(f"PID do Bakta   : {process_pid if process_pid is not None else '?'}")
    print(f"CPUs           : {cpus}")
    print(f"Banco          : {db_path}")
    print(f"Saída          : {job.get('outdir', job.get('output_dir', '?'))}")
    print(f"Log            : {log_path}")
    print("-" * width)
    print(
        f"Genoma atual   : {progress_bar(estimated_percent, 36)} "
        f"{estimated_percent:6.1f}%  (estimado)"
    )
    print(f"Etapa          : {stage}")
    print(f"Tempo decorrido: {format_duration(elapsed_seconds)}")
    print("-" * width)
    print(
        f"Fila estimada  : {progress_bar(q['estimated_percent'], 36)} "
        f"{q['estimated_percent']:6.1f}%"
    )
    print(
        f"Jobs concluídos: {q['completed']}/{q['total']} "
        f"({q['exact_percent']:.1f}% exato por jobs)"
    )
    print(
        "Estados        : "
        f"pending={counts['pending']}  running={counts['running']}  "
        f"done={counts['done']}  failed={counts['failed']}"
    )
    print("-" * width)
    print("Última mensagem do Bakta:")
    print("  " + compact_text(last_message or "Aguardando saída do programa...", width - 4))
    print("=" * width)
    sys.stdout.flush()


def persist_job_progress(
    claimed: str,
    job: dict,
    stage: str,
    estimated_percent: float,
    last_message: str,
    started_monotonic: float,
) -> None:
    """Atualiza o JSON do job para que outros terminais vejam o progresso."""
    job["progress"] = {
        "stage": stage,
        "estimated_percent": round(float(estimated_percent), 1),
        "last_message": " ".join((last_message or "").split())[:1000],
        "elapsed_seconds": int(max(0.0, time.monotonic() - started_monotonic)),
        "updated": now_iso(),
    }
    write_json_atomic(claimed, job)


def run_bakta_streamed(
    cmd: list,
    log_path: str,
    claimed: str,
    job: dict,
    queue_dir: str,
    cpus: int,
    db_path: str,
    dashboard: bool,
    dashboard_interval: float,
) -> int:
    """Executa o Bakta, transmite o log e atualiza o painel/JSON em tempo real."""
    started_monotonic = time.monotonic()
    stage = "Inicializando o Bakta"
    estimated_percent = 1.0
    last_message = "Processo iniciado"
    last_refresh = 0.0
    last_persist = 0.0

    with open(log_path, "wt", encoding="utf-8") as logf:
        logf.write("CMD: " + shlex.join(cmd) + "\n")
        logf.flush()

        process_env = os.environ.copy()
        # O Bakta é um programa Python. Quando a saída é redirecionada para um
        # pipe, o buffering pode atrasar as mensagens usadas pelo painel.
        process_env["PYTHONUNBUFFERED"] = "1"

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
            env=process_env,
        )
        job["bakta_pid"] = process.pid
        persist_job_progress(
            claimed, job, stage, estimated_percent, last_message, started_monotonic
        )

        output_queue: _queue.Queue = _queue.Queue()

        def reader() -> None:
            try:
                assert process.stdout is not None
                for output_line in process.stdout:
                    output_queue.put(output_line)
            finally:
                output_queue.put(None)

        reader_thread = threading.Thread(target=reader, daemon=True)
        reader_thread.start()
        stream_closed = False

        try:
            while True:
                try:
                    item = output_queue.get(timeout=0.25)
                    if item is None:
                        stream_closed = True
                    else:
                        logf.write(item)
                        logf.flush()
                        clean = item.rstrip("\r\n")
                        if clean.strip():
                            last_message = clean
                            estimated_percent, stage = infer_bakta_progress(
                                clean, estimated_percent, stage
                            )
                except _queue.Empty:
                    pass

                now_mono = time.monotonic()
                if now_mono - last_persist >= 1.0:
                    persist_job_progress(
                        claimed,
                        job,
                        stage,
                        estimated_percent,
                        last_message,
                        started_monotonic,
                    )
                    last_persist = now_mono

                if dashboard and now_mono - last_refresh >= dashboard_interval:
                    render_worker_dashboard(
                        job=job,
                        queue_dir=queue_dir,
                        elapsed_seconds=now_mono - started_monotonic,
                        stage=stage,
                        estimated_percent=estimated_percent,
                        last_message=last_message,
                        cpus=cpus,
                        db_path=db_path,
                        log_path=log_path,
                        process_pid=process.pid,
                    )
                    last_refresh = now_mono

                if stream_closed and process.poll() is not None:
                    break
        except KeyboardInterrupt:
            # Interrompe também o processo Bakta para não deixá-lo órfão.
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise

        exit_code = process.wait()
        if exit_code == 0:
            estimated_percent = 100.0
            stage = "Finalizado"
            last_message = "Bakta finalizado; validando arquivos de saída"
        else:
            stage = "Falha na execução"
            last_message = f"Bakta terminou com código {exit_code}"

        persist_job_progress(
            claimed, job, stage, estimated_percent, last_message, started_monotonic
        )
        if dashboard:
            render_worker_dashboard(
                job=job,
                queue_dir=queue_dir,
                elapsed_seconds=time.monotonic() - started_monotonic,
                stage=stage,
                estimated_percent=estimated_percent,
                last_message=last_message,
                cpus=cpus,
                db_path=db_path,
                log_path=log_path,
                process_pid=process.pid,
            )
        return exit_code


def send_completion_email(
    email_to: str,
    smtp_host: str,
    smtp_port: int,
    smtp_user: str,
    smtp_password_env: str,
    output_root: str,
    counts: dict,
    queue_dir: Optional[str] = None,
) -> None:
    """Envia o resumo final via SMTP com STARTTLS, sem armazenar a senha."""
    password = os.environ.get(smtp_password_env, "")
    if not password:
        raise RuntimeError(
            f"variável de ambiente {smtp_password_env} não definida"
        )

    queue_dir = queue_dir or os.path.join(output_root, QUEUE_FOLDER_NAME)
    execution = queue_execution_summary(queue_dir)
    notifier_hostname = socket.gethostname()
    notifier_ip = get_local_ip()
    notifier_stats = execution["per_ip"].get(
        notifier_ip,
        {"total": 0, "done": 0, "failed": 0},
    )
    started_at = execution.get("started_at") or "não registrado"
    finished_at = now_iso()

    failed = int(counts.get("failed", 0))
    status_text = "concluído com falhas" if failed else "concluído com sucesso"
    subject = (
        f"[Bakta] Processamento {status_text}: "
        f"{counts.get('done', 0)} concluído(s), {failed} falha(s)"
    )
    body = (
        "O processamento da fila do Bakta terminou.\n\n"
        f"Status: {status_text}\n"
        f"Concluídos: {counts.get('done', 0)}\n"
        f"Falhas: {failed}\n"
        f"Pendentes: {counts.get('pending', 0)}\n"
        f"Em execução: {counts.get('running', 0)}\n"
        f"Computador notificador: {notifier_hostname}\n"
        f"IP do computador notificador: {notifier_ip}\n"
        f"Jobs processados por este IP: {notifier_stats.get('total', 0)} "
        f"({notifier_stats.get('done', 0)} concluído(s), "
        f"{notifier_stats.get('failed', 0)} falha(s))\n"
        f"Início: {started_at}\n"
        f"Término: {finished_at}\n"
        "\nProcessamento por IP:\n"
        f"{format_ip_work_summary(execution['per_ip'])}"
        f"\nDiretório de resultados: {output_root}\n"
    )

    message = EmailMessage()
    message["From"] = smtp_user
    message["To"] = email_to
    message["Subject"] = subject
    message.set_content(body)

    context = ssl.create_default_context()
    with smtplib.SMTP(smtp_host, smtp_port, timeout=30) as server:
        server.ehlo()
        server.starttls(context=context)
        server.ehlo()
        server.login(smtp_user, password)
        server.send_message(message)

def maybe_notify_completion(
    queue_dir: str,
    output_root: str,
    worker_id: str,
    email_to: str,
    smtp_host: str,
    smtp_port: int,
    smtp_user: str,
    smtp_password_env: str,
) -> bool:
    """Envia uma única notificação quando não há jobs pendentes ou rodando."""
    counts = queue_counts(queue_dir)
    if counts["pending"] or counts["running"]:
        return False
    if counts["done"] + counts["failed"] == 0:
        return False

    sent_marker = os.path.join(queue_dir, EMAIL_SENT_MARKER)
    lock_marker = os.path.join(queue_dir, EMAIL_LOCK_MARKER)
    failure_marker = os.path.join(queue_dir, EMAIL_FAILURE_MARKER)
    if os.path.isfile(sent_marker):
        return False
    if os.path.isfile(failure_marker):
        age = time.time() - os.path.getmtime(failure_marker)
        if age < EMAIL_RETRY_SECONDS:
            return False

    try:
        lock_fd = os.open(lock_marker, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        try:
            if time.time() - os.path.getmtime(lock_marker) > EMAIL_RETRY_SECONDS:
                os.remove(lock_marker)
        except OSError:
            pass
        return False

    try:
        with os.fdopen(lock_fd, "wt", encoding="utf-8") as lock_fh:
            lock_fh.write(f"worker={worker_id}\nstarted={now_iso()}\n")

        counts = queue_counts(queue_dir)
        if counts["pending"] or counts["running"] or os.path.isfile(sent_marker):
            return False

        send_completion_email(
            email_to=email_to,
            smtp_host=smtp_host,
            smtp_port=smtp_port,
            smtp_user=smtp_user,
            smtp_password_env=smtp_password_env,
            output_root=output_root,
            counts=counts,
            queue_dir=queue_dir,
        )
        write_json_atomic(
            sent_marker,
            {
                "status": "sent",
                "sent_at": now_iso(),
                "worker_id": worker_id,
                "email_to": email_to,
                "smtp_host": smtp_host,
                "counts": counts,
                "notifier_hostname": socket.gethostname(),
                "notifier_ip": get_local_ip(),
                "execution_summary": queue_execution_summary(queue_dir),
            },
        )
        if os.path.isfile(failure_marker):
            os.remove(failure_marker)
        print(f"[{worker_id}] E-mail de conclusão enviado para {email_to}")
        return True
    except Exception as exc:
        write_json_atomic(
            failure_marker,
            {
                "status": "failed",
                "failed_at": now_iso(),
                "worker_id": worker_id,
                "email_to": email_to,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "retry_after_seconds": EMAIL_RETRY_SECONDS,
            },
        )
        print(
            f"[{worker_id}] AVISO: falha ao enviar e-mail: {exc}",
            file=sys.stderr,
        )
        return False
    finally:
        try:
            os.remove(lock_marker)
        except FileNotFoundError:
            pass

def bakta_cmd(
    fna: str,
    outdir: str,
    prefix: str,
    locustag: str,
    cpus: int,
    db_path: str,
    genus: str,
    species: str,
    strain: str,
    force: bool,
    complete: bool,
    compliant: bool,
    keep_contig_headers: bool,
    skip_plot: bool,
    gram: str,
    translation_table: int,
    tmp_dir: Optional[str],
) -> list:
    """Monta o comando Bakta sem usar shell=True."""
    cmd = [
        "bakta",
        "--verbose",
        "--db", db_path,
        "--output", outdir,
        "--prefix", prefix,
        "--locus-tag", locustag,
        "--threads", str(cpus),
        "--gram", gram,
        "--translation-table", str(translation_table),
    ]
    if force:
        cmd.append("--force")
    if complete:
        cmd.append("--complete")
    if compliant:
        cmd.append("--compliant")
    if keep_contig_headers:
        cmd.append("--keep-contig-headers")
    if skip_plot:
        cmd.append("--skip-plot")
    if tmp_dir:
        cmd += ["--tmp-dir", tmp_dir]
    if genus:
        cmd += ["--genus", genus]
    if species:
        cmd += ["--species", species]
    if strain:
        cmd += ["--strain", strain]
    cmd.append(fna)
    return cmd

def expected_bakta_outputs(outdir: str, prefix: str) -> list:
    """Arquivos mínimos usados para confirmar que a anotação terminou."""
    return [
        os.path.join(outdir, f"{prefix}.json"),
        os.path.join(outdir, f"{prefix}.gff3"),
        os.path.join(outdir, f"{prefix}.faa"),
        os.path.join(outdir, f"{prefix}.tsv"),
    ]

def work_loop(
    input_root: str,
    output_root: str,
    cpus: int,
    force: bool,
    idle_sleep: float,
    worker_id: str,
    once: bool,
    until_empty: bool,
    notify_email: bool,
    email_to: str,
    smtp_host: str,
    smtp_port: int,
    smtp_user: str,
    smtp_password_env: str,
    db_path: Optional[str],
    complete: bool,
    compliant: bool,
    keep_contig_headers: bool,
    skip_plot: bool,
    gram: str,
    translation_table: int,
    tmp_dir: Optional[str],
    dashboard: bool,
    dashboard_interval: float,
) -> int:
    input_root = os.path.abspath(os.path.expanduser(input_root))
    output_root = os.path.abspath(os.path.expanduser(output_root))
    db_path = resolve_db_path(db_path)
    tmp_dir = (
        os.path.abspath(os.path.expanduser(tmp_dir)) if tmp_dir else None
    )

    runtime_error = validate_runtime(db_path)
    if runtime_error:
        print(f"[worker] ERRO: {runtime_error}", file=sys.stderr)
        return 127 if "executável" in runtime_error else 2
    if tmp_dir:
        os.makedirs(tmp_dir, exist_ok=True)
    queue_dir = os.path.join(output_root, QUEUE_FOLDER_NAME)
    pending = os.path.join(queue_dir, "pending")
    running = os.path.join(queue_dir, "running")
    done = os.path.join(queue_dir, "done")
    failed = os.path.join(queue_dir, "failed")
    logs_dir = os.path.join(queue_dir, "logs")

    for d in (pending, running, done, failed, logs_dir, output_root):
        os.makedirs(d, exist_ok=True)

    ip_auto = get_local_ip()
    owner_id = worker_id
    owner_ip = worker_id if IPV4_RE.match(worker_id) else ip_auto

    def notify_if_requested() -> bool:
        if not notify_email:
            return True
        return maybe_notify_completion(
            queue_dir=queue_dir,
            output_root=output_root,
            worker_id=owner_id,
            email_to=email_to,
            smtp_host=smtp_host,
            smtp_port=smtp_port,
            smtp_user=smtp_user,
            smtp_password_env=smtp_password_env,
        )

    while True:
        claimed = claim_one_job(pending, running)
        if not claimed:
            counts = queue_counts(queue_dir)
            queue_finished = counts["pending"] == 0 and counts["running"] == 0

            if queue_finished:
                notification_ok = notify_if_requested()
                if until_empty:
                    sent_marker = os.path.join(queue_dir, EMAIL_SENT_MARKER)
                    if notify_email and not notification_ok and not os.path.isfile(sent_marker):
                        print(
                            f"[{owner_id}] ERRO: fila terminou, mas o e-mail de "
                            "conclusão não foi enviado. Consulte "
                            f"{os.path.join(queue_dir, EMAIL_FAILURE_MARKER)}",
                            file=sys.stderr,
                        )
                        return 4
                    print(
                        f"[{owner_id}] FILA CONCLUÍDA: "
                        f"done={counts['done']} failed={counts['failed']} "
                        f"output={output_root}"
                    )
                    return 0 if counts["failed"] == 0 else 3

            if once:
                return 0
            time.sleep(idle_sleep)
            continue

        try:
            with open(claimed, "rt") as fh:
                job = json.load(fh)
        except Exception as e:
            bad_dst = os.path.join(failed, os.path.basename(claimed))
            os.replace(claimed, bad_dst)
            print(f"[{owner_id}] Job inválido -> failed: {bad_dst} ({e})", file=sys.stderr)
            if once:
                notify_if_requested()
                return 0
            continue

        folder_name = job.get("source_folder_name") or job.get("sample") or "amostra"
        sample = job.get("sample") or safe_sample_name(folder_name)
        base = sample
        fna = job.get("fna", "")
        gbff = job.get("gbff", "")

        if not os.path.isfile(fna):
            log_path = os.path.join(logs_dir, f"{safe_sample_name(owner_id)}_{sample}.log")
            with open(log_path, "wt") as logf:
                logf.write(f"ERRO: arquivo genomic.fna não encontrado: {fna}\n")
            job.update({
                "end": now_iso(),
                "exit_code": 2,
                "status": "failed",
                "log": log_path,
            })
            with open(claimed, "wt") as fh:
                json.dump(job, fh, indent=2, sort_keys=True)
            final_path = os.path.join(failed, os.path.basename(claimed))
            os.replace(claimed, final_path)
            print(f"[{owner_id}] FAIL {sample}: FNA não encontrado", file=sys.stderr)
            if once:
                notify_if_requested()
                return 0
            continue

        genus = species = acc = strain = ""
        if gbff and os.path.isfile(gbff):
            genus, species, acc, strain = parse_gbff(gbff)

        assembly_accession = job.get("assembly_accession") or extract_assembly_accession(
            folder_name, fna, gbff
        )
        accession_for_locus = assembly_accession or acc

        if (not genus or not species) and "_" in base:
            genus = genus or base.split("_", 1)[0]
            species = species or base.split("_", 1)[1].split("_", 1)[0]

        if strain:
            strain = sanitize(strain.replace(" ", "_"))
            strain = shorten_strain(strain)
            strain = sanitize(strain)

        prefix = safe_sample_name(folder_name)
        locustag = build_locustag(genus, species, accession_for_locus, base)

        outdir = job.get("output_dir") or os.path.join(output_root, folder_name)
        log_owner = safe_sample_name(owner_id)
        log_path = os.path.join(logs_dir, f"{log_owner}_{sample}.log")

        job.update({
            "owner_ip": owner_ip,
            "owner_id": owner_id,
            "owner_hostname": socket.gethostname(),
            "worker_python_pid": os.getpid(),
            "start": now_iso(),
            "prefix": prefix,
            "locustag": locustag,
            "outdir": outdir,
            "genus": genus,
            "species": species.lower() if species else "",
            "strain": strain,
            "assembly_accession": assembly_accession,
            "sequence_accession_first_record": normalize_accession(acc),
            "bakta_db": db_path,
            "bakta_options": {
                "complete": complete,
                "compliant": compliant,
                "keep_contig_headers": keep_contig_headers,
                "skip_plot": skip_plot,
                "gram": gram,
                "translation_table": translation_table,
                "tmp_dir": tmp_dir or "",
                "threads": cpus,
                "verbose": True,
            },
            "log": log_path,
        })
        with open(claimed, "wt") as fh:
            json.dump(job, fh, indent=2, sort_keys=True)

        cmd = bakta_cmd(
            fna=fna,
            outdir=outdir,
            prefix=prefix,
            locustag=locustag,
            cpus=cpus,
            db_path=db_path,
            genus=genus,
            species=species.lower() if species else "",
            strain=strain,
            force=force,
            complete=complete,
            compliant=compliant,
            keep_contig_headers=keep_contig_headers,
            skip_plot=skip_plot,
            gram=gram,
            translation_table=translation_table,
            tmp_dir=tmp_dir,
        )

        print(f"[{owner_id}] RUN {base} -> prefix={prefix} locustag={locustag}")
        exit_code = 1
        try:
            exit_code = run_bakta_streamed(
                cmd=cmd,
                log_path=log_path,
                claimed=claimed,
                job=job,
                queue_dir=queue_dir,
                cpus=cpus,
                db_path=db_path,
                dashboard=dashboard,
                dashboard_interval=dashboard_interval,
            )
            if exit_code == 0:
                missing_outputs = [
                    path for path in expected_bakta_outputs(outdir, prefix)
                    if not os.path.isfile(path)
                ]
                if missing_outputs:
                    with open(log_path, "a", encoding="utf-8") as logf:
                        logf.write(
                            "ERRO: Bakta terminou com código 0, mas faltam "
                            "arquivos obrigatórios:\n"
                        )
                        for missing in missing_outputs:
                            logf.write(f"  - {missing}\n")
                    exit_code = 3
        except FileNotFoundError:
            with open(log_path, "a", encoding="utf-8") as logf:
                logf.write("ERRO: bakta não encontrado no PATH.\n")
            exit_code = 127
        except KeyboardInterrupt:
            with open(log_path, "a", encoding="utf-8") as logf:
                logf.write("ERRO: execução interrompida pelo usuário.\n")
            exit_code = 130
        except Exception as e:
            with open(log_path, "a", encoding="utf-8") as logf:
                logf.write(f"ERRO ao executar Bakta: {e}\n")
            exit_code = 1

        job["end"] = now_iso()
        job["exit_code"] = exit_code
        job["status"] = "done" if exit_code == 0 else "failed"
        final_progress = job.get("progress") or {}
        final_progress["stage"] = (
            "Concluído" if exit_code == 0 else f"Falhou (código {exit_code})"
        )
        final_progress["estimated_percent"] = 100.0 if exit_code == 0 else float(
            final_progress.get("estimated_percent", 0.0)
        )
        final_progress["updated"] = now_iso()
        job["progress"] = final_progress

        with open(claimed, "wt") as fh:
            json.dump(job, fh, indent=2, sort_keys=True)

        dst_folder = done if exit_code == 0 else failed
        final_path = os.path.join(dst_folder, os.path.basename(claimed))
        os.replace(claimed, final_path)

        if exit_code == 0:
            print(f"[{owner_id}] OK  {base} -> {prefix}")
        else:
            print(f"[{owner_id}] FAIL {base} (exit={exit_code}) log={log_path}", file=sys.stderr)

        if once:
            notify_if_requested()
            return 0

def status(output_root: str, watch: bool, interval: float) -> int:
    """Mostra um painel da fila e dos workers ativos."""
    output_root = os.path.abspath(os.path.expanduser(output_root))
    queue_dir = os.path.join(output_root, QUEUE_FOLDER_NAME)
    running_dir = os.path.join(queue_dir, "running")

    while True:
        q = queue_progress(queue_dir)
        counts = q["counts"]
        width = max(80, min(140, shutil.get_terminal_size((110, 30)).columns))
        if watch:
            terminal_clear()

        print("=" * width)
        print(" BAKTA FLEET — STATUS DA FILA")
        print("=" * width)
        print(f"Diretório       : {queue_dir}")
        print(
            f"Progresso       : {progress_bar(q['estimated_percent'], 36)} "
            f"{q['estimated_percent']:6.1f}% estimado"
        )
        print(
            f"Concluídos      : {q['completed']}/{q['total']} "
            f"({q['exact_percent']:.1f}% exato por jobs)"
        )
        print(
            "Estados         : "
            f"pending={counts['pending']}  running={counts['running']}  "
            f"done={counts['done']}  failed={counts['failed']}"
        )

        running_paths = sorted(glob.glob(os.path.join(running_dir, "*.json")))
        if running_paths:
            print("-" * width)
            print("WORKERS EM EXECUÇÃO")
            for path in running_paths[:50]:
                job = read_json_safe(path)
                progress = job.get("progress") or {}
                try:
                    percent = float(progress.get("estimated_percent", 0.0))
                except (TypeError, ValueError):
                    percent = 0.0
                elapsed = progress.get("elapsed_seconds", 0)
                try:
                    elapsed_text = format_duration(float(elapsed))
                except (TypeError, ValueError):
                    elapsed_text = "??:??:??"
                print("-" * width)
                print(
                    f"Amostra : {job.get('sample', os.path.basename(path))} | "
                    f"{progress_bar(percent, 24)} {percent:5.1f}%"
                )
                print(
                    f"Worker  : {job.get('owner_id', '?')} | "
                    f"host={job.get('owner_hostname', '?')} | "
                    f"ip={job.get('owner_ip', '?')} | "
                    f"PID={job.get('bakta_pid', '?')}"
                )
                print(
                    f"Etapa   : {progress.get('stage', 'Aguardando atualização')} | "
                    f"tempo={elapsed_text}"
                )
                print(f"Arquivo : {job.get('fna', '?')}")
                print(
                    "Mensagem: "
                    + compact_text(
                        progress.get("last_message", "Sem mensagem"), width - 10
                    )
                )
        elif q["total"] == 0:
            print("-" * width)
            print("Nenhum job foi encontrado nessa fila.")
            print("Confira o caminho acima ou execute o comando init para criar a fila.")
        else:
            print("-" * width)
            print("Nenhum worker está executando um job neste momento.")

        print("=" * width)
        sys.stdout.flush()
        if not watch:
            return 0
        time.sleep(interval)



def main():
    ap = argparse.ArgumentParser(
        prog=SCRIPT_NAME,
        description=(
            "Etapa anot00/01: percorre as pastas de genomas do NCBI e executa "
            "Bakta por uma fila coordinator–worker."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}",
        help="Exibe o nome canônico e a versão do script",
    )
    ap.add_argument(
        "--input-root",
        default=DEFAULT_INPUT_ROOT,
        help=(
            "Diretório cujas subpastas contêm os arquivos *_genomic.fna e "
            "*_genomic.gbff do NCBI"
        ),
    )
    ap.add_argument(
        "--output-root",
        default=DEFAULT_OUTPUT_ROOT,
        help="Diretório de saída; cada amostra recebe uma subpasta homônima",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_init = sub.add_parser(
        "init",
        help="Cria jobs a partir das subpastas de etapa00_genomas_ncbi",
    )
    ap_init.add_argument("--overwrite", action="store_true", help="Limpa pending/running/done/failed antes de recriar.")

    ap_work = sub.add_parser("work", help="Worker: pega jobs de pending e executa Bakta.")
    ap_work.add_argument("--cpus", type=int, default=4, help="CPUs por execução do Bakta.")
    ap_work.add_argument("--force", action="store_true", help="Passa --force ao Bakta (sobrescreve outdir se existir).")
    ap_work.add_argument(
        "--db",
        default=None,
        help=(
            "Banco do Bakta. Prioridade: --db, variável BAKTA_DB e "
            f"caminho padrão {DEFAULT_DB_PATH}"
        ),
    )
    ap_work.add_argument(
        "--complete",
        action="store_true",
        help=(
            "Passa --complete; use somente quando TODAS as sequências de "
            "cada FASTA forem replicons completos"
        ),
    )
    ap_work.add_argument(
        "--compliant",
        action="store_true",
        help="Passa --compliant para gerar anotação compatível com INSDC",
    )
    ap_work.add_argument(
        "--keep-contig-headers",
        action="store_true",
        help="Preserva os cabeçalhos originais dos contigs",
    )
    ap_work.add_argument(
        "--skip-plot",
        action="store_true",
        help="Não gera os gráficos PNG/SVG do Bakta",
    )
    ap_work.add_argument(
        "--gram",
        choices=("+", "-", "?"),
        default="?",
        help="Tipo de Gram usado na predição de peptídeos-sinal",
    )
    ap_work.add_argument(
        "--translation-table",
        type=int,
        choices=(11, 4, 25),
        default=11,
        help="Tabela de tradução genética",
    )
    ap_work.add_argument(
        "--tmp-dir",
        default=None,
        help="Diretório temporário local opcional para o Bakta",
    )
    ap_work.add_argument("--idle-sleep", type=float, default=2.0, help="Espera quando não há jobs.")
    ap_work.add_argument("--id", required=True, help="ID do worker (ex.: 192.168.0.104).")
    execution_mode = ap_work.add_mutually_exclusive_group()
    execution_mode.add_argument(
        "--once",
        action="store_true",
        help="Executa apenas 1 job e encerra (modo de teste)",
    )
    execution_mode.add_argument(
        "--until-empty",
        action="store_true",
        help=(
            "Processa todos os jobs até pending=0 e running=0; depois envia "
            "a notificação solicitada e encerra automaticamente"
        ),
    )
    ap_work.add_argument(
        "--no-dashboard",
        action="store_true",
        help="Desativa o painel visual no terminal do worker",
    )
    ap_work.add_argument(
        "--dashboard-interval",
        type=float,
        default=1.0,
        help="Intervalo, em segundos, entre atualizações do painel do worker",
    )
    ap_work.add_argument(
        "--notify-email",
        action="store_true",
        help="Este worker envia um único e-mail quando toda a fila terminar",
    )
    ap_work.add_argument(
        "--email-to",
        default=DEFAULT_EMAIL_TO,
        help="Destinatário do aviso de conclusão",
    )
    ap_work.add_argument(
        "--smtp-host",
        default=DEFAULT_SMTP_HOST,
        help="Servidor SMTP",
    )
    ap_work.add_argument(
        "--smtp-port",
        type=int,
        default=DEFAULT_SMTP_PORT,
        help="Porta SMTP com STARTTLS",
    )
    ap_work.add_argument(
        "--smtp-user",
        default=DEFAULT_SMTP_USER,
        help="Conta remetente usada na autenticação SMTP",
    )
    ap_work.add_argument(
        "--smtp-password-env",
        default=DEFAULT_SMTP_PASSWORD_ENV,
        help="Nome da variável de ambiente que contém a senha de aplicativo",
    )

    ap_status = sub.add_parser("status", help="Mostra status e quem está rodando o quê.")
    ap_status.add_argument("--watch", action="store_true", help="Atualiza continuamente.")
    ap_status.add_argument("--interval", type=float, default=5.0, help="Intervalo do --watch.")

    args = ap.parse_args()

    if args.cmd == "work" and args.smtp_port < 1:
        ap.error("--smtp-port deve ser maior que zero")
    if args.cmd == "work" and args.cpus < 1:
        ap.error("--cpus deve ser maior que zero")
    if args.cmd == "work" and args.dashboard_interval <= 0:
        ap.error("--dashboard-interval deve ser maior que zero")
    if args.cmd == "work" and args.notify_email:
        if not args.email_to:
            ap.error("--notify-email exige --email-to ou a variável BAKTA_EMAIL_TO")
        if not args.smtp_user:
            ap.error("--notify-email exige --smtp-user ou a variável BAKTA_SMTP_USER")

    if args.cmd == "init":
        init_queue(args.input_root, args.output_root, overwrite=args.overwrite)
        return 0
    if args.cmd == "work":
        return work_loop(
            args.input_root,
            args.output_root,
            cpus=args.cpus,
            force=args.force,
            idle_sleep=args.idle_sleep,
            worker_id=args.id,
            once=args.once,
            until_empty=args.until_empty,
            notify_email=args.notify_email,
            email_to=args.email_to,
            smtp_host=args.smtp_host,
            smtp_port=args.smtp_port,
            smtp_user=args.smtp_user,
            smtp_password_env=args.smtp_password_env,
            db_path=args.db,
            complete=args.complete,
            compliant=args.compliant,
            keep_contig_headers=args.keep_contig_headers,
            skip_plot=args.skip_plot,
            gram=args.gram,
            translation_table=args.translation_table,
            tmp_dir=args.tmp_dir,
            dashboard=not args.no_dashboard,
            dashboard_interval=args.dashboard_interval,
        )
    if args.cmd == "status":
        return status(args.output_root, watch=args.watch, interval=args.interval)
    return 1

if __name__ == "__main__":
    raise SystemExit(main())