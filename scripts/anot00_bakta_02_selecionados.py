#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
"""Etapa 02 do pipeline Bakta: seleciona alvos para reanotação funcional.

O script lê o JSON e o FAA produzidos pelo Bakta, relaciona as features pelo
``locus``/identificador FASTA e classifica CDS em prioridades ``high``,
``medium`` e ``low``. A saída ``<amostra>.priority_high.faa`` é o contrato de
entrada da etapa Swiss-Prot.

Princípios desta etapa
----------------------
* JSON do Bakta como fonte principal de metadados.
* FAA do Bakta como fonte contratual das sequências proteicas.
* Verificação de identidade entre a sequência no JSON e no FAA.
* Proteínas sORF separadas do fluxo convencional.
* Pseudogenes excluídos do FASTA de busca, mas registrados na auditoria.
* Escrita atômica, hashes SHA-256 e manifesto por amostra.
* Reexecução segura: resultados válidos são ignorados; resultados parciais
  exigem ``--force`` após revisão.
* Linhas de inference.tsv sem locus (por exemplo, ``-``) são ignoradas.
* O arquivo hypotheticals.faa é a autoridade para a categoria hipotética oficial; o JSON é usado como metadado complementar e para auditoria.

Modo em lote::

    python3 anot00_bakta_02_selecionados.py

Modo direto::

    python3 anot00_bakta_02_selecionados.py \
      --bakta-dir /caminho/output00_bakta/AMOSTRA \
      --outdir /caminho/p03_target_priority/bakta/AMOSTRA

A raiz de saída e as subpastas por amostra são criadas automaticamente.

Requisitos: Python 3.9+. Somente biblioteca padrão do Python.
Compatibilidade de referência: Bakta 1.12.x e JSONs recentes do Bakta.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PIPELINE_NAME = "anot00_bakta"
PIPELINE_STAGE = "02_selecionados"
SCRIPT_NAME = "anot00_bakta_02_selecionados.py"
SCRIPT_VERSION = "1.0.4"
NEXT_SCRIPT = "anot01_swissprot_01_search_fleet.py"

DEFAULT_BAKTA_ROOT = Path(
    "/home/bioinfo/anotacao_coordinator_worker/output00_bakta"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/home/bioinfo/anotacao_coordinator_worker/output00_bakta_selecao"
)

QUEUE_DIRECTORY_NAMES = {"bakta_queue", "queue", "logs"}
EMPTY_VALUES = {"", "-", ".", "NA", "N/A", "NONE", "NULL", "NOT_FOUND"}

WEAK_EXACT_PRODUCTS = {
    "hypothetical protein",
    "conserved hypothetical protein",
    "uncharacterized protein",
    "uncharacterised protein",
    "unknown protein",
    "protein of unknown function",
    "putative protein",
    "predicted protein",
}
GENERIC_EXACT_PRODUCTS = {
    "protein",
    "conserved protein",
    "membrane protein",
    "domain-containing protein",
    "family protein",
}
WEAK_PRODUCT_PATTERNS = (
    re.compile(r"\bhypothet(?:ical|h?ethical)\b", re.I),
    re.compile(r"\buncharacteri[sz]ed\b", re.I),
    re.compile(r"\bprotein of unknown function\b", re.I),
    re.compile(r"\bunknown protein\b", re.I),
)
GENERIC_PRODUCT_PATTERNS = (
    re.compile(r"\b(?:DUF|UPF)\d+\b", re.I),
    re.compile(r"\bdomain[- ]containing protein\b", re.I),
    re.compile(r"\bfamily protein\b", re.I),
)
UNCERTAIN_PRODUCT_PATTERNS = (
    re.compile(r"\bputative\b", re.I),
    re.compile(r"\bpredicted\b", re.I),
    re.compile(r"\bprobable\b", re.I),
    re.compile(r"\bpossible\b", re.I),
    re.compile(r"\bpotential\b", re.I),
)

PRIORITY_ORDER = {"high": 0, "medium": 1, "low": 2}

TARGET_COLUMNS = (
    "sample",
    "source_annotation",
    "feature_id",
    "locus_tag",
    "feature_type",
    "contig",
    "start",
    "stop",
    "strand",
    "aa_length",
    "original_product",
    "original_gene",
    "original_genes",
    "original_ec",
    "original_cog",
    "original_go",
    "original_kegg",
    "original_pfam",
    "original_uniref",
    "original_refseq",
    "original_dbxrefs",
    "inference_method",
    "inference_accession",
    "inference_score",
    "inference_evalue",
    "inference_query_cov",
    "inference_subject_cov",
    "inference_identity",
    "truncated",
    "edge",
    "selection_flags",
    "target_priority_refined",
    "priority_rationale",
)

AUDIT_COLUMNS = (
    "sample",
    "locus_tag",
    "feature_type",
    "aa_length",
    "product",
    "gene",
    "decision",
    "priority",
    "reason",
    "flags",
)

REPORT_COLUMNS = ("sample", "category", "name", "count", "description")


class SelectionError(RuntimeError):
    """Erro operacional ou de consistência entre JSON, FAA e auxiliares."""


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    header: str
    sequence: str


@dataclass(frozen=True)
class InferenceRecord:
    method: str = ""
    accession: str = ""
    score: float | None = None
    evalue: float | None = None
    query_cov: float | None = None
    subject_cov: float | None = None
    identity: float | None = None


@dataclass(frozen=True)
class SampleInputs:
    sample: str
    prefix: str
    bakta_dir: Path
    bakta_json: Path
    bakta_faa: Path
    inference_tsv: Path | None
    hypotheticals_faa: Path | None
    hypotheticals_tsv: Path | None
    outdir: Path


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    bakta_dir: Path
    outdir: Path
    explicit_prefix: str | None = None


@dataclass(frozen=True)
class SelectionParameters:
    short_aa_threshold: int
    weak_identity: float
    weak_query_cov: float
    weak_subject_cov: float
    include_sorf_targets: bool
    include_pseudogenes: bool
    allow_id_mismatch: bool
    allow_sequence_mismatch: bool
    allow_hypothetical_mismatch: bool


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def is_empty(value: Any) -> bool:
    return clean(value).upper() in EMPTY_VALUES


def normalized_product(value: str) -> str:
    return " ".join(clean(value).casefold().split())


def bool_text(value: Any) -> str:
    return "true" if bool(value) else "false"


def safe_float(value: Any) -> float | None:
    text = clean(value)
    if not text or text.upper() in EMPTY_VALUES:
        return None
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def unique_clean(values: Iterable[Any]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = clean(value)
        if text and text not in seen:
            seen.add(text)
            result.append(text)
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_metadata(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_file(path),
    }


def optional_file_metadata(path: Path | None) -> dict[str, Any] | None:
    return file_metadata(path) if path is not None and path.is_file() else None


def fingerprint_matches(saved: Mapping[str, Any], path: Path) -> bool:
    if not path.is_file():
        return False
    if saved.get("size_bytes") != path.stat().st_size:
        return False
    expected = clean(saved.get("sha256", ""))
    return not expected or expected == sha256_file(path)


def optional_fingerprint_matches(saved: Any, path: Path | None) -> bool:
    if path is None:
        return saved is None
    return isinstance(saved, Mapping) and fingerprint_matches(saved, path)


def atomic_write_bytes(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "wb", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def atomic_write_text(path: Path, content: str) -> None:
    atomic_write_bytes(path, content.encode("utf-8"))


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    atomic_write_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def atomic_write_tsv(
    path: Path, rows: Iterable[Mapping[str, Any]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "wt",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        try:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fieldnames),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: clean(row.get(field, "")) for field in fieldnames})
            handle.flush()
            os.fsync(handle.fileno())
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def resolve_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise SelectionError(f"{label} não encontrado: {resolved}")
    return resolved


def resolve_directory(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise SelectionError(f"{label} não encontrado: {resolved}")
    return resolved


def path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def open_text_auto(path: Path):
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_json(path: Path) -> dict[str, Any]:
    try:
        with open_text_auto(path) as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SelectionError(f"JSON inválido: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SelectionError(f"JSON precisa conter um objeto: {path}")
    features = value.get("features")
    if not isinstance(features, list):
        raise SelectionError(f"JSON do Bakta sem lista 'features': {path}")
    return value


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    seen: set[str] = set()
    header: str | None = None
    sequence_parts: list[str] = []

    def finish() -> None:
        nonlocal header, sequence_parts
        if header is None:
            return
        identifier = header.split(maxsplit=1)[0]
        sequence = "".join(sequence_parts).replace(" ", "").upper()
        if sequence.endswith("*"):
            sequence = sequence[:-1]
        if not identifier or identifier in seen:
            raise SelectionError(
                f"Identificador FASTA vazio ou duplicado em {path}: {identifier!r}"
            )
        if not sequence:
            raise SelectionError(f"Sequência FASTA vazia em {path}: {identifier}")
        if re.search(r"[^A-Z.-]", sequence):
            raise SelectionError(f"Caracteres inválidos na sequência FASTA: {identifier}")
        seen.add(identifier)
        records.append(FastaRecord(identifier, header, sequence))

    with open_text_auto(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                sequence_parts = []
                if not header:
                    raise SelectionError(f"Cabeçalho FASTA vazio na linha {line_number}: {path}")
            else:
                if header is None:
                    raise SelectionError(
                        f"Sequência antes do cabeçalho na linha {line_number}: {path}"
                    )
                sequence_parts.append(line)
    finish()
    return records


def wrap_sequence(sequence: str, width: int = 60) -> Iterable[str]:
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def clean_target_header(header: str) -> str:
    result = re.sub(r"\s+\[target_priority=[^\]]*\]", "", header)
    result = re.sub(r"\s+\[selection_flags=[^\]]*\]", "", result)
    result = re.sub(r"\s+\[source_annotation=[^\]]*\]", "", result)
    return result.strip()


def write_selected_fasta(
    path: Path,
    records: Sequence[FastaRecord],
    label: str,
    flags_by_locus: Mapping[str, str],
) -> None:
    output: list[str] = []
    for record in records:
        base_header = clean_target_header(record.header)
        flags = flags_by_locus.get(record.identifier, "")
        output.append(
            f">{base_header} [source_annotation=bakta] "
            f"[target_priority={label}] [selection_flags={flags}]\n"
        )
        output.extend(f"{line}\n" for line in wrap_sequence(record.sequence))
    atomic_write_text(path, "".join(output))


def product_class(product: str) -> str:
    normalized = normalized_product(product)
    if not normalized:
        return "missing"
    if normalized in WEAK_EXACT_PRODUCTS:
        return "hypothetical"
    if any(pattern.search(product) for pattern in WEAK_PRODUCT_PATTERNS):
        return "hypothetical"
    if normalized in GENERIC_EXACT_PRODUCTS:
        return "generic"
    if any(pattern.search(product) for pattern in GENERIC_PRODUCT_PATTERNS):
        return "generic"
    if any(pattern.search(product) for pattern in UNCERTAIN_PRODUCT_PATTERNS):
        return "uncertain"
    return "informative"


def feature_locus(feature: Mapping[str, Any]) -> str:
    for key in ("locus", "locus_tag", "id"):
        value = clean(feature.get(key, ""))
        if value:
            return value
    return ""


def feature_type(feature: Mapping[str, Any]) -> str:
    return clean(feature.get("type", "")).casefold()


def feature_sequence_id(feature: Mapping[str, Any]) -> str:
    return clean(feature.get("sequence", feature.get("contig", "")))


def feature_is_pseudogene(feature: Mapping[str, Any]) -> bool:
    return "pseudogene" in feature and bool(feature.get("pseudogene"))


def feature_is_hypothetical(feature: Mapping[str, Any]) -> bool:
    """Identifica somente a classificação hipotética oficial do Bakta.

    O arquivo ``<prefix>.hypotheticals.faa`` contém o subconjunto oficial de
    CDS hipotéticas do Bakta. A mera presença da chave ``hypothetical`` não
    basta, pois alguns JSON podem armazená-la com valor falso. Produtos
    textualmente fracos, como ``uncharacterized protein``, continuam sendo
    selecionados em prioridade alta por ``product_class()``, mas não entram
    nesta validação contratual.
    """
    marker = feature.get("hypothetical")

    if isinstance(marker, bool):
        return marker
    if isinstance(marker, (int, float)):
        return bool(marker)
    if marker is not None:
        marker_text = clean(marker).casefold()
        if marker_text in {"true", "yes", "y", "1"}:
            return True
        if marker_text in {"false", "no", "n", "0", "-", ".", "none", "null"}:
            return False

    return normalized_product(clean(feature.get("product", ""))) == "hypothetical protein"


def feature_is_partial(feature: Mapping[str, Any]) -> bool:
    return bool(clean(feature.get("truncated", ""))) or bool(feature.get("edge", False))


def feature_aa(feature: Mapping[str, Any]) -> str:
    value = clean(feature.get("aa", "")).replace(" ", "").upper()
    if value.endswith("*"):
        value = value[:-1]
    return value


def feature_dbxrefs(feature: Mapping[str, Any]) -> list[str]:
    raw = feature.get("db_xrefs", [])
    if isinstance(raw, str):
        raw_values: Iterable[Any] = re.split(r"\s*,\s*", raw)
    elif isinstance(raw, Sequence):
        raw_values = raw
    else:
        raw_values = []
    return sorted(unique_clean(raw_values))


def split_dbxrefs(dbxrefs: Sequence[str]) -> dict[str, list[str]]:
    buckets: dict[str, list[str]] = {
        "ec": [],
        "cog": [],
        "go": [],
        "kegg": [],
        "pfam": [],
        "uniref": [],
        "refseq": [],
    }
    for xref in dbxrefs:
        prefix, separator, value = xref.partition(":")
        prefix_cf = prefix.casefold()
        canonical = xref if separator else clean(value or prefix)
        if prefix_cf == "ec":
            buckets["ec"].append(canonical)
        elif prefix_cf == "cog":
            buckets["cog"].append(canonical)
        elif prefix_cf == "go":
            buckets["go"].append(canonical)
        elif prefix_cf in {"kegg", "ko", "kofam"}:
            buckets["kegg"].append(canonical)
        elif prefix_cf == "pfam":
            buckets["pfam"].append(canonical)
        elif prefix_cf == "uniref":
            buckets["uniref"].append(canonical)
        elif prefix_cf in {"refseq", "ncbiprotein"}:
            buckets["refseq"].append(canonical)
    return {key: sorted(unique_clean(values)) for key, values in buckets.items()}


def inference_from_feature(feature: Mapping[str, Any]) -> InferenceRecord:
    experts = feature.get("expert", [])
    if isinstance(experts, list) and experts:
        candidates = [hit for hit in experts if isinstance(hit, Mapping)]
        if candidates:
            def expert_key(hit: Mapping[str, Any]) -> tuple[float, float, int]:
                return (
                    safe_float(hit.get("rank")) or 0.0,
                    safe_float(hit.get("score")) or 0.0,
                    int(bool(hit.get("gene"))) + int(bool(hit.get("product"))),
                )

            top = sorted(candidates, key=expert_key, reverse=True)[0]
            accessions = unique_clean(top.get("db_x_refs", top.get("db_xrefs", [])))
            return InferenceRecord(
                method="expert",
                accession=accessions[0] if accessions else "",
                score=safe_float(top.get("score")),
                evalue=safe_float(top.get("evalue")),
                query_cov=safe_float(top.get("query_cov")),
                subject_cov=safe_float(top.get("subject_cov")),
                identity=safe_float(top.get("identity")),
            )

    ips = feature.get("ips")
    if isinstance(ips, Mapping):
        accession = clean(
            ips.get("uniref100_id", ips.get("uniref100", ips.get("uniref", ips.get("accession", ""))))
        )
        return InferenceRecord(
            method="ips", accession=accession, evalue=0.0,
            query_cov=1.0, subject_cov=1.0, identity=1.0,
        )

    ups = feature.get("ups")
    if isinstance(ups, Mapping):
        accession = clean(
            ups.get("uniparc_id", ups.get("uniparc", ups.get("accession", "")))
        )
        return InferenceRecord(
            method="ups", accession=accession, evalue=0.0,
            query_cov=1.0, subject_cov=1.0, identity=1.0,
        )

    for method in ("psc", "pscc"):
        hit = feature.get(method)
        if isinstance(hit, Mapping):
            accession = clean(
                hit.get(
                    "uniref90_id" if method == "psc" else "uniref50_id",
                    hit.get(
                        "uniref90" if method == "psc" else "uniref50",
                        hit.get("uniref", hit.get("accession", "")),
                    ),
                )
            )
            return InferenceRecord(
                method=method,
                accession=accession,
                score=safe_float(hit.get("score")),
                evalue=safe_float(hit.get("evalue")),
                query_cov=safe_float(hit.get("query_cov")),
                subject_cov=safe_float(hit.get("subject_cov")),
                identity=safe_float(hit.get("identity")),
            )
    return InferenceRecord()


def normalize_column_name(value: str) -> str:
    text = clean(value).lstrip("#").casefold()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")



def inference_record_rank(record: InferenceRecord) -> tuple[Any, ...]:
    """Ordena evidências repetidas de forma determinística.

    O Bakta pode registrar mais de uma linha de inferência para um mesmo locus.
    A escolha privilegia uma referência identificada, a classe de evidência,
    a quantidade de métricas disponíveis e, em seguida, métricas mais fortes.
    """
    method_priority = {
        "expert": 6,
        "ips": 5,
        "psc": 4,
        "pscc": 3,
        "ups": 2,
        "inference_tsv": 1,
        "": 0,
    }
    metrics = (
        record.score,
        record.evalue,
        record.query_cov,
        record.subject_cov,
        record.identity,
    )
    metrics_present = sum(value is not None for value in metrics)
    coverages = [
        value for value in (record.query_cov, record.subject_cov)
        if value is not None
    ]
    minimum_coverage = min(coverages) if coverages else -1.0
    identity = record.identity if record.identity is not None else -1.0
    score = record.score if record.score is not None else -1.0
    if record.evalue is None:
        evalue_strength = -1.0
    elif record.evalue <= 0:
        evalue_strength = 999.0
    else:
        evalue_strength = -math.log10(record.evalue)

    return (
        0 if is_empty(record.accession) else 1,
        method_priority.get(record.method, 0),
        metrics_present,
        identity,
        minimum_coverage,
        evalue_strength,
        score,
        record.accession,
    )


def choose_preferred_inference(
    existing: InferenceRecord,
    candidate: InferenceRecord,
) -> InferenceRecord:
    """Seleciona deterministicamente a melhor linha para loci repetidos."""
    return max((existing, candidate), key=inference_record_rank)


def read_inference_tsv(path: Path | None) -> dict[str, InferenceRecord]:
    if path is None:
        return {}
    header: list[str] | None = None
    rows: dict[str, InferenceRecord] = {}
    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("#"):
                candidate = line.lstrip("#")
                if "\t" in candidate and "Locus Tag" in candidate:
                    header = [normalize_column_name(v) for v in candidate.split("\t")]
                continue
            if header is None:
                fields = line.split("\t")
                normalized = [normalize_column_name(v) for v in fields]
                if "locus_tag" in normalized:
                    header = normalized
                    continue
                raise SelectionError(f"Cabeçalho não reconhecido em inference.tsv: {path}")
            values = line.split("\t")
            if len(values) < len(header):
                values.extend([""] * (len(header) - len(values)))
            row = dict(zip(header, values))
            locus = clean(row.get("locus_tag", ""))
            # Bakta usa "-", ".", NA e valores equivalentes em linhas sem locus.
            # Esses valores não são identificadores e precisam ser ignorados.
            if is_empty(locus):
                continue
            accession = clean(row.get("accession", ""))
            metrics_present = any(
                safe_float(row.get(key)) is not None
                for key in ("score", "evalue", "query_cov", "subject_cov", "id", "identity")
            )
            method = "inference_tsv" if accession != "-" or metrics_present else ""
            if accession.startswith("UniRef:UniRef100"):
                method = "ips"
            elif accession.startswith("UniParc:"):
                method = "ups"
            elif accession.startswith("UniRef:UniRef90"):
                method = "psc"
            elif accession.startswith("UniRef:UniRef50"):
                method = "pscc"
            elif accession and accession != "-":
                method = "expert"
            record = InferenceRecord(
                method=method,
                accession="" if is_empty(accession) else accession,
                score=safe_float(row.get("score")),
                evalue=safe_float(row.get("evalue")),
                query_cov=safe_float(row.get("query_cov")),
                subject_cov=safe_float(row.get("subject_cov")),
                identity=safe_float(row.get("id", row.get("identity"))),
            )
            if locus in rows:
                rows[locus] = choose_preferred_inference(rows[locus], record)
            else:
                rows[locus] = record
    return rows


def choose_inference(
    locus: str,
    feature: Mapping[str, Any],
    table: Mapping[str, InferenceRecord],
) -> InferenceRecord:
    return table.get(locus, inference_from_feature(feature))


def classify_target(
    product: str,
    inference: InferenceRecord,
    params: SelectionParameters,
    bakta_hypothetical: bool = False,
) -> tuple[str, list[str], str] | None:
    flags: list[str] = []
    pclass = product_class(product)
    if bakta_hypothetical:
        flags.append("bakta_hypothetical")
        if pclass == "missing":
            flags.append("product_missing")
        else:
            flags.append("product_is_hypothetical")
        return "high", flags, "bakta_hypothetical_feature"
    if pclass == "missing":
        flags.append("product_missing")
        return "high", flags, "missing_product"
    if pclass == "hypothetical":
        flags.append("product_is_hypothetical")
        return "high", flags, "hypothetical_or_uncharacterized_product"
    if pclass == "generic":
        flags.append("generic_product")
        return "medium", flags, "generic_or_domain_level_product"
    if pclass == "uncertain":
        flags.append("uncertain_product_wording")
        return "low", flags, "putative_or_predicted_product"

    weak_metrics: list[str] = []
    if inference.identity is not None and inference.identity < params.weak_identity:
        weak_metrics.append("identity")
        flags.append("weak_identity")
    if inference.query_cov is not None and inference.query_cov < params.weak_query_cov:
        weak_metrics.append("query_cov")
        flags.append("weak_query_coverage")
    if (
        inference.subject_cov is not None
        and inference.subject_cov < params.weak_subject_cov
    ):
        weak_metrics.append("subject_cov")
        flags.append("weak_subject_coverage")
    if weak_metrics:
        return "low", flags, "weak_inference:" + ",".join(weak_metrics)
    return None


def output_paths(outdir: Path, sample: str) -> dict[str, Path]:
    return {
        "targets": outdir / f"{sample}.targets_refined.tsv",
        "audit": outdir / f"{sample}.selection_audit.tsv",
        "high": outdir / f"{sample}.priority_high.faa",
        "medium": outdir / f"{sample}.priority_medium.faa",
        "low": outdir / f"{sample}.priority_low.faa",
        "sorf": outdir / f"{sample}.sorf_separate.faa",
        "pseudogene": outdir / f"{sample}.pseudogene_separate.faa",
        "report": outdir / f"{sample}.target_report.tsv",
        "manifest": outdir / f"{sample}.target_manifest.json",
    }


def parameters_dict(params: SelectionParameters) -> dict[str, Any]:
    return {
        "short_aa_threshold": params.short_aa_threshold,
        "weak_identity": params.weak_identity,
        "weak_query_cov": params.weak_query_cov,
        "weak_subject_cov": params.weak_subject_cov,
        "include_sorf_targets": params.include_sorf_targets,
        "include_pseudogenes": params.include_pseudogenes,
        "allow_id_mismatch": params.allow_id_mismatch,
        "allow_sequence_mismatch": params.allow_sequence_mismatch,
        "allow_hypothetical_mismatch": params.allow_hypothetical_mismatch,
    }


def completion_is_valid(
    inputs: SampleInputs, params: SelectionParameters
) -> tuple[bool, str]:
    paths = output_paths(inputs.outdir, inputs.sample)
    if not paths["manifest"].is_file():
        return False, "manifest_missing"
    try:
        manifest = read_json_object(paths["manifest"])
    except SelectionError as exc:
        return False, str(exc)
    if (
        manifest.get("status") not in {"completed", "completed_with_warnings"}
        or manifest.get("pipeline") != PIPELINE_NAME
        or manifest.get("stage") != PIPELINE_STAGE
        or manifest.get("script_version") != SCRIPT_VERSION
    ):
        return False, "manifest_pipeline_stage_status_or_version_mismatch"
    if manifest.get("parameters") != parameters_dict(params):
        return False, "parameters_changed"

    saved_inputs = manifest.get("inputs", {})
    if not isinstance(saved_inputs, Mapping):
        return False, "manifest_inputs_invalid"
    required = {
        "bakta_json": inputs.bakta_json,
        "bakta_faa": inputs.bakta_faa,
    }
    for key, path in required.items():
        metadata = saved_inputs.get(key)
        if not isinstance(metadata, Mapping) or not fingerprint_matches(metadata, path):
            return False, f"input_fingerprint_mismatch:{key}"
    optional = {
        "inference_tsv": inputs.inference_tsv,
        "hypotheticals_faa": inputs.hypotheticals_faa,
        "hypotheticals_tsv": inputs.hypotheticals_tsv,
    }
    for key, path in optional.items():
        if not optional_fingerprint_matches(saved_inputs.get(key), path):
            return False, f"input_fingerprint_mismatch:{key}"

    saved_outputs = manifest.get("outputs", {})
    if not isinstance(saved_outputs, Mapping):
        return False, "manifest_outputs_invalid"
    for key in ("targets", "audit", "high", "medium", "low", "sorf", "pseudogene", "report"):
        metadata = saved_outputs.get(key)
        if not isinstance(metadata, Mapping) or not fingerprint_matches(metadata, paths[key]):
            return False, f"output_fingerprint_mismatch:{key}"
    return True, "completed_and_valid"


def read_json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SelectionError(f"JSON inválido: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SelectionError(f"JSON precisa conter um objeto: {path}")
    return value


def validate_output_targets(paths: Mapping[str, Path], force: bool) -> None:
    existing = [path for path in paths.values() if path.exists()]
    if existing and not force:
        joined = "\n  ".join(str(path) for path in existing)
        raise SelectionError(
            "Saídas existentes sem conclusão válida. Use --force após revisar:\n  " + joined
        )
    for path in existing:
        if not path.is_file():
            raise SelectionError(f"Destino esperado não é arquivo: {path}")


def discover_primary_json(directory: Path, sample: str, prefix: str | None) -> tuple[Path, str]:
    candidates: list[Path] = []
    if prefix:
        for suffix in (".json", ".json.gz"):
            candidate = directory / f"{prefix}{suffix}"
            if candidate.is_file():
                return candidate.resolve(), prefix
        raise SelectionError(f"JSON do prefixo {prefix!r} não encontrado em {directory}")

    for base in unique_clean([sample, directory.name]):
        for suffix in (".json", ".json.gz"):
            candidate = directory / f"{base}{suffix}"
            if candidate.is_file():
                return candidate.resolve(), base

    for pattern in ("*.json", "*.json.gz"):
        for path in sorted(directory.glob(pattern)):
            name = path.name.casefold()
            if not path.is_file():
                continue
            if name.endswith(".target_manifest.json") or "notification_email" in name:
                continue
            candidates.append(path)
    candidates = sorted(set(candidates))
    valid: list[tuple[Path, str]] = []
    for path in candidates:
        try:
            data = read_json(path)
        except SelectionError:
            continue
        if isinstance(data.get("features"), list):
            stem = path.name[:-8] if path.name.endswith(".json.gz") else path.stem
            valid.append((path.resolve(), stem))
    if not valid:
        raise SelectionError(f"Nenhum JSON principal do Bakta encontrado em {directory}")
    if len(valid) > 1:
        names = ", ".join(path.name for path, _ in valid)
        raise SelectionError(
            f"Mais de um JSON principal do Bakta em {directory}: {names}. Use --prefix."
        )
    return valid[0]


def find_optional(directory: Path, prefix: str, suffix: str) -> Path | None:
    path = directory / f"{prefix}{suffix}"
    return path.resolve() if path.is_file() else None


def make_sample_inputs(
    bakta_dir: Path,
    sample: str,
    outdir: Path,
    explicit_prefix: str | None,
) -> SampleInputs:
    resolved_dir = resolve_directory(bakta_dir, "Diretório Bakta")
    resolved_outdir = outdir.expanduser().resolve()
    if resolved_outdir == resolved_dir or path_is_within(resolved_outdir, resolved_dir):
        raise SelectionError(
            "O diretório de saída não pode ser o diretório Bakta nem ficar dentro dele"
        )
    json_path, prefix = discover_primary_json(resolved_dir, sample, explicit_prefix)
    faa = resolve_file(resolved_dir / f"{prefix}.faa", "FAA principal do Bakta")
    return SampleInputs(
        sample=sample,
        prefix=prefix,
        bakta_dir=resolved_dir,
        bakta_json=json_path,
        bakta_faa=faa,
        inference_tsv=find_optional(resolved_dir, prefix, ".inference.tsv"),
        hypotheticals_faa=find_optional(resolved_dir, prefix, ".hypotheticals.faa"),
        hypotheticals_tsv=find_optional(resolved_dir, prefix, ".hypotheticals.tsv"),
        outdir=resolved_outdir,
    )


def infer_sample(bakta_dir: Path, explicit: str | None = None) -> str:
    sample = clean(explicit) or bakta_dir.name
    if not sample or sample in {".", ".."} or Path(sample).name != sample:
        raise SelectionError(f"Nome de amostra inválido: {sample!r}")
    return sample


def load_selected_samples(args: argparse.Namespace) -> list[str]:
    selected: list[str] = []
    for value in args.sample or []:
        sample = clean(value)
        if sample and sample not in selected:
            selected.append(sample)
    if args.samples_file:
        path = resolve_file(Path(args.samples_file), "Lista de amostras")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            sample = line.strip()
            if sample and not sample.startswith("#") and sample not in selected:
                selected.append(sample)
    for sample in selected:
        if sample in {".", ".."} or Path(sample).name != sample:
            raise SelectionError(f"Nome de amostra inválido: {sample!r}")
    return selected


def directory_looks_like_bakta_output(path: Path) -> bool:
    if not path.is_dir() or path.name in QUEUE_DIRECTORY_NAMES:
        return False
    return any(path.glob("*.json")) or any(path.glob("*.json.gz"))


def resolve_samples(args: argparse.Namespace) -> list[SampleSpec]:
    selected = load_selected_samples(args)
    output_root = Path(args.output_root).expanduser().resolve()

    if args.bakta_dir:
        if len(selected) > 1:
            raise SelectionError("Modo --bakta-dir aceita no máximo uma --sample")
        bakta_dir = Path(args.bakta_dir).expanduser().resolve()
        sample = infer_sample(bakta_dir, selected[0] if selected else None)
        outdir = Path(args.outdir).expanduser().resolve() if args.outdir else output_root / sample
        return [SampleSpec(sample, bakta_dir, outdir, clean(args.prefix) or None)]

    if args.outdir:
        raise SelectionError("--outdir direto exige --bakta-dir")
    if args.prefix:
        raise SelectionError("--prefix direto exige --bakta-dir")

    bakta_root = resolve_directory(Path(args.bakta_root), "Raiz Bakta")
    if selected:
        samples = selected
    else:
        samples = sorted(
            path.name for path in bakta_root.iterdir() if directory_looks_like_bakta_output(path)
        )
    if not samples:
        raise SelectionError(f"Nenhuma amostra Bakta encontrada em {bakta_root}")
    return [SampleSpec(sample, bakta_root / sample, output_root / sample) for sample in samples]


def build_report_rows(sample: str, counts: Counter, flag_counts: Counter) -> list[dict[str, Any]]:
    descriptions = {
        "protein_features_total": "Features CDS e sORF enumeradas no JSON do Bakta",
        "faa_records": "Sequências presentes no FAA principal do Bakta",
        "cds_total": "Features do tipo CDS",
        "sorf_total": "Features do tipo sORF",
        "pseudogene_total": "Features marcadas como pseudogene",
        "targets_total": "Proteínas selecionadas em qualquer prioridade",
        "priority_high": "Produto ausente, hipotético ou não caracterizado",
        "priority_medium": "Produto genérico, DUF, UPF, família ou domínio",
        "priority_low": "Produto incerto ou sustentado por inferência fraca",
        "not_selected": "CDS informativas sem critério de reanotação",
        "sorf_separate": "sORFs preservadas em arquivo separado",
        "pseudogene_separate": "Pseudogenes preservados em arquivo separado",
        "json_without_faa": "Features proteicas do JSON ausentes no FAA",
        "faa_without_json": "Registros FAA ausentes entre CDS/sORF do JSON",
        "sequence_mismatch": "Sequências divergentes entre JSON e FAA",
        "hypothetical_json": "CDS com marcação hipotética reconhecida no JSON",
        "hypothetical_faa": "Registros oficiais no hypotheticals.faa do Bakta",
        "hypothetical_json_only": "Marcações hipotéticas apenas no JSON",
        "hypothetical_faa_only": "Registros oficiais apenas no hypotheticals.faa",
        "hypothetical_missing_main_faa": "Loci hipotéticos ausentes no FAA principal",
        "hypothetical_missing_json": "Loci hipotéticos ausentes nas features do JSON",
    }
    rows: list[dict[str, Any]] = []
    for name, description in descriptions.items():
        rows.append({
            "sample": sample,
            "category": "summary",
            "name": name,
            "count": counts.get(name, 0),
            "description": description,
        })
    for flag, count in sorted(flag_counts.items()):
        rows.append({
            "sample": sample,
            "category": "selection_flag",
            "name": flag,
            "count": count,
            "description": "Ocorrências desta flag entre as proteínas auditadas",
        })
    return rows


def validate_fasta_output(path: Path, expected: Sequence[FastaRecord]) -> None:
    observed = read_fasta(path) if expected else ([] if path.read_text(encoding="utf-8") == "" else read_fasta(path))
    if len(observed) != len(expected):
        raise SelectionError(f"Número de registros alterado no FASTA: {path}")
    for old, new in zip(expected, observed):
        if old.identifier != new.identifier or old.sequence != new.sequence:
            raise SelectionError(
                f"Identificador ou sequência alterado em {path.name}: {old.identifier}"
            )


def validate_disjoint_priorities(priority_records: Mapping[str, Sequence[FastaRecord]]) -> None:
    seen: set[str] = set()
    for label in ("high", "medium", "low"):
        for record in priority_records[label]:
            if record.identifier in seen:
                raise SelectionError(
                    f"Locus presente em mais de uma prioridade: {record.identifier}"
                )
            seen.add(record.identifier)


def bakta_version_info(data: Mapping[str, Any]) -> dict[str, Any]:
    version = data.get("version", {})
    if not isinstance(version, Mapping):
        return {"bakta": "", "db_version": "", "db_type": ""}
    db = version.get("db", {})
    if not isinstance(db, Mapping):
        db = {}
    return {
        "bakta": clean(version.get("bakta", "")),
        "db_version": clean(db.get("version", "")),
        "db_type": clean(db.get("type", "")),
    }


def process_sample(
    inputs: SampleInputs,
    params: SelectionParameters,
    force: bool,
) -> dict[str, Any]:
    valid, reason = completion_is_valid(inputs, params)
    if valid and not force:
        return {"status": "skipped", "reason": reason, "sample": inputs.sample}

    paths = output_paths(inputs.outdir, inputs.sample)
    validate_output_targets(paths, force)

    source_before = {
        "bakta_json": file_metadata(inputs.bakta_json),
        "bakta_faa": file_metadata(inputs.bakta_faa),
        "inference_tsv": optional_file_metadata(inputs.inference_tsv),
        "hypotheticals_faa": optional_file_metadata(inputs.hypotheticals_faa),
        "hypotheticals_tsv": optional_file_metadata(inputs.hypotheticals_tsv),
    }

    data = read_json(inputs.bakta_json)
    features_raw = data.get("features", [])
    features = [feature for feature in features_raw if isinstance(feature, Mapping)]
    fasta_records = read_fasta(inputs.bakta_faa)
    fasta_by_id = {record.identifier: record for record in fasta_records}
    inference_table = read_inference_tsv(inputs.inference_tsv)

    protein_features: dict[str, Mapping[str, Any]] = {}
    warnings: list[str] = []
    for index, feature in enumerate(features):
        ftype = feature_type(feature)
        if ftype not in {"cds", "sorf"}:
            continue
        locus = feature_locus(feature)
        if not locus:
            raise SelectionError(
                f"Feature proteica sem locus no JSON: índice={index}, tipo={ftype}"
            )
        if locus in protein_features:
            raise SelectionError(f"Locus duplicado no JSON do Bakta: {locus}")
        protein_features[locus] = feature

    json_loci = set(protein_features)
    faa_loci = set(fasta_by_id)
    json_without_faa = sorted(json_loci - faa_loci)
    faa_without_json = sorted(faa_loci - json_loci)
    if (json_without_faa or faa_without_json) and not params.allow_id_mismatch:
        raise SelectionError(
            "Inconsistência de loci entre JSON e FAA: "
            f"json_sem_faa={len(json_without_faa)} "
            f"faa_sem_json={len(faa_without_json)}. "
            "Use --allow-id-mismatch somente após revisar."
        )
    if json_without_faa:
        warnings.append(f"json_without_faa={len(json_without_faa)}")
    if faa_without_json:
        warnings.append(f"faa_without_json={len(faa_without_json)}")

    sequence_mismatches: list[str] = []
    for locus in sorted(json_loci & faa_loci):
        json_aa = feature_aa(protein_features[locus])
        if json_aa and json_aa != fasta_by_id[locus].sequence:
            sequence_mismatches.append(locus)
    if sequence_mismatches and not params.allow_sequence_mismatch:
        preview = ", ".join(sequence_mismatches[:10])
        raise SelectionError(
            f"Sequências divergentes entre JSON e FAA: {len(sequence_mismatches)} "
            f"({preview}). Use --allow-sequence-mismatch somente após revisar."
        )
    if sequence_mismatches:
        warnings.append(f"sequence_mismatch={len(sequence_mismatches)}")

    hypothetical_records = (
        read_fasta(inputs.hypotheticals_faa)
        if inputs.hypotheticals_faa
        else []
    )
    hypothetical_faa_loci = {record.identifier for record in hypothetical_records}

    # O arquivo <prefix>.hypotheticals.faa é a fonte oficial do Bakta para
    # proteínas hipotéticas. O JSON contém metadados internos, mas sua
    # representação pode não reproduzir exatamente esse subconjunto.
    #
    # O contrato realmente obrigatório é:
    #   hypotheticals.faa ⊆ FAA principal ∩ features proteicas do JSON
    hypothetical_missing_main_faa = sorted(hypothetical_faa_loci - faa_loci)
    hypothetical_missing_json = sorted(hypothetical_faa_loci - json_loci)
    if hypothetical_missing_main_faa or hypothetical_missing_json:
        message = (
            "hypothetical_contract_mismatch:"
            f"missing_main_faa={len(hypothetical_missing_main_faa)},"
            f"missing_json={len(hypothetical_missing_json)}"
        )
        if not params.allow_hypothetical_mismatch:
            raise SelectionError(
                "Loci do hypotheticals.faa ausentes no FAA principal ou no JSON: "
                f"missing_main_faa={len(hypothetical_missing_main_faa)} "
                f"missing_json={len(hypothetical_missing_json)}. "
                "Use --allow-hypothetical-mismatch somente após revisar."
            )
        warnings.append(message)

    hypothetical_json_loci = {
        locus
        for locus, feature in protein_features.items()
        if (
            feature_type(feature) == "cds"
            and locus in faa_loci
            and feature_is_hypothetical(feature)
        )
    }

    # Diferenças sem perda de locus são registradas para auditoria, mas não
    # interrompem o processamento.
    hypo_only_json = sorted(hypothetical_json_loci - hypothetical_faa_loci)
    hypo_only_faa = sorted(hypothetical_faa_loci - hypothetical_json_loci)
    if inputs.hypotheticals_faa and (hypo_only_json or hypo_only_faa):
        warnings.append(
            f"hypothetical_representation_difference:"
            f"json_only={len(hypo_only_json)},"
            f"faa_only={len(hypo_only_faa)}"
        )

    counts: Counter = Counter()
    counts["protein_features_total"] = len(protein_features)
    counts["faa_records"] = len(fasta_records)
    counts["json_without_faa"] = len(json_without_faa)
    counts["faa_without_json"] = len(faa_without_json)
    counts["sequence_mismatch"] = len(sequence_mismatches)
    counts["hypothetical_json"] = len(hypothetical_json_loci)
    counts["hypothetical_faa"] = len(hypothetical_faa_loci)
    counts["hypothetical_json_only"] = len(hypo_only_json)
    counts["hypothetical_faa_only"] = len(hypo_only_faa)
    counts["hypothetical_missing_main_faa"] = len(hypothetical_missing_main_faa)
    counts["hypothetical_missing_json"] = len(hypothetical_missing_json)

    target_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    flag_counts: Counter = Counter()
    priority_by_locus: dict[str, str] = {}
    flags_by_locus: dict[str, str] = {}
    sorf_loci: set[str] = set()
    pseudogene_loci: set[str] = set()

    for locus, feature in protein_features.items():
        ftype = feature_type(feature)
        if ftype == "cds":
            counts["cds_total"] += 1
        elif ftype == "sorf":
            counts["sorf_total"] += 1

        record = fasta_by_id.get(locus)
        if record is None:
            audit_rows.append({
                "sample": inputs.sample,
                "locus_tag": locus,
                "feature_type": ftype,
                "aa_length": len(feature_aa(feature)),
                "product": clean(feature.get("product", "")),
                "gene": clean(feature.get("gene", "")),
                "decision": "excluded",
                "priority": "",
                "reason": "json_without_faa",
                "flags": "json_without_faa",
            })
            flag_counts["json_without_faa"] += 1
            continue

        product = clean(feature.get("product", ""))
        gene = clean(feature.get("gene", ""))
        genes_raw = feature.get("genes", [])
        genes = unique_clean(genes_raw if isinstance(genes_raw, Sequence) and not isinstance(genes_raw, str) else [genes_raw])
        dbxrefs = feature_dbxrefs(feature)
        dbgroups = split_dbxrefs(dbxrefs)
        inference = choose_inference(locus, feature, inference_table)
        aa_length = len(record.sequence)

        base_flags: list[str] = []
        if aa_length < params.short_aa_threshold:
            base_flags.append("short_protein")
        truncated = clean(feature.get("truncated", ""))
        if truncated:
            base_flags.append(f"truncated_{truncated}")
        if feature.get("edge", False):
            base_flags.append("edge_feature")
        if not gene:
            base_flags.append("gene_missing")
        if not dbgroups["ec"]:
            base_flags.append("ec_missing")
        if not dbgroups["cog"]:
            base_flags.append("cog_missing")
        if not dbgroups["go"]:
            base_flags.append("go_missing")
        if not inference.method:
            base_flags.append("inference_missing")
        if locus in hypothetical_faa_loci:
            base_flags.append("bakta_hypotheticals_faa")
        if feature_is_hypothetical(feature):
            base_flags.append("bakta_hypothetical_json")

        is_pseudo = feature_is_pseudogene(feature)
        if is_pseudo:
            counts["pseudogene_total"] += 1
            base_flags.append("pseudogene")
            if not params.include_pseudogenes:
                pseudogene_loci.add(locus)
                counts["pseudogene_separate"] += 1
                flags_text = ";".join(sorted(set(base_flags)))
                flag_counts.update(set(base_flags))
                flags_by_locus[locus] = flags_text
                audit_rows.append({
                    "sample": inputs.sample,
                    "locus_tag": locus,
                    "feature_type": ftype,
                    "aa_length": aa_length,
                    "product": product,
                    "gene": gene,
                    "decision": "separate",
                    "priority": "pseudogene",
                    "reason": "pseudogene_excluded_from_standard_search",
                    "flags": flags_text,
                })
                continue

        if ftype == "sorf" and not params.include_sorf_targets:
            sorf_loci.add(locus)
            counts["sorf_separate"] += 1
            base_flags.append("sorf")
            flags_text = ";".join(sorted(set(base_flags)))
            flag_counts.update(set(base_flags))
            flags_by_locus[locus] = flags_text
            audit_rows.append({
                "sample": inputs.sample,
                "locus_tag": locus,
                "feature_type": ftype,
                "aa_length": aa_length,
                "product": product,
                "gene": gene,
                "decision": "separate",
                "priority": "sorf",
                "reason": "sorf_excluded_from_standard_search",
                "flags": flags_text,
            })
            continue

        classification = classify_target(
            product,
            inference,
            params,
            bakta_hypothetical=(
                locus in hypothetical_faa_loci
                or feature_is_hypothetical(feature)
            ),
        )
        if classification is None:
            counts["not_selected"] += 1
            flags_text = ";".join(sorted(set(base_flags)))
            flag_counts.update(set(base_flags))
            audit_rows.append({
                "sample": inputs.sample,
                "locus_tag": locus,
                "feature_type": ftype,
                "aa_length": aa_length,
                "product": product,
                "gene": gene,
                "decision": "not_selected",
                "priority": "",
                "reason": "informative_product_with_adequate_evidence",
                "flags": flags_text,
            })
            continue

        priority, class_flags, rationale = classification
        all_flags = sorted(set(base_flags + class_flags))
        flags_text = ";".join(all_flags)
        priority_by_locus[locus] = priority
        flags_by_locus[locus] = flags_text
        counts["targets_total"] += 1
        counts[f"priority_{priority}"] += 1
        flag_counts.update(all_flags)

        target_rows.append({
            "sample": inputs.sample,
            "source_annotation": "bakta",
            "feature_id": clean(feature.get("id", locus)),
            "locus_tag": locus,
            "feature_type": ftype,
            "contig": feature_sequence_id(feature),
            "start": feature.get("start", ""),
            "stop": feature.get("stop", ""),
            "strand": clean(feature.get("strand", "")),
            "aa_length": aa_length,
            "original_product": product,
            "original_gene": gene,
            "original_genes": ";".join(genes),
            "original_ec": ";".join(dbgroups["ec"]),
            "original_cog": ";".join(dbgroups["cog"]),
            "original_go": ";".join(dbgroups["go"]),
            "original_kegg": ";".join(dbgroups["kegg"]),
            "original_pfam": ";".join(dbgroups["pfam"]),
            "original_uniref": ";".join(dbgroups["uniref"]),
            "original_refseq": ";".join(dbgroups["refseq"]),
            "original_dbxrefs": ";".join(dbxrefs),
            "inference_method": inference.method,
            "inference_accession": inference.accession,
            "inference_score": format_float(inference.score),
            "inference_evalue": format_float(inference.evalue),
            "inference_query_cov": format_float(inference.query_cov),
            "inference_subject_cov": format_float(inference.subject_cov),
            "inference_identity": format_float(inference.identity),
            "truncated": truncated,
            "edge": bool_text(feature.get("edge", False)),
            "selection_flags": flags_text,
            "target_priority_refined": priority,
            "priority_rationale": rationale,
        })
        audit_rows.append({
            "sample": inputs.sample,
            "locus_tag": locus,
            "feature_type": ftype,
            "aa_length": aa_length,
            "product": product,
            "gene": gene,
            "decision": "selected",
            "priority": priority,
            "reason": rationale,
            "flags": flags_text,
        })

    target_rows.sort(
        key=lambda row: (
            PRIORITY_ORDER.get(clean(row["target_priority_refined"]), 9),
            clean(row["locus_tag"]),
        )
    )
    audit_rows.sort(key=lambda row: (clean(row["locus_tag"]), clean(row["decision"])))

    priority_records: dict[str, list[FastaRecord]] = {
        label: [record for record in fasta_records if priority_by_locus.get(record.identifier) == label]
        for label in ("high", "medium", "low")
    }
    sorf_records = [record for record in fasta_records if record.identifier in sorf_loci]
    pseudogene_records = [record for record in fasta_records if record.identifier in pseudogene_loci]
    validate_disjoint_priorities(priority_records)

    inputs.outdir.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(paths["targets"], target_rows, TARGET_COLUMNS)
    atomic_write_tsv(paths["audit"], audit_rows, AUDIT_COLUMNS)
    for label in ("high", "medium", "low"):
        write_selected_fasta(paths[label], priority_records[label], label, flags_by_locus)
    write_selected_fasta(paths["sorf"], sorf_records, "sorf_separate", flags_by_locus)
    write_selected_fasta(paths["pseudogene"], pseudogene_records, "pseudogene_separate", flags_by_locus)
    report_rows = build_report_rows(inputs.sample, counts, flag_counts)
    atomic_write_tsv(paths["report"], report_rows, REPORT_COLUMNS)

    for label in ("high", "medium", "low"):
        validate_fasta_output(paths[label], priority_records[label])
    validate_fasta_output(paths["sorf"], sorf_records)
    validate_fasta_output(paths["pseudogene"], pseudogene_records)

    source_after = {
        "bakta_json": file_metadata(inputs.bakta_json),
        "bakta_faa": file_metadata(inputs.bakta_faa),
        "inference_tsv": optional_file_metadata(inputs.inference_tsv),
        "hypotheticals_faa": optional_file_metadata(inputs.hypotheticals_faa),
        "hypotheticals_tsv": optional_file_metadata(inputs.hypotheticals_tsv),
    }
    for key in source_before:
        before = source_before[key]
        after = source_after[key]
        if before is None and after is None:
            continue
        if not isinstance(before, Mapping) or not isinstance(after, Mapping):
            raise SelectionError(f"Arquivo opcional apareceu ou desapareceu durante a execução: {key}")
        if before.get("sha256") != after.get("sha256"):
            raise SelectionError(f"Arquivo original alterado durante a execução: {key}")

    version_info = bakta_version_info(data)
    manifest = {
        "status": "completed_with_warnings" if warnings else "completed",
        "pipeline": PIPELINE_NAME,
        "stage": PIPELINE_STAGE,
        "script": SCRIPT_NAME,
        "script_version": SCRIPT_VERSION,
        "finished_utc": utc_now(),
        "sample": inputs.sample,
        "source_annotation": "bakta",
        "bakta": version_info,
        "inputs": source_after,
        "parameters": parameters_dict(params),
        "counts": dict(sorted(counts.items())),
        "selection_flags": dict(sorted(flag_counts.items())),
        "warnings": warnings,
        "outputs": {
            key: file_metadata(path)
            for key, path in paths.items()
            if key != "manifest"
        } | {"manifest": {"path": str(paths["manifest"])}},
        "next_script": NEXT_SCRIPT,
        "next_input_contract": str(paths["high"]),
    }
    atomic_write_json(paths["manifest"], manifest)
    return {
        "status": manifest["status"],
        "sample": inputs.sample,
        "counts": dict(counts),
        "warnings": warnings,
        "paths": {key: str(path) for key, path in paths.items()},
    }


def validate_fraction(value: float, option: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise SelectionError(f"{option} deve estar entre 0 e 1")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=SCRIPT_NAME,
        description=(
            "Etapa anot00/02: seleciona e prioriza proteínas das saídas do Bakta "
            "para reanotação em bancos funcionais."
        ),
        epilog=f"Próxima etapa: {NEXT_SCRIPT}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {SCRIPT_VERSION}",
        help="Exibe o nome canônico e a versão do script",
    )
    parser.add_argument(
        "--bakta-root",
        default=str(DEFAULT_BAKTA_ROOT),
        help="Raiz com um diretório de saída Bakta por amostra",
    )
    parser.add_argument(
        "--output-root",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Raiz de saída consumida pela etapa Swiss-Prot",
    )
    parser.add_argument("--bakta-dir", help="Modo direto: diretório Bakta de uma amostra")
    parser.add_argument("--outdir", help="Modo direto: diretório de saída da amostra")
    parser.add_argument(
        "--prefix",
        help="Prefixo dos arquivos Bakta no modo direto, quando não puder ser inferido",
    )
    parser.add_argument("--sample", action="append", help="Amostra; opção repetível no lote")
    parser.add_argument("--samples-file", help="Arquivo com uma amostra por linha")
    parser.add_argument(
        "--short-aa-threshold",
        type=int,
        default=90,
        help="Comprimento abaixo do qual a flag short_protein é registrada",
    )
    parser.add_argument(
        "--weak-identity",
        type=float,
        default=0.90,
        help="Identidade abaixo da qual uma anotação informativa vai para prioridade low",
    )
    parser.add_argument(
        "--weak-query-cov",
        type=float,
        default=0.80,
        help="Cobertura da query abaixo da qual a anotação vai para prioridade low",
    )
    parser.add_argument(
        "--weak-subject-cov",
        type=float,
        default=0.80,
        help="Cobertura do alvo abaixo da qual a anotação vai para prioridade low",
    )
    parser.add_argument(
        "--include-sorf-targets",
        action="store_true",
        help="Permite classificar sORFs junto às proteínas convencionais",
    )
    parser.add_argument(
        "--include-pseudogenes",
        action="store_true",
        help="Permite classificar pseudogenes junto às proteínas convencionais",
    )
    parser.add_argument(
        "--allow-id-mismatch",
        action="store_true",
        help="Processa a interseção JSON/FAA quando há loci incompatíveis",
    )
    parser.add_argument(
        "--allow-sequence-mismatch",
        action="store_true",
        help="Aceita divergência entre a sequência do JSON e do FAA",
    )
    parser.add_argument(
        "--allow-hypothetical-mismatch",
        action="store_true",
        help="Aceita loci do hypotheticals.faa ausentes no FAA principal ou no JSON após revisão",
    )
    parser.add_argument("--force", action="store_true", help="Substitui somente saídas conhecidas")
    parser.add_argument("--fail-fast", action="store_true", help="Interrompe o lote no primeiro erro")
    parser.add_argument("--dry-run", action="store_true", help="Exibe amostras sem criar saídas")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.short_aa_threshold < 1:
            raise SelectionError("--short-aa-threshold deve ser maior que zero")
        validate_fraction(args.weak_identity, "--weak-identity")
        validate_fraction(args.weak_query_cov, "--weak-query-cov")
        validate_fraction(args.weak_subject_cov, "--weak-subject-cov")
        sample_specs = resolve_samples(args)
    except (SelectionError, OSError) as exc:
        print(f"[ERRO] {exc}", file=sys.stderr)
        return 2

    params = SelectionParameters(
        short_aa_threshold=args.short_aa_threshold,
        weak_identity=args.weak_identity,
        weak_query_cov=args.weak_query_cov,
        weak_subject_cov=args.weak_subject_cov,
        include_sorf_targets=bool(args.include_sorf_targets),
        include_pseudogenes=bool(args.include_pseudogenes),
        allow_id_mismatch=bool(args.allow_id_mismatch),
        allow_sequence_mismatch=bool(args.allow_sequence_mismatch),
        allow_hypothetical_mismatch=bool(args.allow_hypothetical_mismatch),
    )

    if args.dry_run:
        failed = 0
        for spec in sample_specs:
            try:
                inputs = make_sample_inputs(
                    spec.bakta_dir, spec.sample, spec.outdir, spec.explicit_prefix
                )
                optional = [
                    path.name
                    for path in (
                        inputs.inference_tsv,
                        inputs.hypotheticals_faa,
                        inputs.hypotheticals_tsv,
                    )
                    if path is not None
                ]
                print(
                    f"[DRY-RUN] {inputs.sample}: {inputs.bakta_json.name} + "
                    f"{inputs.bakta_faa.name} -> {inputs.outdir}; "
                    f"auxiliares={','.join(optional) or '-'}"
                )
            except (SelectionError, OSError) as exc:
                failed += 1
                print(f"[FAIL] {spec.sample}: {exc}", file=sys.stderr)
                if args.fail_fast:
                    break
        print(f"samples={len(sample_specs)} failed={failed}")
        return 0 if failed == 0 else 1

    completed = 0
    skipped = 0
    failed = 0
    for spec in sample_specs:
        try:
            inputs = make_sample_inputs(
                spec.bakta_dir, spec.sample, spec.outdir, spec.explicit_prefix
            )
            result = process_sample(inputs=inputs, params=params, force=bool(args.force))
            if result["status"] == "skipped":
                skipped += 1
                print(f"[SKIP] {inputs.sample}: {result['reason']}")
            else:
                completed += 1
                counts = result["counts"]
                warning_text = (
                    f" warnings={len(result.get('warnings', []))}"
                    if result.get("warnings") else ""
                )
                print(
                    f"[OK] {inputs.sample}: high={counts.get('priority_high', 0)} "
                    f"medium={counts.get('priority_medium', 0)} "
                    f"low={counts.get('priority_low', 0)} "
                    f"sorf={counts.get('sorf_separate', 0)} "
                    f"pseudogene={counts.get('pseudogene_separate', 0)}"
                    f"{warning_text}"
                )
        except (SelectionError, OSError) as exc:
            failed += 1
            print(f"[FAIL] {spec.sample}: {exc}", file=sys.stderr)
            if args.fail_fast:
                break
    print(f"completed={completed} skipped={skipped} failed={failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
