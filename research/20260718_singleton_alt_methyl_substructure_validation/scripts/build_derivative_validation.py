#!/usr/bin/env python3
"""Build the read-only positional-singleton methyl-substructure sidecar report."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import re
import sqlite3
import stat
import sys
import traceback
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = ["DejaVu Sans"]


TOPIC_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AUDIT_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "positional_singleton_methyl_multigroup_audit_v1_source_attested"
)
DEFAULT_V10_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
)
DEFAULT_CLAIM_CONTRACT = (
    TOPIC_ROOT.parent
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "claim-contract-v5.md"
)

EXPECTED_SHA256 = {
    "audit_summary": "6d4128be0535f2e16b1382cf038c558054ac42d0d7bde75ab4854b7a5be7aedc",
    "site_audit": "2d5d24790918fca34a32c313b9965f1de8c186c031e31d1a643febba721d90ce",
    "m2_pass_cases": "2fe42d55924d4d73f8a0d7b436cc8c517e9d06834da9fd7a1e882ce24c4f1dd0",
    "v10_site_results": "a8871af3a8c3955bf31aec5eeef0c93aca0683f52cf6d6f1e06fbbb713324f74",
}
EXPECTED_COUNTS = {
    "all_dataset_sites": 469_849,
    "singleton_sites": 50_432,
    "m1_evaluable": 48_347,
    "m1_flagged": 5_961,
    "m2_not_run": 44_471,
    "m2_not_evaluable": 5_913,
    "m2_fail": 18,
    "m2_pass": 30,
}
EXPECTED_HCC1395 = {
    "sites": 8_279,
    "m1_evaluable": 8_074,
    "m1_flagged": 734,
    "m2_not_run": 7_545,
    "m2_not_evaluable": 732,
    "m2_fail": 0,
    "m2_pass": 2,
}
M2_STATUSES = ("NOT_RUN", "NOT_EVALUABLE", "FAIL", "PASS")
TRUTH_LABELS = ("TP", "FP", "UNASSESSED")
DATASET_ORDER = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
REPORT_TITLE = "Positional-singleton focal-ALT 甲基子結構 derivative validation"
SOURCE_CLUSTER_LABEL_RE = re.compile(r"(?<![A-Za-z0-9])1-[12](?![0-9])")


@dataclass(frozen=True)
class Target:
    dataset: str
    chrom: str
    pos: int
    ref: str
    alt: str
    truth: str
    raw_alt_reads: int
    core_reads: int
    group_sizes: tuple[int, int]

    @property
    def key(self) -> tuple[str, str, int, str, str]:
        return self.dataset, self.chrom, self.pos, self.ref, self.alt

    @property
    def short_id(self) -> str:
        return f"{self.chrom}_{self.pos}_{self.ref}_{self.alt}"

    @property
    def display_locus(self) -> str:
        return f"{self.chrom}:{self.pos:,} {self.ref}>{self.alt}"


TARGETS = (
    Target("HCC1395", "chr14", 86_272_476, "A", "T", "TP", 123, 108, (86, 22)),
    Target("HCC1395", "chr22", 47_466_517, "A", "G", "TP", 114, 109, (88, 21)),
)


class ValidationError(RuntimeError):
    """Raised when a frozen source or scientific invariant does not validate."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise ValidationError(f"{label}: expected {expected!r}, observed {actual!r}")


def assert_true(value: bool, label: str) -> None:
    if not value:
        raise ValidationError(label)


def parse_bool_series(series: pd.Series, label: str) -> pd.Series:
    normalized = series.astype(str).str.strip().str.lower()
    allowed = {"true", "false"}
    unexpected = sorted(set(normalized) - allowed)
    if unexpected:
        raise ValidationError(f"{label}: unexpected boolean values {unexpected[:10]}")
    return normalized.eq("true")


def parse_nullable_bool_series(series: pd.Series, label: str) -> pd.Series:
    normalized = series.astype("string").str.strip().str.lower()
    unexpected = sorted(set(normalized.dropna()) - {"true", "false"})
    if unexpected:
        raise ValidationError(f"{label}: unexpected boolean values {unexpected[:10]}")
    return normalized.map({"true": True, "false": False}).astype("boolean")


def contains_source_cluster_label(value: str) -> bool:
    return SOURCE_CLUSTER_LABEL_RE.search(value) is not None


def quote_sql_identifier(value: str) -> str:
    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", value) is None:
        raise ValidationError(f"unsafe SQLite identifier {value!r}")
    return f'"{value}"'


def normalize_sql_scalar(value: Any, template: Any) -> Any:
    if template is None:
        if not pd.isna(value):
            raise ValidationError(f"SQLite round-trip expected null, observed {value!r}")
        return None
    if isinstance(template, (bool, np.bool_)):
        return bool(value)
    if isinstance(template, (int, np.integer)):
        return int(value)
    if isinstance(template, (float, np.floating)):
        return float(value)
    if isinstance(template, str):
        return str(value)
    raise ValidationError(f"unsupported snapshot scalar type {type(template).__name__}")


def roundtrip_snapshot_dataset(
    connection: sqlite3.Connection,
    dataset_id: str,
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], str]:
    assert_true(bool(rows), f"snapshot dataset {dataset_id} must not be empty")
    table_name = f"artifact_{dataset_id}"
    quote_sql_identifier(table_name)
    columns = list(rows[0])
    assert_true(bool(columns), f"snapshot dataset {dataset_id} requires columns")
    assert_true(all(list(row) == columns for row in rows), f"snapshot dataset {dataset_id} column drift")
    for column in columns:
        quote_sql_identifier(column)

    frame = pd.DataFrame([dict(row) for row in rows], columns=columns)
    frame.insert(0, "_row_order", np.arange(len(frame), dtype=np.int64))
    frame.to_sql(table_name, connection, if_exists="fail", index=False)
    selected_columns = ", ".join(quote_sql_identifier(column) for column in columns)
    sql = (
        f"SELECT {selected_columns} FROM {quote_sql_identifier(table_name)} "
        'ORDER BY "_row_order" ASC'
    )
    observed = pd.read_sql_query(sql, connection)
    assert_equal(list(observed.columns), columns, f"SQLite round-trip columns for {dataset_id}")
    assert_equal(len(observed), len(rows), f"SQLite round-trip rows for {dataset_id}")

    normalized: list[dict[str, Any]] = []
    for row_index, template_row in enumerate(rows):
        normalized.append(
            {
                column: normalize_sql_scalar(observed.iloc[row_index][column], template_row[column])
                for column in columns
            }
        )
    assert_equal(normalized, [dict(row) for row in rows], f"SQLite round-trip values for {dataset_id}")
    return normalized, sql


def stage_snapshot_datasets(
    path: Path,
    datasets: Mapping[str, Sequence[Mapping[str, Any]]],
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, str]]:
    if path.exists():
        raise FileExistsError(path)
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    os.close(descriptor)
    connection = sqlite3.connect(path)
    staged: dict[str, list[dict[str, Any]]] = {}
    queries: dict[str, str] = {}
    try:
        connection.execute(
            "CREATE TABLE artifact_metadata "
            "(schema_name TEXT NOT NULL, schema_version TEXT NOT NULL, dataset_count INTEGER NOT NULL)"
        )
        connection.execute(
            "INSERT INTO artifact_metadata VALUES (?, ?, ?)",
            ("intersubmod.singleton_alt_methyl_artifact_staging", "1.0.0", len(datasets)),
        )
        for dataset_id, rows in datasets.items():
            staged[dataset_id], queries[dataset_id] = roundtrip_snapshot_dataset(
                connection,
                dataset_id,
                rows,
            )
        connection.commit()
    finally:
        connection.close()
    return staged, queries


def attach_executed_sql_sources(
    items: Sequence[dict[str, Any]],
    sources: Sequence[Mapping[str, Any]],
    queries: Mapping[str, str],
    staging_sha256: str,
) -> None:
    sources_by_id = {str(source["id"]): source for source in sources}
    for item in items:
        dataset_id = str(item["dataset"])
        upstream_id = str(item["sourceId"])
        upstream = sources_by_id[upstream_id]
        upstream_query = dict(upstream.get("query", {}))
        upstream_tables = [str(value) for value in upstream_query.get("tables_used", [])]
        item["source"] = {
            "id": f"artifact_sql_{dataset_id}",
            "label": f"{upstream['label']} via executed SQLite staging query",
            "path": "artifact_staging.sqlite",
            "sha256": staging_sha256,
            "upstreamSourceId": upstream_id,
            "upstreamSha256": upstream.get("sha256"),
            "query": {
                "sql": queries[dataset_id],
                "engine": f"SQLite {sqlite3.sqlite_version}",
                "language": "SQL",
                "description": (
                    "Actual query executed against the persisted artifact staging database; "
                    "the returned rows were exact-compared with the reviewed snapshot dataset."
                ),
                "tables_used": [f"artifact_{dataset_id}", *upstream_tables],
                "filters": list(upstream_query.get("filters", [])),
                "metric_definitions": list(upstream_query.get("metric_definitions", [])),
            },
        }


def safe_json_dump(path: Path, payload: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def safe_tsv_dump(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def safe_gzip_tsv_dump(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fieldnames: Sequence[str],
) -> None:
    if path.exists():
        raise FileExistsError(path)
    with gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=9) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def ensure_output_dir(output_dir: Path) -> None:
    resolved_topic = TOPIC_ROOT.resolve()
    resolved_output = output_dir.resolve()
    try:
        resolved_output.relative_to(resolved_topic)
    except ValueError as exc:
        raise ValidationError(f"output_dir must stay within {resolved_topic}") from exc
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir()


def readonly(path: Path) -> None:
    path.chmod(0o444)


def source_record(
    path: Path,
    role: str,
    used_for: str,
    expected_sha256: str | None = None,
) -> dict[str, Any]:
    assert_true(path.is_file(), f"missing source file: {path}")
    observed_sha = sha256_path(path)
    if expected_sha256 is not None:
        assert_equal(observed_sha, expected_sha256, f"{role} SHA-256")
    file_stat = path.stat()
    return {
        "role": role,
        "path": str(path.resolve()),
        "size_bytes": file_stat.st_size,
        "mode": oct(stat.S_IMODE(file_stat.st_mode)),
        "sha256": observed_sha,
        "used_for": used_for,
    }


def key_columns(frame: pd.DataFrame) -> pd.Series:
    return pd.Series(
        list(
            zip(
                frame["dataset"].astype(str),
                frame["chrom"].astype(str),
                frame["pos"].astype(int),
                frame["ref"].astype(str),
                frame["alt"].astype(str),
            )
        ),
        index=frame.index,
        dtype=object,
    )


def recompute_positional_singleton_keys(
    frame: pd.DataFrame,
    maximum_adjacent_gap_bp: int = 50_000,
) -> set[tuple[str, str, int, str, str]]:
    required = {"dataset", "chrom", "pos", "ref", "alt"}
    missing = required - set(frame.columns)
    if missing:
        raise ValidationError(f"positional recomputation missing columns: {sorted(missing)}")
    singleton_keys: set[tuple[str, str, int, str, str]] = set()
    for (dataset, chrom), group in frame.groupby(["dataset", "chrom"], sort=False, observed=True):
        ordered = group.sort_values(["pos", "ref", "alt"], kind="mergesort")
        records = list(ordered[["pos", "ref", "alt"]].itertuples(index=False, name=None))
        component_start = 0
        for index in range(1, len(records) + 1):
            component_break = index == len(records)
            if not component_break:
                component_break = int(records[index][0]) - int(records[index - 1][0]) > maximum_adjacent_gap_bp
            if component_break:
                if index - component_start == 1:
                    pos, ref, alt = records[component_start]
                    singleton_keys.add((str(dataset), str(chrom), int(pos), str(ref), str(alt)))
                component_start = index
    return singleton_keys


def load_v10_site_results(path: Path) -> tuple[pd.DataFrame, set[tuple[str, str, int, str, str]]]:
    columns = [
        "dataset",
        "sample",
        "biological_id",
        "truth_label",
        "chrom",
        "pos",
        "ref",
        "alt",
        "ssnv_branch",
        "component_size",
        "region_dir",
        "caller_af",
        "normal_af",
        "normal_dp",
        "n_alt_matrix",
        "n_alt_after_peel",
        "cluster_sizes",
        "stable_null_multigroup",
    ]
    frame = pd.read_csv(path, sep="\t", compression="gzip", usecols=columns, dtype=str)
    assert_equal(len(frame), EXPECTED_COUNTS["all_dataset_sites"], "v10 dataset-site row count")
    frame["pos"] = pd.to_numeric(frame["pos"], errors="raise").astype("int64")
    keys = key_columns(frame)
    assert_equal(int(keys.duplicated().sum()), 0, "v10 duplicate primary keys")
    recomputed = recompute_positional_singleton_keys(frame)
    assert_equal(len(recomputed), EXPECTED_COUNTS["singleton_sites"], "recomputed positional singletons")
    upstream = set(keys[frame["ssnv_branch"].eq("positional_singleton")])
    assert_equal(recomputed, upstream, "recomputed singleton key set vs v10 ssnv_branch")
    singleton_component_sizes = pd.to_numeric(
        frame.loc[frame["ssnv_branch"].eq("positional_singleton"), "component_size"],
        errors="raise",
    )
    assert_true(bool(singleton_component_sizes.eq(1).all()), "v10 singleton rows must all have component_size=1")
    return frame, recomputed


def load_site_audit(
    path: Path,
    singleton_keys: set[tuple[str, str, int, str, str]],
) -> pd.DataFrame:
    columns = [
        "dataset",
        "sample",
        "biological_id",
        "dataset_role",
        "chrom",
        "pos",
        "ref",
        "alt",
        "truth_label",
        "ssnv_branch",
        "component_size",
        "recomputed_component_size",
        "positional_contract_pass",
        "m1_evaluable",
        "m1_status",
        "m1_reason",
        "coarse_ng",
        "unstable",
        "modal_assignment_ari_min",
        "stable_null_multigroup",
        "cluster_sizes",
        "m2_status",
        "m2_bucket",
        "m2_reason_codes",
        "methyl_group_count",
        "core_read_n",
        "m2_aligned_axes",
        "m2_indeterminate_axes",
        "m2_constant_axes",
        "m2_low_power_axes",
        "m2_axis_details_json",
        "g1_status",
        "g2_status",
        "r1_status",
    ]
    frame = pd.read_csv(path, sep="\t", compression="gzip", usecols=columns, dtype=str)
    assert_equal(len(frame), EXPECTED_COUNTS["singleton_sites"], "site-audit row count")
    frame["pos"] = pd.to_numeric(frame["pos"], errors="raise").astype("int64")
    frame["_key"] = key_columns(frame)
    assert_equal(int(frame["_key"].duplicated().sum()), 0, "site-audit duplicate primary keys")
    assert_equal(set(frame["_key"]), singleton_keys, "site-audit singleton key set")
    assert_true(bool(frame["ssnv_branch"].eq("positional_singleton").all()), "site-audit branch identity")
    assert_true(
        bool(pd.to_numeric(frame["component_size"], errors="raise").eq(1).all()),
        "site-audit component_size=1",
    )
    assert_true(
        bool(pd.to_numeric(frame["recomputed_component_size"], errors="raise").eq(1).all()),
        "site-audit recomputed_component_size=1",
    )
    assert_true(parse_bool_series(frame["positional_contract_pass"], "positional_contract_pass").all(), "positional contract")

    frame["_m1_evaluable"] = parse_bool_series(frame["m1_evaluable"], "m1_evaluable")
    unstable = parse_nullable_bool_series(frame["unstable"], "unstable")
    frame["_stable"] = parse_bool_series(frame["stable_null_multigroup"], "stable_null_multigroup")
    coarse = pd.to_numeric(frame["coarse_ng"], errors="coerce")
    ari = pd.to_numeric(frame["modal_assignment_ari_min"], errors="coerce")
    expected_missing = ~frame["_m1_evaluable"]
    assert_true(bool(unstable.isna().eq(expected_missing).all()), "unstable nullability must equal M1 non-evaluable")
    assert_true(bool(ari.isna().eq(expected_missing).all()), "ARI nullability must equal M1 non-evaluable")
    frame["_unstable"] = unstable.fillna(False).astype(bool)
    recomputed_m1 = frame["_m1_evaluable"] & coarse.ge(2) & ~frame["_unstable"] & ari.ge(0.8)
    assert_true(bool(recomputed_m1.eq(frame["_stable"]).all()), "independent M1 formula mismatch")
    assert_true(bool(frame["_stable"].le(frame["_m1_evaluable"]).all()), "M1 flag outside evaluable population")
    assert_true(
        bool(frame["m1_status"].eq("FLAGGED").eq(frame["_stable"]).all()),
        "M1 status and recomputed flag mismatch",
    )
    unexpected_m2 = sorted(set(frame["m2_status"]) - set(M2_STATUSES))
    assert_equal(unexpected_m2, [], "M2 status vocabulary")
    assert_true(
        bool(frame["m2_status"].eq("NOT_RUN").eq(~frame["_stable"]).all()),
        "M2 NOT_RUN must equal M1 not flagged",
    )
    assert_true(
        bool(frame.loc[frame["_stable"], "m2_status"].isin(("NOT_EVALUABLE", "FAIL", "PASS")).all()),
        "M1 flags must have a non-NOT_RUN M2 state",
    )
    assert_true(bool(frame["g1_status"].eq("NOT_RUN").all()), "singleton G1 status")
    assert_true(bool(frame["g2_status"].eq("NOT_RUN").all()), "singleton G2 status")
    assert_true(bool(frame["r1_status"].eq("NOT_RUN").all()), "singleton R1 status")
    return frame


def summarize_group(frame: pd.DataFrame) -> dict[str, Any]:
    sites = len(frame)
    m1_evaluable = int(frame["_m1_evaluable"].sum())
    m1_flagged = int(frame["_stable"].sum())
    m2 = Counter(frame["m2_status"])
    determinate = int(m2["PASS"] + m2["FAIL"])

    def rate(numerator: int, denominator: int) -> float | None:
        return numerator / denominator if denominator else None

    return {
        "sites": sites,
        "m1_evaluable": m1_evaluable,
        "m1_not_evaluable": sites - m1_evaluable,
        "m1_flagged": m1_flagged,
        "m1_not_flagged": sites - m1_flagged,
        "m2_not_run": int(m2["NOT_RUN"]),
        "m2_not_evaluable": int(m2["NOT_EVALUABLE"]),
        "m2_fail": int(m2["FAIL"]),
        "m2_pass": int(m2["PASS"]),
        "m2_determinate": determinate,
        "m1_evaluable_pct_all": rate(m1_evaluable, sites),
        "m1_flagged_pct_all": rate(m1_flagged, sites),
        "m1_flagged_pct_evaluable": rate(m1_flagged, m1_evaluable),
        "m2_determinate_pct_m1": rate(determinate, m1_flagged),
        "m2_pass_pct_all": rate(int(m2["PASS"]), sites),
        "m2_pass_pct_m1": rate(int(m2["PASS"]), m1_flagged),
        "m2_pass_pct_determinate": rate(int(m2["PASS"]), determinate),
    }


def build_census_breakdowns(frame: pd.DataFrame) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    wide: list[dict[str, Any]] = []

    def add(stratum_type: str, stratum: str, subset: pd.DataFrame, dataset: str = "", truth: str = "") -> None:
        wide.append(
            {
                "stratum_type": stratum_type,
                "stratum": stratum,
                "dataset": dataset,
                "truth": truth,
                **summarize_group(subset),
            }
        )

    add("overall", "ALL", frame)
    for dataset in DATASET_ORDER:
        add("dataset", dataset, frame.loc[frame["dataset"].eq(dataset)], dataset=dataset)
    for biological_id in sorted(frame["biological_id"].unique()):
        add(
            "biological_sample",
            biological_id,
            frame.loc[frame["biological_id"].eq(biological_id)],
        )
    for truth in TRUTH_LABELS:
        add("truth", truth, frame.loc[frame["truth_label"].eq(truth)], truth=truth)
    for dataset in DATASET_ORDER:
        for truth in TRUTH_LABELS:
            subset = frame.loc[frame["dataset"].eq(dataset) & frame["truth_label"].eq(truth)]
            add("dataset_truth", f"{dataset}|{truth}", subset, dataset=dataset, truth=truth)

    long_rows: list[dict[str, Any]] = []
    for row in wide:
        status_values = {
            ("M1_EVALUABILITY", "EVALUABLE"): row["m1_evaluable"],
            ("M1_EVALUABILITY", "NOT_EVALUABLE"): row["m1_not_evaluable"],
            ("M1_FLAG", "FLAGGED"): row["m1_flagged"],
            ("M1_FLAG", "NOT_FLAGGED"): row["m1_not_flagged"],
            ("M2", "NOT_RUN"): row["m2_not_run"],
            ("M2", "NOT_EVALUABLE"): row["m2_not_evaluable"],
            ("M2", "FAIL"): row["m2_fail"],
            ("M2", "PASS"): row["m2_pass"],
        }
        for (axis, status_name), count in status_values.items():
            denominator = int(row["sites"])
            long_rows.append(
                {
                    "stratum_type": row["stratum_type"],
                    "stratum": row["stratum"],
                    "dataset": row["dataset"],
                    "truth": row["truth"],
                    "status_axis": axis,
                    "status": status_name,
                    "count": int(count),
                    "denominator": denominator,
                    "percent_of_singleton_stratum": count / denominator if denominator else None,
                }
            )
    return wide, long_rows


def compare_recomputed_to_frozen_summary(
    overall: Mapping[str, Any],
    wide_rows: Sequence[Mapping[str, Any]],
    frozen_summary: Mapping[str, Any],
) -> None:
    expected_overall = {
        "singleton_sites": overall["sites"],
        "m1_evaluable": overall["m1_evaluable"],
        "m1_flagged": overall["m1_flagged"],
        "m1_not_evaluable": overall["m1_not_evaluable"],
        "m2_not_run": overall["m2_not_run"],
        "m2_not_evaluable": overall["m2_not_evaluable"],
        "m2_fail": overall["m2_fail"],
        "m2_pass": overall["m2_pass"],
    }
    assert_equal(expected_overall, frozen_summary["counts"], "recomputed overall vs frozen summary")
    by_dataset = {
        row["stratum"]: {
            key: row[key]
            for key in ("sites", "m1_evaluable", "m1_flagged", "m2_not_run", "m2_not_evaluable", "m2_fail", "m2_pass")
        }
        for row in wide_rows
        if row["stratum_type"] == "dataset"
    }
    assert_equal(by_dataset, frozen_summary["breakdowns"]["dataset"], "dataset census vs frozen summary")
    by_truth = {
        row["stratum"]: {
            key: row[key]
            for key in ("sites", "m1_evaluable", "m1_flagged", "m2_not_run", "m2_not_evaluable", "m2_fail", "m2_pass")
        }
        for row in wide_rows
        if row["stratum_type"] == "truth"
    }
    assert_equal(by_truth, frozen_summary["breakdowns"]["truth"], "truth census vs frozen summary")


def find_target_rows(frame: pd.DataFrame, targets: Sequence[Target]) -> dict[tuple[str, str, int, str, str], pd.Series]:
    result: dict[tuple[str, str, int, str, str], pd.Series] = {}
    for target in targets:
        mask = (
            frame["dataset"].eq(target.dataset)
            & frame["chrom"].eq(target.chrom)
            & frame["pos"].eq(target.pos)
            & frame["ref"].eq(target.ref)
            & frame["alt"].eq(target.alt)
        )
        matches = frame.loc[mask]
        assert_equal(len(matches), 1, f"target row count {target.key}")
        result[target.key] = matches.iloc[0]
    return result


def load_target_assignments(
    path: Path,
    targets: Sequence[Target],
) -> dict[tuple[str, str, int, str, str], dict[str, Any]]:
    by_position = {(target.dataset, target.chrom, target.pos): target for target in targets}
    results: dict[tuple[str, str, int, str, str], dict[str, Any]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if '"dataset":"HCC1395"' not in line:
                continue
            if not any(f'"pos":{target.pos}' in line for target in targets):
                continue
            payload = json.loads(line)
            position_key = (str(payload.get("dataset")), str(payload.get("chrom")), int(payload.get("pos")))
            target = by_position.get(position_key)
            if target is None:
                continue
            posthoc = payload.get("posthoc") or {}
            assert_equal(str(posthoc.get("ref")), target.ref, f"{target.short_id} assignment ref")
            assert_equal(str(posthoc.get("alt")), target.alt, f"{target.short_id} assignment alt")
            assert_equal(str(posthoc.get("truth_label")), target.truth, f"{target.short_id} assignment truth")
            if target.key in results:
                raise ValidationError(f"duplicate stable assignment for {target.key}")
            results[target.key] = payload
    assert_equal(set(results), {target.key for target in targets}, "target stable-assignment key set")
    return results


def natural_read_id(value: str) -> tuple[int, int | str]:
    try:
        return 0, int(value)
    except ValueError:
        return 1, value


def select_medoid_nearest(
    core_ids: Sequence[str],
    group_by_id: Mapping[str, str],
    distance: np.ndarray,
    maximum_per_group: int = 15,
) -> tuple[list[str], list[dict[str, Any]]]:
    if distance.shape != (len(core_ids), len(core_ids)):
        raise ValidationError("distance/core ID dimensions disagree")
    id_to_index = {read_id: index for index, read_id in enumerate(core_ids)}
    if len(id_to_index) != len(core_ids):
        raise ValidationError("duplicate core read IDs")
    groups = sorted(set(group_by_id.values()))
    if groups != ["Group A", "Group B"]:
        raise ValidationError(f"expected Group A/B aliases, observed {groups}")

    selected: list[str] = []
    audit_rows: list[dict[str, Any]] = []
    for group in groups:
        group_ids = sorted(
            [read_id for read_id in core_ids if group_by_id[read_id] == group],
            key=natural_read_id,
        )
        indices = [id_to_index[read_id] for read_id in group_ids]
        within = distance[np.ix_(indices, indices)]
        mean_distances = within.mean(axis=1)
        medoid_id = min(
            zip(group_ids, mean_distances),
            key=lambda item: (float(item[1]), natural_read_id(item[0])),
        )[0]
        medoid_index = id_to_index[medoid_id]
        ranked = sorted(
            group_ids,
            key=lambda read_id: (
                float(distance[medoid_index, id_to_index[read_id]]),
                natural_read_id(read_id),
            ),
        )
        selected_group = ranked[:maximum_per_group]
        selected.extend(selected_group)
        rank_by_id = {read_id: index + 1 for index, read_id in enumerate(ranked)}
        for read_id in ranked:
            audit_rows.append(
                {
                    "group": group,
                    "read_id": read_id,
                    "medoid_read_id": medoid_id,
                    "distance_to_medoid": float(distance[medoid_index, id_to_index[read_id]]),
                    "rank_from_medoid": rank_by_id[read_id],
                    "selected_for_display": read_id in selected_group,
                    "selection_cap_per_group": maximum_per_group,
                }
            )
    return selected, audit_rows


def validate_square_distance(
    distance_frame: pd.DataFrame,
    core_ids: Sequence[str],
    tolerance: float = 1e-12,
) -> tuple[np.ndarray, dict[str, Any]]:
    row_ids = distance_frame["read_id"].astype(str).tolist()
    column_ids = [str(column) for column in distance_frame.columns[1:]]
    assert_equal(len(row_ids), len(set(row_ids)), "distance row IDs unique")
    assert_equal(len(column_ids), len(set(column_ids)), "distance column IDs unique")
    assert_equal(set(row_ids), set(column_ids), "distance row/column identity")
    missing = sorted(set(core_ids) - set(row_ids), key=natural_read_id)
    assert_equal(missing, [], "core reads missing from distance matrix")
    indexed = distance_frame.set_index("read_id")
    core = indexed.loc[list(core_ids), list(core_ids)].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    finite = bool(np.isfinite(core).all())
    symmetric_error = float(np.nanmax(np.abs(core - core.T)))
    diagonal_error = float(np.nanmax(np.abs(np.diag(core))))
    assert_true(finite, "core distance matrix contains non-finite values")
    assert_true(symmetric_error <= tolerance, f"core distance matrix asymmetric: {symmetric_error}")
    assert_true(diagonal_error <= tolerance, f"core distance matrix diagonal nonzero: {diagonal_error}")
    return core, {
        "core_finite": finite,
        "maximum_symmetry_error": symmetric_error,
        "maximum_diagonal_abs": diagonal_error,
        "raw_matrix_rows": len(row_ids),
        "raw_matrix_nonfinite_cells": int(
            distance_frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").isna().sum().sum()
        ),
    }


def validate_methylation(
    methylation: pd.DataFrame,
    reads: pd.DataFrame,
    cpg_sites: pd.DataFrame,
    core_ids: Sequence[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    methyl_ids = methylation["read_id"].astype(str).tolist()
    read_ids = reads["read_id"].astype(str).tolist()
    assert_equal(len(methyl_ids), len(set(methyl_ids)), "methylation row IDs unique")
    assert_equal(set(methyl_ids), set(read_ids), "methylation/read row identity")
    expected_columns = cpg_sites["position"].astype(str).tolist()
    observed_columns = [str(column) for column in methylation.columns[1:]]
    assert_equal(observed_columns, expected_columns, "methylation/CpG column identity")
    assert_equal(len(observed_columns), len(set(observed_columns)), "methylation CpG columns unique")
    missing_core = sorted(set(core_ids) - set(methyl_ids), key=natural_read_id)
    assert_equal(missing_core, [], "core reads missing from methylation matrix")
    indexed = methylation.set_index("read_id").loc[list(core_ids)].apply(pd.to_numeric, errors="coerce")
    values = indexed.to_numpy(float)
    finite_values = values[np.isfinite(values)]
    assert_true(bool(((finite_values >= 0.0) & (finite_values <= 1.0)).all()), "methyl probabilities outside [0,1]")
    return indexed, {
        "raw_methylation_rows": len(methyl_ids),
        "cpg_columns": len(observed_columns),
        "core_called_cells": int(np.isfinite(values).sum()),
        "core_missing_cells": int(np.isnan(values).sum()),
    }


def render_methylation_heatmap(
    output_path: Path,
    target: Target,
    selected_ids: Sequence[str],
    visual_label_by_id: Mapping[str, str],
    group_by_id: Mapping[str, str],
    methyl_core: pd.DataFrame,
) -> None:
    if output_path.exists():
        raise FileExistsError(output_path)
    values = methyl_core.loc[list(selected_ids)].to_numpy(float)
    masked = np.ma.masked_invalid(values)
    cmap = plt.get_cmap("cividis").copy()
    cmap.set_bad("#d7dce2")
    width = max(10.5, min(18.0, 0.14 * values.shape[1] + 6.0))
    height = max(6.8, 0.24 * values.shape[0] + 2.4)
    fig, axis = plt.subplots(figsize=(width, height), constrained_layout=True)
    image = axis.imshow(masked, aspect="auto", interpolation="nearest", vmin=0.0, vmax=1.0, cmap=cmap)
    labels = [visual_label_by_id[read_id] for read_id in selected_ids]
    axis.set_yticks(np.arange(len(labels)))
    axis.set_yticklabels(labels, fontsize=7)
    cpg_positions = [int(column) for column in methyl_core.columns]
    maximum_ticks = 12
    tick_step = max(1, math.ceil(len(cpg_positions) / maximum_ticks))
    x_indices = list(range(0, len(cpg_positions), tick_step))
    if x_indices[-1] != len(cpg_positions) - 1:
        x_indices.append(len(cpg_positions) - 1)
    axis.set_xticks(x_indices)
    axis.set_xticklabels(
        [f"{cpg_positions[index]:,}" for index in x_indices],
        rotation=45,
        ha="right",
        fontsize=7,
    )
    group_a_count = sum(group_by_id[read_id] == "Group A" for read_id in selected_ids)
    axis.axhline(group_a_count - 0.5, color="#171411", linewidth=1.2)
    axis.set_xlabel("CpG genomic position")
    axis.set_ylabel("Deterministic medoid-nearest representative reads")
    axis.set_title(f"{target.display_locus}: focal-ALT read × CpG methylation probability", fontsize=12)
    colorbar = fig.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
    colorbar.set_label("Methylation probability")
    fig.text(
        0.01,
        0.005,
        "Gray = CpG not called. Group A/B are methyl-assignment aliases, not HP labels.",
        fontsize=8,
        color="#4f5660",
    )
    fig.savefig(output_path, dpi=145, facecolor="white")
    plt.close(fig)


def data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def pairwise_off_diagonal_mean(matrix: np.ndarray) -> float:
    if matrix.shape[0] < 2:
        return float("nan")
    mask = ~np.eye(matrix.shape[0], dtype=bool)
    return float(matrix[mask].mean())


def process_target(
    target: Target,
    v10_row: pd.Series,
    audit_row: pd.Series,
    assignment: Mapping[str, Any],
    figures_dir: Path,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    region_dir = Path(str(v10_row["region_dir"]))
    assert_equal(Path(str(assignment["region_dir"])).resolve(), region_dir.resolve(), f"{target.short_id} region_dir")
    reads_path = region_dir / "reads" / "reads.tsv"
    methyl_path = region_dir / "methylation" / "methylation.csv"
    cpg_path = region_dir / "methylation" / "cpg_sites.tsv"
    distance_path = region_dir / "distance" / "BERNOULLI" / "matrix.csv"
    significance_path = region_dir / "clustering" / "significance.json"

    primary = assignment.get("primary_artifacts") or {}
    primary_role_by_path = {
        reads_path.resolve(): "reads",
        methyl_path.resolve(): "methylation_matrix",
        distance_path.resolve(): "distance_matrix",
    }
    for resolved_path, primary_role in primary_role_by_path.items():
        record = primary.get(primary_role) or {}
        assert_equal(Path(str(record.get("path"))).resolve(), resolved_path, f"{target.short_id} {primary_role} path")
        assert_equal(int(record.get("size_bytes")), resolved_path.stat().st_size, f"{target.short_id} {primary_role} size")
        assert_equal(str(record.get("sha256")), sha256_path(resolved_path), f"{target.short_id} {primary_role} hash")

    reads = pd.read_csv(reads_path, sep="\t", dtype={"read_id": str})
    methylation = pd.read_csv(methyl_path, dtype={"read_id": str}, na_values=["NA"])
    cpg_sites = pd.read_csv(cpg_path, sep="\t", dtype={"position": str})
    distance_frame = pd.read_csv(distance_path, dtype={"read_id": str}, na_values=["NA"])
    assert_equal(len(reads), len(set(reads["read_id"])), f"{target.short_id} reads unique")
    assert_equal(set(reads["read_id"]), set(distance_frame["read_id"]), f"{target.short_id} reads/distance rows")

    core_reads = list(assignment.get("core_reads") or [])
    core_ids = [str(item["read_id"]) for item in core_reads]
    assert_equal(len(core_ids), target.core_reads, f"{target.short_id} core read count")
    assert_equal([str(value) for value in assignment.get("read_ids") or []], core_ids, f"{target.short_id} assignment read IDs")
    assert_equal([str(value) for value in assignment.get("labels") or []], [str(item["label"]) for item in core_reads], f"{target.short_id} labels")
    assert_equal(int(audit_row["core_read_n"]), target.core_reads, f"{target.short_id} audit core count")
    assert_equal(int(v10_row["n_alt_matrix"]), target.raw_alt_reads, f"{target.short_id} raw focal-ALT count")
    assert_equal(int(audit_row["m2_status"] == "PASS"), 1, f"{target.short_id} M2 PASS")
    assert_equal(str(audit_row["truth_label"]), target.truth, f"{target.short_id} truth")

    reads_by_id = reads.set_index("read_id")
    raw_alt_ids = set(reads.loc[reads["alt_support"].eq("ALT"), "read_id"])
    assert_equal(len(raw_alt_ids), target.raw_alt_reads, f"{target.short_id} raw ALT rows")
    assert_true(set(core_ids).issubset(raw_alt_ids), f"{target.short_id} core reads must all be focal ALT")
    for item in core_reads:
        read_id = str(item["read_id"])
        assert_equal(str(reads_by_id.loc[read_id, "read_name"]), str(item["read_name"]), f"{target.short_id} read-name join {read_id}")

    distance_core, distance_checks = validate_square_distance(distance_frame, core_ids)
    methyl_core, methyl_checks = validate_methylation(methylation, reads, cpg_sites, core_ids)

    source_labels = sorted(set(str(item["label"]) for item in core_reads))
    assert_equal(len(source_labels), 2, f"{target.short_id} source methyl groups")
    alias = {source_labels[0]: "Group A", source_labels[1]: "Group B"}
    group_by_id = {str(item["read_id"]): alias[str(item["label"])] for item in core_reads}
    group_counts = Counter(group_by_id.values())
    assert_equal(sorted(group_counts.values(), reverse=True), sorted(target.group_sizes, reverse=True), f"{target.short_id} group sizes")

    selected_ids, selection_rows = select_medoid_nearest(core_ids, group_by_id, distance_core, maximum_per_group=15)
    selected_by_group = Counter(group_by_id[read_id] for read_id in selected_ids)
    assert_true(all(count <= 15 for count in selected_by_group.values()), f"{target.short_id} display cap")
    visual_label_by_id: dict[str, str] = {}
    for group in ("Group A", "Group B"):
        for rank, read_id in enumerate(
            [value for value in selected_ids if group_by_id[value] == group],
            start=1,
        ):
            visual_label_by_id[read_id] = f"{group} · {rank:02d}"

    id_to_index = {read_id: index for index, read_id in enumerate(core_ids)}
    group_indices = {
        group: [id_to_index[read_id] for read_id in core_ids if group_by_id[read_id] == group]
        for group in ("Group A", "Group B")
    }
    group_a_distance = distance_core[np.ix_(group_indices["Group A"], group_indices["Group A"])]
    group_b_distance = distance_core[np.ix_(group_indices["Group B"], group_indices["Group B"])]
    between_distance = distance_core[np.ix_(group_indices["Group A"], group_indices["Group B"])]
    within_a = pairwise_off_diagonal_mean(group_a_distance)
    within_b = pairwise_off_diagonal_mean(group_b_distance)
    within_pairs = np.concatenate(
        (
            group_a_distance[~np.eye(group_a_distance.shape[0], dtype=bool)],
            group_b_distance[~np.eye(group_b_distance.shape[0], dtype=bool)],
        )
    )
    pooled_within = float(within_pairs.mean())
    mean_between = float(between_distance.mean())

    methyl_values = methyl_core.to_numpy(float)
    methyl_group_stats: dict[str, dict[str, Any]] = {}
    for group in ("Group A", "Group B"):
        indices = group_indices[group]
        values = methyl_values[indices]
        methyl_group_stats[group] = {
            "mean": float(np.nanmean(values)),
            "called": int(np.isfinite(values).sum()),
            "missing": int(np.isnan(values).sum()),
        }

    selection_output: list[dict[str, Any]] = []
    name_hash_by_id = {
        str(item["read_id"]): sha256_text(str(item["read_name"]))
        for item in core_reads
    }
    for row in selection_rows:
        read_id = str(row["read_id"])
        selection_output.append(
            {
                "locus": target.display_locus,
                "group": row["group"],
                "read_id": read_id,
                "read_name_sha256": name_hash_by_id[read_id],
                "medoid_read_id": row["medoid_read_id"],
                "distance_to_medoid": row["distance_to_medoid"],
                "rank_from_medoid": row["rank_from_medoid"],
                "selected_for_display": row["selected_for_display"],
                "display_label": visual_label_by_id.get(read_id, ""),
                "selection_cap_per_group": row["selection_cap_per_group"],
                "selection_population": "all core focal-ALT reads in the assigned group",
            }
        )

    selected_indices = [id_to_index[read_id] for read_id in selected_ids]
    selected_distance = distance_core[np.ix_(selected_indices, selected_indices)]
    distance_chart_rows: list[dict[str, Any]] = []
    heatmap_cells: list[dict[str, Any]] = []
    for row_index, row_id in enumerate(selected_ids):
        for column_index, column_id in enumerate(selected_ids):
            value = float(selected_distance[row_index, column_index])
            distance_chart_rows.append(
                {
                    "row_read": visual_label_by_id[row_id],
                    "column_read": visual_label_by_id[column_id],
                    "distance": value,
                    "row_group": group_by_id[row_id],
                    "column_group": group_by_id[column_id],
                }
            )
            heatmap_cells.append(
                {
                    "locus": target.display_locus,
                    "heatmap_type": "read_distance",
                    "row_label": visual_label_by_id[row_id],
                    "column_label": visual_label_by_id[column_id],
                    "value": value,
                    "missing": False,
                }
            )

    for read_id in selected_ids:
        for cpg_position, value in methyl_core.loc[read_id].items():
            missing = bool(pd.isna(value))
            heatmap_cells.append(
                {
                    "locus": target.display_locus,
                    "heatmap_type": "read_by_cpg_methylation",
                    "row_label": visual_label_by_id[read_id],
                    "column_label": str(cpg_position),
                    "value": "" if missing else float(value),
                    "missing": missing,
                }
            )

    figure_path = figures_dir / f"{target.short_id}.read_by_cpg_methylation.png"
    render_methylation_heatmap(
        figure_path,
        target,
        selected_ids,
        visual_label_by_id,
        group_by_id,
        methyl_core,
    )

    assignment_axes = assignment.get("associations") or {}
    expected_axes = {"hp_exact", "hp_family", "strand", "start", "end", "length", "mapq", "cpg_called"}
    assert_equal(set(assignment_axes), expected_axes, f"{target.short_id} M2 measured axes")
    assert_true("ps" not in assignment_axes, f"{target.short_id} PS must not be an M2 measured axis")
    assert_true(
        all(not bool(details.get("aligned")) for details in assignment_axes.values()),
        f"{target.short_id} aligned measured axis",
    )
    ps_values = {item.get("latest_ps") for item in core_reads}

    case_stats = {
        "locus": target.display_locus,
        "dataset": target.dataset,
        "truth": target.truth,
        "m2_status": "PASS",
        "caller_af": float(v10_row["caller_af"]),
        "raw_distance_matrix_reads": distance_checks["raw_matrix_rows"],
        "raw_focal_alt_reads": target.raw_alt_reads,
        "core_reads": target.core_reads,
        "group_a_reads": int(group_counts["Group A"]),
        "group_b_reads": int(group_counts["Group B"]),
        "display_group_a_reads": int(selected_by_group["Group A"]),
        "display_group_b_reads": int(selected_by_group["Group B"]),
        "cpg_sites": methyl_checks["cpg_columns"],
        "core_called_methyl_cells": methyl_checks["core_called_cells"],
        "core_missing_methyl_cells": methyl_checks["core_missing_cells"],
        "core_missing_methyl_fraction": methyl_checks["core_missing_cells"]
        / (methyl_checks["core_called_cells"] + methyl_checks["core_missing_cells"]),
        "group_a_mean_methyl_probability": methyl_group_stats["Group A"]["mean"],
        "group_b_mean_methyl_probability": methyl_group_stats["Group B"]["mean"],
        "absolute_group_mean_difference": abs(
            methyl_group_stats["Group A"]["mean"] - methyl_group_stats["Group B"]["mean"]
        ),
        "group_a_within_mean_distance": within_a,
        "group_b_within_mean_distance": within_b,
        "pooled_within_mean_distance": pooled_within,
        "between_group_mean_distance": mean_between,
        "between_to_within_distance_ratio": mean_between / pooled_within,
        "ps_context_unique_values": len(ps_values),
        "ps_is_m2_measured_axis": False,
        "confirmed_cellular_subclone": 0,
        "linear_ancestry": 0,
        "claim_tier": "operational_M2_read_level_residual_epigenetic_partition",
    }
    join_checks = {
        "locus": target.display_locus,
        "region_dir": str(region_dir.resolve()),
        "core_expected": target.core_reads,
        "core_in_reads": len(set(core_ids) & set(reads["read_id"])),
        "core_in_methylation": len(set(core_ids) & set(methylation["read_id"])),
        "core_in_distance_rows": len(set(core_ids) & set(distance_frame["read_id"])),
        "core_in_distance_columns": len(set(core_ids) & set(distance_frame.columns[1:])),
        "read_name_identity_matches": target.core_reads,
        "methyl_cpg_columns": methyl_checks["cpg_columns"],
        "methyl_cpg_identity_pass": True,
        **distance_checks,
        "pass": True,
    }
    target_sources = [
        source_record(reads_path, f"{target.short_id}_reads", "core read identity and focal-ALT membership"),
        source_record(methyl_path, f"{target.short_id}_methylation", "read-by-CpG matrix and static heatmap"),
        source_record(cpg_path, f"{target.short_id}_cpg_sites", "methylation column identity"),
        source_record(distance_path, f"{target.short_id}_distance", "core distance validation, medoids, and distance heatmap"),
        source_record(significance_path, f"{target.short_id}_legacy_clustering_context", "inventory only; not used for M1/M2 decisions"),
    ]
    case_stats["methyl_figure_path"] = str(figure_path.relative_to(figures_dir.parent))
    return case_stats, selection_output, distance_chart_rows, heatmap_cells, [join_checks], target_sources


def percent(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "NA"
    return f"{100.0 * value:.{digits}f}%"


def artifact_sources(source_inventory: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    by_role = {str(row["role"]): row for row in source_inventory}
    return [
        {
            "id": "singleton_audit",
            "label": "Source-attested positional-singleton site audit",
            "path": "positional_singleton_site_audit.tsv.gz",
            "query": {
                "engine": "Python 3 / pandas",
                "language": "Python",
                "description": "Independent count and denominator census over all attested singleton rows.",
                "tables_used": ["positional_singleton_site_audit.tsv.gz"],
                "filters": [
                    "ssnv_branch=positional_singleton",
                    "dataset in seven frozen datasets",
                    "chromosome chr1-22",
                ],
                "metric_definitions": [
                    "M1 flag = coarse_ng>=2 AND unstable=false AND modal_assignment_ari_min>=0.8.",
                    "M2 counts retain NOT_RUN, NOT_EVALUABLE, FAIL, and PASS as disjoint states.",
                ],
            },
            "sha256": by_role["site_audit"]["sha256"],
        },
        {
            "id": "v10_sites",
            "label": "v10 source-locked all-sSNV site results",
            "path": "all_ssnv_site_results.tsv.gz",
            "query": {
                "engine": "Python 3 / pandas",
                "language": "Python",
                "description": "Recomputed dataset/chrom positional components from all 469,849 dataset-sites.",
                "tables_used": ["all_ssnv_site_results.tsv.gz"],
                "filters": [
                    "LongPhase-S recalibrated FILTER=PASS",
                    "autosomal biallelic sSNVs",
                    "adjacent gap <=50,000 bp joins one transitive component",
                ],
                "metric_definitions": [
                    "Positional singleton = transitive component size one within dataset and chromosome.",
                ],
            },
            "sha256": by_role["v10_site_results"]["sha256"],
        },
        {
            "id": "stable_assignments",
            "label": "v10 stable focal-ALT read assignments",
            "path": "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
            "query": {
                "engine": "Python 3 streaming JSONL",
                "language": "Python",
                "description": "Exact two-key extraction followed by read-ID joins to frozen raw matrices.",
                "tables_used": ["all_ssnv_stable_multigroup_read_assignments.jsonl.gz"],
                "filters": [
                    "dataset=HCC1395",
                    "target loci fixed before inspection",
                    "all core reads used for statistics",
                    "display limited to deterministic medoid-nearest representatives",
                ],
            },
            "sha256": by_role["stable_assignments"]["sha256"],
        },
        {
            "id": "claim_contract",
            "label": "Claim contract v5",
            "path": "claim-contract-v5.md",
            "query": {
                "description": "Evidence ladder, M2 measured-axis gate, denominator locks, and claim ceiling.",
                "tables_used": ["claim-contract-v5.md"],
            },
            "sha256": by_role["claim_contract"]["sha256"],
        },
        {
            "id": "chr14_raw",
            "label": "HCC1395 chr14 target raw read/methyl/distance artifacts",
            "path": "chr14_86267476_86277476",
            "query": {
                "engine": "Python 3 / numpy",
                "language": "Python",
                "description": "Exact core-read join, full-core statistics, deterministic representative display.",
                "tables_used": ["reads.tsv", "methylation.csv", "cpg_sites.tsv", "distance/BERNOULLI/matrix.csv"],
                "filters": ["core focal-ALT reads only for calculations", "maximum 15 representatives per Group A/B for display"],
            },
        },
        {
            "id": "chr22_raw",
            "label": "HCC1395 chr22 target raw read/methyl/distance artifacts",
            "path": "chr22_47461517_47471517",
            "query": {
                "engine": "Python 3 / numpy",
                "language": "Python",
                "description": "Exact core-read join, full-core statistics, deterministic representative display.",
                "tables_used": ["reads.tsv", "methylation.csv", "cpg_sites.tsv", "distance/BERNOULLI/matrix.csv"],
                "filters": ["core focal-ALT reads only for calculations", "maximum 15 representatives per Group A/B for display"],
            },
        },
    ]


def build_artifact(
    census_wide: Sequence[Mapping[str, Any]],
    census_long: Sequence[Mapping[str, Any]],
    case_stats: Sequence[Mapping[str, Any]],
    distance_rows_by_target: Mapping[str, Sequence[Mapping[str, Any]]],
    figure_paths_by_target: Mapping[str, Path],
    source_inventory: Sequence[Mapping[str, Any]],
    created_at: str,
    staging_sqlite_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    overall = next(row for row in census_wide if row["stratum_type"] == "overall")
    dataset_rows = [dict(row) for row in census_wide if row["stratum_type"] == "dataset"]
    truth_rows = [dict(row) for row in census_wide if row["stratum_type"] == "truth"]
    hcc = next(row for row in dataset_rows if row["stratum"] == "HCC1395")
    dataset_status = [
        row
        for row in census_long
        if row["stratum_type"] == "dataset" and row["status_axis"] == "M2"
    ]
    hcc_funnel = [
        {"stage": "All HCC1395 sSNVs", "count": 79_687, "denominator_note": "v10 dataset-sites"},
        {"stage": "Positional singleton", "count": 8_279, "denominator_note": "component_size=1"},
        {"stage": "M1 evaluable", "count": 8_074, "denominator_note": "singleton sites"},
        {"stage": "M1 flagged", "count": 734, "denominator_note": "singleton sites"},
        {"stage": "M2 determinate", "count": 2, "denominator_note": "PASS+FAIL"},
        {"stage": "M2 PASS", "count": 2, "denominator_note": "operational yield"},
    ]

    public_source_rows = [
        {
            "role": row["role"],
            "file": Path(str(row["path"])).name,
            "sha256": row["sha256"],
            "size_bytes": row["size_bytes"],
            "used_for": row["used_for"],
        }
        for row in source_inventory
        if row["role"]
        in {
            "audit_summary",
            "site_audit",
            "m2_pass_cases",
            "v10_site_results",
            "stable_assignments",
            "claim_contract",
        }
    ]
    sources = artifact_sources(source_inventory)

    status_chart = {
        "id": "all_dataset_m2_status_chart",
        "title": "Positional-singleton M2 status counts across seven datasets",
        "subtitle": "NOT_RUN, NOT_EVALUABLE, FAIL, and PASS are disjoint; tiny determinate segments require the exact table.",
        "type": "horizontalStackedBar",
        "dataset": "dataset_m2_status",
        "sourceId": "singleton_audit",
        "intent": "composition",
        "question": "How are all positional-singleton sites partitioned across M2 states in each dataset?",
        "rationale": "A stacked count chart preserves the full singleton denominator and keeps non-evaluable states visible.",
        "comparisonContext": {
            "denominator": "all positional-singleton dataset-sites within each dataset",
            "grain": "dataset by M2 status",
            "unit": "dataset-sites",
        },
        "encodings": {
            "x": {"field": "stratum", "type": "nominal", "label": "Dataset"},
            "y": {"field": "count", "type": "quantitative", "label": "Singleton sites"},
            "color": {"field": "status", "type": "nominal", "label": "M2 status"},
        },
        "palette": {"kind": "categorical"},
        "legend": {"position": "bottom", "interactive": True},
        "valueFormat": "number",
        "showDescription": True,
    }
    funnel_chart = {
        "id": "hcc1395_funnel_chart",
        "title": "HCC1395 positional-singleton M1/M2 funnel",
        "subtitle": "M2 determinate=2; NOT_EVALUABLE=732 and NOT_RUN=7,545 remain explicit outside the terminal PASS bar.",
        "type": "funnel",
        "dataset": "hcc1395_funnel",
        "sourceId": "singleton_audit",
        "intent": "funnel",
        "question": "How does the frozen HCC1395 population narrow through positional, M1, and M2 gates?",
        "rationale": "Ordered stage counts show gate attrition without treating NOT_RUN or NOT_EVALUABLE as biological negatives.",
        "comparisonContext": {
            "denominator": "stage-specific frozen populations",
            "grain": "HCC1395 gate stage",
            "unit": "dataset-sites",
        },
        "encodings": {
            "x": {"field": "stage", "type": "ordinal", "label": "Gate stage"},
            "y": {"field": "count", "type": "quantitative", "label": "Sites"},
        },
        "palette": {"kind": "sequential"},
        "valueFormat": "number",
        "showDescription": True,
    }

    charts = [status_chart, funnel_chart]
    blocks: list[dict[str, Any]] = [
        {"id": "title", "type": "markdown", "body": f"# {REPORT_TITLE}"},
        {
            "id": "summary",
            "type": "markdown",
            "sourceId": "singleton_audit",
            "body": (
                "## Summary\n\n"
                f"- 全量重算得到 **{overall['sites']:,}/{EXPECTED_COUNTS['all_dataset_sites']:,}="
                f"{percent(overall['sites']/EXPECTED_COUNTS['all_dataset_sites'], 6)}** positional-singleton dataset-sites。\n"
                f"- M1 evaluable **{overall['m1_evaluable']:,}**，M1 flagged **{overall['m1_flagged']:,}**；"
                f"M2 為 PASS **{overall['m2_pass']:,}**、FAIL **{overall['m2_fail']:,}**、"
                f"NOT_EVALUABLE **{overall['m2_not_evaluable']:,}**、NOT_RUN **{overall['m2_not_run']:,}**。\n"
                f"- `30/50,432={percent(overall['m2_pass_pct_all'], 6)}` 只稱 **observed operational M2-PASS yield**，"
                "不是 biological/subclone prevalence，也不是其 lower bound；`30/48=62.5%` "
                "只是在 M2-determinate 子集中的 conditional rate。\n"
                f"- HCC1395 重算 **{hcc['sites']:,}** singletons、M1 flagged **{hcc['m1_flagged']:,}**、"
                f"M2 PASS **{hcc['m2_pass']:,}**。兩個 PASS 位點的 core join 分別為 108/108 與 109/109。\n"
                "- **confirmed cellular subclone=0、linear ancestry=0 是必要 downstream 尚未完成的未驗證數，"
                "不是真陰性。**"
            ),
        },
        {
            "id": "findings_heading",
            "type": "markdown",
            "body": (
                "## Findings / visual evidence\n\n"
                "圖表先保留完整狀態分母，再展示 HCC1395 唯二 M2 PASS 位點。所有群組只以 "
                "**Group A/B** 表示；這是 methyl-assignment alias，不是 HP。"
            ),
        },
        {"id": "all_dataset_status", "type": "chart", "chartId": status_chart["id"]},
        {
            "id": "all_dataset_status_explanation",
            "type": "markdown",
            "sourceId": "singleton_audit",
            "body": (
                "### 全樣本狀態圖怎麼讀\n\n"
                "NOT_RUN 表示 M1 未標記，NOT_EVALUABLE 表示 M1 已標記但 M2 資訊不足；兩者都不是 FAIL。"
                "PASS/FAIL 很小，精確數字以表格為準，不從像素寬度估算。"
            ),
        },
        {"id": "hcc1395_funnel", "type": "chart", "chartId": funnel_chart["id"]},
        {
            "id": "hcc1395_funnel_explanation",
            "type": "markdown",
            "sourceId": "singleton_audit",
            "body": (
                "### HCC1395 funnel 怎麼讀\n\n"
                "M2 套用於全部 734 個 M1 flags；其中 2 個 determinate 且皆 PASS，另 732 個為 "
                "NOT_EVALUABLE。`2/8,279` 是 HCC1395 operational yield；`2/2` 不能外推成 prevalence。"
            ),
        },
        {"id": "case_statistics", "type": "table", "tableId": "case_statistics_table"},
    ]

    for target in TARGETS:
        target_id = target.short_id
        source_id = "chr14_raw" if target.chrom == "chr14" else "chr22_raw"
        distance_chart = {
            "id": f"{target_id}_distance_heatmap",
            "title": f"{target.display_locus}: representative read-distance matrix",
            "subtitle": "15 deterministic medoid-nearest representatives per Group A/B; medoids and statistics use all core reads.",
            "type": "heatmap",
            "dataset": f"{target_id}_distance_cells",
            "sourceId": source_id,
            "intent": "relationship",
            "question": "Do representative reads retain the full-core methyl-distance block pattern?",
            "rationale": "A square heatmap exposes within- versus between-group distance structure while capping display density.",
            "comparisonContext": {
                "denominator": f"all {target.core_reads} core focal-ALT reads for medoid/statistical calculations",
                "grain": "representative read pair",
                "unit": "Bernoulli methylation distance",
            },
            "encodings": {
                "x": {"field": "row_read", "type": "nominal", "label": "Representative read"},
                "y": {"field": "distance", "type": "quantitative", "label": "Bernoulli distance"},
                "color": {"field": "column_read", "type": "nominal", "label": "Representative read"},
            },
            "palette": {"kind": "sequential"},
            "valueFormat": "number",
            "showDescription": True,
        }
        charts.append(distance_chart)
        image_path = figure_paths_by_target[target_id]
        image_alt = (
            f"{target.display_locus} read by CpG methylation probability heatmap; "
            "rows are deterministic Group A and Group B representatives and gray cells are uncalled CpGs"
        )
        image_body = (
            '<figure style="margin:0;font-family:system-ui,sans-serif;color:#171411">'
            f'<img alt="{image_alt}" src="{data_uri(image_path)}" '
            'style="display:block;width:100%;height:auto;border:0">'
            "<figcaption style=\"margin-top:10px;font-size:13px;line-height:1.5;color:#565d66\">"
            "Read×CpG methylation probability. Gray denotes an uncalled CpG; values were not imputed. "
            "The native portable heatmap maps missing values to the minimum intensity and therefore cannot "
            "represent the required NA state faithfully; this reproducible PNG is embedded in the canonical "
            "artifact without an additional HTML runtime."
            "</figcaption></figure>"
        )
        blocks.extend(
            [
                {
                    "id": f"{target_id}_heading",
                    "type": "markdown",
                    "sourceId": source_id,
                    "body": (
                        f"### {target.display_locus}\n\n"
                        f"TP、M2 PASS；raw focal-ALT={target.raw_alt_reads}，core={target.core_reads}="
                        f"{max(target.group_sizes)}+{min(target.group_sizes)}。PS 只作 provenance context，"
                        "不在八個 M2 measured axes 內。"
                    ),
                },
                {"id": f"{target_id}_distance", "type": "chart", "chartId": distance_chart["id"]},
                {"id": f"{target_id}_methyl", "type": "html", "body": image_body, "sourceId": source_id},
                {
                    "id": f"{target_id}_explanation",
                    "type": "markdown",
                    "sourceId": source_id,
                    "body": (
                        "圖中最多每群 15 reads，只是 deterministic medoid-nearest display；"
                        f"表內 group size、distance 與 methyl summary 均使用全部 {target.core_reads} core reads。"
                        "Block pattern 支持共同 focal ALT 下的 read-level residual epigenetic partition，"
                        "不識別兩個 cellular clones 或 lineage order。"
                    ),
                },
            ]
        )

    blocks.extend(
        [
            {
                "id": "scope_definitions",
                "type": "markdown",
                "sourceId": "claim_contract",
                "body": (
                    "## Scope / definitions\n\n"
                    "- Scope：7 datasets / 6 biological samples / chr1-22 / 469,849 latest LongPhase-S "
                    "recalibrated FILTER=PASS autosomal biallelic sSNV dataset-sites。\n"
                    "- Positional singleton：每個 dataset/chrom 內，相鄰 gap `<=50,000 bp` 連接後的"
                    "傳遞 component，其 `component_size=1`。這不是 read-sharing graph degree=0，兩者不保證等價。\n"
                    "- M1：`coarse_ng>=2 AND unstable=false AND modal_assignment_ari_min>=0.8`。\n"
                    "- M2 measured axes：HP exact、HP family、strand、read start/end/length、MAPQ、called CpGs。"
                    "**PS 不在 M2 measured axes 內。**\n"
                    "- HCC1395 與 HCC1395_DORADO 是同一 biological sample 的 technical datasets，"
                    "不可當獨立生物重現。"
                ),
            },
            {"id": "dataset_denominators", "type": "table", "tableId": "dataset_status_table"},
            {"id": "truth_denominators", "type": "table", "tableId": "truth_status_table"},
            {
                "id": "method",
                "type": "markdown",
                "sourceId": "v10_sites",
                "body": (
                    "## Method\n\n"
                    "1. 從 v10 的 469,849 rows 依 dataset/chrom/position 重新切 50 kb transitive components，"
                    "再與 source-attested singleton key set 做 exact equality。\n"
                    "2. 逐列重算 M1 formula，並獨立彙總 M2 四個互斥狀態的 overall/dataset/biological-sample/"
                    "truth/dataset×truth 分母。\n"
                    "3. 從 stable assignment exact 定位兩個 target，將全部 core read IDs 一對一 join "
                    "reads、methylation 與 distance axes；core distance 必須 finite、symmetric、diag=0。\n"
                    "4. 每群 medoid 由全部 group core reads 的平均 within-group distance決定；依到 medoid "
                    "distance與 read ID tie-break 排序，最多顯示 15 reads。沒有隨機抽樣。"
                ),
            },
            {
                "id": "robustness_limits",
                "type": "markdown",
                "sourceId": "claim_contract",
                "body": (
                    "## Robustness / limits\n\n"
                    "- 兩 target 的 core distance、methyl row/column/read identity 已 exact gate；未呼叫 CpG "
                    "保留為 NA，不做 imputation。\n"
                    "- M2 只處理八個已測軸，**不包含 PS**，也不排除所有技術或生物 confound。\n"
                    "- Candidate-level matched-normal methyl locus、tumor-REF specificity、formal sSNV "
                    "cooccurrence、exact-locus CN、purity、multiplicity 與 CCF 尚未在此 derivative sidecar 執行。\n"
                    "- 單一 focal ALT 的多個 methyl groups 只能提出 ancestral-ALT 下 latent molecular "
                    "substructure 候選；不能證明兩個 clone、clone 數或 linear ancestry。\n"
                    "- `confirmed cellular subclone=0` 與 `linear ancestry=0` 表示未通過／未執行必要"
                    "升級 gate，不是已證明不存在。"
                ),
            },
            {"id": "source_inventory", "type": "table", "tableId": "source_inventory_table"},
            {
                "id": "next_steps",
                "type": "markdown",
                "body": (
                    "## Next steps / questions\n\n"
                    "1. 以至少第二個獨立 sSNV marker 測試 group-specific co-membership，並納入完整 multiple-testing family。\n"
                    "2. 在同 locus 完成 tumor-REF 與 candidate-level matched-normal methyl controls。\n"
                    "3. 取得 authority-reviewed exact-locus allele-specific CN、purity、multiplicity 與 fit-local CCF。\n"
                    "4. 以跨技術同位點、multi-region、single-cell/colony/spatial 等正交資料確認 cellular identity；"
                    "沒有這些證據前維持 M2 read-level claim ceiling。"
                ),
            },
        ]
    )

    tables = [
        {
            "id": "case_statistics_table",
            "title": "HCC1395 M2 PASS case statistics",
            "subtitle": "Full-core calculations; display selection does not change denominators.",
            "dataset": "case_statistics",
            "sourceId": "stable_assignments",
            "columns": [
                {"field": "locus", "label": "Locus", "type": "text"},
                {"field": "truth", "label": "Truth", "type": "text"},
                {"field": "raw_focal_alt_reads", "label": "Raw focal ALT", "type": "number"},
                {"field": "core_reads", "label": "Core reads", "type": "number"},
                {"field": "group_a_reads", "label": "Group A", "type": "number"},
                {"field": "group_b_reads", "label": "Group B", "type": "number"},
                {"field": "cpg_sites", "label": "CpGs", "type": "number"},
                {"field": "between_group_mean_distance", "label": "Between distance", "type": "number"},
                {"field": "pooled_within_mean_distance", "label": "Within distance", "type": "number"},
                {"field": "between_to_within_distance_ratio", "label": "Between/within", "type": "number"},
                {"field": "absolute_group_mean_difference", "label": "|mean methyl Δ|", "type": "number"},
                {"field": "core_missing_methyl_fraction", "label": "Missing methyl cells", "type": "percent"},
            ],
            "defaultSort": {"field": "locus", "direction": "asc"},
        },
        {
            "id": "dataset_status_table",
            "title": "Dataset-level singleton denominators",
            "subtitle": "Every M2 state uses all positional-singleton sites in that dataset as the displayed denominator.",
            "dataset": "dataset_status_table",
            "sourceId": "singleton_audit",
            "columns": [
                {"field": "stratum", "label": "Dataset", "type": "text"},
                {"field": "sites", "label": "Singletons", "type": "number"},
                {"field": "m1_evaluable", "label": "M1 evaluable", "type": "number"},
                {"field": "m1_flagged", "label": "M1 flagged", "type": "number"},
                {"field": "m2_not_run", "label": "M2 NOT_RUN", "type": "number"},
                {"field": "m2_not_evaluable", "label": "M2 NOT_EVALUABLE", "type": "number"},
                {"field": "m2_fail", "label": "M2 FAIL", "type": "number"},
                {"field": "m2_pass", "label": "M2 PASS", "type": "number"},
                {"field": "m2_pass_pct_all", "label": "PASS / all", "type": "percent"},
            ],
            "defaultSort": {"field": "sites", "direction": "desc"},
        },
        {
            "id": "truth_status_table",
            "title": "Truth-level singleton denominators",
            "subtitle": "Truth strata retain NOT_RUN and NOT_EVALUABLE; PASS/all is operational yield.",
            "dataset": "truth_status_table",
            "sourceId": "singleton_audit",
            "columns": [
                {"field": "stratum", "label": "Truth", "type": "text"},
                {"field": "sites", "label": "Singletons", "type": "number"},
                {"field": "m1_evaluable", "label": "M1 evaluable", "type": "number"},
                {"field": "m1_flagged", "label": "M1 flagged", "type": "number"},
                {"field": "m2_not_run", "label": "M2 NOT_RUN", "type": "number"},
                {"field": "m2_not_evaluable", "label": "M2 NOT_EVALUABLE", "type": "number"},
                {"field": "m2_fail", "label": "M2 FAIL", "type": "number"},
                {"field": "m2_pass", "label": "M2 PASS", "type": "number"},
                {"field": "m2_pass_pct_all", "label": "PASS / all", "type": "percent"},
            ],
            "defaultSort": {"field": "sites", "direction": "desc"},
        },
        {
            "id": "source_inventory_table",
            "title": "Frozen source inventory",
            "subtitle": "Reader-safe file identities; exact absolute paths are retained in source_inventory.tsv/json.",
            "dataset": "public_source_inventory",
            "sourceId": "v10_sites",
            "columns": [
                {"field": "role", "label": "Role", "type": "text"},
                {"field": "file", "label": "File", "type": "text"},
                {"field": "sha256", "label": "SHA-256", "type": "text"},
                {"field": "size_bytes", "label": "Bytes", "type": "number"},
                {"field": "used_for", "label": "Used for", "type": "text"},
            ],
            "defaultSort": {"field": "role", "direction": "asc"},
        },
    ]

    snapshot_datasets: dict[str, list[dict[str, Any]]] = {
        "dataset_m2_status": [dict(row) for row in dataset_status],
        "hcc1395_funnel": hcc_funnel,
        "dataset_status_table": dataset_rows,
        "truth_status_table": truth_rows,
        "case_statistics": [dict(row) for row in case_stats],
        "public_source_inventory": public_source_rows,
    }
    for target in TARGETS:
        snapshot_datasets[f"{target.short_id}_distance_cells"] = [
            dict(row) for row in distance_rows_by_target[target.short_id]
        ]
    snapshot_datasets, sql_queries = stage_snapshot_datasets(staging_sqlite_path, snapshot_datasets)
    staging_sha256 = sha256_path(staging_sqlite_path)
    attach_executed_sql_sources(charts, sources, sql_queries, staging_sha256)
    attach_executed_sql_sources(tables, sources, sql_queries, staging_sha256)

    payload = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": REPORT_TITLE,
            "description": "Task Type B derivative validation with full singleton denominators and HCC1395 raw-matrix joins.",
            "generatedAt": created_at,
            "blocks": blocks,
            "charts": charts,
            "tables": tables,
            "cards": [],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "status": "ready",
            "generatedAt": created_at,
            "datasets": snapshot_datasets,
        },
        "sources": sources,
        "package_info": {
            "mode": "portable_html",
            "deliveryMode": "portable_html",
            "readOnly": True,
            "taskType": "B_comprehensive_validation",
        },
    }
    return payload, {
        "schema_name": "intersubmod.singleton_alt_methyl_artifact_sql_provenance",
        "schema_version": "1.0.0",
        "actual_execution": True,
        "database_path": staging_sqlite_path.name,
        "database_sha256": staging_sha256,
        "sqlite_version": sqlite3.sqlite_version,
        "queries": sql_queries,
    }


def validate_artifact_local_contract(payload: Mapping[str, Any]) -> dict[str, Any]:
    manifest = payload["manifest"]
    snapshot = payload["snapshot"]
    assert_equal(payload["surface"], "report", "artifact surface")
    assert_equal(snapshot["status"], "ready", "snapshot status")
    assert_true(bool(manifest["blocks"]), "artifact blocks required")
    assert_equal(manifest["blocks"][0]["type"], "markdown", "first block type")
    assert_true(
        manifest["blocks"][0]["body"].startswith(f"# {manifest['title']}"),
        "first markdown H1 must match manifest title",
    )
    assert_true(any(block["type"] == "chart" for block in manifest["blocks"]), "at least one chart block")
    datasets = snapshot["datasets"]
    assert_true(len(datasets) <= 50, "snapshot dataset count exceeds 50")
    for dataset_id, rows in datasets.items():
        assert_true(isinstance(rows, list), f"snapshot dataset {dataset_id} must be an array")
        assert_true(len(rows) <= 2_000, f"snapshot dataset {dataset_id} exceeds 2,000 rows")
        assert_true(all(isinstance(row, dict) for row in rows), f"snapshot dataset {dataset_id} row shape")
    for table in manifest["tables"]:
        assert_true("defaultSort" in table, f"table {table['id']} missing defaultSort")
        fields = {column["field"] for column in table["columns"]}
        assert_true(table["defaultSort"]["field"] in fields, f"table {table['id']} invalid defaultSort")
    for item in [*manifest.get("charts", []), *manifest.get("tables", [])]:
        sql = item.get("source", {}).get("query", {}).get("sql")
        assert_true(
            isinstance(sql, str) and sql.lstrip().upper().startswith("SELECT "),
            f"native item {item['id']} lacks executed SQL provenance",
        )
    serialized = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    assert_true(len(serialized.encode("utf-8")) <= 3 * 1024 * 1024, "artifact payload exceeds 3 MiB")
    forbidden = ("/big7_disk/", "/bip7_disk/", "file://", "../")
    for token in forbidden:
        assert_true(token not in serialized, f"artifact contains unsafe machine-local token {token}")
    visible = "\n".join(
        block.get("body", "")
        for block in manifest["blocks"]
        if block["type"] in {"markdown", "html"}
    )
    required_claims = (
        "M2 measured axes",
        "PS 不在 M2 measured axes",
        "observed operational M2-PASS yield",
        "不是 biological/subclone prevalence",
        "confirmed cellular subclone=0",
        "linear ancestry=0",
        "不是真陰性",
        "matched-normal",
        "tumor-REF",
        "CN",
        "CCF",
        "Group A/B",
    )
    for claim in required_claims:
        assert_true(claim in visible, f"artifact missing claim guardrail: {claim}")
    assert_true(
        not contains_source_cluster_label(serialized),
        "artifact exposes source labels that can be mistaken for HP",
    )
    return {
        "datasets": len(datasets),
        "rows": sum(len(rows) for rows in datasets.values()),
        "blocks": len(manifest["blocks"]),
        "charts": len(manifest["charts"]),
        "tables": len(manifest["tables"]),
        "html_blocks": sum(block["type"] == "html" for block in manifest["blocks"]),
        "payload_bytes": len(serialized.encode("utf-8")),
        "first_markdown_matches_title": True,
        "all_tables_have_default_sort": True,
        "native_items_have_executed_sql_source": True,
        "machine_local_paths_absent": True,
    }


def chart_map_rows(case_stats: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = [
        {
            "figure_id": "all_dataset_m2_status_chart",
            "section": "Findings / visual evidence",
            "question": "How do all singleton sites partition across M2 states by dataset?",
            "family": "composition",
            "type": "horizontalStackedBar",
            "fields": "dataset,status,count",
            "calculation_population": "all 50,432 positional-singleton dataset-sites",
            "display_population": "all seven dataset aggregates",
            "supported_takeaway": "NOT_RUN/NOT_EVALUABLE/FAIL/PASS counts are conserved by dataset",
            "cannot_support": "biological prevalence or true-negative rate",
            "palette_policy": "categorical status palette plus exact table",
            "delivery": "native artifact chart",
        },
        {
            "figure_id": "hcc1395_funnel_chart",
            "section": "Findings / visual evidence",
            "question": "How does the HCC1395 population narrow through positional, M1, and M2 gates?",
            "family": "funnel",
            "type": "funnel",
            "fields": "stage,count",
            "calculation_population": "all HCC1395 v10 dataset-sites and all 8,279 singletons",
            "display_population": "six frozen stages",
            "supported_takeaway": "8,279 singleton, 734 M1 flags, 2 M2 PASS",
            "cannot_support": "prevalence or biological absence among NOT_RUN/NOT_EVALUABLE",
            "palette_policy": "sequential",
            "delivery": "native artifact chart",
        },
    ]
    stats_by_locus = {row["locus"]: row for row in case_stats}
    for target in TARGETS:
        stats = stats_by_locus[target.display_locus]
        rows.extend(
            [
                {
                    "figure_id": f"{target.short_id}_distance_heatmap",
                    "section": target.display_locus,
                    "question": "Do representative reads retain within/between methyl-distance blocks?",
                    "family": "relationship",
                    "type": "heatmap",
                    "fields": "row_read,column_read,distance",
                    "calculation_population": f"all {target.core_reads} core reads for medoids/statistics",
                    "display_population": (
                        f"{stats['display_group_a_reads']} Group A + "
                        f"{stats['display_group_b_reads']} Group B deterministic representatives"
                    ),
                    "supported_takeaway": "read-overlap molecular distance pattern",
                    "cannot_support": "HP identity, cellular clone, or lineage",
                    "palette_policy": "sequential distance",
                    "delivery": "native artifact heatmap",
                },
                {
                    "figure_id": f"{target.short_id}_methyl_heatmap",
                    "section": target.display_locus,
                    "question": "What CpG-level methylation pattern is visible across deterministic representatives?",
                    "family": "matrix",
                    "type": "static PNG heatmap",
                    "fields": "representative_read,CpG_position,methyl_probability,missing",
                    "calculation_population": f"all {target.core_reads} core reads for summary statistics",
                    "display_population": (
                        f"{stats['display_group_a_reads']} Group A + "
                        f"{stats['display_group_b_reads']} Group B deterministic representatives"
                    ),
                    "supported_takeaway": "read-by-CpG pattern with explicit uncalled cells",
                    "cannot_support": "cellular identity, ancestry, or imputed methyl state",
                    "palette_policy": "cividis sequential; gray=uncalled",
                    "delivery": f"artifact-embedded PNG: {stats['methyl_figure_path']}",
                },
            ]
        )
    return rows


def build(
    audit_root: Path,
    v10_root: Path,
    claim_contract: Path,
    output_dir: Path,
    command: Sequence[str],
) -> dict[str, Any]:
    created_at = utc_now()
    figures_dir = output_dir / "figures"
    figures_dir.mkdir()

    summary_path = audit_root / "positional_singleton_audit_summary.json"
    site_audit_path = audit_root / "positional_singleton_site_audit.tsv.gz"
    pass_cases_path = audit_root / "positional_singleton_m2_pass_cases.tsv"
    success_path = audit_root / "_SUCCESS.json"
    v10_site_path = v10_root / "all_ssnv_site_results.tsv.gz"
    assignments_path = v10_root / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    run_manifest_path = v10_root / "run_manifest.json"

    source_inventory = [
        source_record(summary_path, "audit_summary", "frozen count and source-chain comparison", EXPECTED_SHA256["audit_summary"]),
        source_record(site_audit_path, "site_audit", "M1/M2 census and exact singleton key set", EXPECTED_SHA256["site_audit"]),
        source_record(pass_cases_path, "m2_pass_cases", "frozen M2 PASS case identity", EXPECTED_SHA256["m2_pass_cases"]),
        source_record(success_path, "audit_success", "attestation marker and output digest binding"),
        source_record(v10_site_path, "v10_site_results", "469,849-site positional recomputation and target region lookup", EXPECTED_SHA256["v10_site_results"]),
        source_record(assignments_path, "stable_assignments", "exact target core read assignments"),
        source_record(run_manifest_path, "v10_run_manifest", "v10 source-locked run identity"),
        source_record(claim_contract, "claim_contract", "M1/M2 definitions and claim ceiling"),
    ]

    success = json.loads(success_path.read_text(encoding="utf-8"))
    assert_equal(success.get("pass"), True, "audit _SUCCESS pass")
    assert_equal(success.get("summary", {}).get("sha256"), EXPECTED_SHA256["audit_summary"], "_SUCCESS summary hash")
    assert_equal(success.get("site_audit_sha256"), EXPECTED_SHA256["site_audit"], "_SUCCESS site-audit hash")
    assert_equal(success.get("m2_pass_cases_sha256"), EXPECTED_SHA256["m2_pass_cases"], "_SUCCESS pass-cases hash")

    frozen_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert_equal(frozen_summary.get("pass"), True, "frozen audit pass")
    assert_equal(
        frozen_summary.get("contracts", {}).get("positional_singleton"),
        "same_dataset_chrom_adjacent_gap_gt_50000_component_size_1_v1",
        "positional singleton contract",
    )
    assert_equal(
        frozen_summary.get("contracts", {}).get("claim_ceiling"),
        "M2_read_level_residual_epigenetic_partition",
        "claim ceiling",
    )

    v10_sites, singleton_keys = load_v10_site_results(v10_site_path)
    audit = load_site_audit(site_audit_path, singleton_keys)
    census_wide, census_long = build_census_breakdowns(audit)
    overall = next(row for row in census_wide if row["stratum_type"] == "overall")
    fixed_overall = {
        "all_dataset_sites": len(v10_sites),
        "singleton_sites": overall["sites"],
        "m1_evaluable": overall["m1_evaluable"],
        "m1_flagged": overall["m1_flagged"],
        "m2_not_run": overall["m2_not_run"],
        "m2_not_evaluable": overall["m2_not_evaluable"],
        "m2_fail": overall["m2_fail"],
        "m2_pass": overall["m2_pass"],
    }
    assert_equal(fixed_overall, EXPECTED_COUNTS, "fixed overall counts")
    hcc = next(row for row in census_wide if row["stratum_type"] == "dataset" and row["stratum"] == "HCC1395")
    assert_equal(
        {key: hcc[key] for key in EXPECTED_HCC1395},
        EXPECTED_HCC1395,
        "fixed HCC1395 counts",
    )
    compare_recomputed_to_frozen_summary(overall, census_wide, frozen_summary)

    v10_targets = find_target_rows(v10_sites, TARGETS)
    audit_targets = find_target_rows(audit, TARGETS)
    assignments = load_target_assignments(assignments_path, TARGETS)

    all_case_stats: list[dict[str, Any]] = []
    all_selection_rows: list[dict[str, Any]] = []
    all_heatmap_cells: list[dict[str, Any]] = []
    all_join_checks: list[dict[str, Any]] = []
    distance_rows_by_target: dict[str, list[dict[str, Any]]] = {}
    figure_paths_by_target: dict[str, Path] = {}
    for target in TARGETS:
        (
            case_stats,
            selection_rows,
            distance_rows,
            heatmap_cells,
            join_checks,
            target_sources,
        ) = process_target(
            target,
            v10_targets[target.key],
            audit_targets[target.key],
            assignments[target.key],
            figures_dir,
        )
        all_case_stats.append(case_stats)
        all_selection_rows.extend(selection_rows)
        all_heatmap_cells.extend(heatmap_cells)
        all_join_checks.extend(join_checks)
        distance_rows_by_target[target.short_id] = distance_rows
        figure_paths_by_target[target.short_id] = figures_dir / f"{target.short_id}.read_by_cpg_methylation.png"
        source_inventory.extend(target_sources)

    source_inventory_path = output_dir / "source_inventory.tsv"
    safe_tsv_dump(
        source_inventory_path,
        source_inventory,
        ("role", "path", "size_bytes", "mode", "sha256", "used_for"),
    )
    safe_json_dump(output_dir / "source_inventory.json", source_inventory)

    census_fields = (
        "stratum_type",
        "stratum",
        "dataset",
        "truth",
        "sites",
        "m1_evaluable",
        "m1_not_evaluable",
        "m1_flagged",
        "m1_not_flagged",
        "m2_not_run",
        "m2_not_evaluable",
        "m2_fail",
        "m2_pass",
        "m2_determinate",
        "m1_evaluable_pct_all",
        "m1_flagged_pct_all",
        "m1_flagged_pct_evaluable",
        "m2_determinate_pct_m1",
        "m2_pass_pct_all",
        "m2_pass_pct_m1",
        "m2_pass_pct_determinate",
    )
    safe_tsv_dump(output_dir / "census_breakdown.tsv", census_wide, census_fields)
    safe_tsv_dump(
        output_dir / "status_denominators.tsv",
        census_long,
        (
            "stratum_type",
            "stratum",
            "dataset",
            "truth",
            "status_axis",
            "status",
            "count",
            "denominator",
            "percent_of_singleton_stratum",
        ),
    )
    safe_tsv_dump(
        output_dir / "case_statistics.tsv",
        all_case_stats,
        tuple(all_case_stats[0].keys()),
    )
    safe_tsv_dump(
        output_dir / "target_join_audit.tsv",
        all_join_checks,
        tuple(all_join_checks[0].keys()),
    )
    safe_tsv_dump(
        output_dir / "representative_selection_audit.tsv",
        all_selection_rows,
        (
            "locus",
            "group",
            "read_id",
            "read_name_sha256",
            "medoid_read_id",
            "distance_to_medoid",
            "rank_from_medoid",
            "selected_for_display",
            "display_label",
            "selection_cap_per_group",
            "selection_population",
        ),
    )
    safe_gzip_tsv_dump(
        output_dir / "visual_heatmap_cells.tsv.gz",
        all_heatmap_cells,
        ("locus", "heatmap_type", "row_label", "column_label", "value", "missing"),
    )

    chart_rows = chart_map_rows(all_case_stats)
    safe_tsv_dump(output_dir / "chart_map.tsv", chart_rows, tuple(chart_rows[0].keys()))
    safe_json_dump(output_dir / "chart_map.json", chart_rows)

    machine_summary = {
        "schema_name": "intersubmod.singleton_alt_methyl_substructure_derivative_validation",
        "schema_version": "1.0.0",
        "created_at_utc": created_at,
        "task_type": "B_comprehensive_validation",
        "pass": True,
        "scope": {
            "datasets": list(DATASET_ORDER),
            "biological_samples": sorted(audit["biological_id"].unique()),
            "chromosomes": "chr1-22",
            "all_dataset_sites": len(v10_sites),
            "positional_singleton_definition": (
                "within dataset/chrom, adjacent gap <=50,000 bp joins a transitive component; component_size=1"
            ),
            "read_sharing_degree_zero_equivalent": False,
        },
        "counts": fixed_overall,
        "rates": {
            "singleton_fraction_all_sites": overall["sites"] / len(v10_sites),
            "m1_evaluable_fraction_singletons": overall["m1_evaluable_pct_all"],
            "m1_flagged_fraction_singletons": overall["m1_flagged_pct_all"],
            "m1_flagged_fraction_evaluable": overall["m1_flagged_pct_evaluable"],
            "m2_determinate_fraction_m1": overall["m2_determinate_pct_m1"],
            "m2_pass_operational_yield_all_singletons": overall["m2_pass_pct_all"],
            "m2_pass_fraction_determinate_conditional_only": overall["m2_pass_pct_determinate"],
        },
        "hcc1395": {key: hcc[key] for key in EXPECTED_HCC1395},
        "target_cases": all_case_stats,
        "checks": {
            "source_sha_attestation": True,
            "recomputed_469849_rows": True,
            "recomputed_singleton_key_set_exact": True,
            "m1_formula_exact": True,
            "all_statuses_conserved": True,
            "frozen_summary_breakdowns_exact": True,
            "target_assignment_exact": True,
            "target_core_join_108_of_108": True,
            "target_core_join_109_of_109": True,
            "core_distance_finite": True,
            "core_distance_symmetric": True,
            "core_distance_diagonal_zero": True,
            "methyl_read_and_cpg_identity": True,
            "visual_selection_deterministic": True,
            "groups_exposed_only_as_group_a_b": True,
        },
        "claim_ceiling": "M2_read_level_residual_epigenetic_partition",
        "claim_guardrails": {
            "m2_pass_is_operational_yield_not_prevalence": True,
            "m2_measured_axes_include_ps": False,
            "confirmed_cellular_subclone": 0,
            "confirmed_cellular_subclone_zero_semantics": "unverified_not_true_negative",
            "linear_ancestry": 0,
            "linear_ancestry_zero_semantics": "unverified_not_true_negative",
        },
        "downstream_not_completed": [
            "formal local sSNV cooccurrence and multi-marker association",
            "candidate-level matched-normal methyl control",
            "tumor-REF methyl specificity control",
            "authority-reviewed exact-locus allele-specific CN",
            "purity, multiplicity, and fit-local CCF",
            "orthogonal cellular identity and lineage replication",
        ],
    }
    safe_json_dump(output_dir / "machine_summary.json", machine_summary)

    artifact, sql_provenance = build_artifact(
        census_wide,
        census_long,
        all_case_stats,
        distance_rows_by_target,
        figure_paths_by_target,
        source_inventory,
        created_at,
        output_dir / "artifact_staging.sqlite",
    )
    safe_json_dump(output_dir / "artifact_sql_provenance.json", sql_provenance)
    artifact_profile = validate_artifact_local_contract(artifact)
    artifact_path = output_dir / "artifact.json"
    safe_json_dump(artifact_path, artifact)

    files_before_receipt = sorted(path for path in output_dir.rglob("*") if path.is_file())
    build_receipt = {
        "schema_name": "intersubmod.singleton_alt_methyl_substructure_build_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": created_at,
        "ok": True,
        "command": list(command),
        "inputs": {
            "audit_root": str(audit_root.resolve()),
            "v10_root": str(v10_root.resolve()),
            "claim_contract": str(claim_contract.resolve()),
        },
        "output_dir": str(output_dir.resolve()),
        "counts": fixed_overall,
        "hcc1395": {key: hcc[key] for key in EXPECTED_HCC1395},
        "target_join_checks": all_join_checks,
        "artifact_profile": artifact_profile,
        "output_sha256": {
            str(path.relative_to(output_dir)): sha256_path(path)
            for path in files_before_receipt
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "sqlite": sqlite3.sqlite_version,
        },
        "claim_ceiling": "M2_read_level_residual_epigenetic_partition",
    }
    safe_json_dump(output_dir / "build_receipt.json", build_receipt)

    for path in output_dir.rglob("*"):
        if path.is_file():
            readonly(path)
    return {
        "ok": True,
        "output_dir": str(output_dir.resolve()),
        "artifact": str(artifact_path.resolve()),
        "artifact_sha256": sha256_path(artifact_path),
        "counts": fixed_overall,
        "hcc1395": {key: hcc[key] for key in EXPECTED_HCC1395},
        "target_joins": [
            {
                "locus": row["locus"],
                "core_expected": row["core_expected"],
                "core_in_reads": row["core_in_reads"],
                "core_in_methylation": row["core_in_methylation"],
                "core_in_distance_rows": row["core_in_distance_rows"],
            }
            for row in all_join_checks
        ],
        "artifact_profile": artifact_profile,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-root", type=Path, default=DEFAULT_AUDIT_ROOT)
    parser.add_argument("--v10-root", type=Path, default=DEFAULT_V10_ROOT)
    parser.add_argument("--claim-contract", type=Path, default=DEFAULT_CLAIM_CONTRACT)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    command = [sys.executable, str(Path(__file__).resolve()), *(argv if argv is not None else sys.argv[1:])]
    ensure_output_dir(args.output_dir)
    try:
        result = build(
            args.audit_root.resolve(),
            args.v10_root.resolve(),
            args.claim_contract.resolve(),
            args.output_dir.resolve(),
            command,
        )
    except Exception as exc:
        failure_path = args.output_dir / "failure_receipt.json"
        if not failure_path.exists():
            safe_json_dump(
                failure_path,
                {
                    "schema_name": "intersubmod.singleton_alt_methyl_substructure_failure_receipt",
                    "schema_version": "1.0.0",
                    "created_at_utc": utc_now(),
                    "ok": False,
                    "command": command,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                },
            )
            readonly(failure_path)
        print(
            json.dumps(
                {
                    "ok": False,
                    "output_dir": str(args.output_dir.resolve()),
                    "failure_receipt": str(failure_path.resolve()),
                    "error": str(exc),
                },
                ensure_ascii=False,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 1
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
