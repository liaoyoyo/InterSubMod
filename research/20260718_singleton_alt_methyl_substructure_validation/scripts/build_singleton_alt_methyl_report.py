#!/usr/bin/env python3
"""Build the complete singleton focal-ALT methyl-substructure report dataset.

The script is intentionally deterministic. It recomputes the HCC1395 positional
singleton denominator, verifies the source-attested audit, joins the two HCC1395
M2-PASS loci to exact focal-ALT core reads, and emits a canonical Data Analytics
artifact plus machine-readable receipts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import os
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUTPUT_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = OUTPUT_ROOT / "positional_singleton_methyl_multigroup_audit_v1_source_attested"
SCREEN_ROOT = OUTPUT_ROOT / "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
TOPIC = REPO / "research/20260718_singleton_alt_methyl_substructure_validation"

SUMMARY_PATH = AUDIT_ROOT / "positional_singleton_audit_summary.json"
SUCCESS_PATH = AUDIT_ROOT / "_SUCCESS.json"
AUDIT_TSV = AUDIT_ROOT / "positional_singleton_site_audit.tsv.gz"
M2_PASS_TSV = AUDIT_ROOT / "positional_singleton_m2_pass_cases.tsv"
SITE_RESULTS = SCREEN_ROOT / "all_ssnv_site_results.tsv.gz"
ASSIGNMENTS = SCREEN_ROOT / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
HCC_VCF = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/"
    "HCC1395.longphase_s.recalibrated.pass.vcf.gz"
)

GENCODE_GTF = Path("/big7_disk/liaoyoyo2001/gene_annotation/gencode.v46.basic.annotation.gtf.gz")
CGC_TSV = Path("/big7_disk/liaoyoyo2001/gene_annotation/Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz")
DGIDB_TSV = Path("/big7_disk/liaoyoyo2001/gene_annotation/dgidb_interactions.tsv")

TARGETS = {
    ("HCC1395", "chr14", 86272476, "A", "T"),
    ("HCC1395", "chr22", 47466517, "A", "G"),
}
EXPECTED_DATASET = {
    "COLO829": (7830, 7347, 617, 0, 0, 617, 7213),
    "H1437": (6696, 6298, 813, 2, 2, 809, 5883),
    "H2009": (2853, 2782, 602, 9, 1, 592, 2251),
    "HCC1395": (8279, 8074, 734, 2, 0, 732, 7545),
    "HCC1395_DORADO": (8321, 8016, 962, 2, 1, 959, 7359),
    "HCC1937": (8469, 8069, 843, 14, 9, 820, 7626),
    "HCC1954": (7984, 7761, 1390, 1, 5, 1384, 6594),
}
DISPLAY_PER_GROUP = 4
DISPLAY_CPGS = 8


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def finite_or_none(value: Any, digits: int = 6) -> Any:
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return value
    if not math.isfinite(f):
        return None
    return round(f, digits)


def quantiles(values: Sequence[float]) -> Dict[str, Any]:
    a = np.asarray([x for x in values if math.isfinite(float(x))], dtype=float)
    if not len(a):
        return {"n": 0, "median": None, "q1": None, "q3": None, "iqr": None, "mean": None}
    q1, median, q3 = np.quantile(a, [0.25, 0.5, 0.75], method="linear")
    return {
        "n": int(len(a)),
        "median": finite_or_none(median),
        "q1": finite_or_none(q1),
        "q3": finite_or_none(q3),
        "iqr": finite_or_none(q3 - q1),
        "mean": finite_or_none(np.mean(a)),
    }


def wilson_interval(success: int, total: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if total <= 0:
        return (math.nan, math.nan)
    p = success / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denom
    return center - half, center + half


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def check_display_value(value: Any, limit: int = 500) -> str:
    if isinstance(value, str):
        text = value
    else:
        text = json.dumps(value, ensure_ascii=False, default=str)
    if re.fullmatch(r"[0-9a-f]{64}", text):
        return f"{text[:16]}…{text[-8:]}"
    if len(text) <= limit:
        return text
    count = len(value) if hasattr(value, "__len__") else "n/a"
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return f"{type(value).__name__} n={count}; sha256={digest}"


def assert_equal(actual: Any, expected: Any, label: str, checks: List[Dict[str, Any]]) -> None:
    if isinstance(actual, np.generic):
        actual = actual.item()
    if isinstance(expected, np.generic):
        expected = expected.item()
    passed = actual == expected
    if isinstance(passed, np.generic):
        passed = bool(passed.item())
    checks.append({"check": label, "actual": actual, "expected": expected, "pass": passed})
    if not passed:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def recompute_vcf_components(path: Path) -> Dict[str, int]:
    records: List[Tuple[str, int]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t", 5)
            chrom, pos_text, ref, alt = cols[:4]
            match = re.fullmatch(r"chr([0-9]+)", chrom)
            if not match or not 1 <= int(match.group(1)) <= 22:
                continue
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            records.append((chrom, int(pos_text)))
    singleton = 0
    multi_components = 0
    multi_sites = 0
    group_size = 0
    prior_chrom = ""
    prior_pos = -1

    def finish(size: int) -> Tuple[int, int, int]:
        if size == 1:
            return (1, 0, 0)
        if size > 1:
            return (0, 1, size)
        return (0, 0, 0)

    for chrom, pos in records:
        if group_size == 0:
            prior_chrom, prior_pos, group_size = chrom, pos, 1
            continue
        if chrom != prior_chrom or pos - prior_pos > 50000:
            s, c, n = finish(group_size)
            singleton += s
            multi_components += c
            multi_sites += n
            group_size = 1
        else:
            group_size += 1
        prior_chrom, prior_pos = chrom, pos
    s, c, n = finish(group_size)
    singleton += s
    multi_components += c
    multi_sites += n
    return {
        "autosomal_biallelic": len(records),
        "positional_singleton": singleton,
        "multilocus_components": multi_components,
        "multilocus_sites": multi_sites,
        "all_positional_components": singleton + multi_components,
        "partition_sites": singleton + multi_sites,
    }


def read_target_site_rows() -> pd.DataFrame:
    usecols = [
        "dataset",
        "sample",
        "biological_id",
        "truth_label",
        "chrom",
        "pos",
        "ref",
        "alt",
        "qual",
        "filter",
        "caller_gt",
        "caller_gq",
        "caller_dp",
        "caller_af",
        "caller_ad",
        "normal_af",
        "normal_dp",
        "normal_ad",
        "component_id",
        "component_size",
        "region_dir",
        "n_reads_total",
        "n_alt_raw",
        "n_alt_matrix",
        "n_alt_after_peel",
        "n_alt_peeled",
        "alt_readset_sha256",
        "alt_hp_counts",
        "alt_hp_family_counts",
        "alt_ps_counts",
        "alt_strand_counts",
        "coarse_ng",
        "modal_assignment_ari_min",
        "unstable",
        "cluster_sizes",
        "hp_axis_confound",
        "technical_axis_confound",
        "residual_unexplained_multigroup",
        "strict_methyl_partition_robustness_status",
        "strict_confirm_status",
        "strict_confirm_candidate_is_formal_r1_claim",
        "strand_v",
        "strand_p_perm",
        "start_eta2",
        "start_p_perm",
        "end_eta2",
        "end_p_perm",
        "length_eta2",
        "length_p_perm",
        "mapq_eta2",
        "mapq_p_perm",
        "cpg_called_eta2",
        "cpg_called_p_perm",
    ]
    frames: List[pd.DataFrame] = []
    for chunk in pd.read_csv(
        SITE_RESULTS, sep="\t", usecols=usecols, chunksize=50000, low_memory=False
    ):
        keep = chunk["dataset"].eq("HCC1395") & (
            (
                chunk["chrom"].eq("chr14")
                & chunk["pos"].eq(86272476)
                & chunk["ref"].eq("A")
                & chunk["alt"].eq("T")
            )
            | (
                chunk["chrom"].eq("chr22")
                & chunk["pos"].eq(47466517)
                & chunk["ref"].eq("A")
                & chunk["alt"].eq("G")
            )
        )
        if keep.any():
            frames.append(chunk.loc[keep].copy())
    if not frames:
        raise RuntimeError("No target site rows found")
    result = pd.concat(frames, ignore_index=True)
    if len(result) != 2:
        raise AssertionError(f"Expected exactly two target site rows, got {len(result)}")
    return result


def read_target_assignments() -> Dict[Tuple[str, str, int, str, str], Dict[str, Any]]:
    found: Dict[Tuple[str, str, int, str, str], Dict[str, Any]] = {}
    short_to_full = {(key[0], key[1], key[2]): key for key in TARGETS}
    with gzip.open(ASSIGNMENTS, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            short_key = (
                row.get("dataset"),
                row.get("chrom"),
                int(row.get("pos")),
            )
            if short_key in short_to_full:
                found[short_to_full[short_key]] = row
                if len(found) == len(TARGETS):
                    break
    if set(found) != TARGETS:
        raise AssertionError(f"Target assignment keys mismatch: {sorted(found)}")
    return found


def read_numeric_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"read_id": str}, na_values=["NA", "NaN", "nan", ""])
    frame["read_id"] = frame["read_id"].astype(str)
    return frame.set_index("read_id")


def select_representatives(
    distance: pd.DataFrame, group_for_read: Mapping[str, str], group_names: Sequence[str]
) -> List[str]:
    selected: List[str] = []
    for group in group_names:
        ids = [rid for rid, g in group_for_read.items() if g == group]
        sub = distance.loc[ids, ids].to_numpy(dtype=float)
        np.fill_diagonal(sub, np.nan)
        scores = np.nanmean(sub, axis=1)
        ranked = sorted(zip(ids, scores), key=lambda x: (math.inf if not math.isfinite(x[1]) else x[1], x[0]))
        selected.extend([rid for rid, _ in ranked[: min(DISPLAY_PER_GROUP, len(ranked))]])
    return selected


def gtf_attributes(text: str) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        values[key] = value.strip().strip('"')
    return values


def gene_context(targets: Sequence[Tuple[str, int]]) -> Dict[Tuple[str, int], Dict[str, Any]]:
    wanted_chroms = {chrom for chrom, _ in targets}
    genes: Dict[str, List[Dict[str, Any]]] = {chrom: [] for chrom in wanted_chroms}
    with gzip.open(GENCODE_GTF, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[0] not in wanted_chroms or cols[2] != "gene":
                continue
            attrs = gtf_attributes(cols[8])
            genes[cols[0]].append(
                {
                    "gene_name": attrs.get("gene_name", attrs.get("gene_id", "unknown")),
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_type": attrs.get("gene_type", ""),
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "strand": cols[6],
                }
            )
    cgc = pd.read_csv(CGC_TSV, sep="\t", dtype=str)
    cgc_symbols = set(cgc["GENE_SYMBOL"].fillna("").str.upper())
    dgidb = pd.read_csv(DGIDB_TSV, sep="\t", dtype=str)
    antineoplastic = dgidb["anti_neoplastic"].fillna("").str.upper().eq("TRUE")
    approved = dgidb["approved"].fillna("").str.upper().eq("TRUE")
    dgidb_symbols = set(dgidb.loc[antineoplastic & approved, "gene_name"].fillna("").str.upper())
    result: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for chrom, pos in targets:
        candidates = []
        for gene in genes[chrom]:
            if gene["start"] <= pos <= gene["end"]:
                distance = 0
                relationship = "gene-body overlap"
            elif pos < gene["start"]:
                distance = gene["start"] - pos
                relationship = "nearest downstream"
            else:
                distance = pos - gene["end"]
                relationship = "nearest upstream"
            candidates.append((distance, gene["gene_name"], relationship, gene))
        distance, _, relationship, gene = min(
            candidates,
            key=lambda item: (
                item[0],
                item[1],
                item[3].get("gene_id", ""),
                item[3]["start"],
                item[3]["end"],
            ),
        )
        symbol = str(gene["gene_name"]).upper()
        result[(chrom, pos)] = {
            **gene,
            "relationship": relationship,
            "distance_bp": int(distance),
            "cgc_member": symbol in cgc_symbols,
            "approved_antineoplastic_dgidb": symbol in dgidb_symbols,
        }
    return result


def analyze_locus(
    key: Tuple[str, str, int, str, str],
    assignment: Dict[str, Any],
    site: Mapping[str, Any],
    audit_row: Mapping[str, Any],
    gene: Mapping[str, Any],
    checks: List[Dict[str, Any]],
) -> Dict[str, Any]:
    dataset, chrom, pos, ref, alt = key
    locus = f"{chrom}:{pos:,} {ref}>{alt}"
    artifacts = assignment["primary_artifacts"]
    reads_path = Path(artifacts["reads"]["path"])
    methyl_path = Path(artifacts["methylation_matrix"]["path"])
    distance_path = Path(artifacts["distance_matrix"]["path"])
    for name, path in [("reads", reads_path), ("methyl", methyl_path), ("distance", distance_path)]:
        actual = sha256_file(path)
        expected = artifacts[
            {"reads": "reads", "methyl": "methylation_matrix", "distance": "distance_matrix"}[name]
        ]["sha256"]
        assert_equal(actual, expected, f"{locus} {name} SHA-256", checks)

    reads = pd.read_csv(reads_path, sep="\t", dtype={"read_id": str})
    reads["read_id"] = reads["read_id"].astype(str)
    reads = reads.set_index("read_id")
    methyl = read_numeric_csv(methyl_path)
    distance = read_numeric_csv(distance_path)
    distance.columns = distance.columns.astype(str)
    core = pd.DataFrame(assignment["core_reads"])
    core["read_id"] = core["read_id"].astype(str)
    core_ids = core["read_id"].tolist()
    assert_equal(len(core_ids), int(audit_row["core_read_n"]), f"{locus} core count", checks)
    assert_equal(len(set(core_ids)), len(core_ids), f"{locus} core read IDs unique", checks)
    assert_equal(set(core_ids).issubset(set(reads.index)), True, f"{locus} reads join complete", checks)
    assert_equal(set(core_ids).issubset(set(methyl.index)), True, f"{locus} methyl join complete", checks)
    assert_equal(set(core_ids).issubset(set(distance.index)), True, f"{locus} distance rows join complete", checks)
    assert_equal(set(core_ids).issubset(set(distance.columns)), True, f"{locus} distance cols join complete", checks)
    assert_equal(
        reads.loc[core_ids, "alt_support"].eq("ALT").all(),
        True,
        f"{locus} all core reads focal ALT",
        checks,
    )

    for field in ["read_name", "start", "end", "mapq", "strand"]:
        left = core.set_index("read_id").loc[core_ids, field].astype(str).tolist()
        right = reads.loc[core_ids, field].astype(str).tolist()
        assert_equal(left, right, f"{locus} {field} metadata exact", checks)

    raw_labels = core["label"].astype(str)
    counts = raw_labels.value_counts().to_dict()
    ordered_raw = sorted(counts, key=lambda x: (-counts[x], x))
    name_map = {raw: f"Group {chr(65 + i)}" for i, raw in enumerate(ordered_raw)}
    group_for_read = dict(zip(core_ids, raw_labels.map(name_map)))
    group_names = [name_map[x] for x in ordered_raw]
    group_counts = {name_map[raw]: int(counts[raw]) for raw in ordered_raw}

    core_distance = distance.loc[core_ids, core_ids].astype(float)
    finite_core = np.isfinite(core_distance.to_numpy())
    np.fill_diagonal(finite_core, True)
    assert_equal(bool(finite_core.all()), True, f"{locus} core-core distance finite", checks)
    symmetry = np.nanmax(np.abs(core_distance.to_numpy() - core_distance.to_numpy().T))
    assert_equal(bool(symmetry <= 1e-9), True, f"{locus} distance symmetry", checks)
    assert_equal(
        bool(np.nanmax(np.abs(np.diag(core_distance.to_numpy()))) <= 1e-9),
        True,
        f"{locus} distance diagonal zero",
        checks,
    )

    distance_rows: List[Dict[str, Any]] = []
    within_values: List[float] = []
    for group in group_names:
        ids = [rid for rid in core_ids if group_for_read[rid] == group]
        arr = core_distance.loc[ids, ids].to_numpy(dtype=float)
        vals = arr[np.triu_indices(len(ids), 1)]
        within_values.extend(vals.tolist())
        stats = quantiles(vals.tolist())
        distance_rows.append({"locus": locus, "relation": f"within {group}", "group": group, **stats})
    a_ids = [rid for rid in core_ids if group_for_read[rid] == group_names[0]]
    b_ids = [rid for rid in core_ids if group_for_read[rid] == group_names[1]]
    between_values = core_distance.loc[a_ids, b_ids].to_numpy(dtype=float).ravel().tolist()
    between_stats = quantiles(between_values)
    distance_rows.append({"locus": locus, "relation": "between A–B", "group": "Between", **between_stats})
    within_stats = quantiles(within_values)

    core_methyl = methyl.loc[core_ids].astype(float)
    shared = core_methyl.notna().astype(int).dot(core_methyl.notna().astype(int).T)
    group_methyl: Dict[str, Dict[str, Any]] = {}
    per_cpg: List[Dict[str, Any]] = []
    for group in group_names:
        ids = [rid for rid in core_ids if group_for_read[rid] == group]
        values = core_methyl.loc[ids].to_numpy(dtype=float)
        finite = np.isfinite(values)
        per_read_mean = np.nanmean(values, axis=1)
        per_read_called = finite.sum(axis=1)
        group_methyl[group] = {
            "n_reads": len(ids),
            "cell_weighted_mean": finite_or_none(np.nanmean(values)),
            "coverage": finite_or_none(finite.sum() / finite.size),
            "per_read_mean": quantiles(per_read_mean.tolist()),
            "called_cpg": quantiles(per_read_called.tolist()),
        }
    for cpg in core_methyl.columns:
        a = core_methyl.loc[a_ids, cpg].dropna()
        b = core_methyl.loc[b_ids, cpg].dropna()
        a_cov = len(a) / len(a_ids)
        b_cov = len(b) / len(b_ids)
        delta = float(b.mean() - a.mean()) if len(a) and len(b) else math.nan
        per_cpg.append(
            {
                "cpg": str(cpg),
                "a_mean": finite_or_none(a.mean() if len(a) else math.nan),
                "b_mean": finite_or_none(b.mean() if len(b) else math.nan),
                "delta_b_minus_a": finite_or_none(delta),
                "abs_delta": finite_or_none(abs(delta)),
                "a_coverage": finite_or_none(a_cov),
                "b_coverage": finite_or_none(b_cov),
                "min_coverage": finite_or_none(min(a_cov, b_cov)),
            }
        )
    eligible_cpg = [x for x in per_cpg if x["min_coverage"] is not None and x["min_coverage"] >= 0.5]
    eligible_cpg.sort(key=lambda x: (-float(x["abs_delta"]), x["cpg"]))
    selected_cpg = [x["cpg"] for x in eligible_cpg[:DISPLAY_CPGS]]
    if len(selected_cpg) < DISPLAY_CPGS:
        remaining = [x for x in per_cpg if x["cpg"] not in selected_cpg]
        remaining.sort(key=lambda x: (-float(x["min_coverage"] or 0), -float(x["abs_delta"] or 0), x["cpg"]))
        selected_cpg.extend([x["cpg"] for x in remaining[: DISPLAY_CPGS - len(selected_cpg)]])

    representatives = select_representatives(core_distance, group_for_read, group_names)
    display_labels: Dict[str, str] = {}
    for group in group_names:
        ids = [rid for rid in representatives if group_for_read[rid] == group]
        prefix = group.split()[-1]
        for i, rid in enumerate(ids, start=1):
            display_labels[rid] = f"{prefix}{i:02d}"

    distance_heatmap: List[Dict[str, Any]] = []
    shared_heatmap: List[Dict[str, Any]] = []
    methyl_heatmap: List[Dict[str, Any]] = []
    for rid in representatives:
        base = {
            "read": display_labels[rid],
            "group": group_for_read[rid],
            "read_id": rid,
        }
        distance_row = dict(base)
        shared_row = dict(base)
        for other in representatives:
            field = f"d_{display_labels[other]}"
            distance_row[field] = finite_or_none(core_distance.loc[rid, other])
            shared_row[field] = int(shared.loc[rid, other])
        distance_heatmap.append(distance_row)
        shared_heatmap.append(shared_row)
        methyl_row = dict(base)
        for cpg in selected_cpg:
            methyl_row[f"cpg_{cpg}"] = finite_or_none(core_methyl.loc[rid, cpg])
        methyl_heatmap.append(methyl_row)

    core_indexed = core.set_index("read_id")
    geometry: List[Dict[str, Any]] = []
    focal = int(pos)
    for rid in core_ids:
        row = core_indexed.loc[rid]
        geometry.append(
            {
                "locus": locus,
                "read_id": rid,
                "read_label": display_labels.get(rid, rid),
                "group": group_for_read[rid],
                "start_offset_kb": finite_or_none((float(row["start"]) - focal) / 1000),
                "end_offset_kb": finite_or_none((float(row["end"]) - focal) / 1000),
                "length_kb": finite_or_none((float(row["end"]) - float(row["start"])) / 1000),
                "strand": row["strand"],
                "mapq": int(row["mapq"]),
                "latest_hp": str(row["latest_hp"]),
                "latest_ps": str(row["latest_ps"]),
                "cpg_called": int(core_methyl.loc[rid].notna().sum()),
            }
        )

    latest_hp_counts = Counter(core["latest_hp"].astype(str))
    latest_ps_counts = Counter(core["latest_ps"].astype(str))
    raw_hp_counts = Counter(reads.loc[core_ids, "hp"].astype(str))
    max_cpg = max(per_cpg, key=lambda x: float(x["abs_delta"] or -1))
    global_delta = (
        float(group_methyl[group_names[1]]["cell_weighted_mean"])
        - float(group_methyl[group_names[0]]["cell_weighted_mean"])
    )
    description = (
        "global high/low methylation partition"
        if abs(global_delta) >= 0.15
        else "CpG-pattern partition with similar global methylation means"
    )
    return {
        "key": key,
        "locus": locus,
        "site": dict(site),
        "audit": dict(audit_row),
        "gene": dict(gene),
        "group_names": group_names,
        "group_counts": group_counts,
        "raw_label_map": name_map,
        "core_read_n": len(core_ids),
        "representative_read_n": len(representatives),
        "display_labels": display_labels,
        "distance_rows": distance_rows,
        "within_stats": within_stats,
        "between_stats": between_stats,
        "between_within_delta": finite_or_none(between_stats["median"] - within_stats["median"]),
        "between_within_ratio": finite_or_none(between_stats["median"] / within_stats["median"]),
        "group_methyl": group_methyl,
        "global_methyl_delta_b_minus_a": finite_or_none(global_delta),
        "max_cpg_delta": max_cpg,
        "pattern_description": description,
        "per_cpg": per_cpg,
        "selected_cpg": selected_cpg,
        "distance_heatmap": distance_heatmap,
        "shared_heatmap": shared_heatmap,
        "methyl_heatmap": methyl_heatmap,
        "geometry": geometry,
        "distance_fields": [f"d_{display_labels[rid]}" for rid in representatives],
        "distance_field_labels": {
            f"d_{display_labels[rid]}": display_labels[rid] for rid in representatives
        },
        "cpg_fields": [f"cpg_{cpg}" for cpg in selected_cpg],
        "cpg_field_labels": {f"cpg_{cpg}": cpg for cpg in selected_cpg},
        "latest_hp_counts": dict(latest_hp_counts),
        "latest_ps_counts": dict(latest_ps_counts),
        "raw_hp_counts": dict(raw_hp_counts),
        "full_matrix_n": int(len(distance)),
        "full_distance_finite_fraction": finite_or_none(np.isfinite(distance.to_numpy(dtype=float)).mean()),
        "core_distance_finite_fraction": finite_or_none(np.isfinite(core_distance.to_numpy(dtype=float)).mean()),
        "artifact_hashes": {
            "reads": artifacts["reads"]["sha256"],
            "methyl": artifacts["methylation_matrix"]["sha256"],
            "distance": artifacts["distance_matrix"]["sha256"],
        },
    }


def table_spec(
    table_id: str,
    title: str,
    subtitle: str,
    dataset: str,
    columns: List[Dict[str, Any]],
    sort_field: str,
    direction: str = "desc",
    source_id: str = "derivative_results",
) -> Dict[str, Any]:
    return {
        "id": table_id,
        "title": title,
        "subtitle": subtitle,
        "dataset": dataset,
        "sourceId": source_id,
        "density": "dense",
        "defaultSort": {"field": sort_field, "direction": direction},
        "columns": columns,
    }


def interpolate_hex(start: str, end: str, fraction: float) -> str:
    """Interpolate two six-digit hex colors."""
    fraction = min(1.0, max(0.0, float(fraction)))
    left = tuple(int(start[i : i + 2], 16) for i in (1, 3, 5))
    right = tuple(int(end[i : i + 2], 16) for i in (1, 3, 5))
    rgb = tuple(round(a + (b - a) * fraction) for a, b in zip(left, right))
    return "#" + "".join(f"{value:02x}" for value in rgb)


def readable_text_color(background: str) -> str:
    """Return dark or light text for a hex background."""
    red, green, blue = (int(background[i : i + 2], 16) for i in (1, 3, 5))
    luminance = (0.2126 * red + 0.7152 * green + 0.0722 * blue) / 255
    return "#111827" if luminance >= 0.56 else "#ffffff"


def heatmap_html(
    *,
    block_id: str,
    title: str,
    subtitle: str,
    rows: Sequence[Mapping[str, Any]],
    value_fields: Sequence[str],
    field_labels: Mapping[str, str],
    value_min: float,
    value_max: float,
    start_color: str,
    end_color: str,
    value_format: str,
) -> str:
    """Render an accessible, script-free colored matrix for an HTML artifact block."""
    if not rows or not value_fields:
        raise ValueError(f"{block_id}: heatmap rows and fields must be non-empty")
    span = value_max - value_min

    def format_value(value: Any) -> str:
        if value is None or not math.isfinite(float(value)):
            return "NA"
        numeric = float(value)
        if value_format == "integer":
            return f"{numeric:.0f}"
        if value_format == "percent":
            return f"{numeric * 100:.0f}%"
        return f"{numeric:.2f}"

    def color_for(value: Any) -> Tuple[str, str]:
        if value is None or not math.isfinite(float(value)):
            return "#e5e7eb", "#6b7280"
        numeric = float(value)
        fraction = 0.5 if span <= 0 else (numeric - value_min) / span
        background = interpolate_hex(start_color, end_color, fraction)
        return background, readable_text_color(background)

    group_for_column: Dict[str, str] = {}
    for field in value_fields:
        label = str(field_labels.get(field, field))
        group_for_column[field] = "Group A" if label.startswith("A") else "Group B" if label.startswith("B") else ""

    styles = """
<style>
  .hm-card{max-width:760px;margin:0 auto;padding:18px;border:1px solid #e5e7eb;border-radius:14px;background:#fff;color:#111827;font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
  .hm-card h2{margin:0;font-size:18px;line-height:1.3}.hm-card p{margin:7px 0 12px;color:#4b5563;font-size:13px;line-height:1.45}
  .hm-scroll{max-width:100%;overflow-x:auto;padding-bottom:4px}.hm-grid{margin:0 auto;border-collapse:separate;border-spacing:2px;font-variant-numeric:tabular-nums}
  .hm-grid th{height:32px;padding:3px 6px;background:#f9fafb;color:#374151;font-size:11px;font-weight:650;text-align:center;white-space:nowrap}
  .hm-grid th.hm-row-label{text-align:right}.hm-grid td{min-width:48px;height:40px;padding:2px 4px;border-radius:4px;font-size:10px;font-weight:650;text-align:center}
  .hm-group-a{box-shadow:inset 0 3px 0 #2563eb}.hm-group-b{box-shadow:inset 0 3px 0 #ea580c}
  .hm-row-a{border-right:4px solid #2563eb}.hm-row-b{border-right:4px solid #ea580c}
  .hm-break{border-left:3px solid #111827!important}.hm-legend{display:flex;align-items:center;gap:8px;margin:12px auto 0;max-width:460px;color:#4b5563;font-size:11px}
  .hm-ramp{height:12px;flex:1;border:1px solid #d1d5db;border-radius:999px}
  .hm-groups{display:flex;justify-content:center;gap:16px;margin-top:8px;color:#4b5563;font-size:11px}
  .hm-dot{display:inline-block;width:9px;height:9px;margin-right:5px;border-radius:999px}.hm-a{background:#2563eb}.hm-b{background:#ea580c}
  @media(max-width:520px){.hm-card{padding:12px}.hm-grid td{min-width:42px;height:36px}.hm-card h2{font-size:16px}}
</style>
"""
    header_cells = []
    previous_group = ""
    for field in value_fields:
        label = str(field_labels.get(field, field))
        group = group_for_column[field]
        classes = ["hm-group-a" if group == "Group A" else "hm-group-b" if group == "Group B" else ""]
        if previous_group and group and group != previous_group:
            classes.append("hm-break")
        previous_group = group or previous_group
        header_cells.append(
            f'<th scope="col" class="{html.escape(" ".join(x for x in classes if x))}">'
            f"{html.escape(label)}</th>"
        )

    body_rows = []
    for row_index, row in enumerate(rows):
        read_label = str(row.get("read", f"R{row_index + 1}"))
        group = str(row.get("group", ""))
        row_class = "hm-row-a" if group == "Group A" else "hm-row-b" if group == "Group B" else ""
        cells = []
        previous_group = ""
        for field in value_fields:
            value = row.get(field)
            display = format_value(value)
            background, foreground = color_for(value)
            column_group = group_for_column[field]
            break_class = " hm-break" if previous_group and column_group and column_group != previous_group else ""
            previous_group = column_group or previous_group
            value_attr = "" if value is None else str(value)
            cells.append(
                f'<td class="hm-cell{break_class}" data-heatmap-cell="true" '
                f'data-row="{html.escape(read_label)}" '
                f'data-column="{html.escape(str(field_labels.get(field, field)))}" '
                f'data-value="{html.escape(value_attr)}" '
                f'style="background:{background};color:{foreground}" '
                f'title="{html.escape(read_label)} × {html.escape(str(field_labels.get(field, field)))}: '
                f'{html.escape(display)}">{html.escape(display)}</td>'
            )
        body_rows.append(
            f'<tr><th scope="row" class="hm-row-label {row_class}">{html.escape(read_label)}</th>'
            + "".join(cells)
            + "</tr>"
        )

    gradient = f"linear-gradient(90deg,{start_color},{end_color})"
    return (
        styles
        + f'<figure class="hm-card" data-true-heatmap="{html.escape(block_id)}">'
        + f'<figcaption><h2>{html.escape(title)}</h2><p>{html.escape(subtitle)}</p></figcaption>'
        + '<div class="hm-scroll"><table class="hm-grid"><thead><tr>'
        + '<th scope="col">ALT read</th>'
        + "".join(header_cells)
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table></div>"
        + '<div class="hm-legend" aria-label="Color scale">'
        + f"<span>{html.escape(format_value(value_min))}</span>"
        + f'<span class="hm-ramp" style="background:{gradient}"></span>'
        + f"<span>{html.escape(format_value(value_max))}</span></div>"
        + '<div class="hm-groups"><span><i class="hm-dot hm-a"></i>Group A</span>'
        + '<span><i class="hm-dot hm-b"></i>Group B</span></div>'
        + "</figure>"
    )


def build_artifact(
    generated_at: str,
    datasets: Dict[str, List[Dict[str, Any]]],
    loci: List[Dict[str, Any]],
    headline: Dict[str, Any],
) -> Dict[str, Any]:
    title = "HCC1395 positional-singleton ALT 內甲基子結構完整驗證"
    sources = [
        {
            "id": "singleton_audit",
            "label": "Source-attested positional-singleton methyl audit",
            "path": (
                "output/synthesis/observation_workspaces/"
                "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
                "positional_singleton_methyl_multigroup_audit_v1_source_attested/"
                "positional_singleton_site_audit.tsv.gz"
            ),
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": (
                    "SELECT * FROM read_csv_auto("
                    "'output/synthesis/observation_workspaces/"
                    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
                    "positional_singleton_methyl_multigroup_audit_v1_source_attested/"
                    "positional_singleton_site_audit.tsv.gz', "
                    "delim='\\t', header=true);"
                ),
                "description": "Complete 50,432-site positional-singleton audit across seven dataset rows.",
                "executed_at": generated_at,
                "tables_used": [
                    "positional_singleton_site_audit.tsv.gz",
                    "positional_singleton_audit_summary.json",
                    "_SUCCESS.json",
                ],
                "filters": ["chr1-22", "LongPhase-S recalibrated PASS", "biallelic sSNV"],
                "metric_definitions": [
                    "singleton = same-dataset/chrom 50 kb transitive component size 1",
                    "M1 = coarse_ng>=2, stable, modal assignment ARI min>=0.8",
                    "M2 = all eight measured axes determinate and no aligned confound",
                ],
            },
        },
        {
            "id": "screen_assignments",
            "label": "All-sSNV focal-ALT screen and stable read assignments",
            "path": (
                "output/synthesis/observation_workspaces/"
                "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
                "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/"
                "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
            ),
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": (
                    "SELECT * FROM read_json_auto("
                    "'output/synthesis/observation_workspaces/"
                    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
                    "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/"
                    "all_ssnv_stable_multigroup_read_assignments.jsonl.gz', "
                    "format='newline_delimited') "
                    "WHERE dataset='HCC1395' AND "
                    "((chrom='chr14' AND pos=86272476) OR (chrom='chr22' AND pos=47466517));"
                ),
                "description": "Exact focal-ALT after-peel core read identities and M1 stable group assignments.",
                "executed_at": generated_at,
                "tables_used": [
                    "all_ssnv_site_results.tsv.gz",
                    "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
                ],
                "filters": ["dataset=HCC1395", "M2 PASS target loci"],
                "metric_definitions": [
                    "Core reads are the exact focal-ALT after-peel read set",
                    "Internal labels are relabeled Group A/B and must not be interpreted as HP tags",
                ],
            },
        },
        {
            "id": "region_matrices",
            "label": "Primary read, methylation, and Bernoulli distance matrices",
            "path": (
                "output/synthesis/observation_workspaces/"
                "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
                "intersubmod_all_ssnv_v2_verification_fix/HCC1395/"
            ),
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": (
                    "SELECT * FROM read_csv_auto("
                    "'research/20260718_singleton_alt_methyl_substructure_validation/results/"
                    "locus_*_display_*.tsv', delim='\\t', header=true, "
                    "union_by_name=true, filename=true);"
                ),
                "description": "SHA-verified exact joins of reads.tsv, methylation.csv, and BERNOULLI matrix.csv.",
                "executed_at": generated_at,
                "tables_used": ["reads/reads.tsv", "methylation/methylation.csv", "distance/BERNOULLI/matrix.csv"],
                "filters": ["core focal-ALT read IDs only"],
                "metric_definitions": [
                    "Bernoulli distance is computed over available read-pair CpG evidence",
                    "Shared-CpG heatmaps expose pairwise support for each displayed distance cell",
                ],
            },
        },
        {
            "id": "gene_context",
            "label": "GENCODE v46, COSMIC CGC v104, and DGIdb context",
            "path": "research/20260718_singleton_alt_methyl_substructure_validation/results/gene_context.tsv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": (
                    "SELECT * FROM read_csv_auto("
                    "'research/20260718_singleton_alt_methyl_substructure_validation/results/"
                    "gene_context.tsv', delim='\\t', header=true);"
                ),
                "description": "Gene-body overlap or nearest-gene annotation; context only, not topology or clone truth.",
                "executed_at": generated_at,
                "tables_used": ["GENCODE v46 basic GTF", "COSMIC Cancer Gene Census v104", "DGIdb interactions"],
                "metric_definitions": [
                    "CGC and DGIdb annotations are descriptive context, not locus-specific validation"
                ],
            },
        },
        {
            "id": "derivative_results",
            "label": "Deterministic report build outputs",
            "path": "research/20260718_singleton_alt_methyl_substructure_validation/results/validation_receipt.json",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": (
                    "SELECT * FROM read_json_auto("
                    "'research/20260718_singleton_alt_methyl_substructure_validation/results/"
                    "validation_receipt.json');"
                ),
                "description": "Independent recount, exact joins, summary statistics, and bounded report datasets.",
                "executed_at": generated_at,
                "tables_used": [
                    "hcc1395_singleton_site_audit.tsv.gz",
                    "hcc1395_m2_pass_locus_summary.tsv",
                    "validation_receipt.json",
                ],
                "filters": ["full HCC1395 singleton universe; no site subsampling"],
                "metric_definitions": [
                    "Heatmap display uses up to 12 medoid-nearest reads per group; full group metrics use all core reads",
                    "Operational yield is not biological prevalence",
                ],
            },
        },
    ]

    cards = [
        {
            "id": "singleton_card",
            "description": "完整 HCC1395 positional-singleton dataset-sites；每個 component 恰有一個 sSNV。",
            "dataset": "headline",
            "sourceId": "singleton_audit",
            "metrics": [
                {"label": "Singleton sites / regions", "field": "singleton_sites", "format": "number"},
                {"label": "占 HCC1395 PASS sSNV", "field": "singleton_share_of_hcc_pass", "format": "percent"},
            ],
        },
        {
            "id": "m1_card",
            "description": "M1 是高敏感度、穩定甲基多群 screen；不是 clone prevalence。",
            "dataset": "headline",
            "sourceId": "singleton_audit",
            "metrics": [
                {"label": "M1 stable multigroup", "field": "m1_flagged", "format": "number"},
                {"label": "母體比例", "field": "m1_rate", "format": "percent"},
            ],
        },
        {
            "id": "m2_card",
            "description": "八個已測 HP/strand/read-geometry/CpG-count 軸均可判定且未與群組對齊。",
            "dataset": "headline",
            "sourceId": "singleton_audit",
            "metrics": [
                {"label": "M2 clear residual", "field": "m2_pass", "format": "number"},
                {"label": "每 100k singleton", "field": "m2_per_100k", "format": "number"},
            ],
        },
        {
            "id": "clone_card",
            "description": "G1/G2、formal R1、matched-normal 與 CN/CCF 未完成；0 表示尚未建立，不是真陰性。",
            "dataset": "headline",
            "sourceId": "derivative_results",
            "metrics": [
                {"label": "Established cellular clone/subclone", "field": "established_clone", "format": "number"},
                {"label": "Lineage/order established", "field": "established_lineage", "format": "number"},
            ],
        },
    ]

    charts: List[Dict[str, Any]] = [
        {
            "id": "hcc_screen_chart",
            "title": "HCC1395 singleton 甲基 screen 各階段比例",
            "subtitle": "相對 8,279 singleton sites；這是甲基 screen，不含 G1/G2 遺傳佐證。",
            "type": "bar",
            "dataset": "hcc_screen",
            "sourceId": "singleton_audit",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "stage", "type": "nominal", "label": "Stage"},
                "y": {"field": "share_of_singleton", "type": "quantitative", "label": "Share of singleton sites"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "hcc_status_chart",
            "title": "HCC1395 singleton 的 M2 狀態",
            "subtitle": "NOT_EVALUABLE 是資訊不足；NOT_RUN 是未通過 M1，兩者都不是陰性。",
            "type": "bar",
            "dataset": "hcc_status",
            "sourceId": "singleton_audit",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "status", "type": "nominal", "label": "M2 status"},
                "y": {"field": "count", "type": "quantitative", "label": "Sites"},
            },
            "valueFormat": "number",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "dataset_m1_chart",
            "title": "七個 dataset rows 的 M1 stable multigroup 比例",
            "subtitle": "同一 singleton contract；HCC1395 與 DORADO 是同一生物樣本，不是獨立 replication。",
            "type": "bar",
            "dataset": "dataset_summary_chart",
            "sourceId": "singleton_audit",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "dataset", "type": "nominal", "label": "Dataset"},
                "y": {"field": "m1_rate", "type": "quantitative", "label": "M1 rate"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "dataset_m2_chart",
            "title": "七個 dataset rows 的 M2 PASS operational yield",
            "subtitle": "每 100,000 singleton dataset-sites；不能解讀為 subclone prevalence。",
            "type": "bar",
            "dataset": "dataset_summary_chart",
            "sourceId": "singleton_audit",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "dataset", "type": "nominal", "label": "Dataset"},
                "y": {"field": "m2_per_100k", "type": "quantitative", "label": "M2 PASS per 100k"},
            },
            "valueFormat": "number",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "hcc_group_count_chart",
            "title": "HCC1395 M1 flags 的甲基群數分布",
            "subtitle": "群數來自 M1 operational screen；不等於 cellular clone 數。",
            "type": "bar",
            "dataset": "hcc_group_count",
            "sourceId": "singleton_audit",
            "intent": "distribution",
            "encodings": {
                "x": {"field": "methyl_groups", "type": "ordinal", "label": "M1 methyl groups"},
                "y": {"field": "sites", "type": "quantitative", "label": "Sites"},
            },
            "valueFormat": "number",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "distance_summary_chart",
            "title": "兩個 M2 PASS 位點的群內與群間 Bernoulli 距離",
            "subtitle": "中位數由全部 108/109 core ALT reads 計算；pairwise values 非獨立，不宣告 p-value。",
            "type": "bar",
            "dataset": "distance_summary",
            "sourceId": "region_matrices",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "locus_short", "type": "nominal", "label": "Locus"},
                "y": {"field": "median", "type": "quantitative", "label": "Median distance"},
                "color": {"field": "relation", "type": "nominal", "label": "Pair type"},
            },
            "valueFormat": "number",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "Read-pair relation"},
        },
    ]

    blocks: List[Dict[str, Any]] = [
        {"id": "title", "type": "markdown", "body": f"# {title}"},
        {
            "id": "scope",
            "type": "markdown",
            "body": (
                "**COMPLETE SINGLETON AUDIT / PARTIAL BIOLOGICAL VALIDATION — "
                "完整分析 HCC1395 8,279 個 positional-singleton dataset-sites；"
                "正式 G1/G2 genetic co-segregation、matched-normal、CN/purity/CCF 與 cellular lineage 尚未完成。**"
            ),
            "sourceId": "singleton_audit",
        },
        {
            "id": "executive_summary",
            "type": "markdown",
            "body": (
                "## 重點結論：可找到兩個清楚 ALT-read 甲基子結構，但不能標成已證實 clone\n\n"
                f"- **8,279** 是 HCC1395 的單點區域數，也是 singleton site 數；不是 8,279 個 sSNV 擠在一區。\n"
                f"- M1 stable multigroup 為 **734/8,279 = {headline['m1_rate']*100:.2f}%**，適合當候選 screen。\n"
                f"- 只有 **2/8,279 = {headline['m2_rate']*100:.4f}%** 通過八個 measured-axis guardrails；"
                "這兩點可標為 **M2 residual epigenetic partition candidate**。\n"
                "- **不能標成特殊 cellular clone/subclone**：兩點都是 positional singleton，沒有 local second sSNV；"
                "G1/G2/formal R1 皆 NOT_RUN。VAF 只提供 focal allele burden context，未用來發現甲基群，也不能決定父子關係。"
            ),
            "sourceId": "derivative_results",
        },
        {"id": "headline_metrics", "type": "metric-strip", "cardIds": [x["id"] for x in cards]},
        {
            "id": "denominator_heading",
            "type": "markdown",
            "body": (
                "## 分母與 screen：先分清 79,687、8,279、734 與 2\n\n"
                "HCC1395 有 79,687 個 chr1-22 autosomal biallelic LongPhase-S PASS sSNV。"
                "以同染色體相鄰 gap≤50 kb 形成傳遞 component 後，8,279 個 component 各只有一個 sSNV；"
                "它們才是本報告的完整母體。M1/M2 是甲基軸的兩階段 screen，不是通往 clone 的單一路徑。"
            ),
            "sourceId": "singleton_audit",
        },
        {"id": "hcc_screen", "type": "chart", "chartId": "hcc_screen_chart"},
        {"id": "hcc_status", "type": "chart", "chartId": "hcc_status_chart"},
        {"id": "hcc_status_table", "type": "table", "tableId": "hcc_status_table_spec"},
        {
            "id": "datasets_heading",
            "type": "markdown",
            "body": (
                "## 七個 dataset rows：M1 常見，M2 clear residual 很少\n\n"
                "全資料共有 50,432 個 singleton dataset-sites、5,961 個 M1 flags；"
                "M2 PASS 30、FAIL 18、NOT_EVALUABLE 5,913。30/50,432 是 operational yield，"
                "30/48=62.5% 只適用於高度選擇的 M2-determinate subset，不能外推盛行率。"
            ),
            "sourceId": "singleton_audit",
        },
        {"id": "dataset_m1", "type": "chart", "chartId": "dataset_m1_chart"},
        {"id": "dataset_m2", "type": "chart", "chartId": "dataset_m2_chart"},
        {"id": "dataset_table", "type": "table", "tableId": "dataset_summary_table_spec"},
        {
            "id": "hcc_pair_heading",
            "type": "markdown",
            "body": (
                "## HCC1395 與 DORADO：singleton 層只能做技術一致性描述\n\n"
                "兩者共享同一生物 cell line；exact-site 的 M1 狀態可比，但不是獨立生物重現。"
                "兩側各有 2 個 M2 PASS，卻沒有 exact-locus M2 PASS 重疊；這不否定整體甲基結構，"
                "也不能證明固定 clone，較合理解釋是 singleton eligibility、read sampling 與可評估性差異。"
            ),
            "sourceId": "derivative_results",
        },
        {"id": "hcc_pair_table", "type": "table", "tableId": "hcc_pair_table_spec"},
        {"id": "group_count", "type": "chart", "chartId": "hcc_group_count_chart"},
        {
            "id": "claim_heading",
            "type": "markdown",
            "body": (
                "## Claim ladder 是兩條軸，不是把甲基群直接升級成 clone\n\n"
                "**甲基軸：M1 stable → M2 measured-axis residual。** "
                "**遺傳軸：G1 partner co-segregation → G2 multi-marker molecular haplotype；"
                "formal R1 必須在 G2 之後。** 本報告的兩點只到 M2；"
                "G1/G2/R1 的 NOT_RUN 不是 0 個陰性結果。"
            ),
            "sourceId": "derivative_results",
        },
        {"id": "claim_table", "type": "table", "tableId": "claim_table_spec"},
        {
            "id": "examples_heading",
            "type": "markdown",
            "body": (
                "## 兩個清楚例子：chr14 是整體高／低甲基，chr22 是 CpG pattern 分群\n\n"
                "兩點都用 exact focal-ALT after-peel core reads；internal `1-1/1-2` 只是 cluster label，"
                "報告統一改稱 Group A/B。兩點的 latest HP/PS 都是常數，因此 HP 沒有把兩群切開，"
                "也不能當獨立佐證。"
            ),
            "sourceId": "region_matrices",
        },
        {"id": "locus_table", "type": "table", "tableId": "locus_table_spec"},
        {"id": "distance_summary", "type": "chart", "chartId": "distance_summary_chart"},
    ]

    for i, locus in enumerate(loci, start=1):
        prefix = f"locus_{i}"
        locus_short = locus["locus"]
        core_n = locus["core_read_n"]
        a_n = locus["group_counts"]["Group A"]
        b_n = locus["group_counts"]["Group B"]
        between = locus["between_stats"]["median"]
        within = locus["within_stats"]["median"]
        ratio = locus["between_within_ratio"]
        delta = locus["global_methyl_delta_b_minus_a"]
        distance_values = [
            float(row[field])
            for row in locus["distance_heatmap"]
            for field in locus["distance_fields"]
            if row.get(field) is not None
        ]
        shared_values = [
            float(row[field])
            for row in locus["shared_heatmap"]
            for field in locus["distance_fields"]
            if row.get(field) is not None
        ]
        distance_body = heatmap_html(
            block_id=f"{prefix}_distance_heatmap",
            title=f"{locus_short} representative core-read Bernoulli distance",
            subtitle=(
                "A/B labels mark M1 methyl groups; full core-core distance pairs are 100% finite. "
                "Cells show distance; darker means more dissimilar."
            ),
            rows=locus["distance_heatmap"],
            value_fields=locus["distance_fields"],
            field_labels=locus["distance_field_labels"],
            value_min=min(distance_values),
            value_max=max(distance_values),
            start_color="#f7fbff",
            end_color="#084594",
            value_format="number",
        )
        shared_body = heatmap_html(
            block_id=f"{prefix}_shared_cpg_heatmap",
            title=f"{locus_short} displayed read-pair shared CpG counts",
            subtitle=(
                "Same read order as the distance heatmap. Cells show shared called CpGs; "
                "paler cells warn that the paired distance has less support."
            ),
            rows=locus["shared_heatmap"],
            value_fields=locus["distance_fields"],
            field_labels=locus["distance_field_labels"],
            value_min=min(shared_values),
            value_max=max(shared_values),
            start_color="#fff7ec",
            end_color="#7f0000",
            value_format="integer",
        )
        methyl_body = heatmap_html(
            block_id=f"{prefix}_methyl_probability_heatmap",
            title=f"{locus_short} representative read × CpG methylation probability",
            subtitle=(
                "CpGs are selected by coverage and |Group B−A|. Cells show methylation probability; "
                "gray NA cells had no call."
            ),
            rows=locus["methyl_heatmap"],
            value_fields=locus["cpg_fields"],
            field_labels=locus["cpg_field_labels"],
            value_min=0.0,
            value_max=1.0,
            start_color="#440154",
            end_color="#fde725",
            value_format="percent",
        )
        blocks.extend(
            [
                {
                    "id": f"{prefix}_heading",
                    "type": "markdown",
                    "body": (
                        f"## {locus_short}：{locus['pattern_description']}\n\n"
                        f"全部 core ALT reads n={core_n}；Group A/B={a_n}/{b_n}。"
                        f"群間距離中位數 {between:.3f}，pooled 群內 {within:.3f}，比值 {ratio:.2f}×；"
                        f"Group B−A 全域甲基均值差 {delta:+.3f}。"
                        f"熱圖只顯示每群最多 {DISPLAY_PER_GROUP} 個 deterministic medoid-nearest reads，"
                        "但上述統計均使用全部 core reads。"
                    ),
                    "sourceId": "region_matrices",
                },
                {"id": f"{prefix}_distance", "type": "html", "body": distance_body},
                {"id": f"{prefix}_shared", "type": "html", "body": shared_body},
                {"id": f"{prefix}_methyl", "type": "html", "body": methyl_body},
                {"id": f"{prefix}_geometry", "type": "chart", "chartId": f"{prefix}_geometry_chart"},
            ]
        )
        charts.extend(
            [
                {
                    "id": f"{prefix}_geometry_chart",
                    "title": f"{locus_short} core-read start and end offsets",
                    "subtitle": "Group A/B are not aligned to start/end/length under the formal M2 measured-axis gate.",
                    "type": "scatter",
                    "dataset": f"{prefix}_geometry",
                    "sourceId": "region_matrices",
                    "intent": "relationship",
                    "encodings": {
                        "x": {"field": "start_offset_kb", "type": "quantitative", "label": "Start offset (kb)"},
                        "y": {"field": "end_offset_kb", "type": "quantitative", "label": "End offset (kb)"},
                        "color": {"field": "group", "type": "nominal", "label": "Methyl group"},
                    },
                    "valueFormat": "number",
                    "palette": {"kind": "categorical"},
                    "legend": {"position": "bottom", "interactive": True, "title": "Methyl group"},
                },
            ]
        )

    blocks.extend(
        [
            {
                "id": "vaf_heading",
                "type": "markdown",
                "body": (
                    "## VAF 與基因註解：只能提供 context，不能把群命名為 clone\n\n"
                    "chr14 位點 caller AF=0.827、chr22 AF=1.000；兩者 normal AF 都是 0。"
                    "高 VAF 可由 clonal burden、LOH、copy number、purity 或 sampling 造成。"
                    "甲基群是在 read×CpG 距離中發現，**沒有用 VAF 來切群**；單一 VAF 也不能決定 parent/child。"
                    "GENCODE/CGC/DGIdb 註解只作區域背景，不是甲基子結構或用藥 truth set。"
                ),
                "sourceId": "gene_context",
            },
            {"id": "gene_table", "type": "table", "tableId": "gene_table_spec"},
            {
                "id": "limitations",
                "type": "markdown",
                "body": (
                    "## 限制與可用標記\n\n"
                    "- 可標：`M1 stable ALT-read methyl multigroup`；更嚴格的兩點可標："
                    "`M2 read-level residual epigenetic partition candidate`。\n"
                    "- 不可標：已證實 clone/subclone、clone 數、linear/branching lineage、唯一真實演化樹。\n"
                    "- M2 只檢查八個 measured axes，沒有排除未量測 confound、CN/LOH/purity 或 matched-normal background。\n"
                    "- 兩點都是 positional singleton，無 local partner sSNV；G1/G2/formal R1 因此 NOT_RUN。\n"
                    "- HCC1395 與 HCC1395_DORADO 是同一生物樣本；technical similarity 不能替代獨立生物驗證。\n"
                    "- 代表 read heatmap 是顯示用；全部 effect summaries 使用全 core reads，pairwise distances 不視為獨立樣本。"
                ),
                "sourceId": "derivative_results",
            },
            {"id": "validation_table", "type": "table", "tableId": "validation_table_spec"},
            {"id": "provenance_table", "type": "table", "tableId": "provenance_table_spec"},
            {
                "id": "next_steps",
                "type": "markdown",
                "body": (
                    "## 升級成 clone/subclone claim 的必要下一步\n\n"
                    "1. 完成正式 G1/G2：至少第二個獨立 genetic marker 在同一 molecules 中呈 group-specific co-segregation。\n"
                    "2. 補 tumor-REF 與 matched-normal 同尺度對照，確認不是背景 epiallele 或 allele-independent methyl states。\n"
                    "3. 納入 CN/LOH、purity、multiplicity 與 CCF；VAF 只在此時才可協助 cellular fraction 解釋。\n"
                    "4. 以 single-cell DNA、colony、spatial 或 multi-region truth 建立 cellular identity；之後才能談 clone/subclone 與 lineage。"
                ),
            },
        ]
    )

    tables = [
        table_spec(
            "hcc_status_table_spec",
            "HCC1395 singleton 甲基 screen 完整計數",
            "所有狀態互斥；M2 NOT_EVALUABLE/NOT_RUN 不是 FAIL。",
            "hcc_status",
            [
                {"field": "status", "label": "Status", "type": "text"},
                {"field": "count", "label": "Sites", "format": "number"},
                {"field": "share", "label": "Share of 8,279", "format": "percent"},
                {"field": "interpretation", "label": "Meaning", "type": "text"},
            ],
            "count",
        ),
        table_spec(
            "dataset_summary_table_spec",
            "七個 dataset rows 的 singleton 完整統計",
            "ALL row 是 dataset-site 合計；HCC1395 pair 只算一個 biological sample。",
            "dataset_summary",
            [
                {"field": "dataset", "label": "Dataset", "type": "text"},
                {"field": "sites", "label": "Singleton sites", "format": "number"},
                {"field": "m1_count_rate", "label": "M1 flags (rate)", "type": "text"},
                {"field": "m2_pass", "label": "M2 PASS", "format": "number"},
                {"field": "m2_not_evaluable", "label": "M2 not evaluable", "format": "number"},
                {"field": "m2_per_100k", "label": "M2 / 100k", "format": "number"},
            ],
            "sites",
        ),
        table_spec(
            "hcc_pair_table_spec",
            "HCC1395 technical pair 的 exact singleton-site comparison",
            "只比較兩側都屬 singleton 的 exact variant keys；不是獨立生物 replication。",
            "hcc_pair",
            [
                {"field": "common_exact_singletons", "label": "Common exact singleton loci", "format": "number"},
                {"field": "m1_status_agreement", "label": "M1 state agreement", "format": "percent"},
                {"field": "both_m1_flagged", "label": "Both M1 flagged", "format": "number"},
                {"field": "both_m2_pass", "label": "Exact M2 PASS overlap", "format": "number"},
            ],
            "common_exact_singletons",
        ),
        table_spec(
            "claim_table_spec",
            "本報告可達的證據層級",
            "Status=NOT_RUN 代表未測，不能寫成陰性。",
            "claim_ladder",
            [
                {"field": "evidence_level", "label": "Axis / level", "type": "text"},
                {"field": "status", "label": "Status", "type": "text"},
                {"field": "count", "label": "HCC1395 sites", "format": "number"},
                {"field": "allowed_claim", "label": "Allowed claim", "type": "text"},
            ],
            "count",
        ),
        table_spec(
            "locus_table_spec",
            "HCC1395 兩個 M2 PASS 位點",
            "兩點均 TP、normal AF=0、M2 PASS；G1/G2/R1 NOT_RUN。",
            "locus_summary",
            [
                {"field": "locus", "label": "Locus", "type": "text"},
                {"field": "caller_context", "label": "Caller context", "type": "text"},
                {"field": "group_counts", "label": "Group A/B", "type": "text"},
                {"field": "distance_ratio", "label": "Between/within", "format": "number"},
                {"field": "methyl_delta_b_minus_a", "label": "Methyl Δ B−A", "format": "number"},
                {"field": "interpretation", "label": "Pattern", "type": "text"},
            ],
            "distance_ratio",
        ),
        table_spec(
            "gene_table_spec",
            "兩個 M2 PASS 位點的 gene/cancer/drug context",
            "Annotation is descriptive context only and cannot validate methyl groups or clones.",
            "gene_context",
            [
                {"field": "locus", "label": "Locus", "type": "text"},
                {"field": "gene", "label": "Gene", "type": "text"},
                {"field": "relationship", "label": "Relation", "type": "text"},
                {"field": "distance_bp", "label": "Distance (bp)", "format": "number"},
                {"field": "cancer_drug_context", "label": "CGC / DGIdb", "type": "text"},
            ],
            "distance_bp",
            "asc",
            "gene_context",
        ),
        table_spec(
            "validation_table_spec",
            "獨立重算與 exact-join 驗證",
            "任一核心 fail 都會中止 artifact 生成。",
            "validation",
            [
                {"field": "check", "label": "Check", "type": "text"},
                {"field": "actual", "label": "Actual", "type": "text"},
                {"field": "expected", "label": "Expected", "type": "text"},
                {"field": "status", "label": "Status", "type": "text"},
            ],
            "status",
            "asc",
        ),
        table_spec(
            "provenance_table_spec",
            "來源與 SHA-256",
            "Primary matrices and source-attested audit are immutable inputs; report outputs are derivative.",
            "provenance",
            [
                {"field": "artifact", "label": "Artifact", "type": "text"},
                {"field": "sha256_short", "label": "SHA-256 (short)", "type": "text"},
                {"field": "role", "label": "Role", "type": "text"},
            ],
            "artifact",
            "asc",
        ),
    ]

    return {
        "surface": "report",
        "manifest": {
            "title": title,
            "version": 1,
            "surface": "report",
            "description": (
                "Complete HCC1395 positional-singleton census with two exact focal-ALT read-level "
                "M2-PASS methylation examples, distance/shared-CpG heatmaps, measured-axis guardrails, "
                "and an explicit clone/subclone claim ceiling."
            ),
            "generatedAt": generated_at,
            "blocks": blocks,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "status": "ready",
            "generatedAt": generated_at,
            "datasets": datasets,
        },
        "sources": sources,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=TOPIC / "results")
    parser.add_argument("--artifact", type=Path, default=TOPIC / "artifact.json")
    parser.add_argument(
        "--markdown",
        type=Path,
        default=TOPIC / "20260718_HCC1395_singleton_ALT甲基子結構驗證報告_01.md",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    checks: List[Dict[str, Any]] = []
    generated_at = now_iso()

    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    success = json.loads(SUCCESS_PATH.read_text(encoding="utf-8"))
    assert_equal(summary["pass"], True, "source audit summary pass", checks)
    assert_equal(success["pass"], True, "source atomic marker pass", checks)
    assert_equal(sha256_file(SUMMARY_PATH), success["summary"]["sha256"], "source summary SHA-256", checks)
    assert_equal(sha256_file(AUDIT_TSV), success["site_audit_sha256"], "source site audit SHA-256", checks)
    assert_equal(sha256_file(M2_PASS_TSV), success["m2_pass_cases_sha256"], "source M2 cases SHA-256", checks)

    audit = pd.read_csv(AUDIT_TSV, sep="\t", dtype={"chrom": str, "ref": str, "alt": str})
    assert_equal(len(audit), 50432, "all singleton audit rows", checks)
    assert_equal(int(audit["component_size"].eq(1).sum()), 50432, "all component_size=1", checks)
    assert_equal(int(audit["positional_contract_pass"].map(parse_bool).sum()), 50432, "all positional contract pass", checks)

    dataset_rows: List[Dict[str, Any]] = []
    for dataset, expected in EXPECTED_DATASET.items():
        row = summary["breakdowns"]["dataset"][dataset]
        actual = (
            row["sites"],
            row["m1_evaluable"],
            row["m1_flagged"],
            row["m2_pass"],
            row["m2_fail"],
            row["m2_not_evaluable"],
            row["m2_not_run"],
        )
        assert_equal(actual, expected, f"{dataset} summary counts", checks)
        m1_low, m1_high = wilson_interval(row["m1_flagged"], row["sites"])
        m2_low, m2_high = wilson_interval(row["m2_pass"], row["sites"])
        dataset_rows.append(
            {
                "dataset": dataset,
                **row,
                "m1_rate": row["m1_flagged"] / row["sites"],
                "m1_count_rate": (
                    f"{row['m1_flagged']:,} ({row['m1_flagged'] / row['sites'] * 100:.2f}%)"
                ),
                "m1_evaluable_rate": row["m1_evaluable"] / row["sites"],
                "m2_rate": row["m2_pass"] / row["sites"],
                "m2_per_100k": row["m2_pass"] / row["sites"] * 100000,
                "m1_ci_low": m1_low,
                "m1_ci_high": m1_high,
                "m2_ci_low": m2_low,
                "m2_ci_high": m2_high,
                "m2_determinate": row["m2_pass"] + row["m2_fail"],
            }
        )
    sums = {
        key: sum(int(row[key]) for row in summary["breakdowns"]["dataset"].values())
        for key in ["sites", "m1_evaluable", "m1_flagged", "m2_pass", "m2_fail", "m2_not_evaluable", "m2_not_run"]
    }
    assert_equal(
        tuple(sums[k] for k in ["sites", "m1_evaluable", "m1_flagged", "m2_pass", "m2_fail", "m2_not_evaluable", "m2_not_run"]),
        (50432, 48347, 5961, 30, 18, 5913, 44471),
        "seven-dataset totals",
        checks,
    )
    total_m1_low, total_m1_high = wilson_interval(sums["m1_flagged"], sums["sites"])
    total_m2_low, total_m2_high = wilson_interval(sums["m2_pass"], sums["sites"])
    dataset_rows.append(
        {
            "dataset": "ALL_7_DATASET_ROWS",
            **sums,
            "m1_rate": sums["m1_flagged"] / sums["sites"],
            "m1_count_rate": (
                f"{sums['m1_flagged']:,} ({sums['m1_flagged'] / sums['sites'] * 100:.2f}%)"
            ),
            "m1_evaluable_rate": sums["m1_evaluable"] / sums["sites"],
            "m2_rate": sums["m2_pass"] / sums["sites"],
            "m2_per_100k": sums["m2_pass"] / sums["sites"] * 100000,
            "m1_ci_low": total_m1_low,
            "m1_ci_high": total_m1_high,
            "m2_ci_low": total_m2_low,
            "m2_ci_high": total_m2_high,
            "m2_determinate": sums["m2_pass"] + sums["m2_fail"],
        }
    )

    hcc = audit.loc[audit["dataset"].eq("HCC1395")].copy()
    assert_equal(len(hcc), 8279, "HCC1395 singleton rows", checks)
    keys = hcc[["chrom", "pos", "ref", "alt"]].astype(str).agg("|".join, axis=1)
    assert_equal(int(keys.nunique()), 8279, "HCC1395 unique variant keys", checks)
    truth_counts = hcc["truth_label"].value_counts().to_dict()
    assert_equal(
        (truth_counts.get("TP", 0), truth_counts.get("FP", 0), truth_counts.get("UNASSESSED", 0)),
        (7242, 153, 884),
        "HCC1395 truth breakdown",
        checks,
    )

    vcf_recount = recompute_vcf_components(HCC_VCF)
    assert_equal(vcf_recount["autosomal_biallelic"], 79687, "VCF HCC1395 PASS autosomal biallelic", checks)
    assert_equal(vcf_recount["positional_singleton"], 8279, "VCF independent singleton recount", checks)
    assert_equal(vcf_recount["multilocus_components"], 8222, "VCF multilocus component recount", checks)
    assert_equal(vcf_recount["partition_sites"], 79687, "VCF positional partition conservation", checks)

    hcc_m1 = hcc["m1_status"].eq("FLAGGED")
    hcc_m2_pass = hcc["m2_status"].eq("PASS")
    assert_equal(int(hcc_m1.sum()), 734, "HCC1395 M1 flagged recount", checks)
    assert_equal(int(hcc_m2_pass.sum()), 2, "HCC1395 M2 PASS recount", checks)
    assert_equal(int(hcc["m2_status"].eq("FAIL").sum()), 0, "HCC1395 M2 FAIL recount", checks)
    assert_equal(int(hcc["m2_status"].eq("NOT_EVALUABLE").sum()), 732, "HCC1395 M2 not evaluable", checks)
    assert_equal(int(hcc["m2_status"].eq("NOT_RUN").sum()), 7545, "HCC1395 M2 not run", checks)
    hcc.to_csv(
        args.output_dir / "hcc1395_singleton_site_audit.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    site_rows = read_target_site_rows()
    assignments = read_target_assignments()
    target_audit = hcc.loc[hcc_m2_pass].copy()
    gene_map = gene_context([(key[1], key[2]) for key in sorted(TARGETS)])
    loci: List[Dict[str, Any]] = []
    for key in sorted(TARGETS):
        _, chrom, pos, ref, alt = key
        site_row = site_rows.loc[
            site_rows["dataset"].eq(key[0])
            & site_rows["chrom"].eq(chrom)
            & site_rows["pos"].eq(pos)
            & site_rows["ref"].eq(ref)
            & site_rows["alt"].eq(alt)
        ].iloc[0]
        audit_row = target_audit.loc[
            target_audit["chrom"].eq(chrom)
            & target_audit["pos"].eq(pos)
            & target_audit["ref"].eq(ref)
            & target_audit["alt"].eq(alt)
        ].iloc[0]
        loci.append(
            analyze_locus(
                key,
                assignments[key],
                site_row.to_dict(),
                audit_row.to_dict(),
                gene_map[(chrom, pos)],
                checks,
            )
        )

    hcc_dorado = audit.loc[audit["dataset"].eq("HCC1395_DORADO")].copy()
    join_cols = ["chrom", "pos", "ref", "alt"]
    pair = hcc.merge(hcc_dorado, on=join_cols, suffixes=("_hcc", "_dorado"))
    common = len(pair)
    m1_h = pair["m1_status_hcc"].eq("FLAGGED")
    m1_d = pair["m1_status_dorado"].eq("FLAGGED")
    m2_h = pair["m2_status_hcc"].eq("PASS")
    m2_d = pair["m2_status_dorado"].eq("PASS")
    pair_row = {
        "common_exact_singletons": common,
        "m1_status_agreement": float((m1_h == m1_d).mean()) if common else math.nan,
        "both_m1_flagged": int((m1_h & m1_d).sum()),
        "either_m2_pass": int((m2_h | m2_d).sum()),
        "both_m2_pass": int((m2_h & m2_d).sum()),
        "hcc_only_m2_pass": int((m2_h & ~m2_d).sum()),
        "dorado_only_m2_pass": int((~m2_h & m2_d).sum()),
        "interpretation": "technical singleton-state comparison; same cell line, no exact M2-PASS replication",
    }
    assert_equal(pair_row["both_m2_pass"], 0, "HCC pair exact M2 PASS overlap", checks)

    hcc_m1_low, hcc_m1_high = wilson_interval(734, 8279)
    hcc_m2_low, hcc_m2_high = wilson_interval(2, 8279)
    headline = {
        "hcc_pass_ssnv": 79687,
        "all_positional_components": 16501,
        "singleton_sites": 8279,
        "singleton_share_of_hcc_pass": 8279 / 79687,
        "m1_evaluable": 8074,
        "m1_flagged": 734,
        "m1_rate": 734 / 8279,
        "m1_rate_evaluable": 734 / 8074,
        "m1_ci_low": hcc_m1_low,
        "m1_ci_high": hcc_m1_high,
        "m2_pass": 2,
        "m2_rate": 2 / 8279,
        "m2_per_100k": 2 / 8279 * 100000,
        "m2_ci_low": hcc_m2_low,
        "m2_ci_high": hcc_m2_high,
        "established_clone": 0,
        "established_lineage": 0,
        "established_zero_semantics": "NOT_ESTABLISHED_BECAUSE_REQUIRED_VALIDATION_NOT_RUN",
    }

    hcc_screen = [
        {"stage": "Singleton universe", "count": 8279, "share_of_singleton": 1.0},
        {"stage": "M1 evaluable", "count": 8074, "share_of_singleton": 8074 / 8279},
        {"stage": "M1 stable multigroup", "count": 734, "share_of_singleton": 734 / 8279},
        {"stage": "M2 clear residual", "count": 2, "share_of_singleton": 2 / 8279},
    ]
    status_interpretation = {
        "PASS": "Measured-axis residual partition candidate",
        "FAIL": "Measured confound aligned",
        "NOT_EVALUABLE": "M1 flag, but M2 lacked determinate evidence",
        "NOT_RUN": "M1 not flagged; M2 was not executed",
    }
    hcc_status = []
    for status, count in [("PASS", 2), ("FAIL", 0), ("NOT_EVALUABLE", 732), ("NOT_RUN", 7545)]:
        hcc_status.append(
            {
                "status": status,
                "count": count,
                "share": count / 8279,
                "interpretation": status_interpretation[status],
            }
        )
    group_count_rows = []
    for k, count in sorted(hcc.loc[hcc_m1, "coarse_ng"].astype(int).value_counts().to_dict().items()):
        group_count_rows.append({"methyl_groups": int(k), "sites": int(count), "share_of_m1": int(count) / 734})

    distance_summary: List[Dict[str, Any]] = []
    locus_summary: List[Dict[str, Any]] = []
    gene_rows: List[Dict[str, Any]] = []
    datasets: Dict[str, List[Dict[str, Any]]] = {
        "headline": [headline],
        "hcc_screen": hcc_screen,
        "hcc_status": hcc_status,
        "dataset_summary": dataset_rows,
        "dataset_summary_chart": [x for x in dataset_rows if x["dataset"] != "ALL_7_DATASET_ROWS"],
        "hcc_pair": [pair_row],
        "hcc_group_count": group_count_rows,
    }
    for i, locus in enumerate(loci, start=1):
        short = locus["locus"].split(" ")[0].replace(",", "")
        for row in locus["distance_rows"]:
            distance_summary.append({**row, "locus_short": short})
        site = locus["site"]
        audit_row = locus["audit"]
        locus_summary.append(
            {
                "locus": locus["locus"],
                "truth_label": site["truth_label"],
                "caller_af": finite_or_none(site["caller_af"]),
                "caller_dp": int(site["caller_dp"]),
                "caller_gt": str(site["caller_gt"]),
                "caller_context": (
                    f"AF={float(site['caller_af']):.3f}; DP={int(site['caller_dp'])}; "
                    f"GT={site['caller_gt']}"
                ),
                "normal_af": finite_or_none(site["normal_af"]),
                "normal_dp": int(site["normal_dp"]),
                "core_reads": locus["core_read_n"],
                "group_counts": f"{locus['group_counts']['Group A']}/{locus['group_counts']['Group B']}",
                "minor_group_share": locus["group_counts"]["Group B"] / locus["core_read_n"],
                "between_median": locus["between_stats"]["median"],
                "within_median": locus["within_stats"]["median"],
                "distance_ratio": locus["between_within_ratio"],
                "methyl_delta_b_minus_a": locus["global_methyl_delta_b_minus_a"],
                "max_cpg_delta": locus["max_cpg_delta"]["delta_b_minus_a"],
                "nearest_gap_bp": int(audit_row["nearest_gap_bp"]),
                "m2_status": audit_row["m2_status"],
                "g1_status": audit_row["g1_status"],
                "g2_status": audit_row["g2_status"],
                "r1_status": audit_row["r1_status"],
                "latest_hp": ", ".join(f"{k}:{v}" for k, v in locus["latest_hp_counts"].items()),
                "latest_ps": ", ".join(f"{k}:{v}" for k, v in locus["latest_ps_counts"].items()),
                "interpretation": locus["pattern_description"],
            }
        )
        gene = locus["gene"]
        gene_rows.append(
            {
                "locus": locus["locus"],
                "gene": gene["gene_name"],
                "gene_id": gene["gene_id"],
                "gene_type": gene["gene_type"],
                "relationship": gene["relationship"],
                "distance_bp": gene["distance_bp"],
                "cgc_member": "yes" if gene["cgc_member"] else "no",
                "approved_antineoplastic_dgidb": "yes" if gene["approved_antineoplastic_dgidb"] else "no",
                "cancer_drug_context": (
                    f"CGC={'yes' if gene['cgc_member'] else 'no'}; "
                    f"DGIdb approved antineoplastic="
                    f"{'yes' if gene['approved_antineoplastic_dgidb'] else 'no'}"
                ),
                "interpretation": "context only; not clone/topology truth",
            }
        )
        prefix = f"locus_{i}"
        datasets[f"{prefix}_distance_heatmap"] = locus["distance_heatmap"]
        datasets[f"{prefix}_shared_heatmap"] = locus["shared_heatmap"]
        datasets[f"{prefix}_methyl_heatmap"] = locus["methyl_heatmap"]
        datasets[f"{prefix}_geometry"] = locus["geometry"]
    datasets["distance_summary"] = distance_summary
    datasets["locus_summary"] = locus_summary
    datasets["gene_context"] = gene_rows
    claim_ladder = [
        {
            "axis": "Methylation",
            "level": "M1 stable multigroup",
            "evidence_level": "Methylation / M1 stable",
            "status": "COMPLETE",
            "count": 734,
            "allowed_claim": "Operational stable focal-ALT methyl groups",
            "forbidden_claim": "clone count / prevalence",
        },
        {
            "axis": "Methylation",
            "level": "M2 measured-axis residual",
            "evidence_level": "Methylation / M2 residual",
            "status": "PASS at 2 loci",
            "count": 2,
            "allowed_claim": "Read-level residual epigenetic partition candidate",
            "forbidden_claim": "all confounds excluded / cellular clone",
        },
        {
            "axis": "Genetic",
            "level": "G1/G2 co-segregation",
            "evidence_level": "Genetic / G1–G2",
            "status": "NOT_RUN (singleton)",
            "count": 0,
            "allowed_claim": "No claim; not measured",
            "forbidden_claim": "genetically corroborated subclone",
        },
        {
            "axis": "Cellular/lineage",
            "level": "R1 + orthogonal cellular truth",
            "evidence_level": "Cellular / lineage",
            "status": "NOT_RUN",
            "count": 0,
            "allowed_claim": "No claim; not measured",
            "forbidden_claim": "clone/subclone / ancestry / unique true tree",
        },
    ]
    datasets["claim_ladder"] = claim_ladder

    provenance = [
        {
            "artifact": "positional_singleton_audit_summary.json",
            "sha256": sha256_file(SUMMARY_PATH),
            "role": "source-attested counts and contracts",
        },
        {
            "artifact": "positional_singleton_site_audit.tsv.gz",
            "sha256": sha256_file(AUDIT_TSV),
            "role": "complete 50,432-site audit",
        },
        {
            "artifact": "positional_singleton_m2_pass_cases.tsv",
            "sha256": sha256_file(M2_PASS_TSV),
            "role": "all 30 M2 PASS rows",
        },
        {
            "artifact": "all_ssnv_site_results.tsv.gz",
            "sha256": summary["inputs"]["site_results"]["sha256"],
            "role": "site-level caller/VAF and measured-axis data",
        },
        {
            "artifact": "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
            "sha256": summary["inputs"]["stable_assignments"]["sha256"],
            "role": "exact core-read identities and labels",
        },
    ]
    for locus in loci:
        for role, digest in locus["artifact_hashes"].items():
            provenance.append({"artifact": f"{locus['locus']} {role}", "sha256": digest, "role": "primary region artifact"})
    for row in provenance:
        row["sha256_short"] = f"{row['sha256'][:16]}…{row['sha256'][-8:]}"
    datasets["provenance"] = provenance
    datasets["validation"] = [
        {
            "check": row["check"],
            "actual": check_display_value(row["actual"]),
            "expected": check_display_value(row["expected"]),
            "status": "PASS" if row["pass"] else "FAIL",
        }
        for row in checks
    ]

    artifact = build_artifact(generated_at, datasets, loci, headline)
    json_dump(args.artifact, artifact)

    pd.DataFrame(dataset_rows).to_csv(args.output_dir / "all_dataset_singleton_summary.tsv", sep="\t", index=False)
    pd.DataFrame(locus_summary).to_csv(args.output_dir / "hcc1395_m2_pass_locus_summary.tsv", sep="\t", index=False)
    pd.DataFrame(distance_summary).to_csv(args.output_dir / "hcc1395_m2_pass_distance_summary.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(args.output_dir / "gene_context.tsv", sep="\t", index=False)
    pd.DataFrame(group_count_rows).to_csv(args.output_dir / "hcc1395_m1_group_count_distribution.tsv", sep="\t", index=False)
    for i, locus in enumerate(loci, start=1):
        pd.DataFrame(locus["geometry"]).to_csv(
            args.output_dir / f"locus_{i}_core_read_metadata.tsv", sep="\t", index=False
        )
        pd.DataFrame(locus["distance_heatmap"]).to_csv(
            args.output_dir / f"locus_{i}_display_distance_matrix.tsv", sep="\t", index=False
        )
        pd.DataFrame(locus["shared_heatmap"]).to_csv(
            args.output_dir / f"locus_{i}_display_shared_cpg_matrix.tsv", sep="\t", index=False
        )
        pd.DataFrame(locus["methyl_heatmap"]).to_csv(
            args.output_dir / f"locus_{i}_display_methylation_matrix.tsv", sep="\t", index=False
        )

    receipt = {
        "schema_name": "intersubmod.singleton_alt_methyl_substructure_report_receipt",
        "schema_version": "1.0.0",
        "created_at": generated_at,
        "pass": all(row["pass"] for row in checks),
        "task_type": "B_comprehensive_validation",
        "scope": {
            "all_dataset_singletons": 50432,
            "datasets": 7,
            "hcc1395_pass_ssnv": 79687,
            "hcc1395_singletons": 8279,
            "hcc1395_m1_flags": 734,
            "hcc1395_m2_pass": 2,
            "hcc1395_sites_subsampled": 0,
        },
        "claim_ceiling": "M2_read_level_residual_epigenetic_partition",
        "not_run": ["G1", "G2", "formal_R1", "matched_normal", "CN_purity_CCF", "cellular_lineage"],
        "zero_semantics": {
            "established_cellular_clone_subclone": 0,
            "meaning": "required validation not run; not a true-negative prevalence estimate",
        },
        "vcf_independent_recount": vcf_recount,
        "hcc_pair": pair_row,
        "headline": headline,
        "checks": checks,
        "outputs": {},
        "command": " ".join(sys.argv),
        "python": sys.version,
        "git_head": "0ee2fa1b31fcf6af670efd301251b5b3a24c1a99",
    }
    output_paths = [
        args.artifact,
        args.output_dir / "all_dataset_singleton_summary.tsv",
        args.output_dir / "hcc1395_singleton_site_audit.tsv.gz",
        args.output_dir / "hcc1395_m2_pass_locus_summary.tsv",
        args.output_dir / "hcc1395_m2_pass_distance_summary.tsv",
        args.output_dir / "gene_context.tsv",
        args.output_dir / "hcc1395_m1_group_count_distribution.tsv",
    ]
    for path in output_paths:
        receipt["outputs"][str(path.relative_to(REPO))] = {
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt_path = args.output_dir / "validation_receipt.json"
    json_dump(receipt_path, receipt)

    markdown = f"""<!--
建立時間: {generated_at}
目標: 完整驗證 HCC1395 positional-singleton focal-ALT read 甲基子結構
處理範圍: 7 datasets / 50,432 singleton sites；HCC1395 8,279 sites 全量；兩個 M2 PASS loci read-level exact join
關聯檔案:
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/artifact.json
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/validation_receipt.json
-->

# HCC1395 positional-singleton ALT 內甲基子結構完整驗證

> **結論：** HCC1395 8,279 個 singleton sites 中，M1 stable multigroup 為
> **734（{headline['m1_rate']*100:.2f}%）**；通過八個 measured-axis guardrails 的 M2 PASS 為
> **2（{headline['m2_rate']*100:.4f}%）**。兩點可稱 read-level residual epigenetic partition candidate；
> G1/G2/formal R1、matched-normal 與 CN/CCF 未執行，因此 **不能稱已證實 clone/subclone**。

## 分母

- HCC1395 chr1-22 autosomal biallelic LongPhase-S PASS sSNV：79,687。
- 50 kb positional components：16,501；其中 singleton components/sites：8,279，多位點 components：8,222。
- singleton 真值層：TP 7,242、FP 153、UNASSESSED 884。

## 核心數字

- M1 evaluable：8,074/8,279（{8074/8279*100:.2f}%）。
- M1 stable multigroup：734/8,279（{734/8279*100:.2f}%；95% Wilson CI
  {headline['m1_ci_low']*100:.2f}%–{headline['m1_ci_high']*100:.2f}%）。
- M2 PASS：2/8,279（{2/8279*100:.4f}%；95% Wilson CI
  {headline['m2_ci_low']*100:.4f}%–{headline['m2_ci_high']*100:.4f}%）。
- M2 FAIL：0；NOT_EVALUABLE：732；NOT_RUN：7,545。
- Established cellular clone/subclone：0，語意是必要驗證未跑，不是真陰性。

## 兩個 M2 PASS 例子

| Locus | Core ALT reads | A/B | Between median | Pooled within | Ratio | Methyl Δ B−A | Interpretation |
|---|---:|---:|---:|---:|---:|---:|---|
"""
    for row in locus_summary:
        markdown += (
            f"| {row['locus']} | {row['core_reads']} | {row['group_counts']} | "
            f"{row['between_median']:.3f} | {row['within_median']:.3f} | "
            f"{row['distance_ratio']:.2f}× | {row['methyl_delta_b_minus_a']:+.3f} | "
            f"{row['interpretation']} |\n"
        )
    markdown += """

兩點所有 core reads 的 latest HP/PS 都是同一值，因此 HP 沒有把兩群切開，但也不能提供獨立 clone 佐證。
caller AF（chr14 0.827；chr22 1.000）只作 focal allele burden context，未用於甲基分群。

## Claim ceiling

- 可用：`M1 stable focal-ALT methyl multigroup`。
- 兩個清楚例子可用：`M2 read-level residual epigenetic partition candidate`。
- 不可用：confirmed clone/subclone、clone number、parent-child、linear/branching ancestry、unique true tree。

## 可重現輸出

- 完整 HCC1395 8,279-row audit：`InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/hcc1395_singleton_site_audit.tsv.gz`
- 驗證 receipt：`InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/validation_receipt.json`
- Canonical report artifact：`InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/artifact.json`
"""
    args.markdown.write_text(markdown, encoding="utf-8")
    print(
        json.dumps(
            {
                "pass": receipt["pass"],
                "hcc1395_singletons": 8279,
                "hcc1395_m1_flags": 734,
                "hcc1395_m2_pass": 2,
                "artifact": str(args.artifact),
                "receipt": str(receipt_path),
                "markdown": str(args.markdown),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
