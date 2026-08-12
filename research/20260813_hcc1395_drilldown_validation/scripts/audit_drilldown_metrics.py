#!/usr/bin/env python3
"""Audit one or more drilldown bundles and the comparable topology cohort.

This is deliberately read-only with respect to the input bundles.  It emits
machine-readable tables so every percentage in the accompanying report has an
explicit numerator, denominator, unit, scope, and formula.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


CHROMS = [f"chr{i}" for i in range(1, 23)]
AXES = {
    "hp": ("HPMergedDelta", "delta_beta"),
    "hpfine": ("HPFineF", "stored_HP_fine_F_definition_not_embedded"),
    "allele": ("AlleleDelta", "delta_beta"),
    "cluster": ("ClusterPermanovaF", "pseudo_F_double_dipping_do_not_use_as_evidence"),
    "lineage": ("LineagePermanovaF", "pseudo_F"),
}


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for row in rows for k in row})
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _percent(num: int, den: int) -> float | None:
    return round(100.0 * num / den, 6) if den else None


def _quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    values = sorted(values)
    x = (len(values) - 1) * p
    lo, hi = math.floor(x), math.ceil(x)
    if lo == hi:
        return values[lo]
    return values[lo] * (hi - x) + values[hi] * (x - lo)


def _bh_qvalues(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p values; ties are retained."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [1.0] * n
    running = 1.0
    for rank0 in range(n - 1, -1, -1):
        idx = order[rank0]
        running = min(running, pvalues[idx] * n / (rank0 + 1))
        q[idx] = min(1.0, running)
    return q


def _parse_shard(path: Path, chrom: str) -> tuple[list[dict], dict[str, dict]]:
    text = path.read_text(encoding="utf-8")
    l2_tag = f'window.__DD.L2["{chrom}"]='
    l4_tag = f'window.__DD.L4["{chrom}"]='
    start2 = text.index(l2_tag) + len(l2_tag)
    start4_tag = text.index(";" + l4_tag, start2)
    l2 = json.loads(text[start2:start4_tag])
    start4 = start4_tag + 1 + len(l4_tag)
    l4 = json.loads(text[start4:].rstrip(";\n"))
    return l2, l4


def _selfcheck_status(summary: dict) -> str:
    if summary.get("fail", 0):
        return "BLOCKED"
    if summary.get("skip", 0):
        return "INCOMPLETE"
    return "INTERNAL_QA_PASS"


def _input_integrity(receipt: dict) -> list[dict]:
    rows = []
    for rec in receipt.get("inputs") or []:
        path = Path(rec.get("path", ""))
        kind = "file" if path.is_file() else "directory" if path.is_dir() else "missing"
        actual_size = path.stat().st_size if path.is_file() else None
        expected_sha = rec.get("sha256") or ""
        actual_sha = _sha256(path) if path.is_file() and expected_sha else ""
        if kind == "missing":
            status = "MISSING"
        elif path.is_dir() and not expected_sha:
            status = "UNVERIFIABLE_DIRECTORY"
        elif expected_sha and actual_sha != expected_sha:
            status = "HASH_MISMATCH"
        elif expected_sha:
            status = "HASH_MATCH"
        else:
            status = "NO_HASH"
        rows.append({
            "path": str(path), "kind_now": kind, "status": status,
            "expected_size_bytes": rec.get("size_bytes"), "actual_size_bytes": actual_size,
            "expected_sha256": expected_sha, "actual_sha256": actual_sha,
        })
    return rows


def audit_bundle(bundle: Path, label: str) -> tuple[dict, list[dict], list[dict], list[dict]]:
    receipt = json.loads((bundle / "receipt.json").read_text(encoding="utf-8"))
    regions: list[dict] = []
    methyl: dict[tuple[str, int], dict] = {}
    for chrom in CHROMS:
        l2, l4 = _parse_shard(bundle / "data" / f"L2.{chrom}.js", chrom)
        regions.extend(l2)
        for pos, row in l4.items():
            methyl[(chrom, int(pos))] = row

    site_k: dict[tuple[str, int], int] = {}
    site_regions: Counter = Counter()
    for row in regions:
        for pos in row.get("ap") or []:
            key = (row["c"], int(pos))
            site_k[key] = max(site_k.get(key, 0), int(row.get("k") or 0))
            site_regions[key] += 1

    tree = [r for r in regions if r.get("v")]
    unique = [r for r in tree if r.get("tu")]
    by_status = Counter(r.get("us") or "(missing)" for r in regions)
    by_family = Counter(r.get("fs") or "(missing)" for r in regions)
    by_objective = Counter(r.get("os") or "(missing)" for r in regions)
    summary = receipt.get("selfcheck") or {}

    pngs = list((bundle / "panels").glob("chr*/*.png"))
    igv = list((bundle / "igv").glob("*.js"))
    l2_files = list((bundle / "data").glob("L2.chr*.js"))
    l5_files = list((bundle / "data").glob("L5.chr*.js"))
    asset_rows = []
    for typ, files in (("png_all", pngs), ("png_base", [p for p in pngs if not p.name.endswith(".T.png")]),
                       ("png_tumor_only", [p for p in pngs if p.name.endswith(".T.png")]),
                       ("igv_js", igv), ("L2_plus_L4_shard", l2_files), ("L5_manifest_shard", l5_files)):
        asset_rows.append({"bundle": label, "asset_type": typ, "count": len(files),
                           "bytes": sum(p.stat().st_size for p in files)})

    panel_receipt_bytes = ((receipt.get("panels") or {}).get("bytes"))
    base_png_bytes = next(r["bytes"] for r in asset_rows if r["asset_type"] == "png_base")
    all_png_bytes = next(r["bytes"] for r in asset_rows if r["asset_type"] == "png_all")
    overview = {
        "bundle": label,
        "path": str(bundle),
        "built_at": receipt.get("built_at"),
        "sample": receipt.get("sample"),
        "receipt_schema_name": receipt.get("schema_name"),
        "receipt_schema_version": receipt.get("schema_version"),
        "generator_commit": ((receipt.get("generator") or {}).get("git_commit")),
        "claim_ceiling": receipt.get("claim_ceiling") or "NOT_DECLARED",
        "validation_status_recomputed": _selfcheck_status(summary),
        "selfcheck_pass": summary.get("pass"),
        "selfcheck_fail": summary.get("fail"),
        "selfcheck_skip": summary.get("skip"),
        "regions": len(regions),
        "regions_with_tree": len(tree),
        "tree_coverage_pct_all_regions": _percent(len(tree), len(regions)),
        "unique_best_tree": len(unique),
        "unique_pct_among_tree": _percent(len(unique), len(tree)),
        "unique_pct_all_regions": _percent(len(unique), len(regions)),
        "distinct_ssnv": len(site_k),
        "region_position_pairs": sum(site_regions.values()),
        "multi_region_ssnv": sum(v > 1 for v in site_regions.values()),
        "methyl_summary_rows": len(methyl),
        "methyl_topology_linked": len(set(methyl) & set(site_k)),
        "methyl_topology_linkage_pct": _percent(len(set(methyl) & set(site_k)), len(site_k)),
        "panel_receipt_bytes": panel_receipt_bytes,
        "panel_base_actual_bytes": base_png_bytes,
        "panel_all_actual_bytes": all_png_bytes,
        "panel_receipt_matches_base_only": panel_receipt_bytes == base_png_bytes,
        "panel_unreported_bytes": all_png_bytes - (panel_receipt_bytes or 0),
        "unit_status": dict(sorted(by_status.items())),
        "family_status": dict(sorted(by_family.items())),
        "objective_status": dict(sorted(by_objective.items())),
    }

    coverage_rows = []
    bins = [(str(i), lambda k, i=i: k == i) for i in range(1, 8)] + [("8+", lambda k: k >= 8)]
    meth_keys = set(methyl)
    for name, predicate in bins:
        den = sum(predicate(k) for k in site_k.values())
        num = sum(predicate(k) and site in meth_keys for site, k in site_k.items())
        coverage_rows.append({
            "bundle": label, "metric": "methyl_summary_topology_linkage", "k_bin": name,
            "numerator": num, "denominator": den, "percent": _percent(num, den),
            "unit": "distinct genomic sSNV", "formula": "linked distinct sSNV / topology distinct sSNV",
        })

    axis_rows = []
    for scope, positions in (("all_summary_rows", set(methyl)),
                             ("topology_linked_distinct_ssnv", set(methyl) & set(site_k))):
        for axis, (source_effect, effect_unit) in AXES.items():
            vals = []
            effects = []
            raw_sig = 0
            for pos in positions:
                row = methyl[pos]
                p = row.get(axis + "_p")
                valid = row.get(axis + "_valid")
                ngroup = row.get(axis + "_n")
                if p is None or valid is False or (ngroup is not None and ngroup < 2):
                    continue
                p = float(p)
                vals.append(p)
                raw_sig += p <= 0.05
                if row.get(axis + "_eff") is not None:
                    effects.append(float(row[axis + "_eff"]))
            qvals = _bh_qvalues(vals)
            fdr_sig = sum(q <= 0.05 for q in qvals)
            axis_rows.append({
                "bundle": label, "scope": scope, "axis": axis,
                "source_effect_field": source_effect, "effect_unit": effect_unit,
                "tested_n": len(vals), "raw_p_le_0_05_n": raw_sig,
                "raw_p_le_0_05_pct": _percent(raw_sig, len(vals)),
                "bh_fdr_q_le_0_05_n": fdr_sig,
                "bh_fdr_q_le_0_05_pct": _percent(fdr_sig, len(vals)),
                "multiplicity_family": f"within axis × {scope}",
                "effect_n": len(effects), "effect_median": median(effects) if effects else None,
                "effect_q1": _quantile(effects, 0.25), "effect_q3": _quantile(effects, 0.75),
                "effect_min": min(effects) if effects else None,
                "effect_max": max(effects) if effects else None,
                "interpretation_gate": ("DOUBLE_DIPPING_INVALID_AS_EVIDENCE" if axis == "cluster"
                                        else "EXPLORATORY_NOT_TRUTH_VALIDATED"),
            })

    reads = [int(float(r["NumReads"])) for r in methyl.values() if r.get("NumReads") not in (None, "")]
    cpgs = [int(float(r["NumCpGs"])) for r in methyl.values() if r.get("NumCpGs") not in (None, "")]
    disp = sum(str(r.get("ClusterDispersionWarn", "")).lower() == "true" for r in methyl.values())
    overview.update({
        "methyl_num_reads_median": median(reads) if reads else None,
        "methyl_num_reads_q1": _quantile(reads, 0.25), "methyl_num_reads_q3": _quantile(reads, 0.75),
        "methyl_num_cpgs_median": median(cpgs) if cpgs else None,
        "methyl_num_cpgs_q1": _quantile(cpgs, 0.25), "methyl_num_cpgs_q3": _quantile(cpgs, 0.75),
        "cluster_dispersion_warning_n": disp,
        "cluster_dispersion_warning_pct": _percent(disp, len(methyl)),
    })
    return overview, asset_rows, coverage_rows, axis_rows, _input_integrity(receipt)


def audit_topology_cohort(root: Path) -> list[dict]:
    rows = []
    for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        sample = sample_dir.name
        topo = sample_dir / f"{sample}.topology.jsonl"
        receipt_path = sample_dir / f"{sample}.topology.receipt.json"
        pipeline_path = sample_dir / "pipeline_receipt.json"
        if not topo.is_file() or not receipt_path.is_file():
            continue
        status = Counter()
        family = Counter()
        objective = Counter()
        sites = set()
        schemas = set()
        versions = set()
        samples = set()
        models = set()
        af_basis = set()
        n = tree = unique = pairs = 0
        with topo.open(encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                rec = json.loads(line)
                n += 1
                status[rec.get("unit_status")] += 1
                family[rec.get("family_status")] += 1
                objective[rec.get("objective_status")] += 1
                schemas.add(rec.get("schema_name"))
                versions.add(rec.get("schema_version"))
                samples.add(rec.get("sample"))
                models.add(rec.get("model"))
                af_basis.add(rec.get("af_basis"))
                positions = rec.get("active_positions") or []
                pairs += len(positions)
                sites.update((rec.get("chrom"), int(pos)) for pos in positions)
                if rec.get("representative_best_vertices"):
                    tree += 1
                    unique += bool(rec.get("best_tree_unique"))
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        pipeline = json.loads(pipeline_path.read_text(encoding="utf-8")) if pipeline_path.is_file() else {}
        rows.append({
            "sample": sample, "biological_id": "HCC1395" if sample == "HCC1395_DORADO" else sample,
            "technical_replicate": sample == "HCC1395_DORADO",
            "topology_sha256_actual": _sha256(topo),
            "topology_sha256_receipt": (receipt.get("output") or {}).get("sha256"),
            "hash_match": _sha256(topo) == (receipt.get("output") or {}).get("sha256"),
            "sample_values": "|".join(str(x) for x in sorted(samples, key=str)),
            "sample_identity_all_match": samples == {sample},
            "schema_names": "|".join(str(x) for x in sorted(schemas, key=str)),
            "schema_versions": "|".join(str(x) for x in sorted(versions, key=str)),
            "models": "|".join(str(x) for x in sorted(models, key=str)),
            "af_basis": "|".join(str(x) for x in sorted(af_basis, key=str)),
            "region_n": n, "distinct_ssnv_n": len(sites), "region_position_pair_n": pairs,
            "tree_n": tree, "tree_pct_all_regions": _percent(tree, n),
            "unique_tree_n": unique, "unique_pct_among_tree": _percent(unique, tree),
            "unique_pct_all_regions": _percent(unique, n),
            "ranked_n": status.get("ranked", 0),
            "family_complete_n": family.get("FAMILY_COMPLETE", 0),
            "objective_certified_n": objective.get("OBJECTIVE_CERTIFIED", 0),
            "receipt_all_pass": receipt.get("all_pass"),
            "receipt_all_mutation_bearing_families_complete": receipt.get("all_mutation_bearing_families_complete"),
            "receipt_all_units_objective_certified": receipt.get("all_units_objective_certified"),
            "pipeline_claim_ceiling": pipeline.get("claim_ceiling"),
        })
    return rows


def _record_sample(label: str, rec: dict) -> str | None:
    explicit = rec.get("sample")
    if explicit:
        return str(explicit)
    paths = []
    for key in ("in_bam", "out_bam", "topology", "assignments", "lineage_paths", "sidecar"):
        value = rec.get(key)
        if isinstance(value, str):
            paths.append(value)
        elif isinstance(value, dict):
            paths.extend(str(v) for v in value.values() if isinstance(v, str))
    matches = sorted({m.group(1) for path in paths
                      for m in re.finditer(r"(?:^|[/_.-])(HCC1395_DORADO|HCC1395|HCC1937|HCC1954|H1437|H2009|COLO829)(?:[/_.-]|$)", path)})
    return matches[0] if len(matches) == 1 else None


def audit_lca(pre_dir: Path | None, post_dir: Path | None, expected_sample: str) -> list[dict]:
    if not pre_dir or not post_dir or not pre_dir.is_dir() or not post_dir.is_dir():
        return []

    def load_dir(path: Path) -> dict[str, tuple[Path, dict]]:
        out = {}
        for file in sorted(path.glob("*.json")):
            match = re.search(r"\.(chr[0-9XYM]+)\.", file.name)
            if not match:
                continue
            out[match.group(1)] = (file, json.loads(file.read_text(encoding="utf-8")))
        return out

    pre, post = load_dir(pre_dir), load_dir(post_dir)
    rows = []
    for chrom in sorted(set(pre) & set(post), key=lambda c: (len(c), c)):
        pre_file, a = pre[chrom]
        post_file, b = post[chrom]
        sa, sb = a.get("stats") or {}, b.get("stats") or {}
        non_lca = sorted((set(sa) | set(sb)) - {"lv_written", "lca_resolved", "lca_candidates_sum"})
        rows.append({
            "chrom": chrom,
            "expected_sample": expected_sample,
            "pre_file": str(pre_file), "post_file": str(post_file),
            "pre_sample_inferred": _record_sample("pre", a),
            "post_sample_inferred": _record_sample("post", b),
            "sample_identity_pass": (_record_sample("pre", a) == expected_sample
                                     and _record_sample("post", b) == expected_sample),
            "same_in_bam": a.get("in_bam") == b.get("in_bam"),
            "same_threads": a.get("threads") == b.get("threads"),
            "non_lca_stats_equal_missing_as_zero": all(int(sa.get(k) or 0) == int(sb.get(k) or 0) for k in non_lca),
            "pre_in_bam": a.get("in_bam"), "post_in_bam": b.get("in_bam"),
            "pre_threads": a.get("threads"), "post_threads": b.get("threads"),
            "pre_lv_written": int(sa.get("lv_written") or 0),
            "post_lv_written": int(sb.get("lv_written") or 0),
            "net_new_lv_written": int(sb.get("lv_written") or 0) - int(sa.get("lv_written") or 0),
            "lca_candidates_sum": int(sb.get("lca_candidates_sum") or 0),
            "lca_resolved": int(sb.get("lca_resolved") or 0),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", action="append", required=True,
                    help="LABEL=/absolute/path (repeatable)")
    ap.add_argument("--topology-root", type=Path, required=True)
    ap.add_argument("--lca-pre", type=Path)
    ap.add_argument("--lca-post", type=Path)
    ap.add_argument("--lca-sample", default="HCC1395")
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    overview_rows, asset_rows, coverage_rows, axis_rows, integrity_rows = [], [], [], [], []
    for spec in args.bundle:
        label, raw_path = spec.split("=", 1)
        result = audit_bundle(Path(raw_path), label)
        overview, assets, coverage, axes, integrity = result
        overview_rows.append(overview)
        asset_rows.extend(assets)
        coverage_rows.extend(coverage)
        axis_rows.extend(axes)
        integrity_rows.extend({"bundle": label, **r} for r in integrity)

    cohort = audit_topology_cohort(args.topology_root)
    lca_rows = audit_lca(args.lca_pre, args.lca_post, args.lca_sample)
    _write_csv(args.output_dir / "bundle_overview.csv", overview_rows)
    _write_csv(args.output_dir / "asset_inventory.csv", asset_rows)
    _write_csv(args.output_dir / "input_integrity.csv", integrity_rows)
    _write_csv(args.output_dir / "methylation_coverage_by_k.csv", coverage_rows)
    _write_csv(args.output_dir / "methylation_axis_metrics.csv", axis_rows)
    _write_csv(args.output_dir / "cohort_topology_metrics.csv", cohort)
    if lca_rows:
        _write_csv(args.output_dir / "lca_ab_by_chrom.csv", lca_rows)
    payload = {
        "schema_name": "intersubmod.drilldown.audit.metrics",
        "schema_version": "1.0.0",
        "claim_ceiling": "descriptive observation and internal data-product QA; not truth-set validation",
        "bundle_overview": overview_rows,
        "asset_inventory": asset_rows,
        "input_integrity": integrity_rows,
        "methylation_coverage_by_k": coverage_rows,
        "methylation_axis_metrics": axis_rows,
        "cohort_topology_metrics": cohort,
        "lca_ab_by_chrom": lca_rows,
        "lca_ab_summary": ({
            "pre_file_n": len(list(args.lca_pre.glob("*.json"))),
            "post_file_n": len(list(args.lca_post.glob("*.tag_bam.receipt.json"))),
            "shared_chrom_n": len(lca_rows),
            "same_in_bam_n": sum(r["same_in_bam"] for r in lca_rows),
            "same_threads_n": sum(r["same_threads"] for r in lca_rows),
            "sample_identity_pass_n": sum(r["sample_identity_pass"] for r in lca_rows),
            "non_lca_stats_equal_n": sum(r["non_lca_stats_equal_missing_as_zero"] for r in lca_rows),
            "pre_lv_written": sum(r["pre_lv_written"] for r in lca_rows),
            "post_lv_written": sum(r["post_lv_written"] for r in lca_rows),
            "net_new_lv_written": sum(r["net_new_lv_written"] for r in lca_rows),
            "lca_resolved": sum(r["lca_resolved"] for r in lca_rows),
            "lca_resolved_minus_net_new_lv": (sum(r["lca_resolved"] for r in lca_rows)
                                               - sum(r["net_new_lv_written"] for r in lca_rows)),
            "causal_ab_gate": "FAIL" if any(not r["same_in_bam"] or not r["same_threads"] for r in lca_rows) else "PASS",
        } if lca_rows else None),
        "metric_contract": {
            "distinct_ssnv": "unique (chrom, 1-based position) in topology active_positions",
            "region_position_pair": "one active position membership in one topology region",
            "tree_coverage": "regions with representative_best_vertices / all regions",
            "unique_rate_among_tree": "best_tree_unique regions / regions with a representative tree",
            "bh_fdr": "Benjamini-Hochberg within one axis and one declared scope; q <= 0.05",
            "macro_vs_micro": "cohort CSV is per sample; do not pool until sample-level macro and biological replicate policy are declared",
        },
    }
    (args.output_dir / "metrics_audit.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "bundles": len(overview_rows), "cohort_datasets": len(cohort),
        "outputs": sorted(p.name for p in args.output_dir.iterdir()),
        "validation_status": {r["bundle"]: r["validation_status_recomputed"] for r in overview_rows},
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
