#!/usr/bin/env python3
"""Audit projected HP concordance against every prior latest-tag FP output read.

This regression source usually has no PS column. Passing this audit therefore
provides HP-only evidence and never substitutes for the downstream full HP/PS join.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam

from latest_tag_join import fetch_site_lookup, join_read_rows, normalize_hp


EXPECTED_FP_SITES = 7_745
AUDIT_SCOPE = "HP_ONLY"
PASS_SEMANTICS = "execution_integrity_and_hp_concordance_only_not_ps_completeness"
SITE_PATTERN = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]))_(\d+)$")
DEFAULT_FP_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_single_fp_alt_multicluster_subclone_validation/intersubmod_latest_fp_v1"
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def focal_site(reads_path: Path) -> tuple[str, int]:
    match = SITE_PATTERN.fullmatch(reads_path.parents[2].name)
    if not match:
        raise ValueError(f"Cannot derive focal site from {reads_path}")
    return match.group(1), int(match.group(2))


def reads_tsv_has_ps_column(path: Path) -> bool:
    with path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    return any(column.strip().lower() == "ps" for column in header)


def audit_sample(sample: str, run_root: str, sidecar: str, sidecar_index: str) -> dict[str, Any]:
    sample_root = Path(run_root) / sample
    paths = sorted(sample_root.rglob("reads.tsv"))
    if not paths:
        raise FileNotFoundError(f"No reads.tsv under {sample_root}")
    counts: Counter[str] = Counter()
    hp_counts: Counter[str] = Counter()
    multiplicity_counts: Counter[str] = Counter()
    with pysam.TabixFile(sidecar, index=sidecar_index) as tabix:
        for path in paths:
            chrom, pos1 = focal_site(path)
            if reads_tsv_has_ps_column(path):
                counts["reads_tsv_sites_with_ps_column"] += 1
            else:
                counts["reads_tsv_sites_without_ps_column"] += 1
            lookup = fetch_site_lookup(tabix, chrom, pos1)
            counts["sites"] += 1
            counts["sidecar_rows_fetched"] += lookup.rows_fetched
            counts["sidecar_rows_eligible"] += lookup.rows_eligible
            joined = join_read_rows(path, lookup)
            for row, tags, full_identity_count in joined:
                counts["read_rows"] += 1
                if normalize_hp(row["hp"]) != tags.hp:
                    counts["hp_mismatches"] += 1
                    if counts["hp_mismatches"] <= 10:
                        raise RuntimeError(
                            f"{sample} HP mismatch {chrom}:{pos1} {row['read_name']}: "
                            f"reads.tsv={row['hp']} sidecar={tags.hp}"
                        )
                counts["hp_matches"] += 1
                hp_counts[tags.hp] += 1
                multiplicity_counts[str(full_identity_count)] += 1
    passed = counts["read_rows"] > 0 and counts["hp_matches"] == counts["read_rows"] and not counts["hp_mismatches"]
    return {
        "dataset": sample,
        "sample": sample,
        "audit_scope": AUDIT_SCOPE,
        "pass_semantics": PASS_SEMANTICS,
        "input_reads_root": str(sample_root),
        "input_sidecar": sidecar,
        "input_sidecar_index": sidecar_index,
        "counts": dict(counts),
        "latest_hp_counts": dict(hp_counts),
        "full_identity_multiplicity_per_projected_match": dict(multiplicity_counts),
        "ps_evidence": {
            "status": "NOT_EVALUATED_HP_ONLY_AUDIT",
            "reads_tsv_sites_with_ps_column": counts["reads_tsv_sites_with_ps_column"],
            "reads_tsv_sites_without_ps_column": counts["reads_tsv_sites_without_ps_column"],
            "not_evaluable_reason": (
                "PS concordance is outside this HP-only regression audit; reads.tsv has no PS column "
                f"for {counts['reads_tsv_sites_without_ps_column']}/{counts['sites']} sites."
            ),
            "downstream_full_hp_ps_join_required": True,
        },
        "pass": passed,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fp-root", type=Path, default=DEFAULT_FP_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=7)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    jobs = [
        (
            item["sample"],
            str(args.fp_root),
            item["latest_read_tag_sidecar"]["path"],
            item["latest_read_tag_sidecar_index"]["path"],
        )
        for item in manifest["samples"]
    ]
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
        futures = {pool.submit(audit_sample, *job): job[0] for job in jobs}
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda item: item["sample"])
    totals: Counter[str] = Counter()
    for result in results:
        totals.update(result["counts"])
    passed = (
        len(results) == 7
        and totals["sites"] == EXPECTED_FP_SITES
        and totals["read_rows"] == totals["hp_matches"]
        and totals["hp_mismatches"] == 0
        and all(item["pass"] for item in results)
    )
    payload = {
        "schema_name": "intersubmod.latest_tag_projected_join_audit",
        "schema_version": "1.1.0",
        "created_at_utc": now_utc(),
        "status": "EXECUTION_PASS_HP_CONCORDANCE" if passed else "EXECUTION_FAIL",
        "pass_semantics": PASS_SEMANTICS,
        "contract": {
            "audit_scope": AUDIT_SCOPE,
            "audited_fields": ["HP"],
            "not_audited_fields": ["PS"],
            "source_tag_authority": (
                "same-run LongPhase-S recalibrated sidecar; only HP concordance is audited here"
            ),
            "projection_fields": ["QNAME", "CHROM", "START0", "END0", "MAPQ", "STRAND"],
            "inter_sub_mod_flag_exclusions": ["secondary", "supplementary", "duplicate", "unmapped"],
            "conflicting_exact_identity_tags": "hard_fail",
            "conflicting_projected_identity_tags": "hard_fail",
            "missing_reads_tsv_identity": "hard_fail",
            "truth_or_cooccurrence_used_in_join": False,
            "downstream_full_hp_ps_join_required": True,
            "downstream_responsibility": (
                "analyze_all_ssnv_focal_alt_multigroup.py must exact-join both HP and PS before ALT selection"
            ),
        },
        "evidence_limitations": [
            (
                "A passing result establishes HP concordance only. It does not establish PS presence, "
                "PS concordance, phase-set completeness, or any evidence-ladder claim."
            ),
            (
                f"reads.tsv lacks a PS column for {totals['reads_tsv_sites_without_ps_column']}/"
                f"{totals['sites']} audited sites; PS is formally NOT_EVALUATED in this artifact."
            ),
        ],
        "input_manifest": str(args.manifest.resolve()),
        "input_fp_root": str(args.fp_root.resolve()),
        "totals": dict(totals),
        "datasets": results,
        "samples": results,
        "pass": passed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit: {args.output}")
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output.resolve()), "totals": dict(totals), "pass": passed}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
