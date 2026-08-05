#!/usr/bin/env python3
"""Independent conservation checks for the HCC1395 PS-orientation audit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--neutral-probe", type=Path, required=True)
    parser.add_argument("--panel-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    census = load(args.census)
    neutral = load(args.neutral_probe)
    panel = load(args.panel_summary)
    counts = census["counts"]
    ps_class = counts["ps_class"]
    n_ps = {int(key): int(value) for key, value in counts["n_ps_distribution"].items()}
    mixed = counts["mixed_ps_thresholds"]
    join = counts["read_tag_join"]
    config_space = counts["orientation_configuration_space"]
    neutral_native = neutral["configurations"][0]["region"]
    neutral_flip = neutral["configurations"][1]["region"]

    checks = {
        "census_self_checks_pass": census.get("all_checks_pass") is True,
        "W_tree_ps_partition": sum(int(value) for value in ps_class.values()) == counts["W_tree"] == 8222,
        "W_tree_nps_partition": sum(n_ps.values()) == counts["W_tree"],
        "mixed_ps_count": int(ps_class["multiple"]) == mixed["all"] == 865,
        "mixed_primary_partition": (
            mixed["current_complete_primary"] + mixed["current_incomplete"] == mixed["current_primary"] == 861
        ),
        "mixed_all_partition": (
            mixed["current_complete_primary"]
            + mixed["current_incomplete"]
            + mixed["current_no_primary"]
            == mixed["all"]
        ),
        "orientation_space_including_native": (
            sum(count * (2 ** (blocks - 1)) for blocks, count in n_ps.items() if blocks > 1)
            == config_space["including_native"]
            == 1802
        ),
        "orientation_space_alternatives": (
            sum(count * ((2 ** (blocks - 1)) - 1) for blocks, count in n_ps.items() if blocks > 1)
            == config_space["alternative_flips_only"]
            == 937
        ),
        "exact_read_tag_join": (
            join["alignment_group_exposures"] == join["sidecar_exact_matches"] == 1934226
            and join["sidecar_missing"] == 0
            and join["sidecar_conflicts"] == 0
            and join["alignment_identity_allele_conflicts"] == 0
        ),
        "neutral_region_identity": neutral["scope"]["region"] == "chr20:29077065-29104598",
        "neutral_native_reproduces_canonical": neutral["canonical_reproduction"]["all_pass"] is True,
        "neutral_probe_has_two_orientations": len(neutral["configurations"]) == 2,
        "neutral_T_changes_1_to_3": neutral_native["T"] == 1 and neutral_flip["T"] == 3,
        "neutral_Topo_changes_1_to_2": neutral_native["Topo"] == 1 and neutral_flip["Topo"] == 2,
        "neutral_orientation_sensitive": neutral.get("orientation_sensitive") is True,
        "panel_self_checks_pass": panel.get("all_checks_pass") is True,
        "panel_native_reproduction_12_of_12": panel["counts"]["canonical_reproduction_pass"] == 12,
        "panel_family_unit_sensitive_12_of_12": (
            panel["counts"]["orientation_sensitive"] == 12
            and panel["counts"]["changed_fields"]["family_unit_signature"] == 12
        ),
        "panel_region_T_changes_9": panel["counts"]["changed_fields"]["T"] == 9,
        "panel_region_Topo_changes_8": panel["counts"]["changed_fields"]["Topo"] == 8,
    }
    receipt = {
        "schema_name": "intersubmod.hcc1395_ps_orientation_audit_verification",
        "schema_version": "1.0.0",
        "inputs": {
            "census": {"path": str(args.census.resolve()), "sha256": sha256(args.census)},
            "neutral_probe": {"path": str(args.neutral_probe.resolve()), "sha256": sha256(args.neutral_probe)},
            "panel_summary": {"path": str(args.panel_summary.resolve()), "sha256": sha256(args.panel_summary)},
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
