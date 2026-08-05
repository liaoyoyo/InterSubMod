#!/usr/bin/env python3
"""Quantify reference-only and HP3/HP4 auxiliary denominator effects."""

import argparse
import json
from pathlib import Path


def pct(numerator, denominator):
    return round(100 * numerator / denominator, 3) if denominator else None


def summarize(sample, path):
    doc = json.loads(path.read_text(encoding="utf-8"))
    detail = doc["detail"]
    primary = [u for u in detail if u.get("is_primary_lineage")]
    germline_all = [u for u in detail if u.get("family") in ("1", "2")]
    mutation_plus_h3 = [u for u in detail if u.get("is_primary_lineage")
                        or (u.get("is_h3_auxiliary") and u.get("mutation_bearing"))]
    mutation_plus_aux = [u for u in detail if u.get("is_primary_lineage")
                         or ((u.get("is_h3_auxiliary") or u.get("is_h4_auxiliary"))
                             and u.get("mutation_bearing"))]
    reference = [u for u in detail if u.get("reference_only")]
    h3 = [u for u in detail if u.get("is_h3_auxiliary") and u.get("mutation_bearing")]
    h4 = [u for u in detail if u.get("is_h4_auxiliary") and u.get("mutation_bearing")]
    determined = lambda units: sum(u["L1_class"] == "determined" for u in units)
    by_region = {}
    for unit in detail:
        by_region.setdefault(unit["region"], []).append(unit)
    legacy_multi = 0
    primary_multi = 0
    for units in by_region.values():
        legacy_families = {u["family"] for u in units if u.get("family") in ("1", "2")}
        primary_families = {u["family"] for u in units if u.get("is_primary_lineage")}
        legacy_multi += len(legacy_families) >= 2
        primary_multi += len(primary_families) >= 2
    return {
        "sample": sample,
        "n_regions_with_units": len(by_region),
        "primary_mutation_HP1_HP2": {"n": len(primary), "determined": determined(primary),
                                     "determined_pct": pct(determined(primary), len(primary))},
        "legacy_all_HP1_HP2": {"n": len(germline_all), "determined": determined(germline_all),
                               "determined_pct": pct(determined(germline_all), len(germline_all))},
        "mutation_HP1_HP2_plus_H3_aux": {"n": len(mutation_plus_h3), "determined": determined(mutation_plus_h3),
                                         "determined_pct": pct(determined(mutation_plus_h3), len(mutation_plus_h3))},
        "mutation_HP1_HP2_plus_H3_H4_aux": {
            "n": len(mutation_plus_aux),
            "determined": determined(mutation_plus_aux),
            "determined_pct": pct(determined(mutation_plus_aux), len(mutation_plus_aux)),
        },
        "reference_only_controls": len(reference),
        "H3_auxiliary_mutation_bearing": len(h3),
        "H4_auxiliary_mutation_bearing": len(h4),
        "H4_auxiliary_determined": determined(h4),
        "multiHP_legacy_any_HP1_HP2": legacy_multi,
        "multiHP_primary_mutation_HP1_HP2": primary_multi,
        "multiHP_legacy_pct": pct(legacy_multi, len(by_region)),
        "multiHP_primary_pct": pct(primary_multi, len(by_region)),
        "reference_only_in_primary": sum(u.get("reference_only") and u.get("is_primary_lineage") for u in detail),
        "H3_in_primary": sum(u.get("family") == "3" and u.get("is_primary_lineage") for u in detail),
        "H4_in_primary": sum(u.get("family") == "4" and u.get("is_primary_lineage") for u in detail),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    rows = [summarize(meta["sample"], args.run_root / "samples" / meta["sample"] /
                      f"layered_reconstruction_{meta['sample']}.json") for meta in manifest["samples"]]
    output = {"schema_version": "2.1", "primary_definition": "mutation-bearing HP1/HP2",
              "H3_definition": "unresolved auxiliary",
              "H4_definition": "shared-somatic auxiliary",
              "samples": rows,
              "all_forbidden_counts_zero": all(
                  r["reference_only_in_primary"] == 0
                  and r["H3_in_primary"] == 0
                  and r["H4_in_primary"] == 0
                  for r in rows
              )}
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"LINEAGE ABLATION -> {args.output}")
    for row in rows:
        print(f"  {row['sample']}: primary={row['primary_mutation_HP1_HP2']['n']} "
              f"det={row['primary_mutation_HP1_HP2']['determined_pct']}% "
              f"reference={row['reference_only_controls']} H3?={row['H3_auxiliary_mutation_bearing']} "
              f"H4?={row['H4_auxiliary_mutation_bearing']} "
              f"multiHP={row['multiHP_primary_pct']}%")
    raise SystemExit(0 if output["all_forbidden_counts_zero"] else 1)


if __name__ == "__main__":
    main()
