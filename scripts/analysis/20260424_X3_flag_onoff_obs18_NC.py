#!/usr/bin/env python3
"""
X3 — HCC1395 TO `--germline-hp-only` flag negative control for Thread D.

Compares NG=2 Inner same-hap vs Outer cross-het TP rate gap between:
  - flag=off (normal somatic HP tagging)
  - flag=on (demote somatic HP:i:11/21/33 back to germline HP:i:1/2)

Expected (if Thread D mechanism correct):
  flag=off: gap ~+0.46 (same as obs18 HCC1395)
  flag=on:  gap ~0 (same-hap bucket disappears; cross-het becomes only NG=2 composition)

Data source: /tmp/ism_hp_fix_phase1/{tp,fp}_{off,on}/significance_summary.csv
Reuses obs18 logic; outputs standalone JSON + stdout table.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

BUCKET_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]
NEEDED = ["RegionID", "AlleleDelta", "HPFineNGroups", "Potential_LOH"] + BUCKET_COLS

PATHS = {
    "flag_off": {
        "tp": "/tmp/ism_hp_fix_phase1/tp_off/significance_summary.csv",
        "fp": "/tmp/ism_hp_fix_phase1/fp_off/significance_summary.csv",
    },
    "flag_on": {
        "tp": "/tmp/ism_hp_fix_phase1/tp_on/significance_summary.csv",
        "fp": "/tmp/ism_hp_fix_phase1/fp_on/significance_summary.csv",
    },
}
OUT_JSON = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data/X3_flag_onoff_NC.json")


def categorize_combo(has_hp1, has_hp1s, has_hp2, has_hp2s) -> str:
    if has_hp1 and has_hp1s and not has_hp2 and not has_hp2s:
        return "same_HP1"
    if has_hp2 and has_hp2s and not has_hp1 and not has_hp1s:
        return "same_HP2"
    if has_hp1 and has_hp2s and not has_hp1s and not has_hp2:
        return "cross_het"
    if has_hp1s and has_hp2 and not has_hp1 and not has_hp2s:
        return "cross_het_inv"
    return "other"


def load_combined(tp_path, fp_path):
    dfs = []
    for path, lbl in [(tp_path, 1), (fp_path, 0)]:
        df = pd.read_csv(path, low_memory=False, usecols=NEEDED)
        df["tp_label"] = lbl
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def analyze(df, condition):
    total = len(df)
    ng2 = df[df["HPFineNGroups"] == 2].copy()
    for c in BUCKET_COLS:
        ng2[c] = ng2[c].fillna(0).astype(int)
    ng2["combo"] = ng2.apply(
        lambda r: categorize_combo(
            int(r["HPFineN_HP1"] > 0), int(r["HPFineN_HP1S"] > 0),
            int(r["HPFineN_HP2"] > 0), int(r["HPFineN_HP2S"] > 0),
        ), axis=1,
    )
    ng2["loh_side"] = ng2["Potential_LOH"].apply(
        lambda x: "Inner" if str(x).lower() in ("true", "1") else "Outer"
    )
    ng2["AF_class"] = ng2["AlleleDelta"].clip(0, 1).apply(
        lambda af: "Extreme" if (af < 0.1 or af > 0.9) else
        ("Near-half" if 0.4 <= af <= 0.6 else "Intermediate")
    )

    # Full NG distribution (for flag effect check)
    ng_distribution = df["HPFineNGroups"].fillna(0).astype(int).value_counts().to_dict()
    ng_distribution = {int(k): int(v) for k, v in ng_distribution.items()}

    # NG=2 Extreme subset
    ng2_ext = ng2[ng2["AF_class"] == "Extreme"].copy()
    grp = (ng2_ext.groupby(["loh_side", "combo"])
           .agg(n=("tp_label", "size"), n_tp=("tp_label", "sum"))
           .reset_index())
    grp["tp_rate"] = grp["n_tp"] / grp["n"]

    # Extract canonical obs18 cells
    def cell(side, combo):
        s = grp[(grp["loh_side"] == side) & (grp["combo"] == combo)]
        if s.empty:
            return {"tp_rate": None, "n": 0, "n_tp": 0}
        r = s.iloc[0]
        return {"tp_rate": float(r["tp_rate"]), "n": int(r["n"]), "n_tp": int(r["n_tp"])}

    inner_same = cell("Inner", "same_HP1")
    outer_cross = cell("Outer", "cross_het")
    gap = None
    if inner_same["tp_rate"] is not None and outer_cross["tp_rate"] is not None:
        gap = inner_same["tp_rate"] - outer_cross["tp_rate"]

    return {
        "condition": condition,
        "total_rows": total,
        "ng_distribution_top6": dict(sorted(ng_distribution.items())[:6]),
        "ng2_total": len(ng2),
        "ng2_extreme_n": len(ng2_ext),
        "full_grid": grp.to_dict(orient="records"),
        "canonical": {
            "inner_same_HP1": inner_same,
            "outer_cross_het": outer_cross,
            "gap": gap,
        },
    }


def main():
    results = {}
    for cond, paths in PATHS.items():
        print(f"=== {cond} ===")
        df = load_combined(paths["tp"], paths["fp"])
        r = analyze(df, cond)
        results[cond] = r
        can = r["canonical"]
        print(f"  total={r['total_rows']:,}  NG=2 all={r['ng2_total']:,}  NG=2 Ext={r['ng2_extreme_n']:,}")
        print(f"  NG distribution: {r['ng_distribution_top6']}")
        print(f"  Inner same_HP1: tp_rate={can['inner_same_HP1']['tp_rate']}  n={can['inner_same_HP1']['n']}")
        print(f"  Outer cross_het: tp_rate={can['outer_cross_het']['tp_rate']}  n={can['outer_cross_het']['n']}")
        print(f"  gap = {can['gap']}")
        print()

    # Compare
    gap_off = results["flag_off"]["canonical"]["gap"]
    gap_on = results["flag_on"]["canonical"]["gap"]
    print("=== NEGATIVE CONTROL SUMMARY ===")
    print(f"  gap flag=off = {gap_off}")
    print(f"  gap flag=on  = {gap_on}")
    if gap_off is not None and gap_on is not None:
        collapse = gap_off - gap_on
        print(f"  gap collapse (off - on) = {collapse:.4f}")
        if gap_on is not None and gap_on < 0.1 and gap_off > 0.3:
            print("  VERDICT: Thread D mechanism VERIFIED (gap collapsed under flag=on)")
        elif gap_on is not None and gap_on > 0.3:
            print("  VERDICT: Thread D mechanism NOT clearly verified (gap persists under flag=on)")
        else:
            print("  VERDICT: edge case, manual inspection needed")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump({
            "analysis": "X3_flag_onoff_obs18_negative_control",
            "date": "2026-04-24",
            "sample": "HCC1395",
            "mode": "TO",
            "data": results,
            "summary": {
                "gap_flag_off": gap_off,
                "gap_flag_on": gap_on,
                "gap_collapse_off_minus_on": (gap_off - gap_on) if (gap_off is not None and gap_on is not None) else None,
            },
        }, f, indent=2)
    print(f"\nJSON → {OUT_JSON}")


if __name__ == "__main__":
    main()
