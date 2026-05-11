#!/usr/bin/env python3
"""Phase C G5 prep — build extended master with HP Merged (2-bucket) features.

The primary ``merged_with_vcf.tsv.gz`` lacks G5 columns. Pull them from
per-sample/per-mode canonical + archive ``significance_summary.csv`` files
and left-join on (sample, mode, Chr, Pos).

G5 features attached:
  HP_Ratio, Potential_LOH_ss (sanity), HP1FamilyN, HP2FamilyN,
  HPMergedDelta, HPMergedP, HPMergedSig, GlobalP_HPFamily, CramersV_HPFamily.

Output:
  data/G5/master_g5.tsv.gz

Source map:
  paired_full -> output/canonical/<sample>/paired_full/<latest_run>/
                intersubmod_{tp,fp}/significance_summary.csv
  to_pileup   -> HCC1395  -> output/canonical/HCC1395/to_pileup/
                            20260307_.../intersubmod_{tp,fp}/significance_summary.csv
                 others   -> big7_disk_output/synthesis/research_rounds/archive/
                            202603_early_pilots/<pilot_dir>/step05_intersubmod/
                            intersubmod_{tp,fp}/significance_summary.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BIG7 = Path("/big7_disk/liaoyoyo2001/big7_disk_output")
MASTER_IN = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
OUT_DIR = ROOT / "research/feature_layered_observation/data/G5"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "master_g5.tsv.gz"

G5_COLS = [
    "HP_Ratio", "HP1FamilyN", "HP2FamilyN",
    "HPMergedDelta", "HPMergedP", "HPMergedSig",
    "GlobalP_HPFamily", "CramersV_HPFamily",
]

# ---- Canonical paired_full latest run per sample ---------------------------
CANONICAL_PAIRED = {
    "HCC1395":       "canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO":"canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "HCC1937":       "canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
    "HCC1954":       "canonical/HCC1954/paired_full/20260421_HCC1954_paired_full_full",
    "H2009":         "canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437":         "canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829":       "canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
}
# ---- TO mode sources (HCC1395 in bip8 archive, others in synthesis archive) -
TO_SOURCES = {
    "HCC1395":       "bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod",
    "HCC1395_DORADO":"synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod",
    "HCC1937":       "synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod",
    "HCC1954":       "synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod",
    "H2009":         "synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod",
    "H1437":         "synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod",
    # COLO829 archive step05 is empty; will be NaN for G5 in TO mode
}


def read_ss(sample: str, mode: str, run_path: Path) -> pd.DataFrame:
    frames = []
    for tp_label, subdir in [(1, "intersubmod_tp"), (0, "intersubmod_fp")]:
        ss = run_path / subdir / "significance_summary.csv"
        if not ss.exists():
            # archive layout uses just 'intersubmod/{tp,fp}'
            alt = run_path / subdir / "significance_summary.csv"
            if not alt.exists():
                print(f"[warn] missing: {ss}", flush=True)
                continue
        try:
            d = pd.read_csv(ss, usecols=lambda c: c in (["Chr", "Pos"] + G5_COLS),
                            low_memory=False)
        except Exception as e:
            print(f"[err]  {ss}: {e}", flush=True)
            continue
        d["sample"] = sample
        d["mode"] = mode
        d["_ss_tp_label"] = tp_label
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    print(f"[load] master from {MASTER_IN}", flush=True)
    master = pd.read_csv(MASTER_IN, sep="\t", low_memory=False)
    print(f"[load] master rows={len(master):,}, cols={master.shape[1]}", flush=True)

    frames = []
    for sample, rel in CANONICAL_PAIRED.items():
        p = BIG7 / rel
        if not p.exists():
            p = ROOT / "output" / rel.replace("canonical/", "canonical/", 1)
        if not p.exists():
            print(f"[skip] paired_full {sample}: path missing ({p})", flush=True)
            continue
        df = read_ss(sample, "paired_full", p)
        if not df.empty:
            frames.append(df)
            print(f"[ok]   paired_full {sample}: rows={len(df):,}", flush=True)

    for sample, rel in TO_SOURCES.items():
        p = BIG7 / rel
        if not p.exists():
            print(f"[skip] to_pileup {sample}: path missing ({p})", flush=True)
            continue
        df = read_ss(sample, "to_pileup", p)
        if not df.empty:
            frames.append(df)
            print(f"[ok]   to_pileup {sample}: rows={len(df):,}", flush=True)

    if not frames:
        print("[err]  no G5 data gathered; aborting.", flush=True)
        sys.exit(1)

    g5 = pd.concat(frames, ignore_index=True)
    # Normalize Chr to str for join safety
    g5["Chr"] = g5["Chr"].astype(str)
    master["Chr"] = master["Chr"].astype(str)
    # pos must be int
    g5["Pos"] = pd.to_numeric(g5["Pos"], errors="coerce").astype("Int64")
    master["Pos"] = pd.to_numeric(master["Pos"], errors="coerce").astype("Int64")

    # Drop duplicated keys (should not happen; keep first)
    g5 = g5.drop_duplicates(subset=["sample", "mode", "Chr", "Pos"], keep="first")
    print(f"[g5]   unique joined rows={len(g5):,}", flush=True)

    # Drop internal helper
    g5 = g5.drop(columns=["_ss_tp_label"], errors="ignore")

    merged = master.merge(g5, on=["sample", "mode", "Chr", "Pos"], how="left",
                          suffixes=("", "_ss"))
    # Report join coverage
    for c in G5_COLS:
        if c in merged.columns:
            n_ok = merged[c].notna().sum()
            print(f"[cov]  {c}: {n_ok:,} / {len(merged):,}  "
                  f"({n_ok/len(merged)*100:.1f}%)", flush=True)

    # Write
    print(f"[save] {OUT_TSV}", flush=True)
    merged.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
