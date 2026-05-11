#!/usr/bin/env python3
"""Phase C G9 prep — build extended master with ASM features.

The primary ``merged_with_vcf.tsv.gz`` only carries AlleleDelta (== master AF).
Pull the other 24 G9 features from per-sample/per-mode canonical + archive
``significance_summary.csv`` files and left-join on (sample, mode, Chr, Pos).

G9 features attached (25 total, AlleleDelta already in master -> kept for join):
  - Allele trio      : AlleleDelta, AlleleP, AlleleSig
  - Sample ASM       : SampleASM_Delta, SampleASM_P, SampleASM_Sig,
                       SampleASM_NTumor, SampleASM_NNormal  (paired_full only)
  - HP residual      : HP_Residual_Delta, HP_Residual_P, HP_Residual_Sig
  - Tumor HP         : Tumor_HP_Delta, Tumor_HP_Valid,
                       Tumor_HP1, Tumor_HP2, Tumor_HP_Signed_Delta
  - Normal HP        : Normal_HP_Delta, Normal_HP_Valid,
                       Normal_HP1, Normal_HP2, Normal_HP_Signed_Delta
                       (paired_full only)
  - Normal baseline  : NormalBaseline_Mean, NormalBaseline_Coverage
                       (paired_full only)
  - Combined signed  : HP_Signed_Residual, Combined_HP_Signed_Delta

Note: most Phase B/C/D additions (SampleASM_*, NormalBaseline_*, HP_Residual_*,
Normal_HP_*) exist only in paired_full runs — archive TO significance_summary
predates those columns and will join as NaN.

Output:
  data/G9/master_g9.tsv.gz
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BIG7 = Path("/big7_disk/liaoyoyo2001/big7_disk_output")
MASTER_IN = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
OUT_DIR = ROOT / "research/feature_layered_observation/data/G9"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "master_g9.tsv.gz"

G9_COLS = [
    # Allele trio
    "AlleleDelta", "AlleleP", "AlleleSig",
    # Sample ASM (paired only)
    "SampleASM_Delta", "SampleASM_P", "SampleASM_Sig",
    "SampleASM_NTumor", "SampleASM_NNormal",
    # Normal baseline (paired only)
    "NormalBaseline_Mean", "NormalBaseline_Coverage",
    # HP Residual
    "HP_Residual_Delta", "HP_Residual_P", "HP_Residual_Sig",
    # Tumor HP
    "Tumor_HP_Delta", "Tumor_HP_Valid",
    "Tumor_HP1", "Tumor_HP2",
    "Tumor_HP_Signed_Delta",
    # Normal HP (paired only)
    "Normal_HP_Delta", "Normal_HP_Valid",
    "Normal_HP1", "Normal_HP2",
    "Normal_HP_Signed_Delta",
    # Combined
    "HP_Signed_Residual", "Combined_HP_Signed_Delta",
]

# ---- Canonical paired_full latest run per sample ---------------------------
CANONICAL_PAIRED = {
    "HCC1395":        "canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": "canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "HCC1937":        "canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
    "HCC1954":        "canonical/HCC1954/paired_full/20260421_HCC1954_paired_full_full",
    "H2009":          "canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437":          "canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829":        "canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
}

# ---- TO mode sources (paired with G5 naming convention) --------------------
TO_SOURCES = {
    "HCC1395":        "bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod",
    "HCC1395_DORADO": "synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod",
    "HCC1937":        "synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod",
    "HCC1954":        "synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod",
    "H2009":          "synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod",
    "H1437":          "synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod",
    # COLO829 TO archive step05 is empty
}


def read_ss(sample: str, mode: str, run_path: Path) -> pd.DataFrame:
    frames = []
    for tp_label, subdir in [(1, "intersubmod_tp"), (0, "intersubmod_fp")]:
        ss = run_path / subdir / "significance_summary.csv"
        if not ss.exists():
            print(f"[warn] missing: {ss}", flush=True)
            continue
        try:
            d = pd.read_csv(ss, usecols=lambda c: c in (["Chr", "Pos"] + G9_COLS),
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
    # Drop existing AlleleDelta column so join version is authoritative
    if "AlleleDelta" in master.columns:
        master = master.rename(columns={"AlleleDelta": "AlleleDelta_master"})

    frames = []
    for sample, rel in CANONICAL_PAIRED.items():
        p = BIG7 / rel
        if not p.exists():
            p = ROOT / "output" / rel
        if not p.exists():
            print(f"[skip] paired_full {sample}: path missing ({p})", flush=True)
            continue
        df = read_ss(sample, "paired_full", p)
        if not df.empty:
            frames.append(df)
            print(f"[ok]   paired_full {sample}: rows={len(df):,}  cols={df.shape[1]}", flush=True)

    for sample, rel in TO_SOURCES.items():
        p = BIG7 / rel
        if not p.exists():
            print(f"[skip] to_pileup {sample}: path missing ({p})", flush=True)
            continue
        df = read_ss(sample, "to_pileup", p)
        if not df.empty:
            frames.append(df)
            print(f"[ok]   to_pileup   {sample}: rows={len(df):,}  cols={df.shape[1]}", flush=True)

    if not frames:
        print("[err]  no G9 data gathered; aborting.", flush=True)
        sys.exit(1)

    g9 = pd.concat(frames, ignore_index=True)
    g9["Chr"] = g9["Chr"].astype(str)
    master["Chr"] = master["Chr"].astype(str)
    g9["Pos"] = pd.to_numeric(g9["Pos"], errors="coerce").astype("Int64")
    master["Pos"] = pd.to_numeric(master["Pos"], errors="coerce").astype("Int64")

    g9 = g9.drop_duplicates(subset=["sample", "mode", "Chr", "Pos"], keep="first")
    print(f"[g9]   unique joined rows={len(g9):,}", flush=True)
    g9 = g9.drop(columns=["_ss_tp_label"], errors="ignore")

    merged = master.merge(g9, on=["sample", "mode", "Chr", "Pos"], how="left",
                          suffixes=("", "_ss"))

    # Report join coverage per column
    for c in G9_COLS:
        if c in merged.columns:
            n_ok = merged[c].notna().sum()
            print(f"[cov]  {c}: {n_ok:,} / {len(merged):,}  "
                  f"({n_ok/len(merged)*100:.1f}%)", flush=True)

    # Per mode coverage summary for the key feature NormalBaseline_Mean
    if "NormalBaseline_Mean" in merged.columns:
        for mode in ["paired_full", "to_pileup"]:
            sub = merged[merged["mode"] == mode]
            ok = sub["NormalBaseline_Mean"].notna().sum()
            print(f"[mode-cov] NormalBaseline_Mean {mode}: {ok:,}/{len(sub):,}  "
                  f"({(ok/len(sub)*100 if len(sub) else 0):.1f}%)", flush=True)

    print(f"[save] {OUT_TSV}", flush=True)
    merged.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
