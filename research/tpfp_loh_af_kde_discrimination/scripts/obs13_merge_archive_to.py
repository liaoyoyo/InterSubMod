#!/usr/bin/env python3
"""obs13 | Merge archive TO-mode pilot data into extended master.

Sources (all from ``202603_early_pilots`` archive, post-HP-fix 2026-03-30):
  - HCC1395_DORADO (20260315_hcc1395_dorado_to_pilot)
  - H1437 (20260318_h1437_to_pilot_fastresume)
  - H2009 (20260318_h2009_to_pilot_fastresume)
  - HCC1937 (20260318_hcc1937_to_pilot_fastresume)
  - HCC1954 (20260318_hcc1954_to_pilot_fastresume)
  - COLO829: SKIPPED (step05 is empty)

Per sample-mode:
  1. Read ``intersubmod_tp/label_first_metrics.tsv`` (tp_label=1)
  2. Read ``intersubmod_fp/label_first_metrics.tsv`` (tp_label=0)
  3. bcftools-query AF from ``step02_benchmark_clairs_to/filtered_snv_{tp,fp}.vcf.gz``
  4. Parse ``RegionKey`` (chr:pos:ref:alt) -> ``Chr`` / ``Pos``
  5. Add ``sample`` / ``mode="to_pileup"`` columns
  6. Fill missing master columns with NaN (LOH_Subtype, CN tier, KDE fields)

Output: ``data/master_extended.tsv.gz`` (original master preserved).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import DATA_DIR, MASTER_TSV

ARCHIVE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots"
)

# (sample_name, pilot_dir_name)
ARCHIVE_PILOTS: list[tuple[str, str]] = [
    ("HCC1395_DORADO", "20260315_hcc1395_dorado_to_pilot"),
    ("H1437", "20260318_h1437_to_pilot_fastresume"),
    ("H2009", "20260318_h2009_to_pilot_fastresume"),
    ("HCC1937", "20260318_hcc1937_to_pilot_fastresume"),
    ("HCC1954", "20260318_hcc1954_to_pilot_fastresume"),
]


def parse_region_key(rk: str) -> tuple[str, int, str, str]:
    parts = rk.split(":")
    if len(parts) != 4:
        return ("", -1, "", "")
    chrom, pos, ref, alt = parts
    try:
        return (chrom, int(pos), ref, alt)
    except ValueError:
        return (chrom, -1, ref, alt)


def bcftools_af(vcf_path: Path) -> pd.DataFrame:
    """Return DataFrame with columns [Chr, Pos, Ref, Alt, AF] parsed from VCF."""
    cmd = ["bcftools", "query", "-f",
           "%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n", str(vcf_path)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"[WARN] bcftools failed on {vcf_path}: {exc.stderr[:200]}", file=sys.stderr)
        return pd.DataFrame(columns=["Chr", "Pos", "Ref", "Alt", "AF"])
    rows = []
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 5:
            continue
        try:
            rows.append((parts[0], int(parts[1]), parts[2], parts[3], float(parts[4])))
        except ValueError:
            continue
    return pd.DataFrame(rows, columns=["Chr", "Pos", "Ref", "Alt", "AF"])


def classify_af(af: float) -> str:
    if pd.isna(af):
        return pd.NA
    if af < 0.1:
        return "Extreme"
    if 0.4 <= af <= 0.6:
        return "Near-half"
    return "Intermediate"


def bin_af(af: float) -> str:
    if pd.isna(af):
        return pd.NA
    edges = [-0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for i in range(len(edges) - 1):
        if edges[i] < af <= edges[i + 1]:
            return f"({edges[i]}, {edges[i+1]}]"
    return pd.NA


def load_sample_archive(sample: str, pilot_dir: str) -> pd.DataFrame:
    root = ARCHIVE_ROOT / pilot_dir
    step05 = root / "step05_intersubmod"
    step02 = root / "step02_benchmark_clairs_to"

    frames = []
    for tag, tp_label in [("tp", 1), ("fp", 0)]:
        lfm = step05 / f"intersubmod_{tag}" / "label_first_metrics.tsv"
        vcf = step02 / f"filtered_snv_{tag}.vcf.gz"
        if not lfm.exists():
            print(f"[WARN] missing {lfm}", file=sys.stderr)
            continue
        df = pd.read_csv(lfm, sep="\t", low_memory=False)
        df["tp_label"] = tp_label

        # parse RegionKey -> Chr, Pos, Ref, Alt
        rk_parts = df["RegionKey"].astype(str).apply(parse_region_key)
        df["Chr"] = [p[0] for p in rk_parts]
        df["Pos"] = [p[1] for p in rk_parts]
        df["Ref"] = [p[2] for p in rk_parts]
        df["Alt"] = [p[3] for p in rk_parts]

        # AF from VCF
        if vcf.exists():
            af_df = bcftools_af(vcf)
            df = df.merge(af_df, on=["Chr", "Pos", "Ref", "Alt"], how="left")
        else:
            print(f"[WARN] missing {vcf}", file=sys.stderr)
            df["AF"] = np.nan

        frames.append(df)

    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    combined["sample"] = sample
    combined["mode"] = "to_pileup"
    combined["AF_class"] = combined["AF"].apply(classify_af)
    combined["AF_bin10"] = combined["AF"].apply(bin_af)

    af_ok = combined["AF"].notna().sum()
    total = len(combined)
    print(f"[obs13] {sample}: loaded {total} rows (tp={(combined['tp_label']==1).sum()}, "
          f"fp={(combined['tp_label']==0).sum()}), AF matched {af_ok}/{total} "
          f"({af_ok/max(total,1)*100:.1f}%)")
    return combined


def align_to_master_columns(df: pd.DataFrame, master_cols: list[str]) -> pd.DataFrame:
    """Ensure df has every master column; fill missing with NaN.

    Non-derivable archive columns (filled with NaN):
      NumCpGs, Coverage_Multiple, Coverage_Category, HPFine_NGroups_CF,
      Diploid_Coverage_Used, Potential_LOH, NTumorReads, NNormalReads,
      LOH_Bed_Overlap, LOH_Bed_Annotation, LOH_Subtype, kde_status,
      CovM_used, baseline_x, sample_order, loh_side
    """
    for col in master_cols:
        if col not in df.columns:
            df[col] = np.nan
    return df[master_cols]


def main() -> None:
    master = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip", low_memory=False)
    master_cols = list(master.columns)
    print(f"[obs13] master loaded: {len(master)} rows, {len(master_cols)} cols")

    archive_frames = []
    for sample, pilot_dir in ARCHIVE_PILOTS:
        df = load_sample_archive(sample, pilot_dir)
        if len(df) == 0:
            continue
        df = align_to_master_columns(df, master_cols)
        archive_frames.append(df)

    archive_all = pd.concat(archive_frames, ignore_index=True)
    print(f"[obs13] archive total: {len(archive_all)} rows")

    extended = pd.concat([master, archive_all], ignore_index=True)
    out_path = DATA_DIR / "master_extended.tsv.gz"
    extended.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"[obs13] wrote {out_path}  total rows: {len(extended)}")

    # Summary table
    summary = extended.groupby(["sample", "mode", "tp_label"]).size().reset_index(name="n")
    summary_path = DATA_DIR / "master_extended_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"[obs13] wrote {summary_path}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
