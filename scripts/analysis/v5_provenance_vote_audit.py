#!/usr/bin/env python3
"""T1.2 Vote countMap audit — JOIN baseline/v3f/v5 vote dumps and run 4-path verification.

Inputs:
  - vote_dump_baseline_chr19.tsv.gz  (priority bug 版)
  - vote_dump_v3f_chr19.tsv.gz       (V3F two-layer)
  - vote_dump_v5_chr19.tsv.gz        (HEAD = V5 Layer 1.5 + ploidy fix + threshold 0.9)

Output: read_level_vote_audit summary + figures.

4-path verification (T1.2 plan):
  ① 個案 trace ≥10 條 priority bug 受害 reads
  ② 區域聚集 (chr19 內部 sliding window)
  ③ AF/density 共變代理 (用 cnt_HP1_1 + cnt_HP2_1 高低分組)
  ④ 修正後消失 (baseline 受害 → v3f/v5 hpResult 翻轉比例)
"""
import argparse
import gzip
from pathlib import Path

import pandas as pd


def load_dump(path: Path, ver: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip" if str(path).endswith(".gz") else None)
    df = df.rename(columns={c: f"{c}_{ver}" if c not in ("read_name", "chr", "pos") else c
                            for c in df.columns})
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-dir", required=True,
                    help="dir containing vote_dump_{baseline,v3f,v5}_chr19.tsv.gz")
    ap.add_argument("--out", required=True, help="output summary .md")
    args = ap.parse_args()

    base_dir = Path(args.dump_dir)
    print(f"[1/5] loading 3 dumps from {base_dir}")
    df_b = load_dump(base_dir / "vote_dump_baseline_chr19.tsv.gz", "baseline")
    df_v3f = load_dump(base_dir / "vote_dump_v3f_chr19.tsv.gz", "v3f")
    df_v5 = load_dump(base_dir / "vote_dump_v5_chr19.tsv.gz", "v5")
    print(f"  baseline: {len(df_b)} rows | v3f: {len(df_v3f)} rows | v5: {len(df_v5)} rows")

    print("[2/5] JOIN by (read_name, chr, pos)")
    merged = df_b.merge(df_v3f, on=["read_name", "chr", "pos"], how="inner")
    merged = merged.merge(df_v5, on=["read_name", "chr", "pos"], how="inner")
    print(f"  merged: {len(merged)} rows (3-way intersection)")

    # ── derived columns ──
    # germline majority: 1 if cnt_HP1 > cnt_HP2, 2 if cnt_HP2 > cnt_HP1, 0 if equal
    def majority(a, b):
        return (a > b).astype(int) * 1 + (b > a).astype(int) * 2

    # countMap is identical across versions (only voting logic differs); use baseline cols
    merged["germline_maj"] = majority(merged["cnt_HP1_baseline"], merged["cnt_HP2_baseline"])
    merged["somatic_maj"] = majority(merged["cnt_HP1_1_baseline"], merged["cnt_HP2_1_baseline"])
    merged["has_germline_vote"] = ((merged["cnt_HP1_baseline"] + merged["cnt_HP2_baseline"]) > 0)
    merged["has_somatic_vote"] = (
        (merged["cnt_HP1_1_baseline"] + merged["cnt_HP2_1_baseline"] + merged["cnt_HP3_baseline"]) > 0
    )

    # ── path ① priority bug victims ──
    # Definition: germline majority X, somatic majority Y, X≠Y, both >0
    print("[3/5] path ① priority bug victims")
    victims = merged[
        (merged["has_germline_vote"]) &
        (merged["has_somatic_vote"]) &
        (merged["germline_maj"] > 0) &
        (merged["somatic_maj"] > 0) &
        (merged["germline_maj"] != merged["somatic_maj"])
    ]
    print(f"  raw 雙向矛盾 rows: {len(victims):,}")

    # baseline tagged it the "wrong" direction (matching somatic_maj instead of germline_maj)?
    # baseline hpResult: 1/2 (純 germline) | 11/21 (germline + somatic) | 33 (somatic only) | 0 (untagged)
    def baseline_dir(hp):
        if hp == 1 or hp == 11:
            return 1
        if hp == 2 or hp == 21:
            return 2
        return 0

    def v3f_dir(hp):
        return baseline_dir(hp)  # same encoding

    victims = victims.assign(
        dir_b=victims["hpResult_baseline"].apply(baseline_dir),
        dir_v3f=victims["hpResult_v3f"].apply(v3f_dir),
        dir_v5=victims["hpResult_v5"].apply(v3f_dir),
    )

    # baseline 投到 somatic_maj 的方向（priority bug confirmed victim）
    bug_confirmed = victims[victims["dir_b"] == victims["somatic_maj"]]
    print(f"  baseline 標到 somatic 方向（priority bug confirmed）: {len(bug_confirmed):,}")

    # v3f 把這些 reads 修到 germline_maj 方向？
    if len(bug_confirmed) > 0:
        bug_corrected_v3f = bug_confirmed[bug_confirmed["dir_v3f"] == bug_confirmed["germline_maj"]]
        bug_corrected_v5 = bug_confirmed[bug_confirmed["dir_v5"] == bug_confirmed["germline_maj"]]
        v3f_correction_rate = len(bug_corrected_v3f) / len(bug_confirmed) * 100
        v5_correction_rate = len(bug_corrected_v5) / len(bug_confirmed) * 100
    else:
        v3f_correction_rate = float("nan")
        v5_correction_rate = float("nan")

    # ── path ② regional clustering (sliding window 1Mb on chr19) ──
    print("[4/5] path ② chr19 regional clustering (1Mb windows)")
    bug_confirmed = bug_confirmed.assign(window_1mb=(bug_confirmed["pos"] // 1_000_000).astype(int))
    merged_for_window = merged.assign(window_1mb=(merged["pos"] // 1_000_000).astype(int))
    bug_per_window = bug_confirmed.groupby("window_1mb").size()
    total_per_window = merged_for_window.groupby("window_1mb").size()
    enrichment = (bug_per_window / total_per_window).fillna(0).sort_values(ascending=False)

    # ── path ④ "修正後消失" overall verdict ──
    # path ③ AF/density 共變: 暫用 somatic vote count 作為 density proxy
    print("[5/5] path ③/④ density共變 + 修正後消失")
    high_somatic = bug_confirmed[
        (bug_confirmed["cnt_HP1_1_baseline"] + bug_confirmed["cnt_HP2_1_baseline"]) >= 5
    ]
    low_somatic = bug_confirmed[
        (bug_confirmed["cnt_HP1_1_baseline"] + bug_confirmed["cnt_HP2_1_baseline"]) < 5
    ]

    # ── write summary ──
    summary_md = f"""# T1.2 Read-Level Vote Audit — 4 路徑驗證結果

**Input**: 3 vote dumps from `{args.dump_dir}`
**Region**: chr19 (HCC1395 5kHz, V2b PON-only phased VCF as input)
**Binaries**:
  - baseline = `8b8c1fd` (V2b PON-only flag, getVote priority bug 仍在)
  - v3f      = `380e8d2` (V3F two-layer + INDEL guard)
  - v5       = HEAD `938f0df` (Layer 1.5 + ploidy fix + threshold 0.9)

## Summary

| 指標 | 值 |
|------|---:|
| baseline rows | {len(df_b):,} |
| v3f rows | {len(df_v3f):,} |
| v5 rows | {len(df_v5):,} |
| 3-way merged | {len(merged):,} |
| 雙向矛盾 reads (germline_maj ≠ somatic_maj, both >0) | {len(victims):,} |
| **Priority bug confirmed victims** (baseline 跟 somatic) | **{len(bug_confirmed):,}** |
| V3F 修正比例 (改向 germline_maj) | **{v3f_correction_rate:.2f}%** |
| V5 修正比例 (改向 germline_maj) | **{v5_correction_rate:.2f}%** |

## 4 路徑驗證

### ① 個案 trace（≥10 條）
總數 {len(bug_confirmed):,} 條 — **{'PASS' if len(bug_confirmed) >= 10 else 'FAIL'}**（門檻 ≥10）

前 10 個案例（read_name + pos + countMap + hpResult 三版對比）：

| read_name (前 12 chars) | pos | HP1/HP2 | HP1_1/HP2_1 | germline_maj | somatic_maj | hp_b → hp_v3f → hp_v5 |
|---|---:|---|---|:---:|:---:|---|
"""
    for _, r in bug_confirmed.head(10).iterrows():
        summary_md += (
            f"| {r['read_name'][:12]} | {r['pos']:,} | "
            f"{r['cnt_HP1_baseline']}/{r['cnt_HP2_baseline']} | "
            f"{r['cnt_HP1_1_baseline']}/{r['cnt_HP2_1_baseline']} | "
            f"{int(r['germline_maj'])} | {int(r['somatic_maj'])} | "
            f"{int(r['hpResult_baseline'])} → {int(r['hpResult_v3f'])} → {int(r['hpResult_v5'])} |\n"
        )

    summary_md += f"""

### ② 區域聚集 (chr19 1Mb windows)
Top 10 priority bug 比例最高的 windows：

| window 1Mb | bug confirmed | total reads | enrichment |
|---:|---:|---:|---:|
"""
    for w, frac in enrichment.head(10).items():
        bug_n = int(bug_per_window.get(w, 0))
        tot_n = int(total_per_window.get(w, 0))
        summary_md += f"| chr19:{w}M | {bug_n:,} | {tot_n:,} | {frac*100:.2f}% |\n"

    summary_md += f"""

### ③ Somatic density 共變
分組對比（high somatic vote vs low）：

| 群體 | bug confirmed N | V3F 修正率 |
|---|---:|---:|
| **high** somatic vote ≥5 | {len(high_somatic):,} | {(high_somatic['dir_v3f'] == high_somatic['germline_maj']).mean()*100:.2f}% |
| **low**  somatic vote <5 | {len(low_somatic):,} | {(low_somatic['dir_v3f'] == low_somatic['germline_maj']).mean()*100:.2f}% |

→ **{'PASS' if len(high_somatic) > len(low_somatic) else 'BORDERLINE'}**：High somatic density reads {'是' if len(high_somatic) > len(low_somatic) else '不是'} priority bug 主要受害者

### ④ 修正後消失
- V3F 修正率 = {v3f_correction_rate:.2f}% — **{'PASS (≥80%)' if v3f_correction_rate >= 80 else 'FAIL'}**
- V5 修正率  = {v5_correction_rate:.2f}% — **{'PASS (≥80%)' if v5_correction_rate >= 80 else 'FAIL'}**

## 機制因果結論

如果 4 路徑 ≥3 通過 → priority bug 機制因果**確立**；V5 修對是真實。
如果 ≤2 通過 → 機制詮釋需降級。
"""

    Path(args.out).write_text(summary_md)
    print(f"\n  summary → {args.out}")
    print(f"  bug confirmed: {len(bug_confirmed):,} reads")
    print(f"  V3F correction rate: {v3f_correction_rate:.2f}%")
    print(f"  V5  correction rate: {v5_correction_rate:.2f}%")


if __name__ == "__main__":
    main()
