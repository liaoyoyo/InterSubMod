#!/usr/bin/env python3
"""T1.2 Vote countMap audit — JOIN baseline/v3f/v5 vote dumps and run 4-path verification.

Inputs (suffix-driven):
  - vote_dump_baseline_{suffix}.tsv.gz  (priority bug 版)
  - vote_dump_v3f_{suffix}.tsv.gz       (V3F two-layer)
  - vote_dump_v5_{suffix}.tsv.gz        (HEAD = V5 Layer 1.5 + ploidy fix + threshold 0.9)

  suffix ∈ {chr19, genome}（chr19 = T1.2 pilot；genome = T1.2-F1 全基因組擴展）

Output: read_level_vote_audit summary + per-chromosome stats（genome 模式時）+ Layer 1.5 trigger detection.

4-path verification (T1.2 plan):
  ① 個案 trace ≥10 條 priority bug 受害 reads
  ② 區域聚集 (sliding window 1Mb；genome 模式擴成 per-chr enrichment)
  ③ AF/density 共變代理 (用 cnt_HP1_1 + cnt_HP2_1 高低分組)
  ④ 修正後消失 (baseline 受害 → v3f/v5 hpResult 翻轉比例)

Genome-only extras:
  ⑤ Per-chromosome priority bug enrichment + chr8 hotspot 驗證
  ⑥ Layer 1.5 trigger detection (germline_vote=0 cases V5 hpResult vs V3F)
"""
import argparse
from pathlib import Path

import pandas as pd


def load_dump(path: Path, ver: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip" if str(path).endswith(".gz") else None)
    df = df.rename(columns={c: f"{c}_{ver}" if c not in ("read_name", "chr", "pos") else c
                            for c in df.columns})
    return df


def baseline_dir(hp):
    if hp == 1 or hp == 11:
        return 1
    if hp == 2 or hp == 21:
        return 2
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-dir", required=True,
                    help="dir containing vote_dump_{baseline,v3f,v5}_{suffix}.tsv.gz")
    ap.add_argument("--suffix", default="chr19",
                    help="dump file suffix (chr19 / genome). default chr19")
    ap.add_argument("--out", required=True, help="output summary .md")
    args = ap.parse_args()

    base_dir = Path(args.dump_dir)
    suffix = args.suffix
    is_genome = (suffix == "genome")

    print(f"[1/6] loading 3 dumps from {base_dir} (suffix={suffix})")
    df_b = load_dump(base_dir / f"vote_dump_baseline_{suffix}.tsv.gz", "baseline")
    df_v3f = load_dump(base_dir / f"vote_dump_v3f_{suffix}.tsv.gz", "v3f")
    df_v5 = load_dump(base_dir / f"vote_dump_v5_{suffix}.tsv.gz", "v5")
    print(f"  baseline: {len(df_b):,} rows | v3f: {len(df_v3f):,} rows | v5: {len(df_v5):,} rows")

    print("[2/6] JOIN by (read_name, chr, pos)")
    merged = df_b.merge(df_v3f, on=["read_name", "chr", "pos"], how="inner")
    merged = merged.merge(df_v5, on=["read_name", "chr", "pos"], how="inner")
    print(f"  merged: {len(merged):,} rows (3-way intersection)")

    # ── derived columns ──
    def majority(a, b):
        return (a > b).astype(int) * 1 + (b > a).astype(int) * 2

    merged["germline_maj"] = majority(merged["cnt_HP1_baseline"], merged["cnt_HP2_baseline"])
    merged["somatic_maj"] = majority(merged["cnt_HP1_1_baseline"], merged["cnt_HP2_1_baseline"])
    merged["has_germline_vote"] = ((merged["cnt_HP1_baseline"] + merged["cnt_HP2_baseline"]) > 0)
    merged["has_somatic_vote"] = (
        (merged["cnt_HP1_1_baseline"] + merged["cnt_HP2_1_baseline"] + merged["cnt_HP3_baseline"]) > 0
    )

    # ── path ① priority bug victims ──
    print("[3/6] path ① priority bug victims")
    victims = merged[
        (merged["has_germline_vote"]) &
        (merged["has_somatic_vote"]) &
        (merged["germline_maj"] > 0) &
        (merged["somatic_maj"] > 0) &
        (merged["germline_maj"] != merged["somatic_maj"])
    ]
    print(f"  raw 雙向矛盾 rows: {len(victims):,}")

    victims = victims.assign(
        dir_b=victims["hpResult_baseline"].apply(baseline_dir),
        dir_v3f=victims["hpResult_v3f"].apply(baseline_dir),
        dir_v5=victims["hpResult_v5"].apply(baseline_dir),
    )

    bug_confirmed = victims[victims["dir_b"] == victims["somatic_maj"]]
    print(f"  baseline 標到 somatic 方向（priority bug confirmed）: {len(bug_confirmed):,}")

    if len(bug_confirmed) > 0:
        bug_corrected_v3f = bug_confirmed[bug_confirmed["dir_v3f"] == bug_confirmed["germline_maj"]]
        bug_corrected_v5 = bug_confirmed[bug_confirmed["dir_v5"] == bug_confirmed["germline_maj"]]
        v3f_correction_rate = len(bug_corrected_v3f) / len(bug_confirmed) * 100
        v5_correction_rate = len(bug_corrected_v5) / len(bug_confirmed) * 100
    else:
        v3f_correction_rate = float("nan")
        v5_correction_rate = float("nan")

    # ── path ② regional clustering ──
    print("[4/6] path ② regional clustering")
    bug_confirmed = bug_confirmed.assign(window_1mb=(bug_confirmed["pos"] // 1_000_000).astype(int))
    merged_for_window = merged.assign(window_1mb=(merged["pos"] // 1_000_000).astype(int))

    if is_genome:
        # genome 模式：先做 per-chromosome 統計
        bug_per_chr = bug_confirmed.groupby("chr").size().sort_values(ascending=False)
        total_per_chr = merged.groupby("chr").size()
        chr_enrichment = (bug_per_chr / total_per_chr).fillna(0).sort_values(ascending=False)
        print(f"  per-chr top: {bug_per_chr.head(5).to_dict()}")

        # 區域聚集 1Mb window 用 (chr, window) 複合 key
        bug_confirmed = bug_confirmed.assign(chr_window=bug_confirmed["chr"].astype(str) + ":" + bug_confirmed["window_1mb"].astype(str))
        merged_for_window = merged_for_window.assign(chr_window=merged_for_window["chr"].astype(str) + ":" + merged_for_window["window_1mb"].astype(str))
        bug_per_window = bug_confirmed.groupby("chr_window").size()
        total_per_window = merged_for_window.groupby("chr_window").size()
    else:
        bug_per_window = bug_confirmed.groupby("window_1mb").size()
        total_per_window = merged_for_window.groupby("window_1mb").size()
        bug_per_chr = None
        chr_enrichment = None

    enrichment = (bug_per_window / total_per_window).fillna(0).sort_values(ascending=False)
    # 過濾極端噪音（總 reads <100 的 window 不可信）
    enrichment_filtered = enrichment[total_per_window >= 100].sort_values(ascending=False)

    # ── path ③ density 共變 + ④ 修正後消失 ──
    print("[5/6] path ③ density 共變 + ④ 修正後消失")
    high_somatic = bug_confirmed[
        (bug_confirmed["cnt_HP1_1_baseline"] + bug_confirmed["cnt_HP2_1_baseline"]) >= 5
    ]
    low_somatic = bug_confirmed[
        (bug_confirmed["cnt_HP1_1_baseline"] + bug_confirmed["cnt_HP2_1_baseline"]) < 5
    ]

    # ── ⑥ Layer 1.5 trigger detection (genome 模式 only) ──
    layer15_summary = ""
    if is_genome:
        print("[6/6] path ⑥ Layer 1.5 trigger detection")
        # Layer 1.5 觸發條件：germline_vote=0（無 germline 投票），看 V5 是否標到、V3F 是否標到、差異
        no_germline = merged[~merged["has_germline_vote"]]
        n_no_germline = len(no_germline)
        if n_no_germline > 0:
            no_germline = no_germline.assign(
                dir_b=no_germline["hpResult_baseline"].apply(baseline_dir),
                dir_v3f=no_germline["hpResult_v3f"].apply(baseline_dir),
                dir_v5=no_germline["hpResult_v5"].apply(baseline_dir),
            )
            v3f_tagged = (no_germline["dir_v3f"] != 0).sum()
            v5_tagged = (no_germline["dir_v5"] != 0).sum()
            layer15_diff = v5_tagged - v3f_tagged
            v5_only_tagged = ((no_germline["dir_v5"] != 0) & (no_germline["dir_v3f"] == 0)).sum()
            layer15_summary = (
                f"\n## ⑥ Layer 1.5 觸發偵測（genome only）\n\n"
                f"| 指標 | 值 |\n|---|---:|\n"
                f"| germline_vote=0 reads | {n_no_germline:,} |\n"
                f"| V3F tagged 數 | {v3f_tagged:,} |\n"
                f"| V5 tagged 數 | {v5_tagged:,} |\n"
                f"| **Layer 1.5 額外觸發**（V5 tag 而 V3F 沒 tag） | **{v5_only_tagged:,}** |\n"
                f"| V5 - V3F 差異 | {layer15_diff:+,} |\n\n"
                f"→ {'**Layer 1.5 確實觸發**' if v5_only_tagged > 0 else '**Layer 1.5 在此資料未觸發** — V5 vs V3F 在 germline 缺席區域行為相同'}\n"
            )

    # ── per-chr enrichment table (genome only) ──
    per_chr_md = ""
    if is_genome and bug_per_chr is not None:
        # natural chr ordering: chr1..chr22, chrX, chrY
        chr_order_key = lambda c: (
            int(c.replace("chr", "")) if c.replace("chr", "").isdigit()
            else (100 if c == "chrX" else 101 if c == "chrY" else 999)
        )
        chrs_sorted = sorted(bug_per_chr.index, key=chr_order_key)
        per_chr_md = "\n## ⑤ Per-chromosome priority bug 分布\n\n"
        per_chr_md += "| chr | bug victims | total reads | enrichment ‰ | rank |\n"
        per_chr_md += "|---|---:|---:|---:|---:|\n"
        rank_map = {c: i+1 for i, c in enumerate(bug_per_chr.sort_values(ascending=False).index)}
        for c in chrs_sorted:
            n_bug = int(bug_per_chr.get(c, 0))
            n_total = int(total_per_chr.get(c, 0))
            enrich_permille = (n_bug / n_total * 1000) if n_total > 0 else 0.0
            per_chr_md += f"| {c} | {n_bug:,} | {n_total:,} | {enrich_permille:.3f} | {rank_map.get(c, '-')} |\n"

        # chr8 hotspot 驗證
        chr8_bug = int(bug_per_chr.get("chr8", 0))
        chr8_total = int(total_per_chr.get("chr8", 0))
        chr8_rate = chr8_bug / chr8_total if chr8_total > 0 else 0.0
        all_rate = len(bug_confirmed) / len(merged) if len(merged) > 0 else 0.0
        chr8_fold = chr8_rate / all_rate if all_rate > 0 else 0.0
        per_chr_md += (
            f"\n**chr8 hotspot 驗證**：chr8 enrichment = {chr8_rate*1000:.3f}‰ vs genome-wide {all_rate*1000:.3f}‰ "
            f"→ {chr8_fold:.2f}× {'(高於 genome avg)' if chr8_fold > 1 else '(低於或等於 genome avg)'}\n"
        )

    # ── write summary ──
    region_label = "全基因組 (chr1-22 + chrX/Y)" if is_genome else "chr19"
    summary_md = f"""# T1.2 Read-Level Vote Audit — 4 路徑驗證結果（{region_label}）

**Input**: 3 vote dumps from `{args.dump_dir}` (suffix=`{suffix}`)
**Region**: {region_label} (HCC1395 5kHz, V2b PON-only phased VCF as input)
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

| read_name (前 12 chars) | chr | pos | HP1/HP2 | HP1_1/HP2_1 | germline_maj | somatic_maj | hp_b → hp_v3f → hp_v5 |
|---|---|---:|---|---|:---:|:---:|---|
"""
    for _, r in bug_confirmed.head(10).iterrows():
        summary_md += (
            f"| {r['read_name'][:12]} | {r['chr']} | {r['pos']:,} | "
            f"{r['cnt_HP1_baseline']}/{r['cnt_HP2_baseline']} | "
            f"{r['cnt_HP1_1_baseline']}/{r['cnt_HP2_1_baseline']} | "
            f"{int(r['germline_maj'])} | {int(r['somatic_maj'])} | "
            f"{int(r['hpResult_baseline'])} → {int(r['hpResult_v3f'])} → {int(r['hpResult_v5'])} |\n"
        )

    summary_md += f"""

### ② 區域聚集 (1Mb windows{'，total ≥100 過濾' if is_genome else ''})
Top 10 priority bug 比例最高的 windows：

| window | bug confirmed | total reads | enrichment |
|---|---:|---:|---:|
"""
    target_enrichment = enrichment_filtered if is_genome else enrichment
    for w, frac in target_enrichment.head(10).items():
        bug_n = int(bug_per_window.get(w, 0))
        tot_n = int(total_per_window.get(w, 0))
        label = w if is_genome else f"chr19:{w}M"
        summary_md += f"| {label} | {bug_n:,} | {tot_n:,} | {frac*100:.2f}% |\n"

    summary_md += f"""

### ③ Somatic density 共變
分組對比（high somatic vote vs low）：

| 群體 | bug confirmed N | V3F 修正率 |
|---|---:|---:|
| **high** somatic vote ≥5 | {len(high_somatic):,} | {(high_somatic['dir_v3f'] == high_somatic['germline_maj']).mean()*100 if len(high_somatic) else float('nan'):.2f}% |
| **low**  somatic vote <5 | {len(low_somatic):,} | {(low_somatic['dir_v3f'] == low_somatic['germline_maj']).mean()*100 if len(low_somatic) else float('nan'):.2f}% |

→ **{'PASS' if len(high_somatic) > len(low_somatic) else 'BORDERLINE'}**：High somatic density reads {'是' if len(high_somatic) > len(low_somatic) else '不是'} priority bug 主要受害者

### ④ 修正後消失
- V3F 修正率 = {v3f_correction_rate:.2f}% — **{'PASS (≥80%)' if v3f_correction_rate >= 80 else 'FAIL'}**
- V5 修正率  = {v5_correction_rate:.2f}% — **{'PASS (≥80%)' if v5_correction_rate >= 80 else 'FAIL'}**
{per_chr_md}{layer15_summary}
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
