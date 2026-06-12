#!/usr/bin/env python3
"""取樣門檻觀察分析 — C_min(read-pair 共同 CpG) + min_site_coverage(CpG-column 覆蓋).

全量掃所有 region 的 methylation.csv（read x CpG raw 矩陣, missing=-1.0），
用 mask @ mask.T 高效算 pairwise 共同 CpG 數，累加 histogram（不存所有 pair）。
輸出 baseline tsv + figures，並從 binary run.log 取權威全量 valid-pair-ratio。

用法: python3 analyze_thresholds.py
"""
import sys, os, glob, re, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/big7_disk/liaoyoyo2001/InterSubMod/output/_tmp_sampling_obs_20260612_2201"
EXP = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/experiments/in_progress/2026/06/20260612_2201_sampling_threshold_audit"
BASE = f"{EXP}/baseline"
FIG = f"{EXP}/figures"
os.makedirs(BASE, exist_ok=True); os.makedirs(FIG, exist_ok=True)

CMIN_MAX = 15      # C_min 掃描上限
COV_MAX = 15       # min_site_coverage 掃描上限
PAIR_BINS = 300    # 共同 CpG histogram 上限
COV_BINS = 600     # column 覆蓋 histogram 上限
CLASSES = ["tp", "fp", "fn"]


def analyze_class(cls):
    mcs = glob.glob(f"{OUT}/{cls}/**/methylation/methylation.csv", recursive=True)
    pair_hist = np.zeros(PAIR_BINS, dtype=np.int64)   # 共同 CpG 數分佈
    cov_hist = np.zeros(COV_BINS, dtype=np.int64)     # CpG-column 覆蓋分佈
    n_regions = n_pairs = n_cols = 0
    for mc in mcs:
        try:
            df = pd.read_csv(mc, index_col=0)
        except Exception:
            continue
        if df.shape[0] < 2 or df.shape[1] < 1:
            continue
        M = df.values.astype(np.float64)
        mask = (M != -1.0) & (~np.isnan(M))           # non-missing cell
        mi = mask.astype(np.int32)
        # pairwise 共同 CpG = mask @ mask.T，取上三角（每對 read 一次）
        common = mi @ mi.T
        iu = np.triu_indices(common.shape[0], k=1)
        pc = common[iu]
        pair_hist += np.bincount(np.clip(pc, 0, PAIR_BINS - 1), minlength=PAIR_BINS)
        n_pairs += pc.size
        # CpG-column 覆蓋 = 每欄非 missing read 數
        colcov = mi.sum(axis=0)
        cov_hist += np.bincount(np.clip(colcov, 0, COV_BINS - 1), minlength=COV_BINS)
        n_cols += colcov.size
        n_regions += 1
    return dict(cls=cls, n_regions=n_regions, n_pairs=n_pairs, n_cols=n_cols,
                pair_hist=pair_hist, cov_hist=cov_hist)


def parse_log(cls):
    """從 binary run.log 取權威全量彙總（binary 自己算的 valid pair ratio）。"""
    lp = f"{OUT}/{cls}.run.log"
    out = {}
    if not os.path.exists(lp):
        return out
    txt = open(lp, errors="ignore").read()
    for key, pat in [
        ("total_valid_pairs", r"Total valid read pairs:\s*(\d+)"),
        ("total_invalid_pairs", r"Total invalid pairs.*?:\s*(\d+)"),
        ("valid_pair_ratio", r"Valid pair ratio:\s*([\d.]+)%"),
        ("avg_common_cpg", r"Average common CpG coverage:\s*([\d.]+)"),
        ("regions_with_dist", r"Regions with distance matrices:\s*(\d+)"),
    ]:
        m = re.search(pat, txt)
        if m:
            out[key] = float(m.group(1))
    return out


def cmin_scan(pair_hist):
    """各 C_min 值下「被設 MAX_DIST 的 pair 比例」= P(common < C_min)。"""
    total = pair_hist.sum()
    rows = []
    for c in range(1, CMIN_MAX + 1):
        excluded = pair_hist[:c].sum()  # common in [0, c-1] < c
        rows.append((c, int(excluded), total and 100.0 * excluded / total))
    return total, rows


def cov_scan(cov_hist):
    """各 min_site_coverage 值下「被砍掉的 CpG-column 比例」= P(coverage < k)。"""
    total = cov_hist.sum()
    rows = []
    for k in range(1, COV_MAX + 1):
        excluded = cov_hist[:k].sum()
        rows.append((k, int(excluded), total and 100.0 * excluded / total))
    return total, rows


def main():
    results = {c: analyze_class(c) for c in CLASSES}
    logs = {c: parse_log(c) for c in CLASSES}

    # --- baseline tsv: C_min 掃描 ---
    with open(f"{BASE}/cmin_common_cpg_scan.tsv", "w") as f:
        f.write("class\tC_min\texcluded_pairs\texcluded_pct\ttotal_pairs\tbinary_log_valid_ratio\tbinary_log_avg_common\n")
        for c in CLASSES:
            total, rows = cmin_scan(results[c]["pair_hist"])
            lg = logs[c]
            for (cm, exc, pct) in rows:
                f.write(f"{c}\t{cm}\t{exc}\t{pct:.3f}\t{total}\t"
                        f"{lg.get('valid_pair_ratio','NA')}\t{lg.get('avg_common_cpg','NA')}\n")

    # --- baseline tsv: column 覆蓋掃描 ---
    with open(f"{BASE}/site_coverage_distribution.tsv", "w") as f:
        f.write("class\tmin_site_coverage\texcluded_columns\texcluded_pct\ttotal_columns\n")
        for c in CLASSES:
            total, rows = cov_scan(results[c]["cov_hist"])
            for (k, exc, pct) in rows:
                f.write(f"{c}\t{k}\t{exc}\t{pct:.3f}\t{total}\n")

    # --- 圖 1: 共同 CpG 分佈 + C_min=3 標線 ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, c in zip(axes, CLASSES):
        h = results[c]["pair_hist"]
        x = np.arange(len(h))
        ax.bar(x[:60], h[:60], color="#4C72B0")
        ax.axvline(3, color="red", ls="--", label="C_min=3")
        med = logs[c].get("avg_common_cpg", None)
        if med:
            ax.axvline(med, color="green", ls=":", label=f"avg={med:.1f}")
        ax.set_title(f"{c.upper()}: pairwise common CpG (n_pairs={results[c]['n_pairs']:,})")
        ax.set_xlabel("common CpG count between read-pair"); ax.set_ylabel("pairs"); ax.legend()
    plt.tight_layout(); plt.savefig(f"{FIG}/fig1_common_cpg_distribution.png", dpi=110); plt.close()

    # --- 圖 2: C_min 掃描排除率曲線 ---
    plt.figure(figsize=(7, 5))
    for c, col in zip(CLASSES, ["#4C72B0", "#DD8452", "#55A868"]):
        total, rows = cmin_scan(results[c]["pair_hist"])
        xs = [r[0] for r in rows]; ys = [r[2] for r in rows]
        plt.plot(xs, ys, "o-", color=col, label=f"{c.upper()} (n={total:,})")
    plt.axvline(3, color="red", ls="--", label="current C_min=3")
    plt.xlabel("C_min (min common CpG)"); plt.ylabel("% read-pairs set to MAX_DIST")
    plt.title("C_min scan: how many read-pairs get distance=1.0 penalty")
    plt.legend(); plt.grid(alpha=0.3)
    plt.savefig(f"{FIG}/fig2_cmin_scan.png", dpi=110); plt.close()

    # --- 圖 3: column 覆蓋分佈 + min_site_coverage=5 標線 ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, c in zip(axes, CLASSES):
        h = results[c]["cov_hist"]
        x = np.arange(len(h))
        ax.bar(x[:80], h[:80], color="#8172B3")
        ax.axvline(5, color="red", ls="--", label="min_site_coverage=5")
        ax.set_title(f"{c.upper()}: CpG-column coverage (n_cols={results[c]['n_cols']:,})")
        ax.set_xlabel("reads covering a CpG column"); ax.set_ylabel("columns"); ax.legend()
    plt.tight_layout(); plt.savefig(f"{FIG}/fig3_column_coverage.png", dpi=110); plt.close()

    # --- summary json（供填回 manifest）---
    summ = {}
    for c in CLASSES:
        total_p, rows_c = cmin_scan(results[c]["pair_hist"])
        total_col, rows_k = cov_scan(results[c]["cov_hist"])
        cmin3 = next(r for r in rows_c if r[0] == 3)
        cov5 = next(r for r in rows_k if r[0] == 5)
        summ[c] = {
            "n_regions": results[c]["n_regions"],
            "n_pairs": results[c]["n_pairs"],
            "n_cols": results[c]["n_cols"],
            "cmin3_excluded_pct": round(cmin3[2], 3),
            "min_site_cov5_excluded_col_pct": round(cov5[2], 3),
            "binary_log": logs[c],
        }
    json.dump(summ, open(f"{BASE}/threshold_summary.json", "w"), indent=2)
    print(json.dumps(summ, indent=2))
    print("\nFIGURES:", FIG)
    print("BASELINE:", BASE)


if __name__ == "__main__":
    main()
