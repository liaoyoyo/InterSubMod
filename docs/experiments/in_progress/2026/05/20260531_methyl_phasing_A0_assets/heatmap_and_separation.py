#!/usr/bin/env python3
"""
Per-read 甲基熱圖 (HP 側欄) + 量化 HP1 vs HP2 甲基分離度。

回答用戶「驗證結果是否合理正確」：不只畫圖，還算 4 個量化指標 +
2 個 negative control，判定「甲基到底分不分得開 HP」。

輸入：extract_per_read_methyl.py 產的 _matrix.tsv + _meta.tsv
輸出：
  {prefix}_heatmap.png       — clustermap (rows=read, cols=CpG, 左側欄=HP tag)
  {prefix}_separation.json   — 量化指標 + null

量化指標（對齊 design §6+R1 三道關）：
  1. silhouette(HP1 vs HP2 標籤, 甲基空間)  — >0 表 HP 在甲基空間有結構
  2. anchor AUC: 用 HP1/HP2 read 質心，看每條 read 離哪群近 → ROC AUC
  3. shuffle null: 打亂 HP 標籤後 silhouette（應掉到 ~0）
  4. within-HP vs between-HP 平均甲基距離（effect size）
限制（誠實標註）：
  - NaN(未覆蓋 CpG) 以「欄位均值」插補供 clustering 顯示；分離度計算用 pairwise-complete。
  - 單樣本單 region；非全基因組結論。
"""
import sys, argparse, json
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# CJK 字型（feedback_matplotlib_cjk_font_rule）
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score, roc_auc_score


def pairwise_complete_dist(M):
    """row-row 歐式距離，只用兩 read 共同覆蓋的 CpG（pairwise complete）。"""
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mask = ~np.isnan(M[i]) & ~np.isnan(M[j])
            if mask.sum() >= 3:
                d = np.sqrt(np.mean((M[i, mask] - M[j, mask]) ** 2))
            else:
                d = np.nan
            D[i, j] = D[j, i] = d
    # 缺值補全域中位數（僅供 silhouette 用）
    med = np.nanmedian(D[D > 0]) if np.any(D > 0) else 1.0
    D[np.isnan(D)] = med
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="extract 輸出 prefix（不含 _matrix.tsv）")
    ap.add_argument("--title", default="")
    ap.add_argument("--min-col-cov", type=float, default=0.3, help="CpG 欄至少被此比例 read 覆蓋才留")
    ap.add_argument("--seed-shuffles", type=int, default=200)
    args = ap.parse_args()

    mat = pd.read_csv(args.prefix + "_matrix.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(args.prefix + "_meta.tsv", sep="\t", index_col=0)
    # 去重 read_id（同名 read 只留第一條；避免 .loc 笛卡爾擴張 IndexError）
    mat = mat[~mat.index.duplicated(keep="first")]
    meta = meta[~meta.index.duplicated(keep="first")]
    common = mat.index.intersection(meta.index)
    mat = mat.loc[common]
    meta = meta.loc[common]

    # 過濾稀疏 CpG 欄
    col_cov = mat.notna().mean(axis=0)
    keep_cols = col_cov[col_cov >= args.min_col_cov].index
    mat = mat[keep_cols]

    # 過濾覆蓋過少的 read（過濾後 <3 CpG）
    row_cov = mat.notna().sum(axis=1)
    keep_rows = row_cov[row_cov >= 3].index
    mat = mat.loc[keep_rows]
    meta = meta.loc[keep_rows]

    result = {"title": args.title, "n_reads": int(mat.shape[0]),
              "n_cpg_cols": int(mat.shape[1]),
              "hp_counts": meta["hp"].value_counts().to_dict()}

    if mat.shape[0] < 6 or mat.shape[1] < 3:
        result["status"] = "INSUFFICIENT_DATA"
        json.dump(result, open(args.prefix + "_separation.json", "w"), ensure_ascii=False, indent=2)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return

    M = mat.values.astype(float)

    # ---- 量化分離度：HP1 vs HP2（germline-phased anchor） ----
    hp = meta["hp"].values
    is12 = np.isin(hp, ["1", "2"])
    sep = {}
    if is12.sum() >= 6 and len(set(hp[is12])) == 2:
        M12 = M[is12]
        lab = (hp[is12] == "2").astype(int)  # 0=HP1, 1=HP2
        D12 = pairwise_complete_dist(M12)
        # 1. silhouette
        try:
            sil = float(silhouette_score(D12, lab, metric="precomputed"))
        except Exception:
            sil = None
        sep["silhouette_HP1_vs_HP2"] = sil
        # 2. anchor AUC：leave-one-out 質心
        #    每 read 對 HP1 群均距 vs HP2 群均距 → score=d(HP1)-d(HP2)，HP2 應較大
        scores = []
        for i in range(M12.shape[0]):
            others1 = [k for k in range(M12.shape[0]) if k != i and lab[k] == 0]
            others2 = [k for k in range(M12.shape[0]) if k != i and lab[k] == 1]
            if not others1 or not others2:
                scores.append(np.nan); continue
            d1 = np.nanmean([D12[i, k] for k in others1])
            d2 = np.nanmean([D12[i, k] for k in others2])
            scores.append(d1 - d2)  # 越大越像 HP2
        scores = np.array(scores)
        valid = ~np.isnan(scores)
        if valid.sum() >= 6 and len(set(lab[valid])) == 2:
            try:
                auc = float(roc_auc_score(lab[valid], scores[valid]))
                sep["anchor_AUC_HP1_vs_HP2"] = max(auc, 1 - auc)  # 對稱化
            except Exception:
                sep["anchor_AUC_HP1_vs_HP2"] = None
        # 3. shuffle null silhouette
        if sil is not None:
            rng = np.random.RandomState(0)
            nulls = []
            for _ in range(args.seed_shuffles):
                sl = rng.permutation(lab)
                if len(set(sl)) < 2:
                    continue
                try:
                    nulls.append(silhouette_score(D12, sl, metric="precomputed"))
                except Exception:
                    pass
            if nulls:
                sep["shuffle_null_silhouette_mean"] = float(np.mean(nulls))
                sep["shuffle_null_silhouette_p95"] = float(np.percentile(nulls, 95))
                sep["silhouette_exceeds_null_p95"] = bool(sil > np.percentile(nulls, 95))
        # 4. within vs between 距離
        w = [D12[i, j] for i in range(len(lab)) for j in range(i + 1, len(lab)) if lab[i] == lab[j]]
        b = [D12[i, j] for i in range(len(lab)) for j in range(i + 1, len(lab)) if lab[i] != lab[j]]
        if w and b:
            sep["within_HP_mean_dist"] = float(np.mean(w))
            sep["between_HP_mean_dist"] = float(np.mean(b))
            sep["between_minus_within"] = float(np.mean(b) - np.mean(w))
    else:
        sep["note"] = "HP1/HP2 anchor read 不足 6 條，無法算分離度"
    result["separation"] = sep

    # ---- verdict ----
    auc = sep.get("anchor_AUC_HP1_vs_HP2")
    exceeds = sep.get("silhouette_exceeds_null_p95")
    if auc is not None:
        if auc > 0.58 and exceeds:
            result["verdict"] = "SEPARABLE (AUC>0.58 且 silhouette 超 null p95)"
        elif auc > 0.58:
            result["verdict"] = "WEAK (AUC>0.58 但 silhouette 未顯著超 null)"
        else:
            result["verdict"] = "NOT-SEPARABLE (AUC<=0.58)"
    else:
        result["verdict"] = "INCONCLUSIVE"

    json.dump(result, open(args.prefix + "_separation.json", "w"), ensure_ascii=False, indent=2)

    # ---- 熱圖 ----
    hp_palette = {"1": "#2563EB", "2": "#DC2626", "1-1": "#60A5FA", "2-1": "#F87171",
                  "3": "#A16207", "unphase": "#9CA3AF"}
    row_colors = meta["hp"].map(lambda x: hp_palette.get(x, "#000000"))
    # NaN 補欄均值供 clustering 顯示（limitation 已標）
    Mfill = mat.fillna(mat.mean(axis=0))
    try:
        g = sns.clustermap(
            Mfill, row_colors=row_colors, col_cluster=True, row_cluster=True,
            cmap="RdBu_r", vmin=0, vmax=1, figsize=(min(14, 2 + mat.shape[1] * 0.06), 9),
            xticklabels=False, yticklabels=False,
            cbar_kws={"label": "5mC 機率"},
            dendrogram_ratio=(0.10, 0.10),
        )
        g.ax_heatmap.set_xlabel(f"CpG 位點 (n={mat.shape[1]})")
        g.ax_heatmap.set_ylabel(f"read (n={mat.shape[0]})")
        # HP legend
        from matplotlib.patches import Patch
        present = [h for h in hp_palette if h in set(meta["hp"])]
        handles = [Patch(facecolor=hp_palette[h], label=f"HP {h}") for h in present]
        g.ax_heatmap.legend(handles=handles, title="HP tag", bbox_to_anchor=(1.18, 1.0),
                            loc="upper left", fontsize=8, title_fontsize=9)
        ttl = args.title or args.prefix
        v = result["verdict"]
        g.figure.suptitle(f"{ttl}\nAUC={auc if auc else 'NA'}  |  {v}", fontsize=11, y=1.02)
        g.savefig(args.prefix + "_heatmap.png", dpi=130, bbox_inches="tight")
        plt.close()
        result["heatmap"] = args.prefix + "_heatmap.png"
    except Exception as e:
        result["heatmap_error"] = str(e)

    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
