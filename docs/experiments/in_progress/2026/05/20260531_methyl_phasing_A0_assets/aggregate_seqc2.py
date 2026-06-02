#!/usr/bin/env python3
"""
彙總 seqc2_cn_methyl_chr*.json → 確定性回答 Q1/Q2/Q3。全部統計由本腳本算，落 JSON。
報告須 Read 此 JSON 後才引用數字（feedback_no_fabricated_numbers_in_reports）。

Q1 甲基分群 AUC vs SEQC2 CN/LOH 狀態一致性：
   - 各 status(gain/loss/loh/neutral) 的 anchor_auc 與 auc_minus_p95 分布（median/mean/n）
   - Spearman 相關：seqc2_cn（整數倍體）vs anchor_auc
   - 若 gain/loh 的 auc 顯著 > neutral → 訊號與 CN 狀態相關（confound 嫌疑）
   - 若各 status auc 相近 → 甲基訊號獨立於 CN（非 confound，較好）

Q2 倍體+alignment 造成甲基群差異：
   - gmm_bic_diff（多模態指標）vs status / vs depth：多模態是否集中 gain/高 depth 區
   - depth vs seqc2_cn 一致性（驗證 SEQC2 CN 在本 BAM 成立）

Q3 甲基能否區分幾倍體：
   - depth vs seqc2_cn Spearman（coverage 是否隨 CN 單調 → CN proxy 有效）
   - 各 CN 整數的 anchor_auc / gmm_bic 分布（甲基特徵能否分 CN）
"""
import json, glob, os
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

AD = os.path.dirname(os.path.abspath(__file__))
OUT = AD + "/seqc2_aggregate.json"

def med(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.median(xs)), 4) if xs else None
def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.mean(xs)), 4) if xs else None

def main():
    rows = []
    chroms = []
    for fp in sorted(glob.glob(AD + "/seqc2_cn_methyl_chr*.json")):
        d = json.load(open(fp))
        chroms.append(d["chrom"])
        rows.extend(d["regions"])

    out = {"n_chroms": len(chroms), "chroms": sorted(chroms), "n_regions_total": len(rows)}

    # ---- Q1: status 分層 ----
    by_status = {}
    for s in ("gain", "loss", "loh", "neutral"):
        sub = [r for r in rows if r["seqc2_status"] == s]
        aucs = [r["anchor_auc"] for r in sub]
        deltas = [r["auc_minus_p95"] for r in sub if r["auc_minus_p95"] is not None]
        by_status[s] = {
            "n": len(sub),
            "auc_median": med(aucs), "auc_mean": mean(aucs),
            "delta_median": med(deltas), "delta_mean": mean(deltas),
            "frac_delta_pos": round(float(np.mean([d > 0 for d in deltas])), 3) if deltas else None,
        }
    out["Q1_by_status"] = by_status

    # Q1 顯著性：gain vs neutral, loh vs neutral 的 auc Mann-Whitney
    def mw(a, b):
        a = [x for x in a if x is not None]; b = [x for x in b if x is not None]
        if len(a) < 5 or len(b) < 5:
            return None
        try:
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            return round(float(p), 5)
        except Exception:
            return None
    neu_auc = [r["anchor_auc"] for r in rows if r["seqc2_status"] == "neutral"]
    out["Q1_mannwhitney_auc_vs_neutral"] = {
        "gain_vs_neutral_p": mw([r["anchor_auc"] for r in rows if r["seqc2_status"] == "gain"], neu_auc),
        "loh_vs_neutral_p": mw([r["anchor_auc"] for r in rows if r["seqc2_status"] == "loh"], neu_auc),
        "loss_vs_neutral_p": mw([r["anchor_auc"] for r in rows if r["seqc2_status"] == "loss"], neu_auc),
    }

    # Q1 Spearman: seqc2_cn (整數) vs anchor_auc（只 gain/loss 有 cn）
    cn_pairs = [(r["seqc2_cn"], r["anchor_auc"]) for r in rows
                if isinstance(r["seqc2_cn"], (int, float)) and r["anchor_auc"] is not None]
    if len(cn_pairs) >= 8:
        cns, acs = zip(*cn_pairs)
        rho, p = spearmanr(cns, acs)
        out["Q1_spearman_cn_vs_auc"] = {"rho": round(float(rho), 4), "p": round(float(p), 5), "n": len(cn_pairs)}
    else:
        out["Q1_spearman_cn_vs_auc"] = {"note": "n<8 不足", "n": len(cn_pairs)}

    # ---- Q2: 多模態 vs status / depth ----
    by_status_bic = {}
    for s in ("gain", "loss", "loh", "neutral"):
        bics = [r["gmm_bic_diff"] for r in rows if r["seqc2_status"] == s and r["gmm_bic_diff"] is not None]
        by_status_bic[s] = {
            "n": len(bics), "bic_median": med(bics),
            "frac_bimodal_bic_gt10": round(float(np.mean([b > 10 for b in bics])), 3) if bics else None,
        }
    out["Q2_bimodality_by_status"] = by_status_bic
    # bic vs depth Spearman
    bd = [(r["mean_depth"], r["gmm_bic_diff"]) for r in rows
          if r["mean_depth"] is not None and r["gmm_bic_diff"] is not None]
    if len(bd) >= 8:
        ds, bs = zip(*bd)
        rho, p = spearmanr(ds, bs)
        out["Q2_spearman_depth_vs_bic"] = {"rho": round(float(rho), 4), "p": round(float(p), 5), "n": len(bd)}

    # ---- Q3: depth vs CN（CN proxy 有效性） + 各 CN 的甲基特徵 ----
    dcn = [(r["seqc2_cn"], r["mean_depth"]) for r in rows
           if isinstance(r["seqc2_cn"], (int, float)) and r["mean_depth"] is not None]
    if len(dcn) >= 8:
        cns, ds = zip(*dcn)
        rho, p = spearmanr(cns, ds)
        out["Q3_spearman_cn_vs_depth"] = {"rho": round(float(rho), 4), "p": round(float(p), 5), "n": len(dcn)}
    # 各 CN 整數的 depth + auc
    by_cn = {}
    for r in rows:
        cn = r["seqc2_cn"]
        if isinstance(cn, (int, float)) and float(cn).is_integer():
            cn = int(cn)
            by_cn.setdefault(cn, {"depths": [], "aucs": []})
            if r["mean_depth"] is not None:
                by_cn[cn]["depths"].append(r["mean_depth"])
            if r["anchor_auc"] is not None:
                by_cn[cn]["aucs"].append(r["anchor_auc"])
    out["Q3_by_cn"] = {str(cn): {"n": len(v["depths"]), "depth_median": med(v["depths"]),
                                  "auc_median": med(v["aucs"])}
                       for cn, v in sorted(by_cn.items())}

    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=2)
    # 同時印給 stdout（但報告以 JSON 為準）
    print(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n[aggregate] -> {OUT}")

if __name__ == "__main__":
    main()
