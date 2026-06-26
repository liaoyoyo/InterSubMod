"""[L2 CCF tier — 單變量 VAF→CCF] 三件事:
1. 【最強·CN-independent】祖先≥後代 VAF 梯度 (clonal hierarchy pigeonhole): 每條 nested 邊比 ancestor vs descendant
   VAF — 同區域同局部 CN 故 VAF 可直接比, 不需 CCF 轉換。祖先在更多細胞 → VAF 應 >= 後代。違反率 = data-support。
2. CCF tier 估計 (CN-corrected, 誠實標假設): neutral→min(1,2*VAF); loh→VAF。tier: clonal/major/minor/rare。
3. 跨區 CCF 峰一致性 (CN-clean): 是否離散峰 = 離散 subclone。
輸出 報告夾/data/sm_ccf_tiers.json。
"""
import json
import numpy as np
from collections import Counter, defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
RPT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report"
np.random.seed(11)


def ccf_est(vaf, cn):
    if vaf is None:
        return None
    if cn == "neutral":
        return min(1.0, 2 * vaf)
    if cn == "loh":
        return min(1.0, vaf)
    return None  # gain/loss: multiplicity ambiguous, 不估


def tier(ccf):
    if ccf is None:
        return "undet"
    if ccf >= 0.75:
        return "clonal"
    if ccf >= 0.35:
        return "major_subclone"
    if ccf >= 0.12:
        return "minor_subclone"
    return "rare"


def main():
    reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]

    # 1. 祖先>=後代 VAF 梯度 (nested edges; vaf is per-node in region rec)
    grad_ok = grad_viol = grad_tie = 0
    grad_ok_clean = grad_viol_clean = 0
    deltas = []
    for r in reg:
        vaf = r.get("vaf", {})
        clean = r["cn"] in ("loh", "neutral")
        for anc, des in r["nested_edges"]:
            if anc in vaf and des in vaf:
                d = vaf[anc] - vaf[des]
                deltas.append(d)
                if d > 0.03:
                    grad_ok += 1
                    if clean:
                        grad_ok_clean += 1
                elif d < -0.03:
                    grad_viol += 1
                    if clean:
                        grad_viol_clean += 1
                else:
                    grad_tie += 1
    n_edge = grad_ok + grad_viol + grad_tie
    # null: shuffle ancestor/descendant direction → 50% expected
    grad_support = grad_ok / (grad_ok + grad_viol) if (grad_ok + grad_viol) else None

    # 2. CCF tier (per-node, CN-clean) + 3. peaks
    tier_c = Counter()
    ccf_clean = []
    for r in reg:
        if r["cn"] not in ("loh", "neutral"):
            continue
        for node, v in r.get("vaf", {}).items():
            ce = ccf_est(v, r["cn"])
            if ce is not None:
                tier_c[tier(ce)] += 1
                ccf_clean.append(ce)
    hist, edges = np.histogram(ccf_clean, bins=20, range=(0, 1))
    peaks = sorted([(round(edges[i], 2), int(hist[i])) for i in range(20)], key=lambda x: -x[1])[:6]

    out = {
        "ancestor_ge_descendant_VAF_gradient": {
            "method": "CN-independent: ancestor vs descendant 同區域同局部 CN, VAF 直接比; null=50%",
            "n_nested_edges": n_edge,
            "ancestor_higher": grad_ok, "descendant_higher(violation)": grad_viol, "tie": grad_tie,
            "data_support_rate": round(grad_support, 3) if grad_support else None,
            "clean_only(loh+neutral)": {"ancestor_higher": grad_ok_clean, "violation": grad_viol_clean,
                                        "rate": round(grad_ok_clean / (grad_ok_clean + grad_viol_clean), 3) if (grad_ok_clean + grad_viol_clean) else None},
            "median_delta_vaf": round(float(np.median(deltas)), 3) if deltas else None,
        },
        "ccf_tier_distribution_clean": dict(tier_c),
        "ccf_conversion": "neutral→min(1,2*VAF)（diploid het, purity~1, mult=1）; loh→VAF; gain/loss→undet(multiplicity 歧義)",
        "ccf_peaks_clean(bin,count)": peaks,
        "verdict": "祖先>=後代 VAF 梯度成立率 = clonal hierarchy 的 data-support; CCF 峰離散 = 離散 subclone 層級",
    }
    json.dump(out, open(f"{RPT}/data/sm_ccf_tiers.json", "w"), ensure_ascii=False, indent=1)
    g = out["ancestor_ge_descendant_VAF_gradient"]
    print(f"祖先>=後代 VAF 梯度: ok={g['ancestor_higher']} violation={g['descendant_higher(violation)']} "
          f"tie={g['tie']} → data_support={g['data_support_rate']} (null=0.5)")
    print(f"  clean only: rate={g['clean_only(loh+neutral)']['rate']}")
    print(f"CCF tier (clean): {dict(tier_c)}")
    print(f"CCF peaks: {peaks}")


if __name__ == "__main__":
    main()
