"""[S5b 修正 + 小 summary] 套用對抗稽核 5 修正，產獨立可重算的小 summary JSON。
F1 per-sSNV (非 per-pair) distinct 參與數 + dense-cluster 灌水揭露
F2 偽影量化 (dense-uniform-cluster + CN-gain 疊加; 只量化不清理 — 用戶選 C)
F4 free-margin Monte-Carlo allele-shuffle 對照 Fisher (交付承諾的 MC cross-check)
F5 全標 Tier-R
輸出 sm_summary.json (小檔, 可獨立 recompute)。
"""
import json
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import fisher_exact

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
CNBED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
OUT = f"{A}/sm_summary.json"
np.random.seed(42)


def load_cn():
    iv = defaultdict(list)
    for ln in open(CNBED):
        p = ln.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        iv[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in iv:
        iv[c].sort()
    return iv


def cn_state(cn, chrom, pos):
    for s, e, lab in cn.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"


def mc_pval(rr, ra, ar, aa, K=2000):
    """free-margin allele-shuffle: 固定 n / n_A_alt / n_B_alt, 獨立放置 → 經驗 P(AA_sim<=aa)。"""
    n = rr + ra + ar + aa
    na = ar + aa
    nb = ra + aa
    if n == 0 or na == 0 or nb == 0:
        return None
    cnt = 0
    for _ in range(K):
        sa = np.random.choice(n, na, replace=False)
        mask = np.zeros(n, dtype=bool)
        mask[sa] = True
        sb = np.random.choice(n, nb, replace=False)
        aa_sim = int(mask[sb].sum())
        if aa_sim <= aa:
            cnt += 1
    return cnt / K


def main():
    d = json.load(open(f"{A}/sm_linkage_genomewide.json"))
    led = json.load(open(f"{A}/sm_completeness_ledger.json"))
    nul = json.load(open(f"{A}/sm_verify_null.json"))
    seq = json.load(open(f"{A}/sm_seqc2_overlay.json"))
    cen = d["census"]
    cn = load_cn()

    # ---- F1: per-sSNV distinct participation (powered & somatic both & sameHP) ----
    rel_sets = {"nested": set(), "sibling": set(), "co_linked": set()}
    mut_excl_pairs = []
    for p in d["pairs"]:
        if not (p["powered"] and p["somatic_a"] and p["somatic_b"] and p["same_hp"]):
            continue
        ka, kb = f"{p['chrom']}:{p['a']}", f"{p['chrom']}:{p['b']}"
        if p["rel"] in ("nested_a_in_b", "nested_b_in_a"):
            rel_sets["nested"].update((ka, kb))
        elif p["rel"] == "mutual_excl":
            rel_sets["sibling"].update((ka, kb))
            mut_excl_pairs.append(p)
        elif p["rel"] == "co_linked":
            rel_sets["co_linked"].update((ka, kb))
    distinct_any = rel_sets["nested"] | rel_sets["sibling"] | rel_sets["co_linked"]
    per_sSNV = {k: len(v) for k, v in rel_sets.items()}
    per_sSNV["any_confirmed_sameHP_relation"] = len(distinct_any)

    # ---- F2: dense-uniform-cluster artifact quantification (linked-somatic sSNV) ----
    ls = [(c["chrom"], c["pos"], c.get("vaf")) for c in cen.values()
          if c.get("somatic") and c["n_realized_links"] >= 1]
    ls.sort()
    clusters = []
    cur = [ls[0]]
    for x in ls[1:]:
        if x[0] == cur[-1][0] and x[1] - cur[-1][1] <= 5000:
            cur.append(x)
        else:
            clusters.append(cur); cur = [x]
    clusters.append(cur)
    dense = [c for c in clusters if len(c) >= 5]
    dense_sSNV = sum(len(c) for c in dense)
    cn_counter = Counter(cn_state(cn, c, p) for c, p, _ in ls)
    dense_cn = Counter(cn_state(cn, c[0][0], c[0][1]) for c in dense)

    # ---- F4: Monte-Carlo allele-shuffle cross-check vs Fisher (sample of mutual_excl sameHP) ----
    samp = mut_excl_pairs[:300]
    agree = both_sig = both_ns = disagree = 0
    diffs = []
    for p in samp:
        c = p["cells"]
        rr, ra, ar, aa = c["RR"], c["RA"], c["AR"], c["AA"]
        try:
            _, fp = fisher_exact([[rr, ra], [ar, aa]])
        except ValueError:
            continue
        mp = mc_pval(rr, ra, ar, aa)
        if mp is None:
            continue
        diffs.append(abs(fp - mp))
        fs, ms = fp < 0.05, mp < 0.05
        if fs == ms:
            agree += 1
            both_sig += int(fs)
            both_ns += int(not fs)
        else:
            disagree += 1

    summary = {
        "scope": "HCC1395 single sample ⭐3; Tier-R (same-read ≤50kb) only — same-PS deferred",
        "universe": {"total": d["universe"]["n_union"], "tp": d["universe"]["n_tp"], "fp": d["universe"]["n_fp"]},
        "completeness_ledger_TierR": led["buckets"],
        "ledger_sum_check": led["bucket_sum_check"]["ok"],
        "F1_per_sSNV_distinct": per_sSNV,
        "F1_note": f"distinct linked sSNV={led['buckets']['linked']} 為乾淨分母; pair-count={d['n_pairs_recorded']} "
                   f"被 dense 群灌水 (max {max(c['n_realized_links'] for c in cen.values())} links/sSNV)",
        "F2_artifact_quantification": {
            "linked_somatic_total": len(ls),
            "in_dense_uniform_cluster(>=5 in <=5kb)": dense_sSNV,
            "n_dense_clusters": len(dense),
            "linked_somatic_by_CN_state": dict(cn_counter),
            "dense_clusters_by_CN_state": dict(dense_cn),
            "note": "dense-uniform 群 + CN-gain = segdup/mapping-artifact 或 amplicon 嫌疑; 用戶選 C 故只量化不清理; "
                    "mappability/segdup track 缺 → 殘留偽影為已知限制",
        },
        "F3_gamma_class_reframed": {
            "fp_source_somatic_linked": led["union_recovered_fp_somatic_linked"]["count"],
            "honest_label": "FP-labeled ONT 呼叫 + somatic-like normal(VAF<0.05) + 單分子連鎖 = 候選, 需正交驗證",
            "leading_alternative": "FP-source + clustered-artifact (見 F2) 為首要替代解釋, 非『SEQC2 漏掉的確證 somatic』",
            "convergence_3204_vs_3200_is_tautological": True,
        },
        "F4_fisher_vs_MC": {
            "n_sampled_mutual_excl_sameHP": len(samp),
            "agree": agree, "both_sig": both_sig, "both_ns": both_ns, "disagree": disagree,
            "median_abs_p_diff": round(float(np.median(diffs)), 4) if diffs else None,
            "note": "Fisher(both-margin-fixed) vs free-margin MC: 高一致=Fisher 結論穩健; Fisher-sig 仍不分 subclone/allelic",
            "diffHP_negctrl_sig_frac": round(nul["diffHP_sig_negctrl"] /
                                             (nul["diffHP_sig_negctrl"] + nul["diffHP_ns_negctrl"]), 3),
            "sparse_floor_sig_frac": nul["per_rel"]["sparse"]["frac_sig"],
        },
        "validated_solid": {
            "chr17_48360161_reconstructed": "α→β2 nested + γ∥α sibling, VAF 0.82>0.48>0.18, 0 violations",
            "method_efficient": "group-based single-pass, 5 parallel workers, ~40min genome-wide",
            "evolution": {k: json.load(open(f"{A}/sm_evolution.json"))["aggregate"][k]
                          for k in ("n_components", "n_with_nested", "n_with_sibling", "n_sibling_and_nested", "n_with_cycle")},
        },
        "redlines": "6,146 components = 局部 haplotype-連鎖區域非 subclone; regional≤read-span 非 genome-wide tree; "
                    "甲基 characterize 非偵測; ⭐3 分子證據非 single-cell confirmation",
    }
    json.dump(summary, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("SUMMARY written:")
    print(f"  F1 per-sSNV distinct: {per_sSNV}")
    print(f"  F2 artifact: {dense_sSNV} sSNV in {len(dense)} dense clusters; linked_somatic by CN={dict(cn_counter)}")
    print(f"  F4 Fisher-vs-MC: agree={agree}/{len(samp)} (both_sig={both_sig} both_ns={both_ns} disagree={disagree}) "
          f"median|Δp|={summary['F4_fisher_vs_MC']['median_abs_p_diff']}")
    print(f"  -> {OUT}")


if __name__ == "__main__":
    main()
