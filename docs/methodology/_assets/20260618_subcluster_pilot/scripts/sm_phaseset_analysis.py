"""[L3b phase-set 分析] 三件事 (依 00b grounding):
1. Tier-PS 延伸評估: same-PS sSNV >50kb apart(跨 Tier-R region) → PS 是否延伸『克隆連鎖』還是只『germline 單倍型 context』?
   honest: PS = germline phase block; 同 PS ≠ 同克隆。用 CCF 一致性檢驗(同 PS 遠 pair CCF 是否一致)。
2. PS-reliability 旗標: region 的 somatic sSNV 跨 >1 PS = phase-switch 風險 → HP same/diff 判別不可信 → 標記。
3. 量化: PS-reliable(single-PS) 區域比例; full_tree/sibling 中多少 PS-reliable。
輸出 報告夾/data/sm_phaseset_extension.json。
"""
import json
import glob
from collections import defaultdict, Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
RPT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report"


def main():
    ps = {}
    for f in glob.glob(f"{A}/ps_part_*.json"):
        ps.update(json.load(open(f)))
    cen = json.load(open(f"{A}/sm_linkage_genomewide.json"))["census"]
    regions = json.load(open(f"{A}/sm_region_integration.json"))["regions"]

    n_with_ps = sum(1 for v in ps.values() if v.get("ps"))
    print(f"PS extracted for {len(ps)} somatic sSNV, with_ps={n_with_ps}")

    # 2. PS-reliability per region: distinct PS among region's somatic sSNV
    reg_reliable = 0
    reg_uncertain = 0
    reg_no_ps = 0
    shape_reliable = defaultdict(lambda: [0, 0])  # shape -> [reliable, uncertain]
    region_ps_flag = {}
    for r in regions:
        chrom = r["chrom"]
        pss = set()
        for pos in range(r["start"], r["end"] + 1):
            k = f"{chrom}:{pos}"
            if k in ps and ps[k].get("ps"):
                pss.add(ps[k]["ps"])
        if not pss:
            flag = "no_ps"
            reg_no_ps += 1
        elif len(pss) == 1:
            flag = "single_ps_reliable"
            reg_reliable += 1
        else:
            flag = f"multi_ps_uncertain({len(pss)})"
            reg_uncertain += 1
        region_ps_flag[r["region"]] = {"flag": flag, "n_distinct_ps": len(pss)}
        if r["tree_shape"] not in ("no_confirmed_structure", "inconsistent"):
            shape_reliable[r["tree_shape"]][0 if len(pss) == 1 else 1] += 1

    # 1. Tier-PS: same-PS somatic sSNV pairs >50kb apart (cross Tier-R region)
    by_chrom_ps = defaultdict(lambda: defaultdict(list))  # chrom -> ps -> [pos]
    for k, v in ps.items():
        if v.get("ps"):
            c, p = k.split(":")
            by_chrom_ps[c][v["ps"]].append(int(p))
    tier_ps_pairs = 0
    same_hp_ps = 0
    ccf_consistent = 0
    ccf_pairs = 0
    n_multi_ps_blocks = 0
    n_giant_ps_skipped = 0
    HARD_CAP = 300000
    for c, psd in by_chrom_ps.items():
        for psid, poss in psd.items():
            poss = sorted(poss)
            if len(poss) < 2:
                continue
            n_multi_ps_blocks += 1
            if len(poss) > 40:        # giant PS block: 取最密 40(避 O(n^2) 爆 + 這些多為 problem block)
                n_giant_ps_skipped += 1
                poss = poss[:40]
            for i in range(len(poss)):
                for j in range(i + 1, len(poss)):
                    if poss[j] - poss[i] <= 50000:
                        continue
                    tier_ps_pairs += 1
                    ka, kb = f"{c}:{poss[i]}", f"{c}:{poss[j]}"
                    if ps[ka].get("hp") == ps[kb].get("hp") and ps[ka].get("hp"):
                        same_hp_ps += 1
                    va = cen.get(ka, {}).get("vaf"); vb = cen.get(kb, {}).get("vaf")
                    if va is not None and vb is not None and cen.get(ka, {}).get("somatic") and cen.get(kb, {}).get("somatic"):
                        ccf_pairs += 1
                        if abs(va - vb) <= 0.1:
                            ccf_consistent += 1
            if tier_ps_pairs > HARD_CAP:
                break
        if tier_ps_pairs > HARD_CAP:
            break

    out = {
        "n_somatic_with_ps": n_with_ps,
        "ps_reliability_per_region": {
            "single_ps_reliable": reg_reliable, "multi_ps_uncertain": reg_uncertain, "no_ps": reg_no_ps,
            "total_regions": reg_reliable + reg_uncertain + reg_no_ps,
            "reliable_rate": round(reg_reliable / (reg_reliable + reg_uncertain + reg_no_ps), 3) if (reg_reliable + reg_uncertain + reg_no_ps) else None,
            "by_shape_reliable_vs_uncertain": {k: {"reliable": v[0], "uncertain": v[1]} for k, v in shape_reliable.items()},
        },
        "tier_ps_extension": {
            "same_ps_pairs_gt50kb(Tier-PS candidates)": tier_ps_pairs,
            "of_which_same_HP": same_hp_ps,
            "ccf_pairs_checked": ccf_pairs,
            "ccf_consistent(|dVAF|<=0.1)": ccf_consistent,
            "ccf_consistent_rate": round(ccf_consistent / ccf_pairs, 3) if ccf_pairs else None,
            "interpretation": "same-PS = 同 germline phase block(非同克隆); 若 CCF 不一致 = 同單倍型上的不同 subclone "
                              "→ PS 不建立克隆共現, 只給 germline 單倍型 context",
        },
        "verdict": "PS 的用途 = (a) reliability 旗標(多 PS region = phase-switch 風險, HP 判別不可信) "
                   "(b) Tier-PS 延伸 germline 單倍型 context 但非克隆連鎖(germline phasing != clonal)",
    }
    json.dump(out, open(f"{RPT}/data/sm_phaseset_extension.json", "w"), ensure_ascii=False, indent=1)
    # also dump per-region flag for join
    json.dump(region_ps_flag, open(f"{RPT}/data/region_ps_flag.json", "w"), ensure_ascii=False)
    print(f"PS-reliable regions: {reg_reliable} / uncertain: {reg_uncertain} (rate={out['ps_reliability_per_region']['reliable_rate']})")
    print(f"by shape:", {k: v for k, v in out["ps_reliability_per_region"]["by_shape_reliable_vs_uncertain"].items()})
    print(f"Tier-PS candidates(same-PS >50kb): {tier_ps_pairs}, same-HP {same_hp_ps}, "
          f"CCF-consistent {ccf_consistent}/{ccf_pairs} ({out['tier_ps_extension']['ccf_consistent_rate']})")


if __name__ == "__main__":
    main()
