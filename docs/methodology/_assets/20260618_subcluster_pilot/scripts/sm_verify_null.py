"""[S5 驗證 — null/permutation on 2×2 連鎖]

每對的 2×2 (RR/RA/AR/AA) 已是 read 級等位的聯合分布。Fisher exact test = 「打亂等位標籤」的**精確** permutation
null (對 2×2 無近似)。對每個 confirmed pair (powered & 兩端 somatic) 跑 Fisher：
- 顯著 (p<0.05) + 觀測 AA < 期望 → 真互斥 (sibling/allelic 結構)；觀測 AA > 期望 → 真共現 (nested/co_linked)。
- chance 期望僅 ~5% 對顯著 → 若 confirmed pair 多數顯著 = 連鎖非隨機。
鑑別 subclone vs allelic 靠 HP gate (同 vs 異)，非 Fisher。diff_hp 互斥 = negative control (應同樣顯著但非 subclone)。
另跑小樣本 Monte-Carlo allele-shuffle 對照 Fisher 一致性。輸出 sm_verify_null.json。
"""
import json
from collections import Counter
from scipy.stats import fisher_exact

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
IN = f"{A}/sm_linkage_genomewide.json"
OUT = f"{A}/sm_verify_null.json"
ALPHA = 0.05


def expected_aa(rr, ra, ar, aa):
    n = rr + ra + ar + aa
    if n == 0:
        return 0.0
    n_a_alt = ar + aa
    n_b_alt = ra + aa
    return n_a_alt * n_b_alt / n


def main():
    d = json.load(open(IN))
    res = Counter()           # (group, rel_class, sig?) counts
    per_rel = {}              # rel -> {n, n_sig, frac_sig}
    examples = []
    for p in d["pairs"]:
        if not (p["powered"] and p["somatic_a"] and p["somatic_b"]):
            continue
        c = p["cells"]
        rr, ra, ar, aa = c["RR"], c["RA"], c["AR"], c["AA"]
        try:
            _, pval = fisher_exact([[rr, ra], [ar, aa]])
        except ValueError:
            continue
        sig = pval < ALPHA
        exp = expected_aa(rr, ra, ar, aa)
        direction = "exclusion" if aa < exp else ("co_occur" if aa > exp else "neutral")
        grp = "sameHP" if p["same_hp"] else "diffHP"
        rel = p["rel"]
        per_rel.setdefault(rel, {"n": 0, "n_sig": 0})
        per_rel[rel]["n"] += 1
        per_rel[rel]["n_sig"] += int(sig)
        res[(grp, "sig" if sig else "ns")] += 1
        if p["same_hp"] and sig and aa == 0 and (ar + ra) >= 4 and len(examples) < 40:
            examples.append({"chrom": p["chrom"], "a": p["a"], "b": p["b"], "rel": rel,
                             "cells": c, "p": float(f"{pval:.2e}"), "exp_aa": round(exp, 2),
                             "coread": p["coread"]})
    for rel, v in per_rel.items():
        v["frac_sig"] = round(v["n_sig"] / v["n"], 3) if v["n"] else None

    out = {
        "alpha": ALPHA,
        "method": "Fisher exact = conditional 2x2 test (both margins fixed, hypergeometric); "
                  "free-margin Monte-Carlo allele-shuffle cross-check in sm_summary.json",
        "caveat_F4": "Fisher-sig = 必要非充分: 只證『2x2 非隨機結構』, 不分 subclone vs allelic. "
                     "鑑別 subclone 全靠 HP gate. diffHP negative control ~50% 也顯著即證此點. "
                     "null 真正有效部分 = sparse 類 ~2% 顯著 (≈機率底線, 證 power-gate 不從噪聲造結構).",
        "per_rel": per_rel,
        "sameHP_sig": res[("sameHP", "sig")], "sameHP_ns": res[("sameHP", "ns")],
        "diffHP_sig_negctrl": res[("diffHP", "sig")], "diffHP_ns_negctrl": res[("diffHP", "ns")],
        "interpretation": "Fisher-顯著=非隨機結構(必要非充分); subclone vs allelic 靠 HP gate 非 Fisher; "
                          "diffHP(allelic) ~50% 也顯著=negative control 證 Fisher 無鑑別力; "
                          "sparse 2% ≈ 機率底線 = null 唯一真正有效的驗證",
        "sig_mutual_excl_sameHP_examples": examples,
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"NULL-VERIFY per_rel={ {k: (v['n_sig'], v['n']) for k, v in per_rel.items()} } "
          f"sameHP sig/ns={res[('sameHP','sig')]}/{res[('sameHP','ns')]} "
          f"diffHP(negctrl) sig/ns={res[('diffHP','sig')]}/{res[('diffHP','ns')]} -> {OUT}")


if __name__ == "__main__":
    main()
