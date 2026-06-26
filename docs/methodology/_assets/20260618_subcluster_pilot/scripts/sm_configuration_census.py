"""[2×2 configuration census] 把每個 sSNV 對按 2×2 cell-pattern (RR/RA/AR/AA 哪些非零) 歸類 + 實際計數。
cell: RR=(a-REF,b-REF) RA=(a-REF,b-ALT) AR=(a-ALT,b-REF) AA=(a-ALT,b-ALT)。
配置由 (RA>0, AR>0, AA>0) 決定 (RR 幾乎恆有):
  RA+AR, AA=0     → mutual_excl 互斥 (兩 ALT 從不共現) = sibling(sameHP)/allelic(diffHP)
  RA+AA, AR=0     → nested a⊂b (a-ALT 必伴 b-ALT; b 有額外) → b 祖先→a 後代
  AR+AA, RA=0     → nested b⊂a (b-ALT 必伴 a-ALT; a 有額外) → a 祖先→b 後代
  AA only         → co_linked 共連 (兩 ALT 完美共現)
  RA+AR+AA all    → independent (無乾淨結構)
  其他(缺 AA 且缺一 off-diag) → degenerate/sparse
每配置計: pair 數 + distinct sSNV 數，並 split sameHP/diffHP + CN state。輸出 sm_configuration_census.json。
"""
import json
from collections import Counter, defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
CNBED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
OUT = f"{A}/sm_configuration_census.json"

CONFIG = {
    (1, 1, 0): ("mutual_excl", "互斥 RA+AR (AA=0): 兩 ALT 從不共現 → sibling/allelic"),
    (1, 0, 1): ("nested_a_in_b", "RA+AA (AR=0): a-ALT 必伴 b-ALT → b 祖先→a 後代"),
    (0, 1, 1): ("nested_b_in_a", "AR+AA (RA=0): b-ALT 必伴 a-ALT → a 祖先→b 後代"),
    (0, 0, 1): ("co_linked", "AA only (RA=AR=0): 兩 ALT 完美共現 → 同 lineage"),
    (1, 1, 1): ("independent", "RA+AR+AA 皆有: 無乾淨結構"),
    (1, 0, 0): ("degenerate_RAonly", "只 RA: a 恆 REF (此對無 a-ALT)"),
    (0, 1, 0): ("degenerate_ARonly", "只 AR: b 恆 REF"),
    (0, 0, 0): ("degenerate_RRonly", "只 RR: 無 ALT 共現"),
}


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


def census(pairs, cn, label):
    pair_cnt = Counter()
    ssnv = defaultdict(set)
    cn_by_cfg = defaultdict(Counter)
    for p in pairs:
        c = p["cells"]
        key = (int(c["RA"] > 0), int(c["AR"] > 0), int(c["AA"] > 0))
        name = CONFIG.get(key, ("other", "?"))[0]
        hp = "sameHP" if p["same_hp"] else "diffHP"
        pair_cnt[(name, hp)] += 1
        ka, kb = f"{p['chrom']}:{p['a']}", f"{p['chrom']}:{p['b']}"
        ssnv[(name, hp)].update((ka, kb))
        st = cn_state(cn, p["chrom"], p["a"])
        cn_by_cfg[name][st] += 1
    rows = []
    names = ["mutual_excl", "nested_a_in_b", "nested_b_in_a", "co_linked", "independent",
             "degenerate_RAonly", "degenerate_ARonly", "degenerate_RRonly", "other"]
    for nm in names:
        sp = pair_cnt.get((nm, "sameHP"), 0)
        dp = pair_cnt.get((nm, "diffHP"), 0)
        if sp == 0 and dp == 0:
            continue
        rows.append({"config": nm,
                     "desc": next((v[1] for k, v in CONFIG.items() if v[0] == nm), ""),
                     "pairs_sameHP": sp, "pairs_diffHP": dp, "pairs_total": sp + dp,
                     "distinct_sSNV_sameHP": len(ssnv.get((nm, "sameHP"), set())),
                     "distinct_sSNV_diffHP": len(ssnv.get((nm, "diffHP"), set())),
                     "cn_state_dist": dict(cn_by_cfg[nm])})
    return rows


def main():
    d = json.load(open(f"{A}/sm_linkage_genomewide.json"))
    cn = load_cn()
    all_pairs = d["pairs"]
    powered = [p for p in all_pairs if p["powered"]]
    pow_som = [p for p in powered if p["somatic_a"] and p["somatic_b"]]
    out = {
        "cell_legend": "RR=(a-REF,b-REF) RA=(a-REF,b-ALT) AR=(a-ALT,b-REF) AA=(a-ALT,b-ALT)",
        "filters": {"all_recorded_coread>=2": len(all_pairs),
                    "powered_coread>=6": len(powered),
                    "powered+both_somatic": len(pow_som)},
        "census_all_recorded": census(all_pairs, cn, "all"),
        "census_powered_somatic": census(pow_som, cn, "pow_som"),
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("=== CONFIGURATION CENSUS (powered + both somatic) ===")
    print(f"{'config':18s} {'pairs(same/diff)':18s} {'distinct sSNV(same/diff)':24s} CN(same-axis a)")
    for r in out["census_powered_somatic"]:
        cns = ",".join(f"{k}:{v}" for k, v in sorted(r["cn_state_dist"].items(), key=lambda x: -x[1]))
        print(f"{r['config']:18s} {r['pairs_sameHP']:>7}/{r['pairs_diffHP']:<9} "
              f"{r['distinct_sSNV_sameHP']:>10}/{r['distinct_sSNV_diffHP']:<12} {cns}")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
