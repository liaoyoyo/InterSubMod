#!/usr/bin/env python3
"""[每位點 HP/allele 標籤多樣性] 走保留 region dirs 的 phylo_groups.tsv(peeled tumor reads),
數每位點 distinct hp / distinct HP-family(1/2) / distinct allele(ALT/REF)。
用途: 把「無結構(coarse=1)」拆成「真單群(可測但單一)」vs「無法區分(標籤不足)」。
輸出 hp_allele_diversity.json(key=chrom:pos)。純讀, 零重算。"""
import os, csv, glob, json
from collections import Counter
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"


def work(tsv):
    try:
        region = os.path.dirname(os.path.dirname(tsv))
        parts = os.path.basename(region).split("_")
        chrom, start = parts[0], int(parts[1])
        hps, fams, alls = set(), set(), set()
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            hp = r.get("hp", "")
            hps.add(hp)
            if hp and hp[0] in ("1", "2"):
                fams.add(hp[0])
            al = r.get("alt_support", "")
            if al in ("ALT", "REF"):
                alls.add(al)
        return (f"{chrom}:{start}", {"n_hp": len(hps), "n_hpfam": len(fams), "n_allele": len(alls)})
    except Exception:
        return None


def main():
    tsvs = glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)
    print(f"regions: {len(tsvs)}", flush=True)
    with Pool(24) as p:
        res = [r for r in p.map(work, tsvs, chunksize=60) if r]
    div = dict(res)
    json.dump(div, open(f"{A}/hp_allele_diversity.json", "w"))

    # 直接做無結構拆解
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v4.json"))
    ns = [r for r in rec if r["coarse_ng"] == 1]
    N = len(ns)

    def cls(r):
        d = div.get(f"{r['chrom']}:{int(r['pos'])}", {})
        nf = d.get("n_hpfam", 0)
        na = d.get("n_allele", 0)
        loh = r["cn_state"] == "loh"
        if nf < 2 and na < 2:
            return "無法區分(HP+allele 皆單一)"
        if nf < 2:
            return "僅單 HP-family(無法測 HP 對齊)"
        if na < 2:
            return "僅單 allele(無法測 allele 對齊)"
        return "真單群(雙軸皆≥2 可測但無甲基結構)"

    cat = Counter(cls(r) for r in ns)
    loh_in_singlefam = sum(1 for r in ns if r["cn_state"] == "loh" and div.get(f"{r['chrom']}:{int(r['pos'])}", {}).get("n_hpfam", 0) < 2)
    out = {"no_structure_total": N, "breakdown": dict(cat),
           "LOH_total": sum(1 for r in ns if r["cn_state"] == "loh"),
           "LOH_and_single_HPfam": loh_in_singlefam,
           "matched_diversity": sum(1 for r in ns if f"{r['chrom']}:{int(r['pos'])}" in div)}
    json.dump(out, open(f"{A}/no_structure_breakdown.json", "w"), ensure_ascii=False, indent=1)
    print("DONE " + json.dumps(out, ensure_ascii=False, indent=1), flush=True)


if __name__ == "__main__":
    main()
