#!/usr/bin/env python3
"""[每位點標籤組成 + 修正分類樹] 走 phylo_groups.tsv 抓實際 hp 標籤集合 + allele 集合,
依正確生物學重分類(HP1-1⟺ALT somatic 軸; 單標籤+多結構=subclone候選; LOH+somatic軸=ALT對齊;
無法區分拆 LOH-生物性 vs 低覆蓋-技術性)。輸出 label_composition.json + corrected_tree.json。
比例一律對總數 34736(絕對)。"""
import os, csv, glob, json
from collections import Counter
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"


def work(tsv):
    try:
        region = os.path.dirname(os.path.dirname(tsv)); parts = os.path.basename(region).split("_")
        chrom, start = parts[0], int(parts[1])
        hp, al = set(), set()
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            h = r.get("hp", "")
            if h in ("1", "1-1", "2", "2-1"): hp.add(h)
            a = r.get("alt_support", "")
            if a in ("ALT", "REF"): al.add(a)
        return (f"{chrom}:{start}", {"hp": sorted(hp), "al": sorted(al)})
    except Exception:
        return None


def main():
    tsvs = glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)
    print(f"regions: {len(tsvs)}", flush=True)
    with Pool(24) as p:
        comp = dict(r for r in p.map(work, tsvs, chunksize=60) if r)
    json.dump(comp, open(f"{A}/label_composition.json", "w"))

    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v4.json"))
    N = len(rec)

    # --- verify biology: HP1-1 present => ALT present? ---
    som_loci = [r for r in rec if "1-1" in comp.get(f"{r['chrom']}:{int(r['pos'])}", {}).get("hp", []) or "2-1" in comp.get(f"{r['chrom']}:{int(r['pos'])}", {}).get("hp", [])]
    som_with_alt = sum(1 for r in som_loci if "ALT" in comp.get(f"{r['chrom']}:{int(r['pos'])}", {}).get("al", []))
    bio = {"loci_with_somatic_subhap": len(som_loci), "of_which_have_ALT": som_with_alt,
           "consistency_pct": round(100 * som_with_alt / len(som_loci), 2) if som_loci else None}

    def feat(r):
        c = comp.get(f"{r['chrom']}:{int(r['pos'])}", {}); hp = set(c.get("hp", []))
        f1, s1 = "1" in hp, "1-1" in hp; f2, s2 = "2" in hp, "2-1" in hp
        n_lab = len(hp)
        germ_axis = (f1 or s1) and (f2 or s2)          # 兩單倍型都在 = HP1 vs HP2 軸
        som_axis = (f1 and s1) or (f2 and s2)          # 同單倍型 germline+somatic = REF/ALT 軸
        return n_lab, germ_axis, som_axis, hp

    # --- corrected classification ---
    cat = Counter(); detail = {}
    rows = []
    for r in rec:
        n_lab, germ, som, hp = feat(r)
        coarse = r["coarse_ng"]; loh = r["cn_state"] == "loh"; low = r["n"] < 20
        if n_lab <= 1:
            if coarse >= 2:
                k = "A_subclone候選(單標籤+多結構)"
            else:
                # 真單群: 拆 生物性(LOH/足覆蓋同質) vs 技術性(低覆蓋)
                if loh: k = "B1_真單群-LOH(生物性同質)"
                elif low: k = "B2_無法區分-低覆蓋(技術性)"
                else: k = "B3_真單群-非LOH足覆蓋(生物同質)"
        else:
            if coarse >= 2:
                k = "C_有結構+有軸(" + r["wstate"] + ")"
            else:
                k = "D_可測無結構(有軸但單群)"
        cat[k] += 1
        rows.append((k, r, loh, som, germ))
    # LOH + somatic-axis + structure 高亮(ALT/癌突變對齊好說明)
    loh_som_struct = sum(1 for k, r, loh, som, germ in rows if loh and som and r["coarse_ng"] >= 2)
    out = {"N": N, "biology_check_HP1_1_implies_ALT": bio,
           "categories_abs": {k: {"n": v, "pct_of_total": round(100 * v / N, 2)} for k, v in cat.most_common()},
           "highlight_LOH+somatic_axis+structure": {"n": loh_som_struct, "pct_of_total": round(100 * loh_som_struct / N, 2)},
           "sum_check": sum(cat.values()) == N}
    json.dump(out, open(f"{A}/corrected_tree.json", "w"), ensure_ascii=False, indent=1)
    print("DONE")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
