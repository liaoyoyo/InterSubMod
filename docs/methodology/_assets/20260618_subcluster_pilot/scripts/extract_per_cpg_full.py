#!/usr/bin/env python3
"""[全量逐 CpG 差異 + subclone 驗證] 走全 34,736 region dirs, 每位點 analyze_region(per_cpg_diff v2):
逐 CpG Fisher+MWU 各分組 + cluster-split×label CpG-set overlap → verdict(cis-ASM/subclone_novel/...)。
輸出 phylo_cpp_wg_full_percpg.json(每位點摘要, 不含 per-CpG list 控大小)。並產 label_structure(從 significance.json)。"""
import os, glob, json, sys
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__)); import per_cpg_diff as P

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"


def work(region):
    try:
        base = os.path.basename(region).split("_"); chrom, start = base[0], int(base[1])
        key = f"{chrom}:{start}"
        # label_structure from significance.json
        ls = {}
        sj = f"{region}/clustering/significance.json"
        if os.path.exists(sj):
            d = json.load(open(sj)); L = d.get("label_structure", {})
            ls = {"hp_perm_p": L.get("hp_permanova_p"), "hp_perm_f": L.get("hp_permanova_f"), "hp_disp_warn": L.get("hp_dispersion_warning"),
                  "al_perm_p": L.get("allele_permanova_p"), "al_perm_f": L.get("allele_permanova_f"), "al_disp_warn": L.get("allele_dispersion_warning")}
        res = P.analyze_region(region)
        if res is None:
            return (key, {"label_structure": ls, "verdict": "no_meth"})
        # 壓縮: 只留摘要 + overlap containment
        out = {"label_structure": ls, "verdict": res["verdict"], "best_label_axis": res["best_label_axis"],
               "best_label_nsig": res["best_label_nsig"], "cluster_split_run": res["cluster_split_run"], "cluster_split_nsig": res["cluster_split_nsig"],
               "overlap": {k: v["containment"] for k, v in res["overlap"].items()},
               "axis_nsig_mwu": {k: v["n_sig_mwu"] for k, v in res["groupings"].items()},
               "axis_nsig_fisher": {k: v["n_sig_fisher"] for k, v in res["groupings"].items()}}
        return (key, out)
    except Exception as e:
        return (os.path.basename(region), {"verdict": "ERR:" + str(e)[:40]})


def main():
    regions = [os.path.dirname(os.path.dirname(p)) for p in glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)]
    print(f"regions: {len(regions)}", flush=True)
    with Pool(24) as p:
        res = [r for r in p.map(work, regions, chunksize=20) if r]
    out = dict(res)
    json.dump(out, open(f"{A}/phylo_cpp_wg_full_percpg.json", "w"))
    vc = Counter(v.get("verdict") for v in out.values())
    S = {"n": len(out), "verdict_dist": dict(vc),
         "subclone_novel": sum(1 for v in out.values() if v.get("verdict") == "subclone_novel"),
         "cis_asm_cluster": sum(1 for v in out.values() if str(v.get("verdict", "")).startswith("cis-ASM")),
         "with_label_structure": sum(1 for v in out.values() if v.get("label_structure", {}).get("hp_perm_p") is not None)}
    json.dump(S, open(f"{A}/percpg_summary.json", "w"), ensure_ascii=False, indent=1)
    print("DONE " + json.dumps(S, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
