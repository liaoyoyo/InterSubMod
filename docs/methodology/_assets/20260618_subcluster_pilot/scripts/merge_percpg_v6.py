#!/usr/bin/env python3
"""[v6 records 併逐 CpG verdict + label_structure] records_v5 + phylo_cpp_wg_full_percpg.json → records_v6.json。
加: pc_verdict(cis-ASM/subclone_novel/...) / pc_best_axis / pc_cluster_run / pc_overlap_max /
   label_structure(hp/allele PERMANOVA p+dispersion 旗)。純讀已落檔。"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"


def main():
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v5.json"))
    pc = json.load(open(f"{A}/phylo_cpp_wg_full_percpg.json"))
    for r in rec:
        d = pc.get(f"{r['chrom']}:{int(r['pos'])}", {})
        r["pc_verdict"] = d.get("verdict", "no_meth")
        r["pc_best_axis"] = d.get("best_label_axis")
        r["pc_best_nsig"] = d.get("best_label_nsig", 0)
        r["pc_cluster_run"] = d.get("cluster_split_run")
        r["pc_overlap_max"] = round(max(d.get("overlap", {}).values()), 3) if d.get("overlap") else None
        ls = d.get("label_structure", {})
        r["ls_hp_p"] = ls.get("hp_perm_p"); r["ls_hp_disp"] = ls.get("hp_disp_warn")
        r["ls_al_p"] = ls.get("al_perm_p"); r["ls_al_disp"] = ls.get("al_disp_warn")
    json.dump(rec, open(f"{A}/phylo_cpp_wg_full_records_v6.json", "w"))
    out = {"n": len(rec), "pc_verdict": dict(Counter(r["pc_verdict"] for r in rec)),
           "subclone_novel": sum(1 for r in rec if r["pc_verdict"] == "subclone_novel"),
           "with_ls_hp_p": sum(1 for r in rec if r.get("ls_hp_p") is not None),
           "ls_hp_sig_no_disp": sum(1 for r in rec if (r.get("ls_hp_p") is not None and r["ls_hp_p"] < 0.05 and not r.get("ls_hp_disp")))}
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
