#!/usr/bin/env python3
"""[v4 records 合併甲基] records_v3 + phylo_cpp_wg_full_methylation.json → records_v4.json。
純讀已落檔 join(key=chrom:pos)。每位點加 8 甲基欄(m_*); 無甲基資料(reads 過少/無 CpG)→ None。"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
MFIELDS = ["n_cpg", "mean_beta", "std_beta", "frac_hypo", "frac_hyper", "dbeta_group", "dbeta_tn", "dbeta_hp"]


def main():
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v3.json"))
    meth = json.load(open(f"{A}/phylo_cpp_wg_full_methylation.json"))
    hit = 0
    for r in rec:
        k = f"{r['chrom']}:{int(r['pos'])}"
        m = meth.get(k)
        if m:
            hit += 1
            for f in MFIELDS:
                r["m_" + f] = m.get(f)
        else:
            for f in MFIELDS:
                r["m_" + f] = None
    json.dump(rec, open(f"{A}/phylo_cpp_wg_full_records_v4.json", "w"))
    out = {"n": len(rec), "methylation_matched": hit,
           "with_dbeta_group": sum(1 for r in rec if r.get("m_dbeta_group") is not None),
           "with_dbeta_tn": sum(1 for r in rec if r.get("m_dbeta_tn") is not None),
           "no_methylation": sum(1 for r in rec if r.get("m_mean_beta") is None)}
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
