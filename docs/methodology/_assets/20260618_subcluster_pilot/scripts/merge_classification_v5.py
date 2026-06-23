#!/usr/bin/env python3
"""[v5 records 合併 8類分類 + 列聯型] records_v4 + label_composition + contingency_type → records_v5.json。
cat8 邏輯與 extract_label_composition.py 完全一致(確保儀表板=文件)。純讀已落檔。"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"


def main():
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v4.json"))
    comp = json.load(open(f"{A}/label_composition.json"))
    ct = json.load(open(f"{A}/contingency_type.json"))

    def feat(r):
        c = comp.get(f"{r['chrom']}:{int(r['pos'])}", {}); hp = set(c.get("hp", []))
        f1, s1, f2, s2 = ("1" in hp, "1-1" in hp, "2" in hp, "2-1" in hp)
        n_lab = len(hp)
        germ = (f1 or s1) and (f2 or s2)
        som = (f1 and s1) or (f2 and s2)
        # axis(與 contingency 一致)
        if n_lab == 1:
            axis = "single_label_somatic" if (s1 or s2) else "single_label_germline"
        elif germ:
            axis = "germline_both_families"
        elif som:
            axis = "somatic_one_family"
        else:
            axis = "other"
        return n_lab, axis, som

    for r in rec:
        n_lab, axis, som = feat(r)
        coarse = r["coarse_ng"]; loh = r["cn_state"] == "loh"; low = r["n"] < 20
        if n_lab <= 1:
            cat = "A" if coarse >= 2 else ("B1" if loh else ("B2" if low else "B3"))
        else:
            cat = ("C-" + r["wstate"]) if coarse >= 2 else "D"
        r["cat8"] = cat
        r["axis"] = axis
        r["has_som"] = som  # somatic 子單倍型存在(不論是否也有 germline) — 對齊 doc 的 LOH-somatic=3,553
        r["ctype"] = ct.get(f"{r['chrom']}:{int(r['pos'])}", {}).get("type", "single_cluster")

    json.dump(rec, open(f"{A}/phylo_cpp_wg_full_records_v5.json", "w"))
    out = {"n": len(rec), "cat8": dict(Counter(r["cat8"] for r in rec)),
           "ctype": dict(Counter(r["ctype"] for r in rec)),
           "subclone_candidate(A or many1xsomatic)": sum(1 for r in rec if r["cat8"] == "A" or (r["ctype"] == "many1_結構>標籤" and r["axis"] in ("single_label_somatic", "somatic_one_family"))),
           "loh_somatic_structure": sum(1 for r in rec if r["cn_state"] == "loh" and r["coarse_ng"] >= 2 and r["axis"] in ("single_label_somatic", "somatic_one_family", "germline_both_families")),
           "sumcheck": sum(Counter(r["cat8"] for r in rec).values()) == len(rec)}
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
