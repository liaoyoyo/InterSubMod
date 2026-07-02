#!/usr/bin/env python3
"""[每位點 cluster×label 列聯型] 走 phylo_groups.tsv 建 coarse_label×hp 列聯表,
分類 one2one / many:1(標籤→多群=結構>標籤=subclone-like) / 1:many(群→多標籤=跨標籤) / mixed。
交叉「標籤主軸」(somatic HP1-1主導 / germline 兩家 / 單標籤)。比例對總 34736。
輸出 contingency_type.json + contingency_summary.json。MIN_CELL=3(同 MIN_SZ)。"""
import os, csv, glob, json
from collections import Counter, defaultdict
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"
MIN_CELL = 3


def work(tsv):
    try:
        region = os.path.dirname(os.path.dirname(tsv)); parts = os.path.basename(region).split("_")
        chrom, start = parts[0], int(parts[1])
        T = defaultdict(Counter); hpset = set()
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            cl = r.get("coarse_label", "")
            if cl in ("outlier", "other", ""):
                continue
            hp = r.get("hp", "")
            if hp not in ("1", "1-1", "2", "2-1"):
                continue
            T[cl][hp] += 1; hpset.add(hp)
        clusters = list(T.keys())
        if len(clusters) < 2:
            return (f"{chrom}:{start}", {"type": "single_cluster", "nlab": len(hpset)})
        # many1: some label appears (>=MIN_CELL) in >=2 clusters
        lab_in_clusters = defaultdict(int)
        for cl in clusters:
            for lab, c in T[cl].items():
                if c >= MIN_CELL:
                    lab_in_clusters[lab] += 1
        has_many1 = any(v >= 2 for v in lab_in_clusters.values())
        # 1many: some cluster has >=2 labels (>=MIN_CELL)
        has_1many = any(sum(1 for c in T[cl].values() if c >= MIN_CELL) >= 2 for cl in clusters)
        if has_many1 and has_1many:
            typ = "mixed"
        elif has_many1:
            typ = "many1_結構>標籤"
        elif has_1many:
            typ = "1many_跨標籤"
        else:
            typ = "one2one"
        # 標籤主軸
        f1, s1, f2, s2 = ("1" in hpset, "1-1" in hpset, "2" in hpset, "2-1" in hpset)
        if len(hpset) == 1:
            axis = "single_label_somatic" if (s1 or s2) else "single_label_germline"
        elif (f1 or s1) and (f2 or s2):
            axis = "germline_both_families"
        elif (f1 and s1) or (f2 and s2):
            axis = "somatic_one_family"
        else:
            axis = "other"
        return (f"{chrom}:{start}", {"type": typ, "axis": axis, "nlab": len(hpset)})
    except Exception:
        return None


def main():
    tsvs = glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)
    print(f"regions: {len(tsvs)}", flush=True)
    with Pool(24) as p:
        ct = dict(r for r in p.map(work, tsvs, chunksize=60) if r)
    json.dump(ct, open(f"{A}/contingency_type.json", "w"))

    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v4.json"))
    N = len(rec)
    multi = [r for r in rec if r["coarse_ng"] >= 2]

    def k(r): return f"{r['chrom']}:{int(r['pos'])}"
    # 有結構位點: 列聯型 × 標籤主軸
    by_type = Counter(ct.get(k(r), {}).get("type", "?") for r in multi)
    type_axis = defaultdict(Counter)
    for r in multi:
        d = ct.get(k(r), {})
        type_axis[d.get("type", "?")][d.get("axis", "?")] += 1
    # subclone-like = many1 或 single_label_somatic 多結構
    subclone_like = sum(1 for r in multi if ct.get(k(r), {}).get("type") in ("many1_結構>標籤", "single_cluster") and False)  # placeholder
    out = {"N": N, "n_multi_structured": len(multi),
           "type_abs": {t: {"n": v, "pct_total": round(100 * v / N, 2), "pct_of_multi": round(100 * v / len(multi), 2)} for t, v in by_type.most_common()},
           "type_x_axis": {t: dict(a) for t, a in type_axis.items()}}
    json.dump(out, open(f"{A}/contingency_summary.json", "w"), ensure_ascii=False, indent=1)
    print("DONE"); print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
