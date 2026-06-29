#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLP 突變交叉驗證(2026-06-29):我們的 sSNV 區域 vs COSMIC Cell Lines Project 同 cell line 已知突變。
⚠ 口徑:CLP cell-line 突變多為 "Variant of unknown origin"(無配對 normal,混 germline/somatic)→
   這是「位點存在性 corroborate」非「somatic 一致性」。短讀 WGS/WES vs ONT ClairS,平台/caller 不同 → 預期部分重疊。
指標(per sample): 區域含 ≥1 CLP 突變的比例;癌基因區的 CLP 突變命中;reverse(CLP 突變落我們區域比例)。
用法: python3 clp_crossval.py → 每樣本寫 {workdir}/clp_crossval.json + 印彙總。
"""
import os, gzip, csv, json
from bisect import bisect_left
from collections import defaultdict

ANN = os.environ.get("SM_ANNOT_DIR", "/big7_disk/liaoyoyo2001/gene_annotation")
MUT = f"{ANN}/cosmic_v104/CellLinesProject_GenomeScreensMutant_v104_GRCh38.tsv.gz"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
HCC_FROZEN = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data"
# 我們樣本 → (COSMIC 名, 我們的 data 夾)
SAMPLES = {
    "HCC1395": ("HCC1395", HCC_FROZEN),
    "COLO829": ("COLO-829", f"{MSROOT}/COLO829"),
    "H1437": ("NCI-H1437", f"{MSROOT}/H1437"),
    "H2009": ("NCI-H2009", f"{MSROOT}/H2009"),
    "HCC1937": ("HCC1937", f"{MSROOT}/HCC1937"),
    "HCC1954": ("HCC1954", f"{MSROOT}/HCC1954"),
    "HCC1395_DORADO": ("HCC1395", f"{MSROOT}/HCC1395_DORADO"),
}


def load_clp_mutations():
    """COSMIC 名 → {chrom: sorted [pos]}。只收我們需要的 cell line。"""
    want = set(c for c, _ in SAMPLES.values())
    pos = defaultdict(lambda: defaultdict(list))
    with gzip.open(MUT, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cn = r.get("SAMPLE_NAME", "").strip()
            if cn not in want:
                continue
            chrom = r.get("CHROMOSOME", "").strip()
            try:
                p = int(r.get("GENOME_START", ""))
            except ValueError:
                continue
            if chrom:
                pos[cn][f"chr{chrom}"].append(p)
    for cn in pos:
        for c in pos[cn]:
            pos[cn][c].sort()
    return pos


def region_has_mut(sorted_pos, s, e):
    i = bisect_left(sorted_pos, s)
    return i < len(sorted_pos) and sorted_pos[i] <= e


def main():
    print("loading CLP mutations...", flush=True)
    clp = load_clp_mutations()
    summary = {}
    for ours, (cos, d) in SAMPLES.items():
        tp = os.path.join(d, "topology_per_region.json")
        if not os.path.exists(tp):
            print(f"  {ours}: 無 topology,跳過"); continue
        det = json.load(open(tp))["detail"]
        gp = os.path.join(d, "region_gene_annotation.json")
        gene = json.load(open(gp)).get("regions", {}) if os.path.exists(gp) else {}
        mp = clp.get(cos, {})
        n_mut = sum(len(v) for v in mp.values())
        hit = cancer_hit = cancer_tot = 0
        per_region = {}
        used_mut = set()
        for r in det:
            c = r["chrom"]; s = r["start"]; e = r.get("end") or int(r["region"].split("-")[-1])
            sp = mp.get(c, [])
            h = region_has_mut(sp, s, e) if sp else False
            # 計落在此區的 CLP 突變數(用於 reverse)
            i = bisect_left(sp, s); cnt = 0
            while i < len(sp) and sp[i] <= e:
                used_mut.add((c, sp[i])); cnt += 1; i += 1
            if h:
                hit += 1
                per_region[r["region"]] = cnt
            g = gene.get(r["region"], {})
            if g.get("cancer_genes"):
                cancer_tot += 1
                if h:
                    cancer_hit += 1
        n_region = len(det)
        mut_in_region = len(used_mut)
        summary[ours] = {
            "cosmic_name": cos, "n_clp_mut": n_mut, "n_region": n_region,
            "regions_with_clp_mut": hit, "region_hit_pct": round(100 * hit / max(n_region, 1), 1),
            "cancer_gene_regions": cancer_tot, "cancer_regions_with_clp_mut": cancer_hit,
            "cancer_hit_pct": round(100 * cancer_hit / max(cancer_tot, 1), 1),
            "clp_mut_in_our_regions": mut_in_region, "clp_mut_in_region_pct": round(100 * mut_in_region / max(n_mut, 1), 1),
        }
        json.dump({"summary": summary[ours], "per_region_mut_count": per_region},
                  open(os.path.join(d, "clp_crossval.json"), "w"), ensure_ascii=False)
        s = summary[ours]
        print(f"  {ours:<15} CLP突變={s['n_clp_mut']:>5} 區={n_region:>5} | 區含突變={s['regions_with_clp_mut']:>4}({s['region_hit_pct']}%) "
              f"癌基因區命中={s['cancer_regions_with_clp_mut']}/{cancer_tot}({s['cancer_hit_pct']}%) | CLP突變落我們區={s['clp_mut_in_region_pct']}%", flush=True)
    json.dump(summary, open(os.path.join(HCC_FROZEN, "clp_crossval_summary.json"), "w"), ensure_ascii=False, indent=1)
    print("CROSSVAL DONE")


if __name__ == "__main__":
    main()
