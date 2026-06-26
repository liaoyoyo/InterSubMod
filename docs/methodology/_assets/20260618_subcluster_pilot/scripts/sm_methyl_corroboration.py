"""[L4 ISM 甲基 corroboration — 單變量, genotype-anchored] 全基因組(可用既有 ISM 輸出的)結構區。
對每個 CN-clean 結構區(span<=8kb, 有 covering ISM window): re-pileup region somatic sSNV → per-read 基因型群
(genotype-anchored, 非循環) → join methylation.csv(by read_name) → 兩大 population 間 per-CpG MWU + BH-FDR
→ 甲基是否區分這些遺傳定義群。cis-ASM proxy: 同 CpG 在 HP1 vs HP2 是否也區分(= germline cis 非 subclone)。
輸出 報告夾/data/sm_methyl_corroboration.json。🔴 甲基 = corroborate 非 detect。
"""
import json
import os
import csv
from collections import Counter, defaultdict
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
RPT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report"
OUTROOT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_phylo_wg_full"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
DBETA_MIN = 0.2


def bh(pvals):
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * n / (rank + 1))
        q[i] = prev
    return q


import bisect as _bisect
_CW = None


def _load_cw():
    global _CW
    cc = json.load(open(f"{A}/cis_candidates_resolved.json"))
    cw = defaultdict(list)
    for c in cc:
        try:
            cw[c["chrom"]].append((int(c["pos"]), c.get("region_dir")))
        except Exception:
            pass
    for k in cw:
        cw[k].sort()
    _CW = cw


def find_window(chrom, region_start, region_end, som_pos):
    """用 ISM cis-candidate window([cand,cand+10000]) 覆蓋 region 者。回 (meth_path, reads_path) 或 None。"""
    if _CW is None:
        _load_cw()
    arr = _CW.get(chrom, [])
    poss = [x[0] for x in arr]
    i = _bisect.bisect_right(poss, region_start) - 1
    for j in range(max(0, i - 3), min(len(arr), i + 4)):
        cand, rd = arr[j]
        if rd and cand <= region_start and cand + 10000 >= region_end:
            mp = f"{rd}/methylation/methylation.csv"
            rp = f"{rd}/reads/reads.tsv"
            if os.path.exists(mp) and os.path.exists(rp):
                return mp, rp
    return None


def per_read_geno(tb, chrom, som):
    """som=[(pos,ref,alt)]. 回 {read_name: 'R/A string over som'} (covering 全部)。"""
    calls = defaultdict(dict)
    for pos, ref, alt in som:
        for col in tb.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, stepper="samtools"):
            for pr in col.pileups:
                if pr.is_del or pr.query_position is None or pr.alignment.mapping_quality < MAPQ_MIN:
                    continue
                b = pr.alignment.query_sequence[pr.query_position].upper()
                if b == ref:
                    calls[pr.alignment.query_name][pos] = "R"
                elif b == alt:
                    calls[pr.alignment.query_name][pos] = "A"
    poss = [s[0] for s in som]
    geno = {}
    for rn, c in calls.items():
        if all(p in c for p in poss):
            geno[rn] = "".join(c[p] for p in poss)
    return geno


def main():
    reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]
    cen = json.load(open(f"{A}/sm_linkage_genomewide.json"))["census"]
    tb = pysam.AlignmentFile(TBAM, "rb")
    doable = [r for r in reg if r["cn"] in ("loh", "neutral")
              and r["tree_shape"] not in ("no_confirmed_structure", "inconsistent")
              and r["span"] <= 8000 and (r.get("n_populations") or 0) >= 2]

    n_tested = n_corrob = n_window_missing = n_nopop = 0
    sig_cpg_counts = []
    cis_explained = 0
    examples = []
    for r in doable:
        chrom = r["chrom"]
        som = []
        for pos in range(r["start"], r["end"] + 1):
            k = f"{chrom}:{pos}"
            if k in cen and cen[k].get("somatic"):
                som.append((pos, cen[k]["ref"], cen[k]["alt"]))
        if len(som) < 2:
            continue
        som_pos = [s[0] for s in som]
        w = find_window(chrom, r["start"], r["end"], som_pos)
        if not w:
            n_window_missing += 1
            continue
        mp, rp = w
        # load reads.tsv: read_name -> local_id, hp, is_tumor
        rn2id = {}
        id2hp = {}
        id2tum = {}
        for row in csv.DictReader(open(rp), delimiter="\t"):
            rn2id[row["read_name"]] = row["read_id"]
            id2hp[row["read_id"]] = row["hp"]
            id2tum[row["read_id"]] = row["is_tumor"]
        # load methylation: read_id -> {cpg: beta}
        ml = open(mp).read().strip().split("\n")
        cpgs = ml[0].split(",")[1:]
        meth = {}
        for line in ml[1:]:
            parts = line.split(",")
            meth[parts[0]] = parts[1:]
        # per-read genotype (tumor)
        geno = per_read_geno(tb, chrom, som)
        # group reads by geno (tumor only), via read_name->local_id->meth
        pop_reads = defaultdict(list)  # geno -> [local_id]
        for rn, g in geno.items():
            lid = rn2id.get(rn)
            if lid is None or id2tum.get(lid) != "1" or lid not in meth:
                continue
            pop_reads[g].append(lid)
        pops = [(g, ids) for g, ids in pop_reads.items() if len(ids) >= 3]
        if len(pops) < 2:
            n_nopop += 1
            continue
        pops.sort(key=lambda x: -len(x[1]))
        gA, idsA = pops[0]
        gB, idsB = pops[1]
        n_tested += 1
        # per-CpG MWU between A and B
        pvals = []
        cpg_idx = []
        dbetas = []
        for j, cp in enumerate(cpgs):
            a = [float(meth[i][j]) for i in idsA if meth[i][j] not in ("", "NA")]
            b = [float(meth[i][j]) for i in idsB if meth[i][j] not in ("", "NA")]
            if len(a) >= 3 and len(b) >= 3:
                try:
                    _, p = mannwhitneyu(a, b, alternative="two-sided")
                except ValueError:
                    continue
                pvals.append(p); cpg_idx.append(j); dbetas.append(abs(np.mean(a) - np.mean(b)))
        if not pvals:
            continue
        q = bh(np.array(pvals))
        sig = [(cpg_idx[k], dbetas[k]) for k in range(len(pvals)) if q[k] < 0.05 and dbetas[k] >= DBETA_MIN]
        sig_cpg_counts.append(len(sig))
        if sig:
            n_corrob += 1
            # cis proxy: 這些 sig CpG 是否也在 HP1 vs HP2 區分 (germline cis)
            hp1 = [i for i in (idsA + idsB) if id2hp.get(i) in ("1", "1-1")]
            hp2 = [i for i in (idsA + idsB) if id2hp.get(i) in ("2", "2-1")]
            cis_hit = 0
            if len(hp1) >= 3 and len(hp2) >= 3:
                for j, _ in sig:
                    a = [float(meth[i][j]) for i in hp1 if meth[i][j] not in ("", "NA")]
                    b = [float(meth[i][j]) for i in hp2 if meth[i][j] not in ("", "NA")]
                    if len(a) >= 3 and len(b) >= 3:
                        try:
                            _, p = mannwhitneyu(a, b, alternative="two-sided")
                            if p < 0.05 and abs(np.mean(a) - np.mean(b)) >= DBETA_MIN:
                                cis_hit += 1
                        except ValueError:
                            pass
            if cis_hit >= max(1, len(sig) // 2):
                cis_explained += 1
            if len(examples) < 20:
                examples.append({"region": r["region"], "shape": r["tree_shape"], "n_sig_cpg": len(sig),
                                 "popA": f"{gA}({len(idsA)})", "popB": f"{gB}({len(idsB)})", "cis_hit": cis_hit})
    tb.close()
    out = {
        "scope": "CN-clean(loh/neutral) structured span<=8kb regions with covering ISM window; existing ISM output",
        "n_doable_candidate": len(doable), "n_window_missing": n_window_missing,
        "n_insufficient_populations": n_nopop, "n_tested": n_tested,
        "n_methyl_corroborated(>=1 sig CpG, |db|>=0.2 q<0.05)": n_corrob,
        "corroboration_rate": round(n_corrob / n_tested, 3) if n_tested else None,
        "median_sig_cpg_per_tested": float(np.median(sig_cpg_counts)) if sig_cpg_counts else None,
        "cis_ASM_explained_of_corroborated": cis_explained,
        "cis_explained_rate": round(cis_explained / n_corrob, 3) if n_corrob else None,
        "examples": examples,
        "verdict": "甲基 corroborate genotype-anchored 群的比例 = 輔助程度; cis_explained = 該 corroboration 其實是 "
                   "germline cis-ASM(非 subclone-specific)的比例 → 扣掉才是真 subclone-specific 甲基",
        "redline": "甲基 = corroborate 非 detect; 全基因組大區/缺 window 區需 ISM 重抽整合(後補)",
    }
    json.dump(out, open(f"{RPT}/data/sm_methyl_corroboration.json", "w"), ensure_ascii=False, indent=1)
    print(f"L4 methyl: doable={len(doable)} window_missing={n_window_missing} nopop={n_nopop} tested={n_tested}")
    print(f"  corroborated={n_corrob} ({out['corroboration_rate']}) median sig CpG={out['median_sig_cpg_per_tested']}")
    print(f"  cis-ASM explained of corroborated={cis_explained} ({out['cis_explained_rate']})")


if __name__ == "__main__":
    main()
