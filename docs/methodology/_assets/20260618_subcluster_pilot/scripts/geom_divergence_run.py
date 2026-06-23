#!/usr/bin/env python3
"""[geometry-vs-null divergence 子集重跑] 讓「scipy 幾何平切看得到、但 phylo null 閘擋掉」的弱結構池可觀察。
流程: 分層抽 ~400 位點(主力=S6/S4 n>=12 幾何分歧候選 + S1/S2/S3 對照) → 建 mini-VCF(snv=rec_pos+5000)
→ 重跑 binary(內建 phylo + BERNOULLI matrix) → Python 算 geometry_ng=fcluster(0.7*max, size>=3)
→ 與 coarse_ng(null95)/fine_ng(null90) 並列, divergence = geometry_ng>=2 & coarse_ng<2。
⚠ geometry_ng = 純幾何無 null, 本就會過切單一 clone, 是「閘保守度」診斷, 非 subclone 宣稱。
保留 output 矩陣供 workstation 渲染(不刪)。輸出 geom_records.json + geom_summary.json。"""
import os, csv, glob, json, subprocess, time, shutil
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", TMPDIR="/big7_disk/liaoyoyo2001/tmp")
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from scipy.cluster.hierarchy import fcluster

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
BIN = f"{WT}/build/bin/inter_sub_mod"
TUMOR = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF = "/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR = "/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
OUTROOT = f"{WT}/output/_geomdiv"
WIN = 5000; MINSZ = 3
# 分層目標(n>=12)。主力 = S6/S4 幾何分歧候選; S1/S2/S3 = 幾何應同意的對照。
STRATA = {"S6_NO_STRUCTURE_CLEAN": {"TP": 150, "FP": 50}, "S4_WEAK_FINE_ONLY": {"TP": 90, "FP": 30},
          "S1_CONFIRMED_MULTI_ALIGNED": {"TP": 25, "FP": 15}, "S2_CONFIRMED_MULTI_UNALIGNED": {"TP": 14, "FP": 6},
          "S3_UNSTABLE_MULTI": {"TP": 14, "FP": 6}}


def even_sample(items, k):
    items = sorted(items, key=lambda r: r["n"])
    if len(items) <= k:
        return items
    step = len(items) / k
    return [items[int(i * step)] for i in range(k)]


def chr_key(c):
    c = c.replace("chr", "").replace("X", "23").replace("Y", "24")
    return int(c) if c.isdigit() else 99


def build_minivcf(loci, setlab):
    """loci: list of dict with chrom + snv_pos. 建 mini-VCF(extract 自 filtered_snv_{tp|fp}_chr.vcf.gz)。"""
    pref = "filtered_snv_tp" if setlab == "TP" else "filtered_snv_fp"
    bychr = {}
    for r in loci:
        bychr.setdefault(r["chrom"], set()).add(str(r["snv_pos"]))
    hdr = None; body = []
    for chrom, poss in bychr.items():
        vcf = f"{VCFDIR}/{pref}_{chrom}.vcf.gz"
        if not os.path.exists(vcf):
            continue
        raw = subprocess.run(["zcat", vcf], capture_output=True, text=True).stdout.splitlines()
        if hdr is None:
            hdr = [ln for ln in raw if ln.startswith("#")]
        for ln in raw:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if f[1] in poss:
                body.append(ln)
    mini = f"{A}/geom_mini_{setlab.lower()}.vcf"
    with open(mini, "w") as fo:
        fo.write("\n".join(hdr) + "\n")
        for ln in sorted(set(body), key=lambda l: (chr_key(l.split('\t')[0]), int(l.split('\t')[1]))):
            fo.write(ln + "\n")
    subprocess.run(["bgzip", "-f", mini]); subprocess.run(["tabix", "-f", "-p", "vcf", mini + ".gz"])
    return mini + ".gz", len(set(body))


def geometry_ng(D):
    """scipy 預設樹配色等價: UPGMA Z → fcluster 0.7*max(Z height) flat cut → 群數(size>=MINSZ)。RNG-free。"""
    Z, _ = CR.linkZ(D)
    t = 0.7 * float(Z[:, 2].max())
    lab = fcluster(Z, t=t, criterion="distance")
    sizes = Counter(lab)
    ng = sum(1 for s in sizes.values() if s >= MINSZ)
    return ng, t, len(set(lab))


def tsv_counts(tsv):
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    cl = [r["coarse_label"] for r in rows]; fl = [r["fine_label"] for r in rows]
    cng = len(set(l for l in cl if l not in ("other", "outlier", "-")))
    fng = len(set(l for l in fl if l not in ("other", "outlier", "-")))
    rid = [int(r["read_id"]) for r in rows]
    return cng, fng, len(rows), rid


def main():
    t0 = time.time()
    per = json.load(open(f"{A}/phylo_weak_structure_classified.json"))
    # 分層抽樣
    sample = []
    for st, tgt in STRATA.items():
        for setlab, k in tgt.items():
            pool = [r for r in per if r["wstate"] == st and r["set"] == setlab and r["n"] >= 12]
            pick = even_sample(pool, k)
            for r in pick:
                sample.append({"chrom": r["chrom"], "rec_pos": int(r["pos"]), "snv_pos": int(r["pos"]) + WIN,
                               "set": setlab, "wstate": st, "n_orig": r["n"], "coarse_orig": r["coarse_ng"],
                               "fine_orig": r["fine_ng"], "V_hp": r["V_hp"], "V_allele": r["V_allele"], "aligned": r["aligned"]})
    json.dump(sample, open(f"{A}/geom_sample_loci.json", "w"), ensure_ascii=False, indent=1)
    print(f"sampled {len(sample)} loci: " + str(dict(Counter((s['set'], s['wstate']) for s in sample))), flush=True)

    # 建 mini-VCF + 跑 binary (TP/FP 各一)
    region_map = {}  # (chrom, rec_pos) -> sample rec
    for s in sample:
        region_map[(s["chrom"], s["rec_pos"])] = s
    recs = []
    for setlab in ("TP", "FP"):
        loci = [s for s in sample if s["set"] == setlab]
        if not loci:
            continue
        mini, nv = build_minivcf(loci, setlab)
        print(f"[{setlab}] mini-VCF variants={nv} -> running binary ...", flush=True)
        od = f"{OUTROOT}_{setlab}"; shutil.rmtree(od, ignore_errors=True); os.makedirs(od, exist_ok=True)
        subprocess.run([BIN, "-t", TUMOR, "-n", NORMAL, "-r", REF, "-v", mini, "-w", str(WIN), "-j", "16",
                        "--distance-metric", "BERNOULLI", "--nan-distance-strategy", "SKIP", "-o", od],
                       stdout=open(f"{od}/run.log", "w"), stderr=subprocess.STDOUT)
        regions = [os.path.dirname(os.path.dirname(p)) for p in glob.glob(f"{od}/**/phylo_groups_summary.json", recursive=True)]
        print(f"[{setlab}] binary produced {len(regions)} regions, elapsed={int(time.time()-t0)}s", flush=True)
        for region in regions:
            base = os.path.basename(region); parts = base.split("_")
            if len(parts) < 2:
                continue
            chrom = parts[0]; start = int(parts[1])
            key = (chrom, start)
            smp = region_map.get(key)
            mpath = f"{region}/distance/BERNOULLI/matrix.csv"; tpath = f"{region}/clustering/phylo_groups.tsv"
            if not (os.path.exists(mpath) and os.path.exists(tpath)):
                continue
            try:
                dids, D = CR.loadm(mpath)
                cng, fng, nrr, rid = tsv_counts(tpath)
                idx = [dids.index(str(x)) for x in rid if str(x) in dids]  # C++-peeled complete submatrix
                if len(idx) < 2 * MINSZ:
                    continue
                sub = D[np.ix_(idx, idx)]
                gng, gcut, gflat = geometry_ng(sub)
            except Exception as e:
                print(f"  WARN {base}: {str(e)[:80]}", flush=True); continue
            rec = {"chrom": chrom, "rec_pos": start, "set": setlab,
                   "wstate": smp["wstate"] if smp else "UNMATCHED",
                   "n": nrr, "coarse_ng": cng, "fine_ng": fng, "geometry_ng": gng,
                   "geom_flat_raw": gflat, "geom_cut": round(gcut, 4),
                   "divergence": bool(gng >= 2 and cng < 2),
                   "geom_gt_coarse": bool(gng > cng),
                   "n_orig": smp["n_orig"] if smp else None, "coarse_orig": smp["coarse_orig"] if smp else None,
                   "V_hp": smp["V_hp"] if smp else None, "V_allele": smp["V_allele"] if smp else None,
                   "aligned": smp["aligned"] if smp else None, "region": region}
            recs.append(rec)
    # summary
    div = [r for r in recs if r["divergence"]]
    byw = Counter(r["wstate"] for r in recs)
    div_byw = Counter(r["wstate"] for r in div)
    # geometry_ng vs coarse_ng cross-tab
    xtab = Counter((min(r["coarse_ng"], 3), min(r["geometry_ng"], 5)) for r in recs)
    summary = {
        "source": "subset rerun (geom_divergence_run.py) on phylo_weak_structure_classified sample",
        "n_sampled": len(sample), "n_regions_rerun": len(recs),
        "n_divergence": len(div), "divergence_pct_of_rerun": round(100 * len(div) / len(recs), 2) if recs else 0,
        "divergence_by_wstate": dict(div_byw), "regions_by_wstate": dict(byw),
        "geom_ge2_total": sum(1 for r in recs if r["geometry_ng"] >= 2),
        "geom_gt_coarse_total": sum(1 for r in recs if r["geom_gt_coarse"]),
        "coarse_vs_geometry_xtab": {f"coarse{c}_geom{g}": n for (c, g), n in sorted(xtab.items())},
        "caveat": "geometry_ng = pure scipy 0.7*max flat-cut, NO null, EXPECTED to over-split single clones (gap/silhouette put 1 clone at k=3-4 under epiallele correlation). Diagnostic of gate conservativeness, NOT a subclone/structure claim. Single-sample HCC1395 exploratory.",
        "elapsed_s": int(time.time() - t0),
    }
    for r in recs:
        r.pop("region", None)  # keep records light; region dirs retained on disk for rendering
    json.dump(recs, open(f"{A}/geom_records.json", "w"), ensure_ascii=False)
    json.dump(summary, open(f"{A}/geom_summary.json", "w"), ensure_ascii=False, indent=2)
    print("DONE " + json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
