#!/usr/bin/env python3
"""[全 34736 位點 C++ phylo + tree-aware 群色] 渲染 ALL(含單群,驗分錯) + 保留 matrix(下次改色只 render 免重跑 binary)。
每 chr×TP/FP: binary → 收 summary/tsv → render 全部 region(新 H.cluster_tree_color) → **不刪 od**。
輸出 phylo_cpp_wg_full_records.json + summary.json + figs_cpp_wg_full/。matrices 留 output/_phylo_wg_full/。"""
import os, csv, glob, json, subprocess, shutil, time
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", TMPDIR="/big7_disk/liaoyoyo2001/tmp")
import numpy as np
from collections import Counter
import sys
sys.path.insert(0, os.path.dirname(__file__))
import phylo_cpp_render as PR        # 已接 H.cluster_tree_color
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
BIN = f"{WT}/build/bin/inter_sub_mod"
TUMOR = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF = "/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR = "/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
KEEP = f"{WT}/output/_phylo_wg_full"                      # 保留 matrix(persistent)
FIGS = f"{A}/figs_cpp_wg_full"; PR.FIGS = FIGS
NPROC = 24
CHRS = [f"chr{c}" for c in range(1, 23)]


def cramerV(a, b):
    ca = sorted(set(a)); cb = sorted(set(b)); ia = {c: i for i, c in enumerate(ca)}; ib = {c: i for i, c in enumerate(cb)}
    tab = np.zeros((len(ca), len(cb)))
    for x, y in zip(a, b): tab[ia[x], ib[y]] += 1
    nn = tab.sum()
    if nn == 0 or min(len(ca), len(cb)) < 2: return 0.0
    row = tab.sum(1, keepdims=True); col = tab.sum(0, keepdims=True); exp = row * col / nn
    chi2 = np.nansum((tab - exp) ** 2 / np.where(exp > 0, exp, 1))
    return float(np.sqrt(chi2 / (nn * (min(len(ca), len(cb)) - 1))))


def collect_region(region, setlab):
    sj = f"{region}/clustering/phylo_groups_summary.json"
    if not os.path.exists(sj): return None
    s = json.load(open(sj))
    parts = os.path.basename(region).split("_")
    pos = parts[1] if len(parts) >= 2 else None
    aligned = False; Vhp = Val = 0.0
    if s["coarse_ng"] >= 2:
        rows = [r for r in csv.DictReader(open(f"{region}/clustering/phylo_groups.tsv"), delimiter="\t") if r["coarse_label"] not in ("other", "outlier")]
        if len(rows) >= 2:
            cl = [r["coarse_label"] for r in rows]; hp = [r["hp"] for r in rows]; al = [r["alt_support"] for r in rows]
            Vhp = cramerV(cl, hp); Val = cramerV(cl, al); aligned = max(Vhp, Val) >= 0.3
    return {"chrom": None, "pos": pos, "set": setlab, "n": s["n"], "coarse_ng": s["coarse_ng"], "fine_ng": s["fine_ng"],
            "n_other": s["n_other"], "unstable": s["unstable"], "hidden_het": s["hidden_het"], "aligned": aligned,
            "V_hp": round(Vhp, 2), "V_allele": round(Val, 2)}


def summ(g):
    if not g: return {}
    m = [r for r in g if r["coarse_ng"] >= 2]
    return {"n": len(g), "structure_pct": round(100 * len(m) / len(g), 2),
            "aligned_pct": round(100 * sum(1 for r in m if r["aligned"]) / len(g), 2),
            "no_structure_pct": round(100 * sum(1 for r in g if r["coarse_ng"] < 2) / len(g), 2),
            "fine_multi_pct": round(100 * sum(1 for r in g if r["fine_ng"] > r["coarse_ng"]) / len(g), 2),
            "unstable_pct": round(100 * sum(1 for r in g if r["unstable"]) / len(g), 2),
            "ngroups_dist": dict(Counter(r["coarse_ng"] for r in g))}


def run():
    os.makedirs(KEEP, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
    out = []; t0 = time.time(); nfig = 0
    for chrom in CHRS:
        for setlab, pref in (("TP", "filtered_snv_tp"), ("FP", "filtered_snv_fp")):
            vcf = f"{VCFDIR}/{pref}_{chrom}.vcf.gz"
            if not os.path.exists(vcf): continue
            od = f"{KEEP}/{chrom}_{setlab}"; shutil.rmtree(od, ignore_errors=True); os.makedirs(od, exist_ok=True)
            subprocess.run([BIN, "-t", TUMOR, "-n", NORMAL, "-r", REF, "-v", vcf, "-w", "5000", "-j", "16",
                            "--distance-metric", "BERNOULLI", "--nan-distance-strategy", "SKIP", "-o", od],
                           stdout=open(f"{od}/run.log", "w"), stderr=subprocess.STDOUT)
            regions = [os.path.dirname(os.path.dirname(sj)) for sj in glob.glob(f"{od}/**/phylo_groups_summary.json", recursive=True)]
            for region in regions:
                r = collect_region(region, setlab)
                if r: r["chrom"] = chrom; out.append(r)
            with Pool(NPROC) as pool:                              # render ALL(單群+多群), matrices 留著
                figs = [f for f in pool.map(PR.render_one, regions) if f and f[2]]
            nfig += len(figs)
            json.dump(out, open(f"{A}/phylo_cpp_wg_full_records.json", "w"))
            print(f"[{chrom}/{setlab}] loci={len(regions)} figs={len(figs)} total={len(out)} nfig={nfig} elapsed={int(time.time()-t0)}s", flush=True)
    S = {"TP": summ([r for r in out if r["set"] == "TP"]), "FP": summ([r for r in out if r["set"] == "FP"]),
         "n_records": len(out), "n_figs": nfig, "rendered": "ALL (single+multi, 驗分錯)", "matrices_kept": KEEP,
         "palette": "tree-aware H.cluster_tree_color", "elapsed_s": int(time.time() - t0)}
    json.dump(S, open(f"{A}/phylo_cpp_wg_full_summary.json", "w"), indent=2, ensure_ascii=False)
    print("DONE " + json.dumps(S, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    run()
