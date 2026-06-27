#!/usr/bin/env python3
"""全基因組 C++ 原生 phylo-v4.1 — 一次 binary run 出全位點 phylo_groups.tsv(C++ native), Python 只收集+渲染。
每 chr×TP/FP: binary(內建 compute_phylo_groups) → 收 summary.json+tsv(算對齊) → 渲染多群圖(Pool讀tsv不重算) → 刪矩陣.
輸出 phylo_cpp_wg_records.json + summary.json + figs_cpp_wg/. 用法: python3 phylo_cpp_wg.py [CHRS]."""
import os, csv, glob, json, sys, subprocess, shutil, time
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
import phylo_cpp_render as PR  # render_one (reads C++ tsv, no recompute)
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; BIN = f"{WT}/build/bin/inter_sub_mod"
TUMOR = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"; NORMAL = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF = "/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"; VCFDIR = "/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
os.environ["TMPDIR"] = "/big7_disk/liaoyoyo2001/tmp"
FIGS = f"{A}/figs_cpp_wg"; os.makedirs(FIGS, exist_ok=True); PR.FIGS = FIGS
NPROC = 24
_raw = sys.argv[1].split(",") if len(sys.argv) > 1 else [str(c) for c in range(1, 23)]
CHRS = [c if str(c).startswith("chr") else f"chr{c}" for c in _raw]


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
    base = os.path.basename(region); parts = base.split("_")
    pos = parts[1] if len(parts) >= 2 else None
    aligned = False; Vhp = Val = 0.0
    if s["coarse_ng"] >= 2:
        tsv = f"{region}/clustering/phylo_groups.tsv"
        rows = [r for r in csv.DictReader(open(tsv), delimiter="\t") if r["coarse_label"] not in ("other", "outlier")]
        if len(rows) >= 2:
            cl = [r["coarse_label"] for r in rows]; hp = [r["hp"] for r in rows]; al = [r["alt_support"] for r in rows]
            Vhp = cramerV(cl, hp); Val = cramerV(cl, al); aligned = max(Vhp, Val) >= 0.3
    return {"chrom": region.split("/")[-1].split("_")[0] if False else None, "pos": pos, "set": setlab,
            "n": s["n"], "coarse_ng": s["coarse_ng"], "fine_ng": s["fine_ng"], "n_other": s["n_other"],
            "unstable": s["unstable"], "hidden_het": s["hidden_het"], "aligned": aligned,
            "V_hp": round(Vhp, 2), "V_allele": round(Val, 2), "region": region}


def run():
    out = []; t0 = time.time(); nfig = 0
    for chrom in CHRS:
        for setlab, pref in (("TP", "filtered_snv_tp"), ("FP", "filtered_snv_fp")):
            vcf = f"{VCFDIR}/{pref}_{chrom}.vcf.gz"
            if not os.path.exists(vcf): continue
            od = f"{WT}/output/_cppwg_{chrom}_{setlab}"; shutil.rmtree(od, ignore_errors=True); os.makedirs(od, exist_ok=True)
            subprocess.run([BIN, "-t", TUMOR, "-n", NORMAL, "-r", REF, "-v", vcf, "-w", "5000", "-j", "16",
                            "--distance-metric", "BERNOULLI", "--nan-distance-strategy", "SKIP", "-o", od],
                           stdout=open(f"{od}/run.log", "w"), stderr=subprocess.STDOUT)
            regions = [os.path.dirname(os.path.dirname(sj)) for sj in glob.glob(f"{od}/**/phylo_groups_summary.json", recursive=True)]
            recs = []
            for region in regions:
                r = collect_region(region, setlab)
                if r: r["chrom"] = chrom; recs.append(r)
            out += recs
            # 渲染多群(coarse_ng>=2) — Pool 平行讀 tsv, 不重算
            multi = [r["region"] for r in recs if r["coarse_ng"] >= 2]
            if multi:
                with Pool(NPROC) as pool:
                    figs = [f for f in pool.map(PR.render_one, multi) if f and f[2]]
                nfig += len(figs)
            shutil.rmtree(od, ignore_errors=True)
            for r in out:
                r.pop("region", None)
            json.dump(out, open(f"{A}/phylo_cpp_wg_records.json", "w"))
            print(f"[{chrom}/{setlab}] loci={len(regions)} multi-figs={len(multi)} total={len(out)} figs={nfig} elapsed={int(time.time()-t0)}s", flush=True)
    TP = [r for r in out if r["set"] == "TP"]; FP = [r for r in out if r["set"] == "FP"]
    def summ(g):
        if not g: return {}
        m = [r for r in g if r["coarse_ng"] >= 2]
        return {"n": len(g), "structure_pct": round(100 * len(m) / len(g), 2),
                "aligned_pct": round(100 * sum(1 for r in m if r["aligned"]) / len(g), 2),
                "unaligned_pct": round(100 * sum(1 for r in m if not r["aligned"]) / len(g), 2),
                "no_structure_pct": round(100 * sum(1 for r in g if r["coarse_ng"] < 2) / len(g), 2),
                "fine_multi_pct": round(100 * sum(1 for r in g if r["fine_ng"] > r["coarse_ng"]) / len(g), 2),
                "other_pct": round(100 * sum(1 for r in g if r["n_other"] > 0) / len(g), 2),
                "unstable_pct": round(100 * sum(1 for r in g if r["unstable"]) / len(g), 2),
                "hidden_het_pct": round(100 * sum(1 for r in g if r["hidden_het"]) / len(g), 2),
                "ngroups_dist": dict(Counter(r["coarse_ng"] for r in g))}
    S = {"TP": summ(TP), "FP": summ(FP), "n_figs": nfig, "elapsed_s": int(time.time() - t0), "source": "C++ native phylo_groups.tsv"}
    json.dump(S, open(f"{A}/phylo_cpp_wg_summary.json", "w"), indent=2)
    print("DONE", json.dumps(S, ensure_ascii=False))


if __name__ == "__main__":
    run()
