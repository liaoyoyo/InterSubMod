#!/usr/bin/env python3
"""
34_control_loci_cohesion_cistest.py — per-tag cohesion + normal-anchored cis-test
across BRCA2 + 9 control loci (TP hi/lo cohesion, TP-LOH, germline-het, FP).

(1) Per-tag cohesion: silhouette (within vs between) for tumor HP1/HP1-1/HP2/HP2-1.
(2) Normal-anchored cis-test (HP1 axis): per shared-CpG beta of
       A = normal HP1 (pre-mutation baseline)
       B = tumor  HP1 (germline-haplotype reads w/o somatic allele)
       C = tumor  HP1-1 (somatic-allele reads)
    Deltas: d_somatic = C-B (within-tumor ASM); d_cis = C-A; d_drift = B-A.
    cis-specific candidate  <=>  |d_cis| sig AND |d_cis| >> |d_drift|.
    (Caveat: assumes tumor/normal HP1 = same parental haplotype; flagged + checked.)

Read-only input: genome_survey_v2/control_msa_out/.../level1_raw_methylation_details.tsv.gz
                 + msa_brca2/.../level1...  (BRCA2 may already be in control set)
Outputs: genome_survey_v2/control_cohesion_cistest.json + figures/control_cohesion_cistest.png
"""
import gzip, json, os, glob
import numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import sys; sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod")
    from scripts.lib.plot_setup import setup_plot_style
    setup_plot_style(base_size=10, dpi=150)
except Exception:
    pass

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
L1 = glob.glob(os.path.join(ROOT, "genome_survey_v2/control_msa_out/*/level1_raw_methylation_details.tsv.gz"))[0]
MAN = json.load(open(os.path.join(ROOT, "genome_survey_v2/control_loci_manifest.json")))
man = {m["locus"].split(":")[1]: m for m in MAN}
THR = 0.5; MIN_N = 3; MIN_PAIR = 4

# load: per somatic_pos -> per (bam,hp,allele) -> per cpg -> {read: max meth}
D = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
f = gzip.open(L1, "rt"); h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
for line in f:
    p = line.rstrip("\n").split("\t")
    try: mc = float(p[ix["meth_call"]])
    except ValueError: continue
    spos = p[ix["somatic_pos"]]; cpg = p[ix["methyl_pos"]]; rid = p[ix["read_id"]]
    bam = p[ix["bam_source_id"]]; hp = p[ix["haplotype_tag"]]; al = p[ix["somatic_allele_type"]]
    cell = D[spos][(bam, hp, al)][cpg]
    if rid not in cell or mc > cell[rid]:
        cell[rid] = mc

def read_matrix(spos, bam):
    """return {read_id: {cpg: beta_call}}, {read_id:(hp,allele)} for one bam (max over allele/hp dup)."""
    mat = defaultdict(dict); info = {}
    for (b, hp, al), cpgd in D[spos].items():
        if b != bam: continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items():
                if cpg not in mat[rid] or v > mat[rid][cpg]:
                    mat[rid][cpg] = v
                info[rid] = (hp, al)
    return mat, info

def beta_per_cpg(spos, bam, hp_filter):
    """mean meth-frac per cpg over reads whose hp in hp_filter (set), min MIN_N reads."""
    agg = defaultdict(list)
    for (b, hp, al), cpgd in D[spos].items():
        if b != bam or hp not in hp_filter: continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items():
                agg[(cpg, rid)] = v  # dedup read per cpg
    # collapse to per-cpg fraction
    bycpg = defaultdict(list)
    for (cpg, rid), v in agg.items():
        bycpg[cpg].append(v)
    return {c: (np.mean([1 if x >= THR else 0 for x in vs]), len(vs)) for c, vs in bycpg.items()}

def silhouette_by_tag(spos):
    mat, info = read_matrix(spos, "tumor")
    reads = [r for r in mat if len(mat[r]) >= 5]
    tags = {r: info[r][0] for r in reads}
    n = len(reads)
    if n < 6: return {}
    dist = np.full((n, n), np.nan)
    for i in range(n):
        dist[i, i] = 0
        for j in range(i + 1, n):
            sh = set(mat[reads[i]]) & set(mat[reads[j]])
            if len(sh) >= 5:
                d = np.mean([abs(mat[reads[i]][c] - mat[reads[j]][c]) for c in sh])
                dist[i, j] = dist[j, i] = d
    res = {}
    labs = sorted(set(tags.values()))
    for li in labs:
        idx_in = [k for k in range(n) if tags[reads[k]] == li]
        if len(idx_in) < 3: continue
        sils = []
        for k in idx_in:
            a = np.nanmean([dist[k, m] for m in idx_in if m != k])
            bs = [np.nanmean([dist[k, m] for m in range(n) if tags[reads[m]] == lj])
                  for lj in labs if lj != li]
            bs = [b for b in bs if not np.isnan(b)]
            if bs and not np.isnan(a):
                sils.append((min(bs) - a) / max(a, min(bs)))
        if sils: res[li] = {"n": len(idx_in), "silhouette": round(float(np.mean(sils)), 3)}
    return res

def cis_test(spos):
    """HP1-axis: A=normal HP1, B=tumor HP1, C=tumor HP1-1 per-CpG paired."""
    A = beta_per_cpg(spos, "normal", {"1"})
    B = beta_per_cpg(spos, "tumor", {"1"})
    C = beta_per_cpg(spos, "tumor", {"1-1"})
    shared = set(A) & set(B) & set(C)
    shared = [c for c in shared if A[c][1] >= MIN_N and B[c][1] >= MIN_N and C[c][1] >= MIN_N]
    if len(shared) < MIN_PAIR:
        return {"n_shared_cpg": len(shared), "testable": False}
    a = np.array([A[c][0] for c in shared]); b = np.array([B[c][0] for c in shared]); c_ = np.array([C[c][0] for c in shared])
    def wp(x, y):
        try: return round(float(stats.wilcoxon(x, y).pvalue), 4)
        except Exception: return None
    return {
        "n_shared_cpg": len(shared), "testable": True,
        "d_somatic_C_minus_B": round(float(np.mean(c_ - b)), 4), "p_somatic": wp(c_, b),
        "d_cis_C_minus_normalA": round(float(np.mean(c_ - a)), 4), "p_cis": wp(c_, a),
        "d_drift_B_minus_normalA": round(float(np.mean(b - a)), 4), "p_drift": wp(b, a),
        "mean_beta": {"normalHP1": round(float(a.mean()), 3), "tumorHP1": round(float(b.mean()), 3), "tumorHP1-1": round(float(c_.mean()), 3)},
    }

results = {}
for spos in D:
    m = man.get(spos, {"tag": "?", "cls": "?", "loh": None})
    results[spos] = {"locus": f"chr?:{spos}", "tag": m["tag"], "cls": m["cls"], "loh": m.get("loh"),
                     "ari_blind_allele": m.get("ari_blind_allele"),
                     "cohesion": silhouette_by_tag(spos), "cis": cis_test(spos)}
    # fix locus chrom
    results[spos]["locus"] = m.get("locus", f"?:{spos}")

json.dump(results, open(os.path.join(ROOT, "genome_survey_v2/control_cohesion_cistest.json"), "w"), indent=1, ensure_ascii=False)

# ---- print comparison table ----
print(f"{'tag':10s}{'locus':18s}{'cls':6s}{'loh':4s}{'HP1':>6}{'HP1-1':>7}{'HP2':>6} | {'d_som':>7}{'d_cis':>7}{'d_drift':>8}{'p_cis':>8} | verdict")
order = sorted(results, key=lambda s: (results[s]["cls"], -(results[s].get("ari_blind_allele") or -9)))
for s in order:
    r = results[s]; co = r["cohesion"]; ci = r["cis"]
    hp1 = co.get("1", {}).get("silhouette"); hp11 = co.get("1-1", {}).get("silhouette"); hp2 = co.get("2", {}).get("silhouette")
    def fmt(x): return f"{x:.3f}" if isinstance(x, float) else "  -  "
    if ci.get("testable"):
        dcis = ci["d_cis_C_minus_normalA"]; ddr = ci["d_drift_B_minus_normalA"]; pcis = ci["p_cis"]
        cis_spec = (pcis is not None and pcis < 0.05 and abs(dcis) > 1.8 * abs(ddr) + 0.02)
        verdict = "CIS-cand" if cis_spec else ("drift" if abs(ddr) >= abs(dcis) else "weak/ns")
        print(f"{r['tag']:10s}{r['locus']:18s}{r['cls']:6s}{str(r['loh'])[:3]:4s}{fmt(hp1):>6}{fmt(hp11):>7}{fmt(hp2):>6} | "
              f"{ci['d_somatic_C_minus_B']:>7.3f}{dcis:>7.3f}{ddr:>8.3f}{(pcis if pcis is not None else -1):>8} | {verdict}")
    else:
        print(f"{r['tag']:10s}{r['locus']:18s}{r['cls']:6s}{str(r['loh'])[:3]:4s}{fmt(hp1):>6}{fmt(hp11):>7}{fmt(hp2):>6} | cis untestable (shared CpG={ci['n_shared_cpg']}, need normal HP1 + tumor HP1 + HP1-1)")

# ---- figure: cohesion (HP1-1 vs HP1) + cis vs drift, colored by cls ----
CLSCOL = {"TP": "#d62728", "het-NULL": "#1f77b4", "FP": "#7f7f7f"}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
for s in results:
    r = results[s]; co = r["cohesion"]; ci = r["cis"]
    col = CLSCOL.get(r["cls"], "#999")
    hp11 = co.get("1-1", {}).get("silhouette"); hp1 = co.get("1", {}).get("silhouette")
    if hp11 is not None and hp1 is not None:
        ax1.scatter(hp1, hp11, c=col, s=90, alpha=0.8, edgecolors="k", linewidths=0.6)
        ax1.annotate(r["tag"], (hp1, hp11), fontsize=7.5, xytext=(3, 3), textcoords="offset points")
    if ci.get("testable"):
        ax2.scatter(abs(ci["d_drift_B_minus_normalA"]), abs(ci["d_cis_C_minus_normalA"]),
                    c=col, s=90, alpha=0.8, edgecolors="k", linewidths=0.6)
        ax2.annotate(r["tag"], (abs(ci["d_drift_B_minus_normalA"]), abs(ci["d_cis_C_minus_normalA"])),
                     fontsize=7.5, xytext=(3, 3), textcoords="offset points")
ax1.plot([-.2, .5], [-.2, .5], ls="--", color="#aaa", lw=1)
ax1.set_xlabel("HP1 germline silhouette"); ax1.set_ylabel("HP1-1 somatic silhouette")
ax1.set_title("Per-tag cohesion: somatic(HP1-1) vs germline(HP1)\n(above diagonal = somatic reads more cohesive)")
ax1.grid(True, alpha=0.2)
lim = max(0.25, *[abs(results[s]["cis"].get("d_cis_C_minus_normalA", 0)) for s in results if results[s]["cis"].get("testable")] or [0.25])
ax2.plot([0, lim], [0, lim], ls="--", color="#aaa", lw=1)
ax2.set_xlabel("|d_drift|  tumorHP1 − normalHP1 (germline drift)")
ax2.set_ylabel("|d_cis|  tumorHP1-1 − normalHP1 (somatic vs baseline)")
ax2.set_title("Normal-anchored cis test\n(above diagonal = somatic-allele reads diverge MORE than drift = cis-candidate)")
ax2.grid(True, alpha=0.2)
from matplotlib.lines import Line2D
ax2.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=v, label=k, markersize=9) for k, v in CLSCOL.items()], fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(ROOT, "figures/control_cohesion_cistest.png"), dpi=150, bbox_inches="tight")
print("\nwrote figures/control_cohesion_cistest.png")
