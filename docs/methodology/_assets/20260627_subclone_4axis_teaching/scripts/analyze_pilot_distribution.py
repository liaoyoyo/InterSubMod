#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
觀察 pilot cis-control 分布 + 驗證「tumor genotype 軸 vs normal germline-HP 軸」是否對齊。
corr(tumor_Δβ, normal_Δβ)≈0 有兩解:(a)無 cis-confound (b)兩軸不對齊→normal-HP 無 control 力。
本分析用 tumor reads 的 HP tag 組成判斷 popA/popB 是否沿 germline HP 分(對齊) 或同 HP 內(不對齊)。
輸出:scatter PNG + axis_alignment.json。compute batch(§13.0,跑完另批讀+寫報告)。
"""
import json, os
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
ASSETS = os.path.normpath(os.path.join(HERE, ".."))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3

pilot = json.load(open(f"{DATA}/pilot_cis_control.json"))
pairs = pilot["all_cpg_pairs"]
det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}

# ---- (1) 散點圖: signed tumor_Δβ vs signed normal_Δβ + 四象限 ----
tx = np.array([p["normal_dbeta"] for p in pairs])  # x = germline-ASM 軸
ty = np.array([p["tumor_dbeta"] for p in pairs])    # y = subclone-候選 軸
fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter(tx, ty, s=8, alpha=0.35, c="#3366cc", edgecolors="none")
for v in (0.2, -0.2):
    ax.axhline(v, color="#cc3333", ls="--", lw=0.8)
    ax.axvline(v, color="#888", ls="--", lw=0.8)
ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
ax.set_xlabel("normal HP1-HP2 dbeta  (germline-ASM baseline)")
ax.set_ylabel("tumor clusterA-clusterB dbeta  (subclone candidate)")
corr_val = np.corrcoef(tx, ty)[0, 1]
ax.set_title(f"Pilot cis-control: tumor(genotype-axis) vs normal(HP-axis)\nn={len(pairs)} CpG  corr={corr_val:.3f}")
# 象限標註(tumor 高訊號為主)
ax.text(0.55, 0.92, "tumor-hi & normal-hi\n=germline-ASM", transform=ax.transAxes, fontsize=8, color="#aa6600", ha="center")
ax.text(0.12, 0.92, "tumor-hi & normal-lo\n=subclone? (if axes aligned)", transform=ax.transAxes, fontsize=8, color="#117733", ha="center")
ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
fig.tight_layout()
png = f"{ASSETS}/pilot_cis_control_scatter.png"
fig.savefig(png, dpi=130); plt.close(fig)

# ---- (2) 軸對齊: 重算 tumor genotype cluster, 取每群 reads 的 HP tag 組成 ----
def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        p = f"{VCFD}/filtered_snv_{src}_{chrom}.vcf.gz"
        if not os.path.exists(p): continue
        try:
            tb = pysam.TabixFile(p)
            for ln in tb.fetch(chrom, max(0, s - 1), e + 1):
                f = ln.split("\t"); pos = int(f[1])
                if pos not in m: m[pos] = (f[3].upper(), f[4].strip().upper())
        except Exception: pass
    return m

def geno_str(a, som):
    rq = {rr: q for q, rr in a.get_aligned_pairs(matches_only=True)}
    seq = a.query_sequence; g = []
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None or seq is None: g.append("-"); continue
        b = seq[q].upper()
        g.append("A" if b == alt else ("R" if b == ref else "-"))
    return "".join(g)

def hp_norm(a):
    # tumor BAM HP tag = str("1-1"/"1"/"2"/"3"); normal = int. 取 "-" 前的數字為 haplotype。
    if not a.has_tag("HP"): return 0
    v = a.get_tag("HP")
    try: return int(str(v).split("-")[0])
    except Exception: return 0

tb = pysam.AlignmentFile(TBAM, "rb")
align = []
for region in {p["region"] for p in pairs}:
    r = det.get(region)
    if not r: continue
    chrom, s = r["chrom"], r["start"]; e = int(region.split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: continue
    pop_hp = defaultdict(Counter)   # genotype -> HP Counter
    pop_n = Counter()
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        g = geno_str(a, som)
        if "-" in g: continue
        hp = hp_norm(a)
        pop_hp[g][hp] += 1; pop_n[g] += 1
    pops = [g for g, n in pop_n.most_common() if n >= MINREAD]
    if len(pops) < 2: continue
    A, B = pops[0], pops[1]
    def hp_frac(c):
        tot = sum(c.values()) or 1
        h1 = c.get(1, 0); h2 = c.get(2, 0); h3 = c.get(3, 0); h0 = c.get(0, 0)
        return {"H1": round(h1 / tot, 2), "H2": round(h2 / tot, 2), "H3": round(h3 / tot, 2), "unphased": round(h0 / tot, 2), "n": tot}
    fa, fb = hp_frac(pop_hp[A]), hp_frac(pop_hp[B])
    # 對齊判定: popA 與 popB 是否各自被不同 HP 主導
    def dom(f): return "H1" if f["H1"] >= 0.6 else ("H2" if f["H2"] >= 0.6 else ("unphased" if f["unphased"] >= 0.6 else "mixed"))
    da, db = dom(fa), dom(fb)
    if da in ("H1", "H2") and db in ("H1", "H2") and da != db:
        verdict = "ALIGNED(normal-HP 可 control)"
    elif da == db and da in ("H1", "H2"):
        verdict = "SAME-HP(subclone 在同一 germline HP 內→normal-HP 不對齊)"
    else:
        verdict = "MIXED/UNPHASED(無法用 HP 對齊)"
    align.append({"region": region, "popA": A, "popA_hp": fa, "popB": B, "popB_hp": fb,
                  "domA": da, "domB": db, "alignment": verdict})

vc = Counter(a["alignment"].split("(")[0] for a in align)
out = {"scatter_png": os.path.relpath(png, ASSETS),
       "corr_tumor_vs_normal_dbeta": round(float(corr_val), 3),
       "n_regions_checked": len(align),
       "alignment_dist": dict(vc),
       "regions": align,
       "interpretation": "corr≈0 + SAME-HP/MIXED 主導 → normal-HP 基線與 tumor genotype 軸不對齊,故 corr≈0 非『無 cis-confound』而是『比錯軸』;subclone-specificity 對 genotype 軸仍 UNDETERMINED。"}
json.dump(out, open(f"{DATA}/pilot_axis_alignment.json", "w"), ensure_ascii=False, indent=1)
print("AXIS-ALIGN DONE", json.dumps({"corr": round(float(corr_val),3), "n_regions": len(align), "dist": dict(vc)}, ensure_ascii=False))
for a in align:
    print(f"  {a['region']} popA={a['popA']}({a['domA']},n{a['popA_hp']['n']}) popB={a['popB']}({a['domB']},n{a['popB_hp']['n']}) → {a['alignment']}")
print("scatter →", png)
