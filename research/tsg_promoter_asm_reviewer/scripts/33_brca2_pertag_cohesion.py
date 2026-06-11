#!/usr/bin/env python3
"""
33_brca2_pertag_cohesion.py — per-HP-tag within/between-group cohesion at BRCA2.

Answers: 「HP1-1 / HP1 / HP2 (/HP2-1/HP3) 的群內聚集度與群外聚集度」+ mixed-pattern
(which tags cluster tightly vs spread). Uses surviving MSA Level 1 (no re-run).

Method (validated M0 observed-only):
  - read x CpG methylation matrix (tumor reads), max-collapse the 2 unlabeled mod rows per (read,CpG)
  - read-read distance = mean |meth diff| over OBSERVED-shared CpGs (no imputation), min 5 shared
  - per-tag silhouette = (b-a)/max(a,b): a=mean within-tag dist, b=mean nearest-other-tag dist
    -> high silhouette = 群內緊 群外遠 (tight cohesive cluster); ~0 = 與其他 tag 混雜
  - also report per-tag mean within / between distance + MDS embedding colored by tag

Read-only input: msa_brca2/output/brca2_cloned/level1_raw_methylation_details.tsv.gz
Outputs: genome_survey_v2/brca2_pertag_cohesion.json + figures/brca2_pertag_cohesion.png
"""
import gzip, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
try:
    import sys; sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod")
    from scripts.lib.plot_setup import setup_plot_style
    setup_plot_style(base_size=11, dpi=150)
except Exception:
    pass

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
L1 = os.path.join(ROOT, "msa_brca2/output/brca2_cloned/level1_raw_methylation_details.tsv.gz")
OUTJSON = os.path.join(ROOT, "genome_survey_v2/brca2_pertag_cohesion.json")
OUTFIG = os.path.join(ROOT, "figures/brca2_pertag_cohesion.png")
THR = 0.5; MIN_SHARED = 5

# ---- load Level1, tumor only, max-collapse mods per (read,CpG) ----
f = gzip.open(L1, "rt"); hdr = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hdr)}
# read -> {cpg: max meth_call}; read -> (hp, allele)
mat = defaultdict(dict); rinfo = {}
for line in f:
    p = line.rstrip("\n").split("\t")
    if p[ix["bam_source_id"]] != "tumor":
        continue
    rid = p[ix["read_id"]]; cpg = p[ix["methyl_pos"]]
    try: mc = float(p[ix["meth_call"]])
    except ValueError: continue
    if cpg not in mat[rid] or mc > mat[rid][cpg]:
        mat[rid][cpg] = mc
    rinfo[rid] = (p[ix["haplotype_tag"]], p[ix["somatic_allele_type"]])

reads = [r for r in mat if len(mat[r]) >= MIN_SHARED]
tags = [rinfo[r][0] for r in reads]
alleles = [rinfo[r][1] for r in reads]
from collections import Counter
print("tumor reads (>=5 CpG):", len(reads), "| tag dist:", dict(Counter(tags)), "| allele:", dict(Counter(alleles)))

# ---- observed-only read-read distance ----
n = len(reads)
D = np.full((n, n), np.nan)
for i in range(n):
    D[i, i] = 0.0
    ci = mat[reads[i]]
    for j in range(i + 1, n):
        cj = mat[reads[j]]
        shared = set(ci) & set(cj)
        if len(shared) >= MIN_SHARED:
            d = np.mean([abs(ci[c] - cj[c]) for c in shared])
            D[i, j] = D[j, i] = d

# ---- per-group cohesion (silhouette on observed pairs) ----
def group_cohesion(labels):
    labs = sorted(set(labels))
    res = {}
    sil_all = {}
    for li in labs:
        idx_in = [k for k in range(n) if labels[k] == li]
        if len(idx_in) < 2:
            res[li] = {"n": len(idx_in), "within": None, "between": None, "silhouette": None}
            continue
        sils = []; wins = []; btws = []
        for k in idx_in:
            din = [D[k, m] for m in idx_in if m != k and not np.isnan(D[k, m])]
            a = np.mean(din) if din else np.nan
            # nearest other group
            best = np.inf
            for lj in labs:
                if lj == li: continue
                idx_o = [m for m in range(n) if labels[m] == lj]
                do = [D[k, m] for m in idx_o if not np.isnan(D[k, m])]
                if do: best = min(best, np.mean(do))
            b = best if best != np.inf else np.nan
            if not np.isnan(a) and not np.isnan(b):
                sils.append((b - a) / max(a, b)); wins.append(a); btws.append(b)
        res[li] = {"n": len(idx_in),
                   "within_mean": round(float(np.mean(wins)), 4) if wins else None,
                   "between_mean": round(float(np.mean(btws)), 4) if btws else None,
                   "silhouette": round(float(np.mean(sils)), 4) if sils else None,
                   "n_scored": len(sils)}
        sil_all[li] = sils
    return res, sil_all

hp_res, hp_sil = group_cohesion(tags)
al_res, al_sil = group_cohesion(alleles)
# combined HP x allele label (e.g. 1-1/alt)
combo = [f"{rinfo[r][0]}/{rinfo[r][1]}" for r in reads]
co_res, co_sil = group_cohesion(combo)

out = {
    "locus": "chr13:32315128 BRCA2/ZAR1L", "n_tumor_reads": n, "min_shared_cpg": MIN_SHARED,
    "loh_status": "nonLOH (pipeline) + CN-gain (SEQC2 CNV bed chr13:18.4-49.0Mb)",
    "per_HP_tag_cohesion": hp_res,
    "per_allele_cohesion": al_res,
    "per_HPxallele_cohesion": co_res,
    "interpretation": "silhouette>0.1 => 群內緊群外遠(內聚); ~0 => 與其他 tag 混雜; <0 => 比其他群還散",
}
json.dump(out, open(OUTJSON, "w"), indent=1, ensure_ascii=False)
print(json.dumps(out, ensure_ascii=False, indent=1))

# ---- figure: MDS colored by HP tag + per-tag silhouette bars ----
from sklearn.manifold import MDS
Dfill = np.where(np.isnan(D), np.nanmax(D), D)  # for embedding only
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0, n_init=4, normalized_stress="auto")
emb = mds.fit_transform(Dfill)
TAGCOL = {"1": "#1f77b4", "1-1": "#d62728", "2": "#2ca02c", "2-1": "#ff7f0e", "3": "#9467bd"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
for t in sorted(set(tags)):
    pts = [emb[k] for k in range(n) if tags[k] == t]
    pts = np.array(pts)
    lab = {"1": "HP1 germline", "1-1": "HP1-1 somatic", "2": "HP2", "2-1": "HP2-1", "3": "HP3"}.get(t, t)
    ax1.scatter(pts[:, 0], pts[:, 1], s=42, c=TAGCOL.get(t, "#888"), alpha=0.7, edgecolors="white", linewidths=0.5,
                label=f"{lab} (n={len(pts)})")
ax1.set_title("BRCA2 read embedding by HP tag (MDS of observed-only meth distance)")
ax1.set_xlabel("MDS-1"); ax1.set_ylabel("MDS-2"); ax1.legend(fontsize=9); ax1.grid(True, alpha=0.2)

bars = [(t, hp_res[t]["silhouette"]) for t in sorted(hp_res) if hp_res[t]["silhouette"] is not None]
ax2.bar([b[0] for b in bars], [b[1] for b in bars],
        color=[TAGCOL.get(b[0], "#888") for b in bars], alpha=0.85)
ax2.axhline(0, color="#444", lw=1); ax2.axhline(0.1, ls="--", color="#aaa", lw=1)
for i, b in enumerate(bars):
    ax2.text(i, b[1] + (0.01 if b[1] >= 0 else -0.03), f"{b[1]:.3f}", ha="center", fontsize=10)
ax2.set_ylabel("per-tag silhouette (群內緊←高 / 混雜←0)")
ax2.set_title("Per-HP-tag cohesion (within vs between-group)")
ax2.set_xticks(range(len(bars))); ax2.set_xticklabels(
    [{"1": "HP1\ngermline", "1-1": "HP1-1\nsomatic", "2": "HP2"}.get(b[0], b[0]) for b in bars], fontsize=9)
ax2.grid(True, axis="y", alpha=0.2)
plt.tight_layout(); os.makedirs(os.path.dirname(OUTFIG), exist_ok=True)
plt.savefig(OUTFIG, dpi=150, bbox_inches="tight")
print("wrote", OUTFIG)
