#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part C — stratified sampling of "verdict-ignored significant PERMANOVA" loci from the
Phase-2 CORRECTED ISM scan, + 2-panel inspection figures, for the human 4-way labelling
workstation (真結構 / germline-only / 離散度假象 / sample-confound).

Target set = loci the OLD verdict ignored despite HP-aligned PERMANOVA structure:
  legacy VerificationClass ∈ {Weak, Noise} (ignored)  AND
  LabelHPPermanovaValid & LabelHPPermanovaP <= 0.05    (HP-aligned PERMANOVA significant)

4 strata (user's axes, mapped to corrected columns):
  high_label_F   : highest LabelHPPermanovaF (strongest HP separation the verdict missed)
  near_boundary  : HP_AUC_All closest to the 0.70 structure threshold (ambiguous)
  sample_only    : cluster structure present (ClusterPermanovaP<=0.05) but HP-misaligned
                   (HP_AUC_All <= 0.60) — looks like sample/methyl split, not haplotype
  germline_asm   : GermlineAsmDbeta_Sig=true (germline ASM)
  candidate_som  : (SubcloneDbeta/SomaticResidual)_Sig & not GermlineAsmDbeta_Sig

All numbers read from the Phase-2 significance_summary (no fabrication). N per stratum = TOPN.
Figure render adapts develop gen_inspection_figures.py to any chromosome + NHD distance.

Usage: python3 select_and_render.py
"""
import csv, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P2 = "/big7_disk/liaoyoyo2001/ism_phase2_scan"
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)
FIGDIR = os.path.join(OUTDIR, "figs")
os.makedirs(FIGDIR, exist_ok=True)
TOPN = 8

# ---------------------------------------------------------------- load summaries
def num(r, k):
    try:
        v = float(r.get(k, "") or "nan"); return None if v != v else v
    except (ValueError, TypeError):
        return None
def tb(r, k):
    return str(r.get(k, "")).lower() == "true"

ROWS = []
for cls in ("tp", "fp"):
    p = f"{P2}/HCC1395_{cls}/significance_summary.csv"
    for r in csv.DictReader(open(p)):
        r["_cls"] = cls
        ROWS.append(r)
print(f"[select] loaded {len(ROWS)} loci (TP+FP) from Phase-2 corrected scan")

# ---------------------------------------------------------------- target set
def in_target(r):
    legacy = r.get("VerificationClass_Legacy", "")
    return (legacy in ("Weak", "Noise")
            and tb(r, "LabelHPPermanovaValid")
            and (num(r, "LabelHPPermanovaP") or 1) <= 0.05)

target = [r for r in ROWS if in_target(r)]
print(f"[select] target (verdict-ignored significant PERMANOVA) = {len(target)}")

# ---------------------------------------------------------------- strata
def pick(cands, keyfn, reverse=True, n=TOPN):
    cc = [r for r in cands if keyfn(r) is not None]
    cc.sort(key=keyfn, reverse=reverse)
    return cc[:n]

strata = {}
strata["high_label_F"] = pick(target, lambda r: num(r, "LabelHPPermanovaF"))
strata["near_boundary"] = pick(target, lambda r: -abs((num(r, "HP_AUC_All") or 0) - 0.70), reverse=True)
strata["sample_only"] = pick(
    [r for r in target if (num(r, "ClusterPermanovaP") or 1) <= 0.05 and (num(r, "HP_AUC_All") or 1) <= 0.60],
    lambda r: num(r, "ClusterPermanovaF"))
strata["germline_asm"] = pick([r for r in target if tb(r, "GermlineAsmDbeta_Sig")],
                              lambda r: abs(num(r, "GermlineAsmDbeta") or 0))
strata["candidate_som"] = pick(
    [r for r in target if (tb(r, "SubcloneDbeta_HP1_Sig") or tb(r, "SubcloneDbeta_HP2_Sig")
                           or tb(r, "SomaticResidualDbeta_Sig")) and not tb(r, "GermlineAsmDbeta_Sig")],
    lambda r: max(abs(num(r, "SubcloneDbeta_HP1") or 0), abs(num(r, "SubcloneDbeta_HP2") or 0),
                  abs(num(r, "SomaticResidualDbeta") or 0)))
for k, v in strata.items():
    print(f"[select]   {k}: {len(v)}")

# ---------------------------------------------------------------- figure render (NHD, any chr)
HP_ORDER = {"1": 0, "HP1": 0, "1-1": 1, "2": 2, "HP2": 2, "2-1": 3}
HP_NAME = {0: "HP1(germ)", 1: "HP1-1(carrier)", 2: "HP2(germ)", 3: "HP2-1(carrier)", 4: "unlabeled"}
HP_COLOR = {0: "#1f5fd0", 1: "#7fb0f0", 2: "#d02020", 3: "#f08080", 4: "#bbbbbb"}
TN_COLOR = {True: "#e8820e", False: "#2ca02c"}

def region_inner(cls, chrom, pos):
    g = glob.glob(f"{P2}/HCC1395_{cls}/filtered_snv_{cls}/{chrom}/{chrom}_{pos}/{chrom}_*")
    return g[0] if g else None

def load_reads(inner):
    meta = {}
    with open(os.path.join(inner, "reads", "reads.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            meta[row["read_id"]] = (row.get("hp", ""), row.get("is_tumor", "1") == "1")
    return meta

def load_mat(path):
    ids, rows = [], []
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            ids.append(p[0])
            rows.append([np.nan if v in ("NA", "") else float(v) for v in p[1:]])
    return ids, np.array(rows, dtype=float)

def sb(ax, ids, order, meta, kind):
    arr = np.zeros((len(order), 1, 3))
    for row, i in enumerate(order):
        hp, ist = meta.get(ids[i], ("", True))
        c = HP_COLOR[HP_ORDER.get(hp, 4)] if kind == "hp" else TN_COLOR[ist]
        arr[row, 0] = matplotlib.colors.to_rgb(c)
    ax.imshow(arr, aspect="auto"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("HP" if kind == "hp" else "T/N", fontsize=7)

def render(cls, chrom, pos, label, r):
    inner = region_inner(cls, chrom, pos)
    if not inner:
        return None
    meta = load_reads(inner)
    mids, M = load_mat(os.path.join(inner, "methylation", "methylation.csv"))
    if M.size == 0:
        return None
    dpath = os.path.join(inner, "distance", "NHD", "matrix.csv")
    dids, D = load_mat(dpath) if os.path.exists(dpath) else (None, None)
    order = sorted(range(len(mids)), key=lambda i: (HP_ORDER.get(meta.get(mids[i], ("", 0))[0], 4),
                                                    0 if meta.get(mids[i], ("", True))[1] else 1))
    Ms = M[order, :]
    Ds = D[np.ix_(np.array(order), np.array(order))] if (D is not None and dids == mids) else None
    n = len(mids)
    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(1, 7, width_ratios=[0.08, 0.08, 2.2, 0.18, 0.08, 0.08, 2.2], wspace=0.04)
    sb(fig.add_subplot(gs[0, 0]), mids, order, meta, "hp")
    sb(fig.add_subplot(gs[0, 1]), mids, order, meta, "tn")
    axm = fig.add_subplot(gs[0, 2])
    cmap = plt.cm.RdYlBu_r.copy(); cmap.set_bad("#dddddd")
    im = axm.imshow(Ms, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    axm.set_title(f"Methylation β (read×CpG), HP-grouped  [n={n}, CpG={M.shape[1]}]", fontsize=9)
    axm.set_xlabel("CpG", fontsize=8); axm.set_ylabel("reads (HP-grouped)", fontsize=8)
    grp = [HP_ORDER.get(meta.get(mids[i], ("", 0))[0], 4) for i in order]
    for rr in range(1, n):
        if grp[rr] != grp[rr - 1]:
            axm.axhline(rr - 0.5, color="black", lw=0.8)
    plt.colorbar(im, ax=axm, fraction=0.025, pad=0.01, label="β")
    if Ds is not None:
        sb(fig.add_subplot(gs[0, 4]), mids, order, meta, "hp")
        sb(fig.add_subplot(gs[0, 5]), mids, order, meta, "tn")
        axd = fig.add_subplot(gs[0, 6])
        imd = axd.imshow(Ds, aspect="auto", cmap="magma_r", interpolation="nearest")
        axd.set_title("Read×Read distance (NHD), same HP order", fontsize=9)
        axd.set_xlabel("reads (HP-grouped)", fontsize=8)
        for rr in range(1, n):
            if grp[rr] != grp[rr - 1]:
                axd.axhline(rr - 0.5, color="white", lw=0.6); axd.axvline(rr - 0.5, color="white", lw=0.6)
        plt.colorbar(imd, ax=axd, fraction=0.025, pad=0.01, label="dist")
    g = lambda k: r.get(k, "NA") or "NA"
    title = (f"[{label}] {chrom}:{pos}  NumReads={g('NumReads')} NReadsValid={g('NReadsValid')}  "
             f"HP-AUC(all/tum/norm)={g('HP_AUC_All')}/{g('HP_AUC_Tumor')}/{g('HP_AUC_Normal')}  "
             f"CramersV={g('CramersV')} LabelF={g('LabelHPPermanovaF')} p={g('LabelHPPermanovaP')}\n"
             f"VC={g('VerificationClass')} (legacy {g('VerificationClass_Legacy')})  Strength={g('StrengthScore')}/{g('StrengthGrade')}  "
             f"germΔβ={g('GermlineAsmDbeta')}({g('GermlineAsmDbeta_Sig')}) "
             f"subΔβ={g('SubcloneDbeta_HP1')}/{g('SubcloneDbeta_HP2')} somΔβ={g('SomaticResidualDbeta')}({g('SomaticResidualDbeta_Sig')})")
    fig.suptitle(title, fontsize=8, y=1.0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=HP_COLOR[k]) for k in range(5)] + \
              [plt.Rectangle((0, 0), 1, 1, color=TN_COLOR[True]), plt.Rectangle((0, 0), 1, 1, color=TN_COLOR[False])]
    fig.legend(handles, [HP_NAME[k] for k in range(5)] + ["tumor", "normal"],
               loc="lower center", ncol=7, fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.03))
    out = os.path.join(FIGDIR, f"{label}_{chrom}_{pos}.png")
    fig.savefig(out, dpi=100, bbox_inches="tight"); plt.close(fig)
    return os.path.basename(out)

# ---------------------------------------------------------------- build cases.json + render
cases = []
seen = set()
for label, items in strata.items():
    for r in items:
        chrom, pos, cls = r["Chr"], r["Pos"], r["_cls"]
        key = f"{chrom}:{pos}"
        if key in seen:
            continue
        seen.add(key)
        fig = render(cls, chrom, pos, label, r)
        cases.append({
            "id": f"{chrom}_{pos}", "stratum": label, "chrom": chrom, "pos": pos, "cls": cls.upper(),
            "fig": fig,
            "metrics": {k: r.get(k, "") for k in (
                "NumReads", "NReadsValid", "HP_AUC_All", "HP_AUC_Tumor", "HP_AUC_Normal",
                "CramersV", "LabelHPPermanovaF", "LabelHPPermanovaP", "ClusterPermanovaF", "ClusterPermanovaP",
                "VerificationClass", "VerificationClass_Legacy", "StrengthScore", "StrengthGrade",
                "GermlineAsmDbeta", "GermlineAsmDbeta_Sig", "SubcloneDbeta_HP1", "SubcloneDbeta_HP2",
                "SomaticResidualDbeta", "SomaticResidualDbeta_Sig", "Potential_LOH", "Significant")},
        })
        print(f"  rendered {label} {key} -> {fig}")

json.dump({"source": f"{P2} (Phase-2 corrected, ±5000 SKIP, develop binary)",
           "target_def": "legacy VC in {Weak,Noise} AND LabelHPPermanovaValid AND LabelHPPermanovaP<=0.05",
           "n_target": len(target), "n_cases": len(cases), "strata_counts": {k: len(v) for k, v in strata.items()},
           "cases": cases},
          open(os.path.join(OUTDIR, "cases.json"), "w"), ensure_ascii=False, indent=1)
print(f"[select] wrote cases.json ({len(cases)} cases) + {len(os.listdir(FIGDIR))} figs")
