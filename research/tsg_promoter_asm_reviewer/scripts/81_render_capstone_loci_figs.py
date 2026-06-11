#!/usr/bin/env python3
"""
81 - Render read×CpG methylation heatmaps for curated capstone loci (HCC1395),
spanning the CramersV × Δβ 2×2 (consensus / Δβ-only-TP / FP-prone / BRCA2 / excluded),
so the interactive explorer can show a figure per locus on click.

Output: genome_survey_v2/cn_confound/cross_sample/capstone_loci_figs.json
"""
import os, csv, json, base64, io
os.environ.setdefault("TMPDIR", "/big7_disk/liaoyoyo2001/tmp")
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
BAM_TP = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
          "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
WINDOW = 600
ML_HIGH, ML_LOW = 200, 50
MAXG = 45


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def hp(r):
    try:
        return str(r.get_tag("HP"))
    except KeyError:
        return None


def mcalls(r, s, e):
    o = {}
    try:
        mod = r.modified_bases
    except Exception:
        return o
    if not mod:
        return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False) if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None and s <= rf <= e:
                o[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)
    return o


def render(bam, chrom, pos, title):
    var0 = pos - 1
    s, e = max(0, var0 - WINDOW), var0 + WINDOW
    groups = {}
    for r in bam.fetch(chrom, s, e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h is None:
            continue
        m = mcalls(r, s, e)
        if len(m) >= 3:
            groups.setdefault(h, []).append(m)
    order = [g for g in ["1", "1-1", "2", "2-1"] if g in groups]
    if not order or sum(len(groups[g]) for g in order) < 6:
        return None
    rng = np.random.RandomState(42)
    for g in order:
        if len(groups[g]) > MAXG:
            sel = sorted(rng.choice(len(groups[g]), MAXG, replace=False))
            groups[g] = [groups[g][i] for i in sel]
    cols = sorted(set(c for g in order for rd in groups[g] for c in rd))
    cidx = {c: i for i, c in enumerate(cols)}
    rows, rowg = [], []
    for g in order:
        for rd in sorted(groups[g], key=lambda d: np.mean([v for v in d.values() if v in (0, 1)] or [0.5]), reverse=True):
            row = np.full(len(cols), np.nan)
            for c, v in rd.items():
                row[cidx[c]] = v
            rows.append(row); rowg.append(g)
    M = np.where(np.isnan(rows), 3, np.where(np.array(rows) == -1, 2, rows)).astype(int)
    cmap = ListedColormap(["#2563eb", "#dc2626", "#d1d5db", "#ffffff"])
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    ax.imshow(M, aspect="auto", cmap=cmap, norm=BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], 4), interpolation="nearest")
    y = 0
    lab_color = {"1": "#1d4ed8", "1-1": "#16a34a", "2": "#9333ea", "2-1": "#ca8a04"}
    for g in order:
        n = sum(1 for x in rowg if x == g)
        ax.axhline(y - 0.5, color="black", lw=0.7)
        ax.text(-0.015 * len(cols), y + n / 2, f"HP{g}", va="center", ha="right", fontsize=7,
                color=lab_color.get(g, "#333"), fontweight="bold", rotation=90)
        y += n
    ax.axhline(y - 0.5, color="black", lw=0.7)
    vcol = min(range(len(cols)), key=lambda i: abs(cols[i] - var0))
    ax.axvline(vcol, color="#f59e0b", lw=1, ls="--")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(f"{len(cols)} CpG ±{WINDOW}bp", fontsize=7)
    ax.set_title(title, fontsize=7.5)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=108, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def main():
    tp = list(csv.DictReader(open(f"{EX}/HCC1395_tp/significance_summary.csv")))
    fp = list(csv.DictReader(open(f"{EX}/HCC1395_fp/significance_summary.csv")))
    for rws in (tp, fp):
        for r in rws:
            r["_cv"] = num(r, "CramersV") or 0.0
            r["_db"] = abs(num(r, "HPMergedDelta") or 0.0)
            r["_loh"] = tb(r, "Potential_LOH")
            r["_sig"] = tb(r, "Significant")

    curated = []
    # consensus (both strong)
    cons = sorted([r for r in tp if r["_cv"] >= 0.5 and r["_db"] >= 0.15 and r["_sig"]],
                  key=lambda r: -(r["_cv"] + r["_db"]))[:35]
    for r in cons:
        curated.append((r, "consensus"))
    # Δβ-only TP (no clustering)
    do = sorted([r for r in tp if r["_cv"] < 0.1 and r["_db"] >= 0.2],
                key=lambda r: -r["_db"])[:18]
    for r in do:
        curated.append((r, "dbeta_only_TP"))
    # FP-prone
    fpp = sorted([r for r in fp if r["_db"] >= 0.15 and r["_cv"] < 0.1 and r["_loh"]],
                 key=lambda r: -r["_db"])[:18]
    for r in fpp:
        curated.append((r, "FP_prone"))
    # BRCA2 + excluded
    for r in tp:
        if r.get("Chr") == "chr13" and r.get("Pos") == "32315128":
            curated.append((r, "clustering_only_BRCA2"))
    print(f"[81] curated {len(curated)} loci")

    bam = pysam.AlignmentFile(BAM_TP, "rb")
    figs = {}
    for i, (r, cat) in enumerate(curated):
        ch, ps = r.get("Chr"), int(r.get("Pos"))
        title = f"{ch}:{ps}\nCV={r['_cv']:.2f} Δβ={num(r,'HPMergedDelta'):+.2f} {cat}"
        fb = render(bam, ch, ps, title)
        if fb:
            figs[f"{ch}:{ps}"] = dict(category=cat, cv=round(r["_cv"], 3),
                                      db=round(num(r, "HPMergedDelta"), 3), loh=r["_loh"],
                                      sig=r["_sig"], fig=fb)
        if (i + 1) % 20 == 0:
            print(f"   ...{i+1}/{len(curated)} rendered={len(figs)}")
    bam.close()
    outp = f"{CS}/capstone_loci_figs.json"
    with open(outp, "w") as f:
        json.dump(figs, f)
    print(f"[81] wrote {outp} ({os.path.getsize(outp)//1024} KB, {len(figs)} figs)")


if __name__ == "__main__":
    main()
