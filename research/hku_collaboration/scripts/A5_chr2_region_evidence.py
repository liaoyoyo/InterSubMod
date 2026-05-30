#!/usr/bin/env python3
"""
A5 / Figure F7 — chr2:18,068,480-18,100,683 (32 kb) IGV-style region evidence
=============================================================================
Maps to PI presentation slide #1 example "往前推三點" (15 SNV positions over 32 kb LOH window).

Demonstrates:
  - 15 candidate positions in 32 kb chr2 window, all inside SEQC2 LOH
    (chr2:16,146,119-22,100,000, ~16 Mb)
  - 4 SEQC2 TP ∩ 5 ClairS-TO ssrs PASS → 4 true somatics + 1 FP candidate
    (chr2:18,086,020)
  - HP tag distribution (HP2-related 86.7% over the 32 kb window) → mono-haplotype collapse
  - chr2:18,086,020 (DP=77, A=49%, G=51%) looks 50:50 but is an FP because LOH eliminated one
    haplotype; the apparent 50:50 comes from ~50% normal contamination.

Figure layout (3 rows + side callout):
  Row 1 : chr2 ideogram + LOH bar (16.1M-22.1M) + yellow triangle at 18,086,020 + zoom box.
  Row 2 : 15-site lollipop: x = position, y = DP. Color = TP+ClairS / ClairS-only / no-call.
  Row 3 : Per-position HP tag stacked bars (HP1 / HP2 / HP1-1 / HP2-1 / HP33 / no-HP).

Outputs:
  figures/F7_chr2_region_evidence.png
  data/F7_14sites_status.tsv
"""

from __future__ import annotations

import gzip
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle

# -------- CJK font injection (per project rule) --------
for f in ["/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
          "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
          "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if os.path.exists(f):
        font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = ["DejaVu Sans", "Droid Sans Fallback", "WenQuanYi Micro Hei",
                                "Noto Sans CJK SC", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42

# -------- Paths --------
BAM = "/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
TP_VCF = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
SSRS_VCF = ("/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/"
            "ONT_5khz_simplex_5mCG_5hmCG/pileup/HCC1395_methyl_PASS.vcf")
LOH_BED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"

OUT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration")
OUT_FIG = OUT_ROOT / "figures" / "F7_chr2_region_evidence.png"
OUT_TSV = OUT_ROOT / "data" / "F7_14sites_status.tsv"
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# -------- Window & sites --------
CHROM = "chr2"
WIN_START = 18_068_480
WIN_END = 18_100_683
FP_POS = 18_086_020

# 15 candidate positions from plan v4.
PLAN_SITES = [
    18068480, 18072546, 18073987, 18076409, 18086020,
    18090947, 18091584, 18093414, 18094114, 18096033,
    18096269, 18096341, 18099697, 18100544, 18100683,
]

# hg38 chr2 length
CHR2_LEN = 242_193_529

# Color palette (PI report design)
COL_BOTH    = "#7B1FA2"   # TP+ClairS → true somatic (purple)
COL_CLAIRS  = "#F9A825"   # ClairS-only → FP candidate (yellow/amber)
COL_TPONLY  = "#9E9E9E"   # TP-only, not called (grey)
COL_FP      = "#D32F2F"   # highlight chr2:18086020 (red)

# HP integer keys (after normalization from string tags):
#   "1"   -> 1   (HP1)
#   "2"   -> 2   (HP2)
#   "1-1" -> 11  (HP1-1 somatic)
#   "2-1" -> 21  (HP2-1 somatic)
#   "3-3" / "33" -> 33 (HP33 ambiguous)
#   None / "none" -> 0 (untagged)
HP_COLORS = {
    1:   "#90CAF9",   # HP1   light blue
    2:   "#1565C0",   # HP2   dark  blue
    11:  "#FFAB91",   # HP1-1 light red
    21:  "#C62828",   # HP2-1 dark red
    33:  "#9E9E9E",   # HP33  grey (ambiguous)
    0:   "#FFFFFF",   # no-HP white
}
HP_LABELS = {1: "HP1", 2: "HP2", 11: "HP1-1", 21: "HP2-1", 33: "HP33", 0: "no-HP"}
HP_ORDER  = [1, 2, 11, 21, 33, 0]


def normalize_hp_tag(raw) -> int:
    """Map BAM HP tag (string or int) to canonical integer key."""
    if raw is None:
        return 0
    s = str(raw).strip().lower()
    if s in ("", "none", "."):
        return 0
    if s == "1":
        return 1
    if s == "2":
        return 2
    if s in ("1-1", "11"):
        return 11
    if s in ("2-1", "21"):
        return 21
    if s in ("3-3", "33"):
        return 33
    # Fallback: try int parse
    try:
        i = int(s)
        return i if i in (1, 2, 11, 21, 33, 0) else 33
    except ValueError:
        return 33  # any unrecognized non-empty tag → ambiguous bucket


# ---------------------------------------------------------------- helpers
def load_tp_positions(vcf_path: str, chrom: str, start: int, end: int) -> set[int]:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    mode = "rt" if vcf_path.endswith(".gz") else "r"
    out: set[int] = set()
    with opener(vcf_path, mode) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom and start <= int(f[1]) <= end:
                out.add(int(f[1]))
    return out


def load_ssrs_calls(vcf_path: str, chrom: str, start: int, end: int) -> dict[int, dict]:
    """Return {pos: {ref, alt, dp, af, au, cu, gu, tu, filter}}."""
    out: dict[int, dict] = {}
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            pos = int(f[1])
            if pos < start or pos > end:
                continue
            fmt = f[8].split(":")
            sam = f[9].split(":")
            d = dict(zip(fmt, sam))
            out[pos] = {
                "ref": f[3], "alt": f[4], "filter": f[6],
                "dp": int(d.get("DP", 0)),
                "af": float(d.get("AF", 0.0)),
                "au": int(d.get("AU", 0)), "cu": int(d.get("CU", 0)),
                "gu": int(d.get("GU", 0)), "tu": int(d.get("TU", 0)),
            }
    return out


def load_loh_intervals(bed_path: str, chrom: str) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    with open(bed_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] != chrom:
                continue
            if f[3].strip().lower() == "loh":
                out.append((int(f[1]), int(f[2])))
    return out


def site_in_loh(pos: int, intervals: list[tuple[int, int]]) -> tuple[bool, tuple[int, int] | None]:
    for s, e in intervals:
        if s <= pos <= e:
            return True, (s, e)
    return False, None


def fetch_hp_counts_at_pos(bam: pysam.AlignmentFile, chrom: str, pos: int,
                            half_window: int = 50) -> Counter:
    """Count HP tags among reads overlapping `pos`."""
    c: Counter = Counter()
    for r in bam.fetch(chrom, max(0, pos - half_window), pos + half_window):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        if r.reference_start > pos or (r.reference_end or 0) < pos:
            continue
        try:
            raw = r.get_tag("HP")
        except KeyError:
            raw = None
        c[normalize_hp_tag(raw)] += 1
    return c


def fetch_region_hp_counts(bam: pysam.AlignmentFile, chrom: str, start: int, end: int) -> Counter:
    """Count unique reads in region by normalized HP tag."""
    c: Counter = Counter()
    seen: set[str] = set()
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        if r.query_name in seen:
            continue
        seen.add(r.query_name)
        try:
            raw = r.get_tag("HP")
        except KeyError:
            raw = None
        c[normalize_hp_tag(raw)] += 1
    return c


def dominant_hp_label(counts: Counter) -> str:
    if not counts:
        return "none"
    total = sum(counts.values())
    if total == 0:
        return "none"
    top, n = counts.most_common(1)[0]
    return f"{HP_LABELS.get(top, str(top))}={n}/{total} ({n / total * 100:.0f}%)"


# ---------------------------------------------------------------- data prep
print("[A5] Loading TP VCF...", flush=True)
tp_pos = load_tp_positions(TP_VCF, CHROM, WIN_START, WIN_END)
print(f"      {len(tp_pos)} TP sites in window: {sorted(tp_pos)}", flush=True)

print("[A5] Loading ClairS-TO ssrs VCF...", flush=True)
ssrs = load_ssrs_calls(SSRS_VCF, CHROM, WIN_START, WIN_END)
print(f"      {len(ssrs)} ssrs PASS sites in window: {sorted(ssrs)}", flush=True)

print("[A5] Loading LOH BED...", flush=True)
loh_chr2 = load_loh_intervals(LOH_BED, CHROM)
print(f"      {len(loh_chr2)} LOH intervals on chr2", flush=True)

print("[A5] Opening BAM + extracting HP counts per position...", flush=True)
bam = pysam.AlignmentFile(BAM, "rb")
region_hp_counts = fetch_region_hp_counts(bam, CHROM, WIN_START, WIN_END)
print(f"      32 kb region HP profile: {dict(region_hp_counts)}", flush=True)
per_site_hp: dict[int, Counter] = {}

# Union all positions: plan list + ssrs calls + TP calls (deduped)
all_positions = sorted(set(PLAN_SITES) | set(tp_pos) | set(ssrs.keys()))
for p in all_positions:
    per_site_hp[p] = fetch_hp_counts_at_pos(bam, CHROM, p)
bam.close()

# ---------------------------------------------------------------- TSV output
@dataclass
class Site:
    pos: int
    in_tp: bool
    in_ssrs: bool
    in_loh: bool
    vcf_dp: int
    vcf_af: float
    hp_dominant: str


sites: list[Site] = []
for p in all_positions:
    in_tp = p in tp_pos
    in_ssrs = p in ssrs
    in_loh, _ = site_in_loh(p, loh_chr2)
    dp = ssrs[p]["dp"] if in_ssrs else 0
    af = ssrs[p]["af"] if in_ssrs else 0.0
    sites.append(Site(p, in_tp, in_ssrs, in_loh, dp, af, dominant_hp_label(per_site_hp[p])))

with open(OUT_TSV, "w") as fh:
    fh.write("chrom\tpos\tin_tp_vcf\tin_clairs_ssrs\tin_loh\tvcf_dp\tvcf_af\thp_dominant\n")
    for s in sites:
        fh.write(f"{CHROM}\t{s.pos}\t{int(s.in_tp)}\t{int(s.in_ssrs)}\t{int(s.in_loh)}\t"
                 f"{s.vcf_dp}\t{s.vcf_af:.4f}\t{s.hp_dominant}\n")
print(f"[A5] TSV written: {OUT_TSV} ({len(sites)} rows)", flush=True)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(16, 9), dpi=150)
# Reserve right ~26% for callout box
gs = fig.add_gridspec(nrows=3, ncols=2,
                       height_ratios=[1.0, 2.6, 1.9],
                       width_ratios=[3.6, 1.0],
                       hspace=0.75, wspace=0.10,
                       left=0.06, right=0.985, top=0.93, bottom=0.07)

ax_chr = fig.add_subplot(gs[0, 0])
ax_lol = fig.add_subplot(gs[1, 0])
ax_hp = fig.add_subplot(gs[2, 0])
ax_callout = fig.add_subplot(gs[:, 1])

# ----- Row 1: chr2 ideogram -----
ax_chr.set_xlim(0, CHR2_LEN)
ax_chr.set_ylim(-1.4, 1.4)
ax_chr.add_patch(Rectangle((0, -0.35), CHR2_LEN, 0.7,
                            facecolor="#ECEFF1", edgecolor="#37474F", linewidth=0.9))
# All LOH intervals (light) and target LOH (highlight)
for s, e in loh_chr2:
    ax_chr.add_patch(Rectangle((s, -0.35), e - s, 0.7,
                                facecolor="#FFCDD2", edgecolor="none", alpha=0.55))
# Highlight the LOH containing the window
hit, target_loh = site_in_loh(FP_POS, loh_chr2)
if hit:
    ts, te = target_loh
    ax_chr.add_patch(Rectangle((ts, -0.35), te - ts, 0.7,
                                facecolor="#E53935", edgecolor="#B71C1C",
                                linewidth=1.0, alpha=0.85))
# Yellow triangle for FP candidate
ax_chr.plot([FP_POS], [0.85], marker="v", color="#FBC02D",
            markersize=14, markeredgecolor="#5D4037", markeredgewidth=0.8, zorder=5)
ax_chr.annotate(f"chr2:{FP_POS:,}",
                xy=(FP_POS, 0.85), xytext=(FP_POS, 1.25),
                ha="center", fontsize=8.5, color="#5D4037")
# Zoom box (window bounds projected on chr2)
zoom_w = WIN_END - WIN_START
ax_chr.add_patch(Rectangle((WIN_START - zoom_w * 30, -0.55),
                            zoom_w * 61, 1.1,
                            facecolor="none", edgecolor="#000000", linewidth=1.4))
ax_chr.text(WIN_START - zoom_w * 30, -0.95,
            f"zoom: chr2:{WIN_START / 1e6:.3f}-{WIN_END / 1e6:.3f} Mb (~{zoom_w / 1000:.1f} kb)",
            fontsize=7.5, color="#000000", va="top")

ax_chr.set_yticks([])
ax_chr.set_xticks([0, 50e6, 100e6, 150e6, 200e6, CHR2_LEN])
ax_chr.set_xticklabels(["0", "50M", "100M", "150M", "200M", f"{CHR2_LEN / 1e6:.0f}M"], fontsize=9)
ax_chr.set_xlabel("chr2 position (hg38)", fontsize=10)
ax_chr.set_title(f"Row 1 — chr2 ideogram with all LOH segments (pink); target LOH "
                  f"{target_loh[0] / 1e6:.1f}-{target_loh[1] / 1e6:.1f} Mb (red, "
                  f"{(target_loh[1] - target_loh[0]) / 1e6:.1f} Mb) contains 32 kb zoom window",
                  fontsize=10, loc="left", pad=4)
for s in ("top", "right", "left"):
    ax_chr.spines[s].set_visible(False)

# ----- Row 2: lollipop chart -----
ax_lol.set_xlim(WIN_START - 200, WIN_END + 200)

# Background LOH shading across whole window (all sites are in LOH)
ax_lol.axvspan(WIN_START - 200, WIN_END + 200,
                facecolor="#FFEBEE", alpha=0.45, zorder=0)

# Compute max DP for y-axis
max_dp = max([s.vcf_dp for s in sites if s.vcf_dp > 0] + [80])
y_top = max_dp * 1.55

for s in sites:
    if s.in_tp and s.in_ssrs:
        col = COL_BOTH
    elif s.in_ssrs and not s.in_tp:
        col = COL_FP if s.pos == FP_POS else COL_CLAIRS
    elif s.in_tp and not s.in_ssrs:
        col = COL_TPONLY
    else:
        col = "#BDBDBD"
    height = s.vcf_dp if s.vcf_dp > 0 else y_top * 0.05  # short stub for no-call sites
    ax_lol.vlines(s.pos, 0, height, color=col,
                  linewidth=2.6 if s.pos == FP_POS else 1.6,
                  alpha=0.95, zorder=2)
    ax_lol.scatter([s.pos], [height], s=160 if s.pos == FP_POS else 75,
                    color=col, edgecolor="black", linewidth=0.7, zorder=3)
    # Annotate DP / AF on the lollipop head (only for called sites)
    if s.vcf_dp > 0:
        ax_lol.text(s.pos, height + y_top * 0.025,
                    f"DP={s.vcf_dp}\nAF={s.vcf_af:.2f}",
                    ha="center", va="bottom",
                    fontsize=6.8, color="#212121")

# FP red callout box around chr2:18086020
ax_lol.annotate("", xy=(FP_POS, ssrs[FP_POS]["dp"] + y_top * 0.16),
                xytext=(FP_POS + 6800, y_top * 0.88),
                arrowprops=dict(arrowstyle="->", color=COL_FP, lw=1.4))
ax_lol.text(FP_POS + 6800, y_top * 0.92,
             f"FP candidate\nchr2:{FP_POS:,}\nDP=77, A=49%, G=51%",
             color=COL_FP, fontsize=8.4, ha="left", va="bottom",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                       edgecolor=COL_FP, linewidth=1.0))

# Legend for Row 2
leg_handles = [
    mpatches.Patch(color=COL_BOTH, label="TP+ClairS PASS (true somatic)"),
    mpatches.Patch(color=COL_CLAIRS, label="ClairS-only PASS (FP candidate)"),
    mpatches.Patch(color=COL_FP, label="FP focus: chr2:18,086,020"),
    mpatches.Patch(color=COL_TPONLY, label="TP-only (not called by ClairS)"),
    mpatches.Patch(color="#BDBDBD", label="No call (plan-listed stub)"),
    mpatches.Patch(facecolor="#FFEBEE", edgecolor="#E53935",
                   label="LOH region chr2:16.1M-22.1M"),
]
ax_lol.legend(handles=leg_handles, loc="lower center", fontsize=7.4,
              ncol=3, framealpha=0.95, handlelength=1.4, columnspacing=1.2,
              bbox_to_anchor=(0.5, -0.18))

ax_lol.set_ylim(0, y_top)
ax_lol.set_xlim(WIN_START - 800, WIN_END + 800)
ax_lol.set_ylabel("VCF DP (read depth)", fontsize=10)
ax_lol.set_title(f"Row 2 — 15-site lollipop: 4 TP∩ClairS (purple), 1 ClairS-only FP "
                  f"(red @ {FP_POS:,}), 10 plan-listed sites with no call (stubs)",
                  fontsize=10, loc="left", pad=4)
ax_lol.grid(True, axis="y", linestyle=":", alpha=0.4)
for s in ("top", "right"):
    ax_lol.spines[s].set_visible(False)

# Make X-axis ticks at site positions for alignment with Row 3
ax_lol.set_xticks([s.pos for s in sites])
ax_lol.set_xticklabels([])  # blank here; labelled in Row 3

# ----- Row 3: HP stacked bars -----
ax_hp.set_xlim(WIN_START - 800, WIN_END + 800)
bar_width = 1500  # ~1.5 kb wide per bar (visible on 32 kb axis)

for s in sites:
    counts = per_site_hp[s.pos]
    total = sum(counts.values())
    bottom = 0.0
    if total == 0:
        ax_hp.add_patch(Rectangle((s.pos - bar_width / 2, 0), bar_width, 1.0,
                                   facecolor="#FAFAFA", edgecolor="#BDBDBD", linewidth=0.6))
        continue
    for hp in HP_ORDER:
        v = counts.get(hp, 0)
        if v == 0:
            continue
        frac = v / total
        ax_hp.add_patch(Rectangle((s.pos - bar_width / 2, bottom), bar_width, frac,
                                   facecolor=HP_COLORS[hp], edgecolor="black", linewidth=0.3))
        bottom += frac
    # Read total on top
    ax_hp.text(s.pos, 1.04, f"n={total}", ha="center", fontsize=6.2, color="#37474F")
    # Red border for FP focus
    if s.pos == FP_POS:
        ax_hp.add_patch(Rectangle((s.pos - bar_width / 2 - 80, -0.04),
                                   bar_width + 160, 1.10,
                                   facecolor="none", edgecolor=COL_FP, linewidth=1.8))

# Caveat annotation at FP position
hp_region_total = sum(region_hp_counts.values())
hp2_related = region_hp_counts.get(2, 0) + region_hp_counts.get(21, 0)
hp2_frac = hp2_related / hp_region_total * 100 if hp_region_total else 0.0
ax_hp.annotate(f"HP2-related {hp2_frac:.1f}%\n(region-wide haplotype collapse)",
                xy=(FP_POS, 1.05), xytext=(FP_POS + 6800, 1.30),
                ha="left", fontsize=8.0, color=COL_FP,
                arrowprops=dict(arrowstyle="->", color=COL_FP, lw=1.2),
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=COL_FP, linewidth=1.0))

# Legend for HP
hp_handles = [mpatches.Patch(facecolor=HP_COLORS[h], label=HP_LABELS[h],
                              edgecolor="black", linewidth=0.3)
              for h in HP_ORDER]
ax_hp.legend(handles=hp_handles, loc="lower left", fontsize=7.6, ncol=6,
              bbox_to_anchor=(0.0, -0.95), framealpha=0.95, handlelength=1.2,
              columnspacing=0.6)

ax_hp.set_ylim(-0.06, 1.45)
ax_hp.set_ylabel("HP tag fraction", fontsize=10)
ax_hp.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_hp.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=8.5)
ax_hp.set_xlabel(f"chr2 position (zoom: {WIN_START:,}-{WIN_END:,})", fontsize=10)
ax_hp.set_title("Row 3 — Per-site HP tag composition (bar = HP fraction, n = covered reads)",
                 fontsize=10, loc="left", pad=4)
ax_hp.grid(True, axis="y", linestyle=":", alpha=0.35)
for s in ("top", "right"):
    ax_hp.spines[s].set_visible(False)

# X-tick labels = position
xtick_pos = [s.pos for s in sites]
xtick_labels = [f"{p:,}".replace(",", " ") for p in xtick_pos]
ax_hp.set_xticks(xtick_pos)
ax_hp.set_xticklabels(xtick_labels, rotation=45, ha="right", fontsize=7.2)

# ----- Right column: callout box -----
ax_callout.axis("off")
ax_callout.set_xlim(0, 1)
ax_callout.set_ylim(0, 1)

ax_callout.add_patch(FancyBboxPatch((0.02, 0.04), 0.96, 0.92,
                                     boxstyle="round,pad=0.02",
                                     facecolor="#FFF8E1",
                                     edgecolor=COL_FP, linewidth=1.4))
callout_lines = [
    (f"chr2:{FP_POS:,} — FP 候選", 11.5, "bold", COL_FP),
    ("─" * 30, 8.5, "normal", "#5D4037"),
    ("• ClairS-TO ssrs PASS:", 9.5, "bold", "#212121"),
    (f"    DP=77, A=49%, G=51%", 9, "normal", "#212121"),
    (f"    AF={ssrs[FP_POS]['af']:.4f} (看似 50:50)", 9, "normal", "#212121"),
    ("", 5, "normal", "#000000"),
    ("• SEQC2 TP VCF: 不存在", 9.5, "bold", "#212121"),
    ("    (high-confidence sSNV", 8.5, "normal", "#424242"),
    ("     未列入 truth set)", 8.5, "normal", "#424242"),
    ("", 5, "normal", "#000000"),
    ("• LOH 區段內:", 9.5, "bold", "#212121"),
    (f"    chr2:{target_loh[0]:,}", 8.6, "normal", "#212121"),
    (f"    -{target_loh[1]:,}", 8.6, "normal", "#212121"),
    (f"    ({(target_loh[1] - target_loh[0]) / 1e6:.1f} Mb)", 8.6, "normal", "#212121"),
    ("", 5, "normal", "#000000"),
    ("• HP2-related (32 kb 區域):", 9.5, "bold", "#212121"),
    (f"    HP2:{region_hp_counts.get(2, 0)}"
     f" + HP2-1:{region_hp_counts.get(21, 0)}", 8.6, "normal", "#212121"),
    (f"    = {hp2_related}/{hp_region_total}"
     f" ({hp2_frac:.1f}%)", 8.6, "normal", "#212121"),
    (f"    HP1:{region_hp_counts.get(1, 0)}, "
     f"HP1-1:{region_hp_counts.get(11, 0)},", 8.0, "normal", "#424242"),
    (f"    HP33:{region_hp_counts.get(33, 0)}, "
     f"no-HP:{region_hp_counts.get(0, 0)}", 8.0, "normal", "#424242"),
    ("", 5, "normal", "#000000"),
    ("• 機制:", 9.5, "bold", "#5D4037"),
    ("    LOH 崩潰單倍型 →", 8.6, "normal", "#212121"),
    ("    一邊變 100% 高，", 8.6, "normal", "#212121"),
    ("    看似 50:50 是因", 8.6, "normal", "#212121"),
    ("    ~50% normal contam.", 8.6, "normal", "#212121"),
]
y0 = 0.93
for txt, sz, weight, col in callout_lines:
    ax_callout.text(0.06, y0, txt, fontsize=sz, color=col, weight=weight,
                    transform=ax_callout.transAxes, ha="left", va="top",
                    family=plt.rcParams["font.family"])
    y0 -= 0.0315 + (sz - 8.5) * 0.0030

# ----- Super title -----
fig.suptitle(f"F7 — chr2:{WIN_START:,}-{WIN_END:,} (~{(WIN_END - WIN_START) / 1000:.1f} kb) "
              f"region evidence: 15 sites × LOH × HP collapse → why chr2:{FP_POS:,} is FP",
              fontsize=12.5, fontweight="bold", x=0.06, y=0.985, ha="left")

# ----- Save -----
fig.savefig(OUT_FIG, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"[A5] Figure written: {OUT_FIG}", flush=True)

# ---------------------------------------------------------------- summary
print("\n========================= SUMMARY =========================")
print(f"Window: {CHROM}:{WIN_START:,}-{WIN_END:,} ({(WIN_END - WIN_START) / 1000:.1f} kb)")
print(f"Target LOH: {CHROM}:{target_loh[0]:,}-{target_loh[1]:,} "
      f"({(target_loh[1] - target_loh[0]) / 1e6:.1f} Mb)")
print(f"Sites: {len(sites)} (plan {len(PLAN_SITES)} ∪ TP {len(tp_pos)} ∪ ssrs {len(ssrs)})")
print(f"  TP ∩ ssrs (true somatic): {sum(1 for s in sites if s.in_tp and s.in_ssrs)}")
print(f"  ssrs-only (FP candidates): "
      f"{sum(1 for s in sites if s.in_ssrs and not s.in_tp)}")
print(f"  TP-only (missed by ClairS): "
      f"{sum(1 for s in sites if s.in_tp and not s.in_ssrs)}")
print(f"  no-call (plan listed): "
      f"{sum(1 for s in sites if not s.in_tp and not s.in_ssrs)}")
print(f"  all sites inside target LOH: "
      f"{sum(1 for s in sites if s.in_loh)} / {len(sites)}")
print(f"Region HP profile (32 kb): {dict(region_hp_counts)}")
print(f"  HP2-related (HP2 + HP2-1): {hp2_related}/{hp_region_total} ({hp2_frac:.1f}%)")
print(f"Outputs:")
print(f"  Figure: {OUT_FIG}")
print(f"  TSV   : {OUT_TSV}")
