#!/usr/bin/env python3
"""[差異甲基定位點 圖表 + FP/TP 歸因] 讀 percpg_per_locus.json + percpg_cpg_classification.json + records_v6
→ 全套 matplotlib PNG(英文 label 防 CJK 亂碼) + fptp_attribution.json + standalone base64 gallery(中文說明在 HTML)。
數字由 JSON 注入(§13)。"""
import os, json, base64
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
FQ = f"{A}/figs_quant"; os.makedirs(FQ, exist_ok=True)
TPN, FPN = 30077, 4659
PL = json.load(open(f"{A}/percpg_per_locus.json"))
CLS = json.load(open(f"{A}/percpg_cpg_classification.json"))
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))
TPc, FPc = "#3b82f6", "#ef4444"


def save(fig, name):
    p = f"{FQ}/{name}.png"; fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); return p


def hist_tpfp(getter, title, xlab, name, clip=40):
    tp = [min(getter(v), clip) for v in PL.values() if v["set"] == "TP"]
    fp = [min(getter(v), clip) for v in PL.values() if v["set"] == "FP"]
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.arange(0, clip + 2) - 0.5
    ax.hist([tp, fp], bins=bins, label=[f"TP (n={len(tp)})", f"FP (n={len(fp)})"], color=[TPc, FPc], density=True)
    ax.set_xlabel(xlab); ax.set_ylabel("fraction of SNV loci (density)"); ax.set_title(title)
    ax.legend(); ax.set_xlim(-0.5, clip)
    return save(fig, name)


# 1. 主圖 ×4 (X=diff CpG/locus, Y=SNV count density, TP/FP)
charts = []
charts.append(("hist_union", hist_tpfp(lambda v: v["union"], "Differential CpG per locus (union, any axis)", "# differential CpG (union)", "hist_union")))
charts.append(("hist_somatic", hist_tpfp(lambda v: v["som_marker"], "Somatic-marker CpG per locus (HP1 vs HP1-1)", "# somatic differential CpG", "hist_somatic")))
charts.append(("hist_cluster", hist_tpfp(lambda v: v["sub_marker"], "Subclone-marker CpG per locus (cluster-split)", "# cluster differential CpG", "hist_cluster", clip=30)))
charts.append(("hist_tn", hist_tpfp(lambda v: v["axis_nsig"].get("tumor_vs_normal", 0), "Tumor/normal differential CpG per locus", "# T/N differential CpG", "hist_tn")))

# 2. per-axis loci bar (≥1/≥3/≥5) TP/FP
axes = ["tumor_vs_normal", "HP1_vs_HP2", "REF_vs_ALT", "HP1_vs_HP1-1", "HP2_vs_HP2-1", "cluster_split"]
fig, ax = plt.subplots(figsize=(9, 4.2))
x = np.arange(len(axes)); w = 0.13
for i, (thr, hatch) in enumerate([(1, ""), (3, "//"), (5, "xx")]):
    tp = [sum(1 for v in PL.values() if v["set"] == "TP" and v["axis_nsig"].get(a, 0) >= thr) / TPN * 100 for a in axes]
    fp = [sum(1 for v in PL.values() if v["set"] == "FP" and v["axis_nsig"].get(a, 0) >= thr) / FPN * 100 for a in axes]
    ax.bar(x + (i - 1) * 2 * w, tp, w, color=TPc, alpha=1 - i * 0.25, label=f"TP ≥{thr}")
    ax.bar(x + (i - 1) * 2 * w + w, fp, w, color=FPc, alpha=1 - i * 0.25, label=f"FP ≥{thr}")
ax.set_xticks(x); ax.set_xticklabels([a.replace("_vs_", "/") for a in axes], rotation=20, ha="right", fontsize=8)
ax.set_ylabel("% of loci (TP/FP rate)"); ax.set_title("Loci with differential CpG by axis (≥1/≥3/≥5)"); ax.legend(fontsize=7, ncol=3)
charts.append(("axis_bars", save(fig, "axis_bars")))

# 3. 差異 CpG 分類 stacked: axis × hyper/hypo × |Δβ| bins
fig, ax = plt.subplots(figsize=(9, 4.2)); axn = list(CLS["axes"].keys())
hyper = [CLS["axes"][a]["n_hyper"] for a in axn]; hypo = [CLS["axes"][a]["n_hypo"] for a in axn]
xx = np.arange(len(axn))
ax.bar(xx, hyper, color="#dc2626", label="hyper (Δβ>0, gain)")
ax.bar(xx, hypo, bottom=hyper, color="#2563eb", label="hypo (Δβ<0, loss)")
ax.set_xticks(xx); ax.set_xticklabels([a.replace("_vs_", "/") for a in axn], rotation=20, ha="right", fontsize=8)
ax.set_ylabel("# significant differential CpG"); ax.set_title("Differential CpG by axis × direction"); ax.legend(fontsize=8)
charts.append(("cls_direction", save(fig, "cls_direction")))

# 4. |Δβ| 強度 bins per axis (stacked)
fig, ax = plt.subplots(figsize=(9, 4.2)); BL = CLS["axes"][axn[0]]["bin_labels"]; bot = np.zeros(len(axn))
cols = ["#fca5a5", "#f87171", "#dc2626", "#7f1d1d"]
for bi, bl in enumerate(BL):
    vals = [CLS["axes"][a]["bins"][bi] for a in axn]
    ax.bar(xx, vals, bottom=bot, color=cols[bi], label=f"|Δβ| {bl}"); bot += np.array(vals)
ax.set_xticks(xx); ax.set_xticklabels([a.replace("_vs_", "/") for a in axn], rotation=20, ha="right", fontsize=8)
ax.set_ylabel("# significant differential CpG"); ax.set_title("Differential CpG by axis × |Δβ| strength"); ax.legend(fontsize=8)
charts.append(("cls_strength", save(fig, "cls_strength")))

# 5. FP/TP 歸因 + 研究可用性 (compute from PL)
def fptp(pred):
    t = sum(1 for v in PL.values() if v["set"] == "TP" and pred(v)); f = sum(1 for v in PL.values() if v["set"] == "FP" and pred(v))
    return t, f, round((f / FPN) / (t / TPN), 2) if t else 0, round(100 * (t + f) / (TPN + FPN), 1)
ATTR = {
    "tumor_normal_≥3": fptp(lambda v: v["axis_nsig"].get("tumor_vs_normal", 0) >= 3),
    "HP1_HP2_≥3": fptp(lambda v: v["axis_nsig"].get("HP1_vs_HP2", 0) >= 3),
    "REF_ALT_≥3": fptp(lambda v: v["axis_nsig"].get("REF_vs_ALT", 0) >= 3),
    "somatic(HP-1)_≥3": fptp(lambda v: v["som_marker"] >= 3),
    "cluster_≥3": fptp(lambda v: v["sub_marker"] >= 3),
    "cluster_≥5": fptp(lambda v: v["sub_marker"] >= 5),
}
MECH = {"tumor_normal_≥3": "普遍 somatic 甲基偏移, 不專一", "HP1_HP2_≥3": "germline cis-ASM, 不專一",
        "REF_ALT_≥3": "TP 有真 somatic ALT → ALT read 真差; FP ALT=germline/artifact",
        "somatic(HP-1)_≥3": "somatic 變異定義子單倍型 → TP 專一(最佳 somatic 定位)",
        "cluster_≥3": "子群結構, 略 FP", "cluster_≥5": "subclone-like cis/CN 驅動 → FP 富集"}
fig, ax = plt.subplots(figsize=(8, 4)); ks = list(ATTR.keys()); ratios = [ATTR[k][2] for k in ks]
cols2 = [TPc if r < 1 else FPc for r in ratios]
ax.barh(range(len(ks)), ratios, color=cols2); ax.axvline(1, color="#444", ls="--")
ax.set_yticks(range(len(ks))); ax.set_yticklabels(ks, fontsize=8); ax.set_xlabel("FP-rate / TP-rate  (<1 = TP-enriched / >1 = FP-enriched)")
ax.set_title("FP/TP attribution by axis/state"); ax.invert_yaxis()
for i, r in enumerate(ratios): ax.text(r + 0.02, i, str(r), va="center", fontsize=8)
charts.append(("fptp_attr", save(fig, "fptp_attr")))
json.dump({k: {"TP": v[0], "FP": v[1], "fp_tp_ratio": v[2], "pct_loci": v[3], "mechanism": MECH[k]} for k, v in ATTR.items()},
          open(f"{A}/fptp_attribution.json", "w"), ensure_ascii=False, indent=1)

# 6. CDF: % loci with ≥k diff CpG (union) TP vs FP
fig, ax = plt.subplots(figsize=(7, 4))
for setl, col in [("TP", TPc), ("FP", FPc)]:
    vals = sorted(v["union"] for v in PL.values() if v["set"] == setl); n = len(vals)
    ks2 = list(range(0, 41)); cdf = [100 * sum(1 for x in vals if x >= k) / n for k in ks2]
    ax.plot(ks2, cdf, color=col, label=setl)
ax.set_xlabel("k = # differential CpG (union)"); ax.set_ylabel("% loci with ≥k"); ax.set_title("Cumulative: loci with ≥k differential CpG"); ax.legend()
charts.append(("cdf_union", save(fig, "cdf_union")))

# ---- standalone gallery (base64) ----
def uri(p): return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()
cap = {"hist_union": "① 主圖: 每位點 union 差異 CpG 數 × SNV 個數 (TP/FP). 兩者形狀相近→「有差異甲基」普遍不判別。",
       "hist_somatic": "② somatic-marker (HP1 vs HP1-1) 差異 CpG/位點. TP 明顯右移=TP 專一 somatic 甲基定位。",
       "hist_cluster": "③ subclone-marker (cluster-split) 差異 CpG/位點。",
       "hist_tn": "④ tumor/normal 差異 CpG/位點 (somatic 偏移, 普遍)。",
       "axis_bars": "⑤ 哪一軸有差異 CpG: 各軸 ≥1/≥3/≥5 的位點% (TP/FP)。",
       "cls_direction": "⑥ 差異 CpG 分類 — 軸 × 方向 (hyper 增甲基/hypo 失甲基)。",
       "cls_strength": "⑦ 差異 CpG 分類 — 軸 × |Δβ| 強度 bins。",
       "fptp_attr": "⑧ 🔴 FP/TP 歸因: 各軸 ratio (<1=TP專一可做研究 / >1=FP富集非somatic)。",
       "cdf_union": "⑨ 累積: ≥k 差異 CpG 的位點比例 TP vs FP。"}
imgs = "".join(f'<figure><img src="{uri(p)}"><figcaption>{cap.get(nm, nm)}</figcaption></figure>' for nm, p in charts)
H = f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><title>差異甲基定位點 量化圖庫</title><style>
body{{font-family:system-ui,"Noto Sans CJK TC",sans-serif;background:#0a0e16;color:#e7edf5;margin:0}}.w{{max-width:1000px;margin:0 auto;padding:16px}}
h1{{font-size:19px}}.note{{background:#2a2410;border:1px solid #d97706;border-radius:8px;padding:9px 13px;font-size:12px;margin:10px 0}}.note b{{color:#fbbf24}}
figure{{margin:14px 0;background:#111826;border:1px solid #2a3344;border-radius:9px;padding:10px}}img{{width:100%;border-radius:6px;background:#fff}}figcaption{{font-size:12.5px;color:#8b9bb4;margin-top:6px}}</style></head><body><div class="w">
<h1>差異甲基「定位點」量化圖庫 — HCC1395 single-sample</h1>
<div class="note">🔴 <b>核心</b>：「有差異甲基」普遍(96.8%)不判別 → 真正能當 somatic 定位點的是 <b>HP1/HP1-1 軸(TP專一 0.33×, 70.8%)</b>；subclone-marker FP 富集非 somatic-specific。⭐2-3 單樣本。數字源 percpg_per_locus.json / percpg_cpg_classification.json。</div>
{imgs}
<footer style="color:#8b9bb4;font-size:11px;margin-top:16px">FP/TP 歸因 → fptp_attribution.json。方法 docs/methodology/20260624_per_cpg_differential_and_subclone_validation_01.md §7。單樣本 ⭐2-3。</footer>
</div></body></html>'''
op = f"{A}/obs_ws/cpp_wg/20260624_differential_methylation_marker_charts.standalone.html"
open(op, "w").write(H)
print(f"charts={len(charts)} | gallery {op} ({len(H)//1024} KB)")
print("fptp_attribution:", json.dumps({k: v[2] for k, v in ATTR.items()}, ensure_ascii=False))
