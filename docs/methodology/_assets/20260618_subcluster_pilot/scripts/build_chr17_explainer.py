#!/usr/bin/env python3
"""[chr17 subclone worked-example 圖文解釋頁] §0-§7 敘述 + 完整克隆樹/SVG read 關聯/甲基熱圖/per-CpG 歸因 內嵌。
§13-A JSON 注入。輸出 ../../20260625_chr17_subclone_worked_example_01.standalone.html。"""
import json, io, base64
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_chr17_subclone_worked_example_01.standalone.html"
D = json.load(open(f"{A}/chr17_complete_data.json")); T = json.load(open(f"{A}/chr17_tree_data.json"))
reads = D["reads"]; cpgs = D["cpgs"]; lc = D["lineage_counts"]
TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; AC = "#D97757"; BLU = "#4A6E8A"; PUR = "#8A6FA0"; GRN = "#5B8A5B"; RED = "#C0563F"
LABS = ["ancestral(no_somatic)", "γ_subclone(sibling)", "L1_α_only(ancestor)", "L2_αβ(descendant)"]
LS = {"ancestral(no_somatic)": "ancestral", "γ_subclone(sibling)": "γ subclone", "L1_α_only(ancestor)": "L1 α-only", "L2_αβ(descendant)": "L2 α+β"}
LC = {"ancestral(no_somatic)": MUT, "γ_subclone(sibling)": PUR, "L1_α_only(ancestor)": BLU, "L2_αβ(descendant)": AC}
plt.rcParams.update({"font.size": 10, "text.color": TX, "axes.edgecolor": BD, "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=115, bbox_inches="tight"); plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def fig_tree():
    fig, ax = plt.subplots(figsize=(8.5, 4)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 8)
    N = {"anc": (6, 7, MUT, f"germline hap1\nancestral · {lc.get('ancestral(no_somatic)',0)}"),
         "γ": (2.5, 4, PUR, f"γ subclone\n+48357368 · {lc.get('γ_subclone(sibling)',0)}"),
         "α": (8.5, 4.5, BLU, f"α subclone\n+48365089 · {lc.get('L1_α_only(ancestor)',0)+lc.get('L2_αβ(descendant)',0)}"),
         "L1": (7, 1.5, BLU, f"L1 α-only · {lc.get('L1_α_only(ancestor)',0)}"), "L2": (10.5, 1.5, AC, f"L2 α+β · {lc.get('L2_αβ(descendant)',0)}")}
    for a, b in [("anc", "γ"), ("anc", "α"), ("α", "L1"), ("α", "L2")]:
        ax.annotate("", xy=(N[b][0], N[b][1] + 0.6), xytext=(N[a][0], N[a][1] - 0.6), arrowprops=dict(arrowstyle="-|>", color=TX, lw=1.5))
    for k, (x, y, c, l) in N.items():
        ax.add_patch(plt.Circle((x, y), 0.92, color=c, alpha=0.85)); ax.text(x, y, l, ha="center", va="center", fontsize=7, color="white")
    ax.text(3.6, 5.8, "sibling 互斥", fontsize=8, color=PUR); ax.text(10.7, 3.1, "nested", fontsize=8, color=AC)
    return b64(fig)


def svg_bar():
    grp = {g: [r for r in reads if r["lineage"] == g] for g in LABS}
    rh = 9; x0 = 88; cw = 40; cols = ["γ", "α", "β1", "β2"]
    H = sum(len(grp[g]) * rh + 16 for g in LABS if grp[g]) + 30
    s = [f'<svg viewBox="0 0 540 {H}" xmlns="http://www.w3.org/2000/svg" font-size="10" font-family="system-ui">']
    for k, nm in enumerate(cols):
        s.append(f'<text x="{x0+k*cw}" y="13" font-size="10">{nm}</text>')
    s.append(f'<text x="{x0+4*cw+10}" y="13" font-size="9">methyl</text>')
    y = 22
    for g in LABS:
        if not grp[g]:
            continue
        s.append(f'<text x="6" y="{y+11}" fill="{LC[g]}" font-weight="bold" font-size="10">{LS[g]} ({len(grp[g])})</text>')
        for r in sorted(grp[g], key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0):
            for k, nm in enumerate(cols):
                al = r["geno"].get(nm); col = "#cfcabf" if al not in ("REF", "ALT") else (BLU if al == "REF" else RED)
                s.append(f'<rect x="{x0+k*cw}" y="{y}" width="{cw-3}" height="{rh-1.5}" fill="{col}"/>')
            mv = np.nanmean([np.nan if v is None else v for v in r["meth"]]); mv = 0 if np.isnan(mv) else mv
            s.append(f'<rect x="{x0+4*cw+10}" y="{y}" width="{mv*100:.0f}" height="{rh-1.5}" fill="#7a9a7a"/>')
            y += rh
        y += 16
    s.append(f'<rect x="{x0}" y="{H-14}" width="9" height="9" fill="{BLU}"/><text x="{x0+12}" y="{H-6}">REF</text><rect x="{x0+48}" y="{H-14}" width="9" height="9" fill="{RED}"/><text x="{x0+60}" y="{H-6}">ALT</text></svg>')
    return "".join(s)


def fig_meth():
    ordered = []
    for g in LABS:
        gg = sorted([r for r in reads if r["lineage"] == g], key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0)
        ordered += [(g, r) for r in gg]
    mat = np.array([[np.nan if v is None else v for v in r["meth"]] for _, r in ordered])
    fig, ax = plt.subplots(figsize=(11, 4.2)); cm = plt.cm.RdBu_r.copy(); cm.set_bad("#e8e4dc")
    ax.imshow(mat, aspect="auto", cmap=cm, vmin=0, vmax=1); y = 0
    for g in LABS:
        c = sum(1 for l, _ in ordered if l == g)
        if c == 0:
            continue
        if y > 0:
            ax.axhline(y - 0.5, color=TX, lw=1.1)
        ax.text(-2.5, y + c / 2, LS[g], va="center", ha="right", fontsize=9.5, color=LC[g], fontweight="bold"); y += c
    ax.set_xlabel(f"CpG (n={len(cpgs)})"); ax.set_yticks([]); ax.set_title("甲基熱圖（按基因型 lineage 分組）", fontsize=10)
    return b64(fig)


def fig_attrib():
    bc = T["best_axis_counts"]; fig, ax = plt.subplots(figsize=(7, 2.8))
    axn = list(bc.keys()); vals = [bc[k] for k in axn]
    cmap = {"ALT@48365089(α)": BLU, "lineage_L1_vs_L2": AC, "ALT@48365161(β2)": GRN, "ALT@48362515(β1)": "#C98A5B"}
    ax.barh(range(len(axn)), vals, color=[cmap.get(k, MUT) for k in axn])
    ax.set_yticks(range(len(axn))); ax.set_yticklabels([k.replace("ALT@", "").replace("lineage_", "") for k in axn], fontsize=8.5)
    for i, v in enumerate(vals):
        ax.text(v + 0.2, i, str(v), va="center", fontsize=8.5)
    ax.set_title("每差異 CpG 最佳歸因軸", fontsize=10); ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


C = {"tree": fig_tree(), "svg": svg_bar(), "meth": fig_meth(), "attrib": fig_attrib()}
sn = D["snv_norm"]
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>chr17 subclone worked example</title><style>
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:{TX};background:#FAF9F5;line-height:1.75}}
.wrap{{max-width:880px;margin:0 auto;padding:0 24px 80px}}
header{{background:{TX};color:#FAF9F5;padding:30px 24px}} header h1{{margin:0;font-size:23px}} header .sub{{color:#cfc9bf;font-size:14px;margin-top:5px}}
nav{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid {BD};padding:9px 24px;font-size:12.5px;z-index:9}} nav a{{color:{MUT};text-decoration:none;margin-right:12px}}
h2{{font-size:18px;margin:30px 0 10px;color:{TX}}} h2 .n{{color:{AC};margin-right:8px}}
p{{margin:9px 0}} img,svg{{max-width:100%;border:1px solid {BD};border-radius:6px;margin:10px 0;background:white;display:block}} svg{{padding:10px}}
.fig{{text-align:center}} .cap{{font-size:12.5px;color:{MUT};margin-top:-4px;margin-bottom:14px}}
.box{{background:white;border:1px solid {BD};border-left:4px solid {AC};border-radius:0 6px 6px 0;padding:12px 16px;margin:14px 0}}
.red{{border-left-color:{RED};background:#fbf0ec}} .ok{{border-left-color:{GRN};background:#eef4ee}}
pre{{background:#f3efe7;border:1px solid {BD};border-radius:6px;padding:12px;font-size:12.5px;overflow-x:auto;line-height:1.5}}
table{{border-collapse:collapse;font-size:13px;margin:8px 0}} th,td{{border:1px solid {BD};padding:4px 9px;text-align:center}} th{{background:#efe9df}}
b.k{{color:{AC}}} footer{{margin-top:34px;padding-top:14px;border-top:1px solid {BD};font-size:12px;color:{MUT}}}
.tag{{display:inline-block;background:{AC};color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>chr17:48360161 — 一個 subclone 從「候選」到「完整重建」<span class="tag">⭐3 HCC1395</span></h1>
<div class="sub">用一個位點把整套 subclonal reconstruction 的方法、陷阱、洞見串成可重述的故事</div></div></header>
<nav class="wrap" style="max-width:880px"><a href="#s1">①候選</a><a href="#s2">②連鎖</a><a href="#s3">③矛盾</a><a href="#s4">④完整樹</a><a href="#s5">⑤甲基</a><a href="#s6">⑥教訓</a></nav>
<div class="wrap">

<p style="font-size:15px"><b>為什麼選這個位點</b>：它把每一環都示範了——甲基如何「提出候選」、為何甲基單獨不能定 subclone、sSNV 連鎖如何「非循環確認」、一個矛盾如何揭露漏掉的 somatic、完整克隆樹如何重建、甲基最終的正確角色。讀懂這一個，就讀懂整個 pipeline。</p>

<h2 id="s1"><span class="n">①</span>起點：甲基 cis-test 篩出的「候選」</h2>
<p>34,736 位點 → 甲基結構 1,139 候選 → normal-anchored cis-test 洗掉循環假象 → 9 個「tumor-specific、非循環」候選，<b class="k">chr17:48360161 是其一</b>。此時只知「這裡甲基不是 germline cis、像 tumor 特有」——<b>還不知是不是真 subclone</b>（甲基切群有 double-dip 循環風險，單獨永遠不能回答）。</p>

<h2 id="s2"><span class="n">②</span>第一層遺傳錨：sSNV read-level 連鎖</h2>
<p>要非循環確認需<b>遺傳證據</b>：同一條 long read 上多個 somatic SNV 的共現。此窗 filtered-TP 有 3 個：β1(48362515)、α(48365089)、β2(48365161)。逐 read 共現 2×2 顯示：<b>α×β2 的 REF-ALT=0 → β2 嵌套在 α 內</b>（α 祖先、β2 衍生）；β1×β2 完美共連。→ α 較早，其上累積 β 形成<b>後代 subclone</b>（nested 線性演化，有方向）。</p>

<h2 id="s3"><span class="n">③</span>一個矛盾：HP1-1 但三個 sSNV 全 REF</h2>
<p>按 α/β 基因型分群後，有一群 read 在我們 3 個 sSNV <b>全 REF</b>（看似「無 somatic 的祖先」），但其中 <b class="k">5 條被 LongPhase-S 標 HP1-1</b>。</p>
<div class="box"><b>🔑 關鍵推理</b>：HP1-1 定義 = 「hap1 上、被 phase-block 內<b>任一</b> somatic 標記的 read」。所以<b>「HP1-1 但已知 sSNV 全 REF」⟹ 這條 read 必然帶著一個我們沒看到的 somatic</b>。HP tag 比手上 VCF「知道更多」。</div>

<h2 id="s4"><span class="n">④</span>找出漏掉的 γ → 完整克隆樹</h2>
<p>順矛盾追：這 5 條 read 的 phase-set PS=46496608（~1.9Mb），span 48326-48381kb。filtered-TP 內沒有第 4 個，但 longphase 輸入集 + FP 集裡有 <b class="k">γ=48357368(C→T)</b>。逐 read 定型：<b>那 5 條 HP1-1 read 在 γ 是 ALT</b>，γ normal=35/0 REF + tumor VAF 0.18 = <b>ONT 偵測候選</b>（SEQC2 未收 HighConf/superSet＝缺金標準確認）；<b>γ×α 互斥(AA=0) → γ、α 是兩 sibling 分支</b>。</p>
<div class="box red">🔴 <b>γ 不在 SEQC2 HighConf/superSet</b>（落 callable HC region 內＝缺確認，非主動判 FP；落到我們 pipeline filtered-FP 集），所以只用 filtered-TP 時漏掉它、把 γ subclone 誤標「ancestral root」。補回 γ → 完整局部克隆樹：</div>
<div class="fig"><img src="{C['tree']}" alt="tree" style="max-width:680px"></div>
<p class="cap">完整局部克隆樹：germline → γ ∥ α 兩 sibling 分支；α 內再分 L1(α-only) + L2(α+β nested 後代)。4 細胞群，全 4 位點 normal=REF（α/β SEQC2 HighConf 確認；γ ONT 候選）。</p>
<table><tr><th>somatic</th><th>pos</th><th>tumor R/A</th><th>normal R/A</th></tr>
<tr><td>γ</td><td>48357368*</td><td>{sn['γ']['tumor_REF']}/{sn['γ']['tumor_ALT']}</td><td>{sn['γ']['normal_REF']}/{sn['γ']['normal_ALT']}</td></tr>
<tr><td>α</td><td>48365089</td><td>{sn['α']['tumor_REF']}/{sn['α']['tumor_ALT']}</td><td>{sn['α']['normal_REF']}/{sn['α']['normal_ALT']}</td></tr>
<tr><td>β1</td><td>48362515</td><td>{sn['β1']['tumor_REF']}/{sn['β1']['tumor_ALT']}</td><td>{sn['β1']['normal_REF']}/{sn['β1']['normal_ALT']}</td></tr>
<tr><td>β2</td><td>48365161</td><td>{sn['β2']['tumor_REF']}/{sn['β2']['tumor_ALT']}</td><td>{sn['β2']['normal_REF']}/{sn['β2']['normal_ALT']}</td></tr></table>
<p class="cap">*γ=SEQC2 未收（callable HC region 內 ONT 候選），normal=REF + VAF 0.18 somatic-like。每條 read 的基因型（下，藍REF 紅ALT）直接看出 4 個克隆群：</p>
<div class="fig">{C['svg']}</div>

<h2 id="s5"><span class="n">⑤</span>甲基的正確角色：characterize 已驗 lineage（非循環）</h2>
<p>lineage 一旦由 <b>sSNV（遺傳）定義好</b>，再看甲基就<b>不循環</b>。實測 <b class="k">L1 祖先 vs L2 後代 16 個 CpG 甲基顯著不同</b>；per-CpG 可歸因：<b>23 個關連 α（祖先 ASM）、6 個關連 L1→L2 轉變</b>。</p>
<div class="fig"><img src="{C['meth']}" alt="meth"></div>
<p class="cap">甲基熱圖（reads 按基因型 lineage 分組）：甲基模式隨遺傳定義的克隆群變化 = characterize。</p>
<div class="fig"><img src="{C['attrib']}" alt="attrib" style="max-width:620px"></div>
<p class="cap">每差異 CpG 歸因到基因驅動者（α 祖先 / L1→L2 後代 / β）。</p>
<div class="box red">🔴 <b>對照——甲基單獨做不到</b>：ISM 無監督甲基 clustering 只切出 1-2(L0/γ 為主) 與 1-1(L1+L2 <b>混</b>)；BERNOULLI 距離 UPGMA 樹也一樣——<b>只切得出「有無 somatic」，切不出 L1↔L2 subclone</b>。→ 甲基<b>提出候選 + 刻畫已驗 lineage</b>，subclone 細分必須靠遺傳錨。</div>

<h2 id="s6"><span class="n">⑥</span>三個方法學教訓</h2>
<div class="box"><b>1. 甲基 ≠ subclone caller</b>：genotype-同質群內甲基結構永遠 double-dip；甲基=characterize/corroborate，subclone 非循環確認來自<b>同 germline-HP 上 ≥2 somatic SNV 的 read 共現</b>。</div>
<div class="box"><b>2.「HP-tag 但已知 sSNV 全 REF」= 偵測漏掉 somatic 的訊號</b>；完整重建應用 <b>longphase 輸入集(methyl_PASS) + normal 確認</b>，非只 SEQC2-filtered TP（會移除如 γ 的 clonally-informative 候選——ONT 偵測但 SEQC2 未收）。</div>
<div class="box"><b>3. 兩個尺度</b>：HP tag 在 phase-block(~Mb) 標 lineage（靠 phasing）；read 共現只在 ≤read-長(~10-50kb) 直接連。互補，但跨 phase-block 無法 read 連（reconstruction GAP）。</div>
<div class="box ok"><b>對論文</b>：這體現「Subclonal reconstruction using somatic haplotagging and methylation profiles」—— <b>somatic haplotagging(sSNV 連鎖+HP tag)=重建骨幹</b>（定 γ/α/L1/L2 克隆群）、<b>methylation profiles=characterize</b>（刻畫表觀，非偵測驅動）。🔴 邊界：⭐3 單樣本、局部克隆階層非 genome-wide tree。</div>

<footer>chr17:48360161 · HCC1395 ⭐3 · §13-A 由 chr17_complete_data.json + chr17_tree_data.json 注入 · 敘述源 InterSubMod/docs/methodology/20260625_chr17_subclone_worked_example_explained_01.md</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes")
