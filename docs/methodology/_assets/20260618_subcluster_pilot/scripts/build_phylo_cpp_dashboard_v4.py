#!/usr/bin/env python3
"""[全位點觀察分類儀表板 v4] v3 + ⓪原始數據盤點 + 每 aggregate「SVG 且 數字+%」 + 摘要一行可展開
 + 更多參數交叉表 + 甲基(β)觀察(per-locus kv-grid + aggregate + by 狀況)。
數字全由 records_v4 / dashboard_stats_v4 注入(§13-A: generator 不手打)。圖名 JS 決定式(無 FIGMAP)。"""
import json, os

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUTD = f"{A}/obs_ws/cpp_wg"; os.makedirs(OUTD, exist_ok=True)
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))
def pcv_idx(v):
    return 0 if v == "no_cluster" else 1 if str(v).startswith("cis-ASM") else 2 if v == "subclone_novel" else 3 if v == "cluster_no_sigCpG" else 4 if v == "weak_scattered" else 5
def ls_idx(p, disp):
    return 0 if (p is None or p >= 0.05) else (2 if disp else 1)
PCV = ["no_cluster", "cis-ASM", "⭐subclone_novel", "no_sigCpG", "weak_scattered", "no_meth"]
S = json.load(open(f"{A}/dashboard_stats_v4.json"))
# 差異甲基定位點 量化(④b, aggregate only, 不加 per-locus 欄控 size)
QPL = json.load(open(f"{A}/percpg_per_locus.json")) if os.path.exists(f"{A}/percpg_per_locus.json") else {}
QCLS = json.load(open(f"{A}/percpg_cpg_classification.json")) if os.path.exists(f"{A}/percpg_cpg_classification.json") else {"axes": {}}
QATTR = json.load(open(f"{A}/fptp_attribution.json")) if os.path.exists(f"{A}/fptp_attribution.json") else {}
CHRS = [f"chr{c}" for c in range(1, 23)]; ck = {c: i for i, c in enumerate(CHRS)}
CNS = ["gain", "loss", "loh", "neutral"]; cni = {s: i for i, s in enumerate(CNS)}
WST = ["S1", "S2", "S3", "S4", "S6"]; wi = {s: i for i, s in enumerate(WST)}
ALS = ["aligned", "unaligned", "NA"]; ali = {s: i for i, s in enumerate(ALS)}
# 8 類分類(corrected tree) + cluster×label 列聯型
CAT8 = ["A", "B1", "B2", "B3", "C-S1", "C-S2", "C-S3", "D"]; cat8i = {c: i for i, c in enumerate(CAT8)}
CTYPE = ["one2one", "many1_結構>標籤", "1many_跨標籤", "mixed", "single_cluster"]; cti = {c: i for i, c in enumerate(CTYPE)}
from collections import Counter as _C
CAT8_CNT = _C(r.get("cat8", "D") for r in rec); CTYPE_CNT = _C(r.get("ctype", "single_cluster") for r in rec)

# figs symlink
figlink = f"{OUTD}/figs"
if not (os.path.islink(figlink) and os.path.realpath(figlink) == f"{A}/figs_cpp_wg_full"):
    try:
        if os.path.islink(figlink) or os.path.exists(figlink): os.remove(figlink)
        os.symlink(f"{A}/figs_cpp_wg_full", figlink)
    except Exception: pass
figset = set(os.listdir(f"{A}/figs_cpp_wg_full")) if os.path.isdir(f"{A}/figs_cpp_wg_full") else set()
# T/N 對照圖 symlink (tumor-vs-normal 甲基)
figlink_tn = f"{OUTD}/figs_tn"; TNDIR = f"{A}/figs_cpp_wg_full_tn"
if not (os.path.islink(figlink_tn) and os.path.realpath(figlink_tn) == TNDIR):
    try:
        if os.path.islink(figlink_tn) or os.path.exists(figlink_tn): os.remove(figlink_tn)
        os.symlink(TNDIR, figlink_tn)
    except Exception: pass
figset_tn = set(os.listdir(TNDIR)) if os.path.isdir(TNDIR) else set()


def r2(x):  # round float to 2dp, -9 sentinel for None
    return round(x, 2) if isinstance(x, (int, float)) else -9


# ---- 緊湊 D (v3 22 欄 + 7 甲基欄) ----
# 22 mb / 23 dbg / 24 dtn / 25 dhp / 26 ncpg  (hypo/hyper 移到 aggregate 控檔案大小 <3.27MB 顯示門檻)
D = []
for r in rec:
    p = int(r["pos"]); snv = p + 5000
    fn = f"cpp_{r['chrom']}_{p}_{p+10000}.png"
    cv = r.get("cn_value")
    D.append([ck.get(r["chrom"], 0), snv, 0 if r["set"] == "TP" else 1, r["n"], r["coarse_ng"], r["fine_ng"],
              r["n_other"], 1 if r["unstable"] else 0, 1 if r["hidden_het"] else 0, r.get("V_hp", 0), r.get("V_allele", 0),
              cni.get(r.get("cn_state", "neutral"), 3), (cv if isinstance(cv, (int, float)) else -1), 1 if r.get("is_loh") else 0,
              1 if fn in figset else 0, r.get("n_tumor", -1) if r.get("n_tumor") is not None else -1,
              r.get("n_normal", -1) if r.get("n_normal") is not None else -1, r.get("geometry_ng", -1) if r.get("geometry_ng") is not None else -1,
              wi.get(r.get("wstate", "S6"), 4), r.get("structure_tier", 3), r.get("geom_divergence", 0), ali.get(r.get("align_status", "NA"), 2),
              r2(r.get("m_mean_beta")), r2(r.get("m_dbeta_group")), r2(r.get("m_dbeta_tn")), r2(r.get("m_dbeta_hp")),
              r.get("m_n_cpg") if r.get("m_n_cpg") is not None else -1,
              cat8i.get(r.get("cat8", "D"), 7), cti.get(r.get("ctype", "single_cluster"), 4),
              1 if (r.get("cat8") == "A" or (str(r.get("ctype", "")).startswith("many1") and r.get("axis") in ("single_label_somatic", "somatic_one_family"))) else 0,
              1 if (r.get("cn_state") == "loh" and r["coarse_ng"] >= 2 and r.get("has_som")) else 0,
              pcv_idx(r.get("pc_verdict")), ls_idx(r.get("ls_hp_p"), r.get("ls_hp_disp"))])


# ================= SVG 生成器 (數字由 stats 注入) =================
def svg_stacked(segs, w=440, h=24):
    tot = sum(v for _, v, _ in segs) or 1; x = 0.0; p = []
    for lab, v, c in segs:
        ww = w * v / tot
        p.append(f'<rect x="{x:.1f}" y="0" width="{ww:.1f}" height="{h}" fill="{c}"/>')
        if ww > 44: p.append(f'<text x="{x+ww/2:.1f}" y="{h/2+4:.0f}" font-size="10" fill="#fff" text-anchor="middle">{100*v/tot:.0f}%</text>')
        x += ww
    return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}" role="img">{"".join(p)}</svg>'


def svg_hist(dist, color="#0891b2", w=440, h=64):
    try:  # numeric keys → 數值排序; 帶標籤桶(如 '0-.05') → 維持插入序
        keys = sorted(dist, key=lambda k: float(k))
    except ValueError:
        keys = list(dist.keys())
    mx = max(dist.values()) or 1; bw = w / max(1, len(keys)); p = []
    for i, k in enumerate(keys):
        bh = h * 0.74 * dist[k] / mx; x = i * bw
        p.append(f'<rect x="{x+1:.1f}" y="{h-bh-12:.1f}" width="{bw-2:.1f}" height="{bh:.1f}" fill="{color}"/>')
        p.append(f'<text x="{x+bw/2:.1f}" y="{h-2:.0f}" font-size="8" fill="#9aa3b2" text-anchor="middle">{k}</text>')
    return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}">{"".join(p)}</svg>'


def svg_funnel(stages, w=440):
    mx = stages[0]["n"] or 1; rh = 24; gap = 5; p = []
    for i, s in enumerate(stages):
        ww = w * s["n"] / mx; y = i * (rh + gap); x = (w - ww) / 2
        p.append(f'<rect x="{x:.1f}" y="{y}" width="{ww:.1f}" height="{rh}" rx="3" fill="#0891b2" opacity="{1-i*0.13:.2f}"/>')
        p.append(f'<text x="{w/2:.0f}" y="{y+rh/2+4:.0f}" font-size="10" fill="#fff" text-anchor="middle">{s["stage"]} {s["n"]:,} ({100*s["n"]/mx:.0f}%)</text>')
    h = len(stages) * (rh + gap); return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}">{"".join(p)}</svg>'


def svg_compare(pairs, w=300, h=86):
    mx = max(v for _, v, _ in pairs) or 1; slot = w / len(pairs); bw = slot * 0.5; p = []
    for i, (lab, v, c) in enumerate(pairs):
        bh = h * 0.62 * v / mx; x = i * slot + (slot - bw) / 2
        p.append(f'<rect x="{x:.1f}" y="{h-bh-16:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="{c}"/>')
        p.append(f'<text x="{x+bw/2:.1f}" y="{h-bh-20:.0f}" font-size="11" fill="#e6e8ec" text-anchor="middle">{v}%</text>')
        p.append(f'<text x="{x+bw/2:.1f}" y="{h-2:.0f}" font-size="9.5" fill="#9aa3b2" text-anchor="middle">{lab}</text>')
    return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}">{"".join(p)}</svg>'


def numrows(pairs):  # pairs = [(label, count, color_or_None)] → 數字+% 逐列(monospace)
    tot = sum(c for _, c, _ in pairs) or 1; rr = []
    for lab, c, col in pairs:
        dot = f'<span style="color:{col}">●</span> ' if col else ''
        rr.append(f'<div class="nr"><span class="nrl">{dot}{lab}</span><span class="nrv">{c:,}</span><span class="nrp">{100*c/tot:.1f}%</span></div>')
    return '<div class="nrows">' + "".join(rr) + '</div>'


def distrows(dist):  # 直方對應的 數字+% 列
    tot = sum(dist.values()) or 1
    keys = list(dist.keys())
    return '<div class="nrows">' + "".join(f'<div class="nr"><span class="nrl">{k}</span><span class="nrv">{v:,}</span><span class="nrp">{100*v/tot:.1f}%</span></div>' for k, v in dist.items()) + '</div>'


CNC = {"gain": "#ef4444", "loss": "#3b82f6", "loh": "#a855f7", "neutral": "#64748b"}
WC = {"S1": "#db2777", "S2": "#9333ea", "S3": "#d97706", "S4": "#64748b", "S6": "#334155"}


def stat_card(title, svg, rows_html, note=""):
    return f'<div class="stcard"><div class="stt">{title}</div>{svg}{rows_html}<div class="sx">{note}</div></div>'


# ================= ⓪ 原始數據盤點 =================
raw = S["raw"]
inv_cards = "".join([
    f'<div class="bignum"><div class="bn">{S["n"]:,}</div><div class="bl">全位點 (HCC1395 single-sample)</div></div>',
    f'<div class="bignum"><div class="bn" style="color:#7cc4ff">{S["tp"]:,}</div><div class="bl">TP <b>{S["tp_pct"]}%</b>（結構率 {raw["tp_structured_pct"]}%）</div></div>',
    f'<div class="bignum"><div class="bn" style="color:#ffa3a3">{S["fp"]:,}</div><div class="bl">FP <b>{S["fp_pct"]}%</b>（結構率 {raw["fp_structured_pct"]}%）</div></div>',
    f'<div class="bignum"><div class="bn">{raw["median_n"]}</div><div class="bl">中位 reads/位點（TP {raw["median_n_tp"]} / FP {raw["median_n_fp"]}）</div></div>',
    f'<div class="bignum"><div class="bn">{S["tn_total"]["tumor"]:,}<span style="font-size:13px;color:var(--mut)">/{S["tn_total"]["normal"]:,}</span></div><div class="bl">tumor/normal reads 總（中位 {S["median_nt"]}/{S["median_nn"]}）</div></div>',
    f'<div class="bignum"><div class="bn">{S["meth_overall"]["mean_beta_median"]}</div><div class="bl">中位 mean β（{S["meth_overall"]["n_with_group"]:,} 位點有群間 Δβ）</div></div>',
])
pc = S["per_chrom"]
chrom_tbl = '<table class="ct"><tr><th>chr</th>' + "".join(f'<th>{c[3:]}</th>' for c in CHRS) + '</tr>'
chrom_tbl += '<tr><td class="l">TP</td>' + "".join(f'<td>{pc[c]["tp"]:,}</td>' for c in CHRS) + '</tr>'
chrom_tbl += '<tr><td class="l">FP</td>' + "".join(f'<td>{pc[c]["fp"]:,}</td>' for c in CHRS) + '</tr>'
chrom_tbl += '<tr><td class="l">結構</td>' + "".join(f'<td>{pc[c]["struct"]:,}</td>' for c in CHRS) + '</tr></table>'
L0 = f'<div class="biggrid">{inv_cards}</div><div class="sx" style="margin:8px 0 3px">每染色體 TP / FP / 有結構(coarse≥2) 位點數：</div><div style="overflow-x:auto">{chrom_tbl}</div>'

# ================= ① 摘要(一行可展開) =================
def sumitem(head, body):
    return f'<details class="sumrow"><summary>{head}</summary><div class="sb">{body}</div></details>'
L1 = "".join([
    sumitem("① 這是<b>觀察/刻畫層</b>，不是 subclone caller、不是 TP/FP filter",
            "單樣本 HCC1395（tumor-only + SEQC2 TP/FP 標籤）探索，⭐2-3，無外部單細胞真值。所有切群/對齊/甲基結果是「可觀察的描述」，不可當偵測器或過濾器宣稱。"),
    sumitem("② <b>對齊（V_hp/V_allele）只敘述 = cis-ASM 跡象，非 subclone</b>",
            f"對齊 = 切群結果與既有 a-priori 軸（haplotype/allele）關聯（CramérV≥0.3）。這是 cis 等位特異甲基（cis-ASM）的跡象，<b>非</b>體細胞亞克隆。{S['align_status'].get('aligned',0):,} 位點對齊、{S['align_status'].get('unaligned',0):,} 未對齊。意義 baseline-dependent，扣 cis 前不可當獨立 subclone 證據。"),
    sumitem(f"③ <b>confident-multi 反判別</b>：FP {S['confident_multi_fp_pct']}% &gt; TP {S['confident_multi_tp_pct']}%",
            f"「確認多群（S1+S2）」在 FP（{S['confident_multi_fp_pct']}%）比 TP（{S['confident_multi_tp_pct']}%）更常見 → 多群結構<b>不</b>偏向真體細胞位點，反而 FP 更多結構（FP 結構率 {raw['fp_structured_pct']}% vs TP {raw['tp_structured_pct']}%）。確認多群 = cis-ASM/CN 驅動，非 subclone。"),
    sumitem("④ <b>幾何切群 = 純過切診斷</b>，非結構宣稱",
            f"幾何切群數 = scipy 0.7×max 平切（無 null 閘），會把單群位點過切。{S['geom_divergence']:,} 個位點「幾何≥2 但嚴閘判單群」= 可疑漏切的肉眼對照，不可當 subclone。"),
    sumitem("⑤ <b>S4 次閾值 = chance-level（ratio 1.15）</b>，sub-threshold candidate",
            f"S4（coarse&lt;2 但 fine≥2，僅 null90 寬閘過）TP 率與 FP 率比 ≈1.15，與噪音統計不可分 → 是 sub-threshold candidates pending validation，<b>非</b>已驗證真弱結構。S5（殘留弱結構）結構性不可達 = 0。"),
    sumitem("⑥ <b>甲基 Δβ = baseline-dependent</b>（本次新增觀察）",
            f"群間 Δβ（驅動切群的甲基差）中位 {S['meth_overall']['dbeta_group_median']}；但 <b>FP（{S['meth_by_set']['FP']['dbeta_group']}）比 TP（{S['meth_by_set']['TP']['dbeta_group']}）更大</b> → 甲基異質性同樣<b>不</b>偏向真位點，扣 cis-ASM 前非獨立 subclone 證據。T−N Δβ 中位 {S['meth_overall']['dbeta_tn_median']}（somatic 甲基偏移，全位點皆有）。"),
])

# ================= ② 整體統計(SVG + 數字+%) =================
mo = S["meth_overall"]; mc_ = S["meth_coverage"]
L2 = '<div class="statgrid">' + "".join([
    stat_card("SEQC2 CN 狀態", svg_stacked([(s, S["cn_state"][s], CNC[s]) for s in CNS]),
              numrows([(s, S["cn_state"][s], CNC[s]) for s in CNS])),
    stat_card("structure 5 態 (S1–S6)", svg_stacked([(s, S["wstate"][s], WC[s]) for s in WST]),
              numrows([(s, S["wstate"][s], WC[s]) for s in WST]), "S1 確認對齊/S2 不對齊/S3 不穩定/S4 次閾值(chance)/S6 乾淨單群"),
    stat_card("structure tier", svg_stacked([("確認", S["tier"]["1"], "#16a34a"), ("邊緣", S["tier"]["2"], "#d97706"), ("無結構", S["tier"]["3"], "#475569")]),
              numrows([("確認(tier1)", S["tier"]["1"], "#16a34a"), ("邊緣(tier2)", S["tier"]["2"], "#d97706"), ("無結構(tier3)", S["tier"]["3"], "#475569")])),
    stat_card("對齊狀態 (有結構位點)", svg_stacked([("aligned", S["align_status"].get("aligned", 0), "#22c55e"), ("unaligned", S["align_status"].get("unaligned", 0), "#f59e0b"), ("NA單群", S["align_status"].get("NA", 0), "#475569")]),
              numrows([("aligned cis-ASM", S["align_status"].get("aligned", 0), "#22c55e"), ("unaligned", S["align_status"].get("unaligned", 0), "#f59e0b"), ("NA(單群)", S["align_status"].get("NA", 0), "#475569")])),
    stat_card("tumor/normal read 組成（總）", svg_stacked([("tumor", S["tn_total"]["tumor"], "#f97316"), ("normal", S["tn_total"]["normal"], "#22c55e")]),
              numrows([("tumor", S["tn_total"]["tumor"], "#f97316"), ("normal", S["tn_total"]["normal"], "#22c55e")]), f"中位 tumor {S['median_nt']} / normal {S['median_nn']}"),
    stat_card("coarse 切群數 (null95 嚴閘)", svg_hist(S["coarse_dist"], "#db2777"), distrows(S["coarse_dist"])),
    stat_card("fine 切群數 (null90 寬閘)", svg_hist(S["fine_dist"], "#9333ea"), distrows(S["fine_dist"])),
    stat_card("幾何切群數 (0.7cut)", svg_hist(S["geom_dist"], "#0891b2"), distrows(S["geom_dist"]), "純幾何過切, 僅肉眼/null 對照, 非 subclone"),
    stat_card("判別 funnel", svg_funnel(S["funnel"]), numrows([(s["stage"], s["n"], None) for s in S["funnel"]])),
    stat_card("confident-multi: TP vs FP（反判別）", svg_compare([("TP", S["confident_multi_tp_pct"], "#60a5fa"), ("FP", S["confident_multi_fp_pct"], "#f87171")]),
              numrows([("TP confident-multi", S["wstate_x_set"]["S1"]["tp"] + S["wstate_x_set"]["S2"]["tp"], "#60a5fa"), ("FP confident-multi", S["wstate_x_set"]["S1"]["fp"] + S["wstate_x_set"]["S2"]["fp"], "#f87171")], ), "FP%&gt;TP% = 反判別, confident-multi=cis-ASM 非 subclone"),
    stat_card("甲基 mean β 分佈", svg_hist(S["mean_beta_hist"], "#0ea5e9"), distrows(S["mean_beta_hist"]), f"中位 {mo['mean_beta_median']}（β=0 未甲基→1 全甲基）"),
    stat_card("群間 |Δβ| 分佈（多群位點）", svg_hist(S["dbeta_group_hist"], "#14b8a6"), distrows(S["dbeta_group_hist"]), f"中位 {mo['dbeta_group_median']}（驅動切群的甲基差; {mo['n_with_group']:,} 位點）"),
    stat_card("CpG（甲基位點）數 / 位點", svg_hist(S["ncpg_hist"], "#38bdf8"), distrows(S["ncpg_hist"]), f"中位 {mc_['ncpg_median']} CpG/位點（{mc_['ncpg_min']}–{mc_['ncpg_max']}）｜全基因組共 {mc_['ncpg_total']:,} CpG 觀測"),
    stat_card("甲基資料覆蓋（位點數量）", "", numrows([("有甲基觀測", mc_["n_with_methylation"], "#0ea5e9"), ("有群間 Δβ（多群）", mc_["n_with_dbeta_group"], "#14b8a6"), ("有 Δβ(T−N)", mc_["n_with_dbeta_tn"], "#f97316"), ("有 Δβ(HP)", mc_["n_with_dbeta_hp"], "#c084fc")]), f"全 {S['n']:,} 位點皆有甲基；{mc_['n_with_dbeta_group']:,} 多群位點可算群間 Δβ"),
    stat_card("8 類分類（修正樹）", svg_stacked([("A", CAT8_CNT["A"], "#eab308"), ("B", CAT8_CNT["B1"] + CAT8_CNT["B2"] + CAT8_CNT["B3"], "#6b7280"), ("C", CAT8_CNT["C-S1"] + CAT8_CNT["C-S2"] + CAT8_CNT["C-S3"], "#3b82f6"), ("D", CAT8_CNT["D"], "#475569")]),
              numrows([("A subclone候選⭐", CAT8_CNT["A"], "#eab308"), ("B 無法區分(LOH主因)", CAT8_CNT["B1"] + CAT8_CNT["B2"] + CAT8_CNT["B3"], "#6b7280"), ("C-S1 cis-ASM", CAT8_CNT["C-S1"], "#3b82f6"), ("C-S2 未對齊", CAT8_CNT["C-S2"], "#60a5fa"), ("C-S3 不穩定", CAT8_CNT["C-S3"], "#93c5fd"), ("D 可測無結構", CAT8_CNT["D"], "#475569")]), "A=唯一 subclone 訊號軸; B=測不了; D=有軸但單群"),
    stat_card("cluster×label 列聯型（有結構位點）", svg_stacked([("①many:1", CTYPE_CNT["many1_結構>標籤"], "#a855f7"), ("②1:many", CTYPE_CNT["1many_跨標籤"], "#f59e0b"), ("1:1", CTYPE_CNT["one2one"], "#22c55e"), ("mixed", CTYPE_CNT["mixed"], "#6b7280")]),
              numrows([("①many:1 subclone-like", CTYPE_CNT["many1_結構>標籤"], "#a855f7"), ("②1:many 跨標籤(無ASM)", CTYPE_CNT["1many_跨標籤"], "#f59e0b"), ("1:1 對齊 cis-ASM", CTYPE_CNT["one2one"], "#22c55e"), ("mixed 複雜", CTYPE_CNT["mixed"], "#6b7280")]), "①結構>標籤=subclone訊號; ②跨標籤=無ASM/trans"),
], ) + '</div>'

# ================= ③ 分類交叉表(更多參數 + 甲基 by 狀況) =================
def xt_wstate_set():
    h = '<table><tr><th>state</th><th>TP</th><th>FP</th><th>TP%</th><th>FP%</th><th>FP/TP比</th></tr>'
    for w in WST:
        tp = S["wstate_x_set"][w]["tp"]; fp = S["wstate_x_set"][w]["fp"]
        tpp = 100 * tp / S["tp"]; fpp = 100 * fp / S["fp"]; rr = round(fpp / tpp, 2) if tpp else 0
        cls = ' class="hl"' if rr > 1.1 and w in ("S1", "S2") else ''
        h += f'<tr{cls}><td class="l"><span style="color:{WC[w]}">●</span>{w}</td><td>{tp:,}</td><td>{fp:,}</td><td>{tpp:.1f}</td><td>{fpp:.1f}</td><td>{rr}</td></tr>'
    return h + '</table>'


def xt_cn_struct():
    h = '<table><tr><th>CN狀態</th><th>多群</th><th>單群</th><th>多群%</th></tr>'
    for s in CNS:
        m = S["cn_x_structure"][s]["multi"]; sg = S["cn_x_structure"][s]["single"]; t = m + sg
        h += f'<tr><td class="l"><span style="color:{CNC[s]}">●</span>{s}</td><td>{m:,}</td><td>{sg:,}</td><td>{round(100*m/t,1) if t else 0}</td></tr>'
    return h + '</table>'


def xt_align_cn():
    h = '<table><tr><th>CN狀態</th><th>aligned</th><th>unaligned</th><th>NA</th><th>aligned%</th></tr>'
    for s in CNS:
        a = S["align_x_cn"][s]; t = sum(a.values())
        h += f'<tr><td class="l"><span style="color:{CNC[s]}">●</span>{s}</td><td>{a["aligned"]:,}</td><td>{a["unaligned"]:,}</td><td>{a["NA"]:,}</td><td>{round(100*a["aligned"]/t,1) if t else 0}</td></tr>'
    return h + '</table>'


def xt_coarse_fine():
    h = '<table><tr><th>coarse＼fine</th>' + "".join(f'<th>{f}{"+" if f==5 else ""}</th>' for f in range(1, 6)) + '</tr>'
    for c in range(1, 6):
        h += f'<tr><td class="l">{c}{"+" if c==5 else ""}</td>' + "".join(f'<td>{S["coarse_x_fine"][str(c)][str(f)]:,}</td>' for f in range(1, 6)) + '</tr>'
    return h + '</table>'


def xt_flags():
    h = '<table><tr><th>flag</th><th>TP%</th><th>FP%</th><th>FP/TP</th></tr>'
    rows = [("unstable(seed分歧)", S["flag_x_set"]["unstable"]), ("hidden_het", S["flag_x_set"]["hidden_het"]),
            ("可疑漏切(geom_div)", S["flag_x_set"]["geom_divergence"]), ("有 other 殘群", S["other_rate"])]
    for lab, d in rows:
        rr = round(d["fp"] / d["tp"], 2) if d["tp"] else 0
        cls = ' class="hl"' if rr > 1.3 else ''
        h += f'<tr{cls}><td class="l">{lab}</td><td>{d["tp"]}</td><td>{d["fp"]}</td><td>{rr}</td></tr>'
    return h + '</table>'


def xt_vaxis():
    h = '<table><tr><th>主導軸(有結構)</th><th>TP</th><th>FP</th></tr>'
    for k, lab in [("hp", "HP 主導"), ("allele", "ALT/REF 主導"), ("weak", "皆弱(<0.3)")]:
        h += f'<tr><td class="l">{lab}</td><td>{S["vaxis_x_set"]["tp"][k]:,}</td><td>{S["vaxis_x_set"]["fp"][k]:,}</td></tr>'
    return h + '</table>'


def xt_meth_situation():
    h = '<table><tr><th>狀況</th><th>中位 mean β</th><th>中位 Δβ群</th><th>中位 Δβ(T−N)</th></tr>'
    h += '<tr><td class="l" colspan="4" style="background:#0b1222;color:var(--mut)">— by structure 態 —</td></tr>'
    for w in WST:
        d = S["meth_by_wstate"][w]
        h += f'<tr><td class="l"><span style="color:{WC[w]}">●</span>{w}</td><td>{d["mean_beta"]}</td><td>{d["dbeta_group"] if d["dbeta_group"] is not None else "—"}</td><td>{d["dbeta_tn"]}</td></tr>'
    h += '<tr><td class="l" colspan="4" style="background:#0b1222;color:var(--mut)">— by CN 狀態 —</td></tr>'
    for s in CNS:
        d = S["meth_by_cn"][s]
        h += f'<tr><td class="l"><span style="color:{CNC[s]}">●</span>{s}</td><td>{d["mean_beta"]}</td><td>{d["dbeta_group"] if d["dbeta_group"] is not None else "—"}</td><td>{d["dbeta_tn"]}</td></tr>'
    h += '<tr><td class="l" colspan="4" style="background:#0b1222;color:var(--mut)">— by set —</td></tr>'
    for st in ("TP", "FP"):
        d = S["meth_by_set"][st]
        h += f'<tr class="{"hl" if st=="FP" else ""}"><td class="l">{st}</td><td>{d["mean_beta"]}</td><td>{d["dbeta_group"]}</td><td>{d["dbeta_tn"]}</td></tr>'
    return h + '</table>'


def xbox(title, html, note=""):
    return f'<div class="xbox"><b>{title}</b>{html}{f"<div class=sx>{note}</div>" if note else ""}</div>'


L3 = '<div class="xgrid">' + "".join([
    xbox("判別主表: structure 態 × TP/FP", xt_wstate_set(), "高亮=confident-multi 在 FP 更多(反判別)"),
    xbox("CN 狀態 × 結構", xt_cn_struct()),
    xbox("對齊 × CN 狀態", xt_align_cn(), "LOH-unmask: loh 高 aligned% 是 confound 通道"),
    xbox("coarse × fine 切群數(過切診斷)", xt_coarse_fine(), "對角線=一致; 右上=fine 過切"),
    xbox("噪音 flag × TP/FP", xt_flags(), "FP/TP>1.3 高亮=FP 偏向(噪音標記)"),
    xbox("主導軸(V_hp vs V_allele) × set", xt_vaxis()),
    xbox("甲基(β) by 狀況", xt_meth_situation(), "Δβ群: FP>TP=甲基異質性非判別; baseline-dependent"),
], ) + '</div>'

# ====== ④b 差異甲基定位點（aggregate） ======
def _qax():
    h = '<table><tr><th>軸</th><th>有差異位點</th><th>差異CpG總數</th><th>hyper%(增甲基)</th></tr>'
    for ax, d in QCLS.get("axes", {}).items():
        h += f'<tr><td class="l">{ax.replace("_vs_","/")}</td><td>{d["n_loci"]:,}</td><td>{d["n_sig"]:,}</td><td>{d["hyper_pct"]}</td></tr>'
    return h + '</table>'
def _qattr():
    h = '<table><tr><th>軸/狀況</th><th>TP</th><th>FP</th><th>FP/TP</th><th>%位點</th><th>研究可用</th></tr>'
    for k, d in QATTR.items():
        use = '✅TP專一' if d["fp_tp_ratio"] < 0.6 else ('⚠不專一' if d["fp_tp_ratio"] <= 1.2 else '✗FP富集')
        cls = ' class="hl"' if d["fp_tp_ratio"] < 0.6 else ''
        h += f'<tr{cls}><td class="l">{k}</td><td>{d["TP"]:,}</td><td>{d["FP"]:,}</td><td>{d["fp_tp_ratio"]}</td><td>{d["pct_loci"]}</td><td>{use}</td></tr>'
    return h + '</table>'
_ut = sum(1 for v in QPL.values() if v.get("set") == "TP" and v.get("union", 0) >= 3); _uf = sum(1 for v in QPL.values() if v.get("set") == "FP" and v.get("union", 0) >= 3)
L4b = (f'<div class="note">⚠ <b>「有差異甲基」普遍({100*(_ut+_uf)/max(1,len(QPL)):.0f}% 有 ≥3 差異 CpG)→ 不判別</b>；真正 somatic 定位點 = <b>HP1/HP1-1 軸(TP 專一)</b>；subclone-marker FP 富集非 somatic-specific。完整 9 張圖（直方/分類/CDF/歸因）: '
       f'<a href="20260624_differential_methylation_marker_charts.standalone.html" style="color:#7fb4e8">→ 差異甲基定位點圖庫</a></div>'
       f'<div class="xgrid"><div class="xbox"><b>各軸差異 CpG（數量 + 方向）</b>{_qax()}<div class="sx">差異 CpG 總 {QCLS.get("total_sig_cpg",0):,}；somatic 軸(T/N・HP-1) hyper 偏高=somatic 多為增甲基</div></div>'
       f'<div class="xbox"><b>🔴 FP/TP 歸因 + 研究可用性</b>{_qattr()}<div class="sx">FP/TP&lt;0.6=TP 專一(可做研究)；&gt;1.2=FP 富集(非 somatic)。機制見方法文件 §7</div></div></div>')

SECTIONS = f"""
<details open><summary>⓪ 原始數據盤點（TP/FP 總量・覆蓋・per-chrom）</summary>{L0}</details>
<details open><summary>① 摘要 + 紅線（每行可展開）</summary>{L1}</details>
<details open><summary>② 整體統計觀察（SVG ＋ 數字＋比例）</summary>{L2}</details>
<details><summary>③ 分類交叉表（切群/分類各參數 ＋ 甲基 by 狀況）</summary>{L3}</details>
<details><summary>④b 差異甲基定位點（site-first 量化 ＋ FP/TP 歸因 ＋ 研究可用性）</summary>{L4b}</details>
"""

TPL = r"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>全位點觀察分類儀表板 v4</title><style>
:root{--bg:#0a0e16;--card:#111826;--fg:#e7edf5;--mut:#8b9bb4;--line:#2a3344;--ac:#D97757}
*{box-sizing:border-box}body{font-family:system-ui,"Noto Sans CJK TC",sans-serif;background:var(--bg);color:var(--fg);margin:0;line-height:1.5}
.wrap{max-width:1480px;margin:0 auto;padding:14px}
h1{font-size:19px;margin:2px 0}.sub{color:var(--mut);font-size:12.5px}
details{background:var(--card);border:1px solid var(--line);border-radius:8px;margin:7px 0;padding:3px 12px}
summary{cursor:pointer;font-weight:600;font-size:13.5px;padding:6px 0}
.sumrow{background:#0d1422;margin:5px 0;padding:2px 10px}.sumrow summary{font-weight:500;font-size:12.5px}.sumrow .sb{color:var(--mut);font-size:11.5px;padding:4px 0 7px;line-height:1.6}.sumrow b{color:#fbbf24}
.biggrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(210px,1fr));gap:9px;margin:6px 0}
.bignum{background:#0b1222;border:1px solid var(--line);border-radius:8px;padding:9px 11px}.bn{font-size:22px;font-weight:700;font-variant-numeric:tabular-nums}.bl{color:var(--mut);font-size:11px;margin-top:2px}.bl b{color:var(--fg)}
.statgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(290px,1fr));gap:10px;margin:8px 0}
.stcard{background:#0b1222;border:1px solid var(--line);border-radius:7px;padding:8px 10px}.stt{font-size:12px;font-weight:600;margin-bottom:5px}.sx{color:var(--mut);font-size:10.5px;margin-top:4px}
.nrows{margin-top:5px;font-variant-numeric:tabular-nums}.nr{display:flex;justify-content:space-between;gap:8px;font-size:11px;padding:1px 0;border-bottom:1px dotted #1c2536}.nrl{color:var(--fg)}.nrv{color:var(--mut);text-align:right;min-width:54px}.nrp{color:#7fb4e8;text-align:right;min-width:46px}
.xgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(330px,1fr));gap:12px}.xbox{background:#0b1222;border:1px solid var(--line);border-radius:7px;padding:8px 10px}.xbox b{font-size:12px}
table{border-collapse:collapse;font-size:11.5px;margin:4px 0;width:100%;font-variant-numeric:tabular-nums}td,th{border:1px solid var(--line);padding:3px 7px;text-align:right}th{background:#1a2130}td.l{text-align:left}tr.hl td{background:#10241a;font-weight:600}
table.ct{font-size:10px}table.ct td,table.ct th{padding:2px 4px}
.bar{display:flex;flex-wrap:wrap;gap:7px;align-items:center;background:var(--card);border:1px solid var(--line);border-radius:8px;padding:9px;margin:8px 0;position:sticky;top:0;z-index:9}
.bar label{font-size:11.5px;color:var(--mut)}input,select{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:5px;padding:4px 6px;font-size:12px}
input[type=search]{width:150px}.btn{background:#1a2130;border:1px solid var(--line);border-radius:5px;color:var(--fg);padding:5px 9px;cursor:pointer;font-size:12px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:10px;margin:9px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden}.card.j-in{border-color:#16a34a}.card.j-mb{border-color:#d97706}.card.j-ex{border-color:#dc2626}
.hd{padding:8px 10px}.ttl{font-weight:700;font-size:13px;cursor:pointer}.badges{margin-top:3px;display:flex;flex-wrap:wrap;gap:4px}
.badge{font-size:10.5px;padding:1px 6px;border-radius:9px;border:1px solid var(--line)}
.b-tp{background:#10243a;color:#7cc4ff}.b-fp{background:#3a1a1a;color:#ffa3a3}.b-gain{background:#3a1414;color:#ffb3b3;border-color:#ef4444}.b-loss{background:#10233a;color:#9fc6ff;border-color:#3b82f6}.b-loh{background:#241038;color:#d6b3ff;border-color:#a855f7}.b-neu{background:#1a2130;color:var(--mut)}
.b-al{background:#13261c;color:#9ff0c4}.b-tn{background:#0f2a24;color:#7fe3c4}.b-div{background:#3a2418;color:#ffba80;border-color:#d97706}.b-meth{background:#06262e;color:#5fd6c8;border-color:#14b8a6}
.b-t1{background:#13261c;color:#86efac;border-color:#16a34a}.b-t2{background:#2a2410;color:#fcd34d;border-color:#d97706}.b-t3{background:#1a2130;color:var(--mut)}
.b-sc{background:#3a2e0a;color:#fcd34d;border-color:#eab308;font-weight:700}.b-un{background:#16161e;color:#6b7280;border-color:#33384a}.b-cstr{background:#101e2e;color:#8fb8e0;border-color:#2a4a6a}
.desc{color:var(--mut);font-size:11px;margin-top:3px}
.figs{cursor:pointer;background:#0b1222;border-top:1px solid var(--line)}.figs img{width:100%;display:block}.nofig{padding:18px;text-align:center;color:var(--mut);font-size:11px}
.exp{padding:0;border-top:1px solid var(--line)}.exp summary{font-size:11.5px;padding:6px 10px;color:var(--mut)}
.kvg{display:grid;grid-template-columns:auto 1fr;gap:2px 10px;padding:4px 12px 9px;font-size:11.5px}.kvg .k{color:var(--mut)}.kvg .v{font-weight:600;font-variant-numeric:tabular-nums}.kvg .h{grid-column:1/3;color:#7fb4e8;font-size:10.5px;font-weight:700;margin-top:5px;border-bottom:1px solid var(--line);padding-bottom:1px}
.jrow{display:flex;gap:5px;padding:7px 10px;border-top:1px solid var(--line)}.jb{flex:1;padding:4px;border:1px solid var(--line);border-radius:5px;background:#0b1222;color:var(--fg);cursor:pointer;font-size:11px}.jb.on-in{background:#14532d}.jb.on-mb{background:#3a2f10}.jb.on-ex{background:#3a1414}
.modal{display:none;position:fixed;inset:0;background:#000a;z-index:50;overflow:auto;padding:20px}.mc{max-width:1180px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:10px;padding:14px}.mc img{width:100%;border-radius:6px;background:#fff}
.pager{text-align:center;margin:10px 0}.pager button{padding:5px 12px}
</style></head><body><div class="wrap">
<h1>全位點觀察分類儀表板 v4 — 原始盤點 + 數字比例 + 完整位點資訊 + 甲基觀察</h1>
<div class="sub">本 session full run（演化樹 tree-aware）+ SEQC2 CN + tumor/normal + 切群數(null95/null90/幾何) + 甲基(β)。路徑式, 需與同層 figs/ 一起開。__STAMP__</div>
__SECTIONS__
<div class="bar">
<label>🔍 <input type="search" id="q" placeholder="chr1:39078224 或 chr8"></label>
<label>判別 <select id="fws"><option value="">全</option><option value="0">S1 確認對齊</option><option value="1">S2 不對齊</option><option value="2">S3 不穩定</option><option value="3">S4 次閾值</option><option value="4">S6 乾淨單群</option></select></label>
<label>tier <select id="ftier"><option value="">全</option><option value="1">1 確認</option><option value="2">2 邊緣</option><option value="3">3 無</option></select></label>
<label>set <select id="fset"><option value="">全</option><option value="0">TP</option><option value="1">FP</option></select></label>
<label>分類 <select id="fcat"><option value="">全8類</option><option value="0">A subclone候選</option><option value="B">B 無法區分(全)</option><option value="1">B1 LOH</option><option value="2">B2 低覆蓋</option><option value="3">B3 非LOH</option><option value="4">C-S1 cis-ASM</option><option value="5">C-S2 未對齊</option><option value="6">C-S3 不穩定</option><option value="7">D 可測無結構</option></select></label>
<label>列聯型 <select id="fct"><option value="">全</option><option value="1">①many:1(subclone-like)</option><option value="2">②1:many(跨標籤)</option><option value="0">1:1對齊</option><option value="3">mixed</option></select></label>
<label><input type="checkbox" id="fsc"> ⭐subclone候選</label>
<label><input type="checkbox" id="flhs"> LOH+somatic結構</label>
<label>逐CpG驗 <select id="fpcv"><option value="">全</option><option value="2">⭐subclone_novel</option><option value="1">cis-ASM(cluster=label)</option><option value="3">no_sigCpG</option><option value="4">weak</option><option value="0">no_cluster</option></select></label>
<label>CN <select id="fcn"><option value="">全</option><option value="0">gain</option><option value="1">loss</option><option value="2">loh</option><option value="3">neutral</option></select></label>
<label>CN值 <input type="number" id="cnmin" placeholder="min" style="width:50px" step="0.5">-<input type="number" id="cnmax" placeholder="max" style="width:50px" step="0.5"></label>
<label><input type="checkbox" id="floh"> LOH</label>
<label>coarse <input type="number" id="gmin" value="1" style="width:40px">-<input type="number" id="gmax" value="9" style="width:40px"></label>
<label>n≥ <input type="number" id="nmin" placeholder="0" style="width:50px"></label>
<label>Δβ群≥ <input type="number" id="dbmin" placeholder="0" style="width:50px" step="0.05"></label>
<label><input type="checkbox" id="funs"> unstable</label>
<label><input type="checkbox" id="ffc"> fine&gt;coarse</label>
<label><input type="checkbox" id="fdiv"> 可疑漏切</label>
<label>排序 <select id="sort"><option value="1">位置</option><option value="3">n</option><option value="4">coarse</option><option value="17">幾何</option><option value="12">CN值</option><option value="22">mean β</option><option value="23">Δβ群</option><option value="vmax">對齊V</option></select><button class="btn" id="sdir">▼</button></label>
<span class="sub" id="cnt"></span><button class="btn" id="ej">匯出判讀</button>
</div>
<div id="grid" class="grid"></div><div class="pager" id="pager"></div>
<details><summary>⑤ 參考 / 方法 / provenance</summary><div class="sub" style="padding:6px 0">
每位點兩張圖（點卡片開 modal）：① <b>tumor phylo</b> = 左 UPGMA 樹(tree-aware 群色) ｜ 中 甲基 read×CpG(RdBu_r) ｜ 右 read×read 距離(magma, 暗=近)，側欄 HP|ALT（phylo 切群 tumor-only, 故無 T/N 側欄）；② <b>tumor-vs-normal 甲基對照</b> = 上 tumor / 下 normal 甲基熱圖(RdBu_r)，側欄 T/N(橘 tumor/綠 normal)|HP，看 tumor 甲基是否異於 normal。<br>
判別 = weak-structure 5 態(S1 確認多群對齊=cis-ASM 候選 / S2 不對齊 / S3 不穩定 / S4 次閾值 chance-level / S6 乾淨單群)。對齊=cis-ASM 跡象只敘述不歸類。切群數三法: null95(嚴主判)/null90(寬候選)/幾何0.7cut(肉眼過切)。<br>
甲基欄: mean β=全 read×CpG 平均(0 未甲基→1 全甲基); Δβ群=coarse leaf 群間 mean β max−min(驅動切群的甲基差); Δβ(T−N)=tumor−normal; Δβ(HP)=HP1−HP2; hypo β&lt;0.3 / hyper β&gt;0.7 佔比。<b>Δβ baseline-dependent, 扣 cis-ASM 前非獨立 subclone 證據</b>。<br>
8 類分類(修正決策樹): A subclone候選(單標籤+多結構)⭐ / B 無法區分(LOH主因, 測不了) / C-S1 cis-ASM / C-S2 未對齊 / C-S3 不穩定 / D 可測無結構。列聯型: ①many:1(結構>標籤=subclone訊號) / ②1:many(跨標籤=無ASM/trans) / 1:1(對齊) / mixed。⭐subclone候選=A或①many:1×somatic; LOH-som=LOH+somatic軸+結構。完整定義+各類數量比例: docs/methodology/20260624_subclone_classification_decision_tree_corrected_01.md。<br>
逐CpG 判別(site-first MWU + cluster×label CpG-set overlap): ⭐subclone_novel(523, coherent CpG 不被標籤解釋, 但 FP 富集 1.78× 非 somatic-specific) / cis-ASM(cluster=label CpG) / no_sigCpG / weak。HP label PERMANOVA(significance.json, ✓sig/⚠sig+disp)。完整方法: docs/methodology/20260624_per_cpg_differential_and_subclone_validation_01.md。<br>
來源: phylo_cpp_wg_full_records_v6.json + dashboard_stats_v4.json + phylo_cpp_wg_full_percpg.json + label_composition.json + contingency_type.json + phylo_cpp_wg_full_methylation.json + significance.json(label_structure) + SEQC2 CNV(hg38) + reads.tsv。§13 數字注入不手打。單樣本 HCC1395 ⭐2-3。
</div></details>
</div>
<div class="modal" id="modal" onclick="if(event.target===this)closeM()"><div class="mc" id="mc"></div></div>
<script>
const D=__D__, CHRS=__CHRS__, CNS=['gain','loss','loh','neutral'], WST=['S1','S2','S3','S4','S6'];
const C={ck:0,ps:1,set:2,n:3,cg:4,fg:5,oth:6,uns:7,hh:8,vhp:9,val:10,cns:11,cnv:12,loh:13,fig:14,nt:15,nn:16,geom:17,wst:18,tier:19,gdiv:20,als:21,mb:22,dbg:23,dtn:24,dhp:25,ncpg:26,cat:27,ct:28,sc:29,lhs:30,pcv:31,lshp:32};
const CAT8=['A','B1','B2','B3','C-S1','C-S2','C-S3','D'],CTYPE=['1:1','①many:1','②1:many','mixed','單群'];
const PCV=['no_cluster','cis-ASM','⭐subclone_novel','no_sigCpG','weak_scattered','no_meth'],LSL=['—','✓sig','⚠sig+disp'];
const $=id=>document.getElementById(id);
const figName=r=>'cpp_'+CHRS[r[C.ck]]+'_'+(r[C.ps]-5000)+'_'+(r[C.ps]+5000)+'.png';
const tnFigName=r=>'cpp_'+CHRS[r[C.ck]]+'_'+(r[C.ps]-5000)+'_'+(r[C.ps]+5000)+'_tn.png';
const key=r=>CHRS[r[C.ck]]+':'+r[C.ps];
const nv=v=>v<=-9?'—':v;  // -9 sentinel
let J=JSON.parse(localStorage.getItem('phylo_v4_judge')||'{}'),page=0,sdir=-1,PER=24;
function filt(){let a=D.slice();
 const q=$('q').value.trim().toLowerCase();
 if(q)a=a.filter(r=>(CHRS[r[C.ck]]+':'+r[C.ps]).toLowerCase().includes(q)||CHRS[r[C.ck]].toLowerCase()===q);
 const fw=$('fws').value;if(fw!=='')a=a.filter(r=>r[C.wst]==fw);
 const ft=$('ftier').value;if(ft!=='')a=a.filter(r=>r[C.tier]==ft);
 const fs=$('fset').value;if(fs!=='')a=a.filter(r=>r[C.set]==fs);
 const fcat=$('fcat').value;if(fcat==='B')a=a.filter(r=>r[C.cat]>=1&&r[C.cat]<=3);else if(fcat!=='')a=a.filter(r=>r[C.cat]==fcat);
 const fct=$('fct').value;if(fct!=='')a=a.filter(r=>r[C.ct]==fct);
 if($('fsc').checked)a=a.filter(r=>r[C.sc]);if($('flhs').checked)a=a.filter(r=>r[C.lhs]);
 const fpcv=$('fpcv').value;if(fpcv!=='')a=a.filter(r=>r[C.pcv]==fpcv);
 const fc=$('fcn').value;if(fc!=='')a=a.filter(r=>r[C.cns]==fc);
 if($('cnmin').value!=='')a=a.filter(r=>r[C.cnv]>=0&&r[C.cnv]>=+$('cnmin').value);if($('cnmax').value!=='')a=a.filter(r=>r[C.cnv]>=0&&r[C.cnv]<=+$('cnmax').value);
 if($('floh').checked)a=a.filter(r=>r[C.loh]);
 a=a.filter(r=>r[C.cg]>=+$('gmin').value&&r[C.cg]<=+$('gmax').value);
 if($('nmin').value!=='')a=a.filter(r=>r[C.n]>=+$('nmin').value);
 if($('dbmin').value!=='')a=a.filter(r=>r[C.dbg]>=0&&r[C.dbg]>=+$('dbmin').value);
 if($('funs').checked)a=a.filter(r=>r[C.uns]);if($('ffc').checked)a=a.filter(r=>r[C.fg]>r[C.cg]);if($('fdiv').checked)a=a.filter(r=>r[C.gdiv]);
 const sk=$('sort').value;a.sort((x,y)=>{const f=v=>sk==='vmax'?Math.max(v[C.vhp],v[C.val]):(sk==1?v[C.ck]*1e10+v[C.ps]:v[+sk]);return sdir*(f(x)-f(y));});
 return a;}
function cnBadge(r){const st=CNS[r[C.cns]];const cls={0:'b-gain',1:'b-loss',2:'b-loh',3:'b-neu'}[r[C.cns]];return `<span class="badge ${cls}">${st==='loh'?'LOH':st+(r[C.cnv]>=0?' CN'+r[C.cnv]:'')}</span>`;}
function tierBadge(r){const t=r[C.tier];return `<span class="badge b-t${t}">${WST[r[C.wst]]}·tier${t}</span>`;}
function methBadge(r){return r[C.dbg]>=0?`<span class="badge b-meth">Δβ群 ${r[C.dbg]}</span>`:'';}
function catBadge(r){const c=CAT8[r[C.cat]];const cls=c==='A'?'b-sc':(c[0]==='B'?'b-un':(c[0]==='C'?'b-cstr':'b-neu'));const ct=r[C.ct]<4?' '+CTYPE[r[C.ct]]:'';return `<span class="badge ${cls}">${c}${r[C.sc]?'⭐':''}${ct}</span>`+(r[C.lhs]?'<span class="badge b-loh">LOH-som</span>':'');}
function pcvBadge(r){const v=r[C.pcv];if(v===0||v===5)return '';return `<span class="badge ${v===2?'b-sc':(v===1?'b-cstr':'b-un')}">逐CpG ${PCV[v]}</span>`;}
function alignDesc(r){const vh=r[C.vhp],va=r[C.val];return Math.max(vh,va)>=0.3?`對齊 ${va>=vh?'allele':'hp'} 軸 (V_hp ${vh}/V_al ${va}) — cis-ASM 跡象，非 subclone`:`未對齊 germline 軸 (V_hp ${vh}/V_al ${va})`;}
function structText(r){return {0:'確認多群·對齊 germline（cis-ASM 候選，非 subclone）',1:'確認多群·不對齊（結構但無標籤對應）',2:'不穩定多群（seed 分歧 borderline）',3:'次閾值候選（僅 null90 過，chance-level）',4:'乾淨單群（無監督切不出≥2）'}[r[C.wst]];}
function clusterText(r){const cg=r[C.cg],fg=r[C.fg],g=r[C.geom];const v=(g>=2&&cg<2)?`<span style="color:#ffba80">幾何 ${g} 群但嚴閘判單群=可疑漏切</span>`:(fg>cg?`寬閘多切（${cg}→${fg}）`:'各法一致');return `null95 <b>${cg}</b> ｜ null90 <b>${fg}</b> ｜ 幾何 <b>${g>=0?g:'?'}</b> → ${v}`;}
function kvg(r){return `<div class="kvg">
  <div class="h">分類（修正決策樹）</div>
  <span class="k">8 類</span><span class="v">${CAT8[r[C.cat]]}${r[C.sc]?' ⭐subclone候選':''}${r[C.lhs]?' · LOH-somatic':''}</span><span class="k">列聯型</span><span class="v">${r[C.ct]<4?CTYPE[r[C.ct]]:'單群'}</span>
  <div class="h">甲基 (β)</div>
  <span class="k">mean β</span><span class="v">${nv(r[C.mb])}</span><span class="k">n_CpG (甲基位點數)</span><span class="v">${r[C.ncpg]>=0?r[C.ncpg]:'—'}</span>
  <span class="k">Δβ 群間(切群驅動)</span><span class="v">${nv(r[C.dbg])}</span><span class="k">Δβ (T−N)</span><span class="v">${nv(r[C.dtn])}</span>
  <span class="k">Δβ (HP1−HP2)</span><span class="v">${nv(r[C.dhp])}</span><span class="k" style="color:#6b7a90">註</span><span class="v" style="font-weight:400;color:var(--mut)">Δβ baseline-dep·hypo/hyper見熱圖</span>
  <div class="h">結構 / 切群</div>
  <span class="k">歸類</span><span class="v" style="font-weight:500">${structText(r)}</span><span class="k">coarse/fine/幾何</span><span class="v">${r[C.cg]}/${r[C.fg]}/${r[C.geom]>=0?r[C.geom]:'?'}</span>
  <span class="k">unstable</span><span class="v">${r[C.uns]?'是':'否'}</span><span class="k">hidden-het / other</span><span class="v">${r[C.hh]?'是':'否'} / ${r[C.oth]}</span>
  <div class="h">對齊 / label PERMANOVA / 逐CpG (cis-ASM 跡象, 非 subclone)</div>
  <span class="k">V_hp / V_allele</span><span class="v">${r[C.vhp]} / ${r[C.val]}</span><span class="k">主導軸</span><span class="v">${Math.max(r[C.vhp],r[C.val])>=0.3?(r[C.val]>=r[C.vhp]?'allele':'hp'):'皆弱'}</span>
  <span class="k">HP label PERMANOVA</span><span class="v">${LSL[r[C.lshp]]}</span><span class="k">逐CpG 判別</span><span class="v" style="font-weight:500">${PCV[r[C.pcv]]}</span>
  <div class="h">覆蓋 / CN</div>
  <span class="k">reads (T/N)</span><span class="v">${r[C.n]} (T${r[C.nt]>=0?r[C.nt]:'?'}/N${r[C.nn]>=0?r[C.nn]:'?'})</span><span class="k">SEQC2 CN</span><span class="v">${CNS[r[C.cns]]}${r[C.cnv]>=0?' '+r[C.cnv]:''}${r[C.loh]?' LOH':''}</span>
 </div>`;}
function card(r){const k=key(r),fn=r[C.fig]?figName(r):null;const j=(J[k]||{}).c;
 return `<div class="card${j?' j-'+j:''}"><div class="hd"><div class="ttl" onclick="openM('${k}')">${k}</div>
  <div class="badges">${r[C.set]===0?'<span class="badge b-tp">TP</span>':'<span class="badge b-fp">FP</span>'}${catBadge(r)}${pcvBadge(r)}${tierBadge(r)}${cnBadge(r)}<span class="badge b-tn">T${r[C.nt]>=0?r[C.nt]:'?'}/N${r[C.nn]>=0?r[C.nn]:'?'}</span><span class="badge ${r[C.gdiv]?'b-div':'b-al'}">切群 ${r[C.cg]}/${r[C.fg]}/${r[C.geom]>=0?r[C.geom]:'?'}${r[C.gdiv]?' ⚠':''}</span>${methBadge(r)}${r[C.loh]?'<span class="badge b-loh">LOH</span>':''}</div>
  <div class="desc">${alignDesc(r)}</div></div>
  ${fn?`<div class="figs" onclick="openM('${k}')"><img loading="lazy" src="figs/${fn}"></div>`:'<div class="nofig">無圖（reads 過少）</div>'}
  <details class="exp"><summary>完整觀察數據（甲基・結構・對齊・覆蓋・CN）</summary>${kvg(r)}</details>
  <div class="jrow"><button class="jb${j==='in'?' on-in':''}" onclick="setJ('${k}','in')">應含</button><button class="jb${j==='mb'?' on-mb':''}" onclick="setJ('${k}','mb')">可能漏</button><button class="jb${j==='ex'?' on-ex':''}" onclick="setJ('${k}','ex')">排除</button></div></div>`;}
function draw(){const a=filt();const pages=Math.max(1,Math.ceil(a.length/PER));page=Math.min(page,pages-1);
 $('cnt').textContent=a.length.toLocaleString()+' 個';
 $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
 $('pager').innerHTML=`<button class="btn" ${page<=0?'disabled':''} onclick="page--;draw()">‹</button> 第 ${page+1}/${pages} 頁 <button class="btn" ${page>=pages-1?'disabled':''} onclick="page++;draw()">›</button>`;}
['q','fws','ftier','fset','fcat','fct','fsc','flhs','fpcv','fcn','cnmin','cnmax','floh','gmin','gmax','nmin','dbmin','funs','ffc','fdiv','sort'].forEach(id=>$(id).addEventListener('input',()=>{page=0;draw();}));
$('sdir').onclick=()=>{sdir*=-1;$('sdir').textContent=sdir===-1?'▼':'▲';page=0;draw();};
window.setJ=(k,c)=>{J[k]=J[k]||{};J[k].c=(J[k].c===c?null:c);localStorage.setItem('phylo_v4_judge',JSON.stringify(J));draw();};
window.openM=k=>{const r=D.find(x=>key(x)===k);if(!r)return;const fn=r[C.fig]?figName(r):null;
 $('mc').innerHTML=`<div style="display:flex;justify-content:space-between"><h2 style="margin:0">${k} ${catBadge(r)} ${tierBadge(r)} ${methBadge(r)}</h2><button class="btn" onclick="closeM()">關閉 ✕</button></div>
  <div class="desc" style="margin:6px 0">${alignDesc(r)}</div>
  <div class="sub" style="margin:6px 0 2px">① tumor phylo（UPGMA 樹·甲基·距離, tumor-only）：</div>${fn?`<img src="figs/${fn}">`:'<div class="nofig">無圖</div>'}
  ${r[C.dtn]>-9?`<div class="sub" style="margin:11px 0 2px">② tumor-vs-normal 甲基對照（上 tumor / 下 normal, RdBu_r）：Δβ(T−N) ${r[C.dtn]}</div><img src="figs_tn/${tnFigName(r)}" loading="lazy" onerror="this.outerHTML='<div class=nofig>T/N 圖未生成（reads 過少）</div>'">`:''}
  ${kvg(r)}
  <div class="sub" style="margin-top:8px">距離對角塊乾淨+coarse 側欄對齊 HP/ALT=真分群。對齊/甲基 Δβ=cis-ASM 跡象非 subclone（baseline-dependent）。T/N 對照看 tumor 甲基是否異於 normal。</div>`;$('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});
$('ej').onclick=()=>{const o=Object.entries(J).filter(([k,v])=>v.c).map(([k,v])=>({locus:k,choice:v.c}));const b=new Blob([JSON.stringify({n:o.length,judgments:o},null,1)],{type:'application/json'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='phylo_v4_judgments.json';a.click();};
draw();
</script></body></html>"""

out = (TPL.replace("__D__", json.dumps(D, separators=(",", ":")))
          .replace("__CHRS__", json.dumps(CHRS))
          .replace("__SECTIONS__", SECTIONS)
          .replace("__STAMP__", "2026-06-24 v4 · 8類分類+列聯型 · SEQC2 CN + T/N + 甲基β · single-sample ⭐2-3"))
outp = f"{OUTD}/20260623_phylo_cpp_observation_dashboard_v4.html"
open(outp, "w").write(out)
print(f"WROTE {outp} ({len(out)//1024} KB) | {len(D)} loci")
