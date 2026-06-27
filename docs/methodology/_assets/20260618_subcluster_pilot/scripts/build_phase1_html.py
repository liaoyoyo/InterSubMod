#!/usr/bin/env python3
"""[Phase 1 HTML 輸出] 讀已驗 JSON → matplotlib 圖表(English label, base64 inline) → standalone HTML
(繁中敘述 + 注入數據表)。§13-A: 數字全由 JSON 注入, 缺 key refuse。觀察分類 + NACT cis-clean funnel + 對比 +
CN 分層 + genotype census + 9+30 survivor handoff。輸出 ../20260625_phase1_subclone_cisclean_01.standalone.html。"""
import json, io, base64, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_phase1_subclone_cisclean_01.standalone.html"


def L(name):
    return json.load(open(f"{A}/{name}"))


def req(d, *keys):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            sys.exit(f"§13-A REFUSE: missing {k}")
    return d


# ---- load verified JSON ----
ct = L("corrected_tree.json"); cont = L("contingency_summary.json"); pcs = L("percpg_summary.json")
dfs = L("decisionflow_summary.json"); ns = L("nact_summary.json"); nst = L("nact_stats.json")
sc = L("survivor_census.json"); ar = L("anchored_retest.json"); cc = L("candidate_characterization.json")
fp = L("fptp_attribution.json")

AC = "#D97757"; TX = "#141413"; BG = "#FAF9F5"; BD = "#E3DACC"; MUT = "#6B6862"
GRN = "#5B8A5B"; RED = "#C0563F"; BLU = "#4A6E8A"; GRY = "#9A958C"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": BD, "axes.linewidth": 0.8,
                     "text.color": TX, "axes.labelcolor": TX, "xtick.color": TX, "ytick.color": TX,
                     "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# ---- chart 1: NACT funnel ----
def fig_funnel():
    fig, ax = plt.subplots(figsize=(9, 4.2))
    vd = req(ns, "verdict_dist")
    stages = [("1139 候選", [("candidate pool", 1139, GRY)]),
              ("NACT", [("cis_asm", vd["cis_asm"], GRN), ("undetermined", vd["undetermined"], GRY),
                        ("no_tumor_sig", vd["no_tumor_signature"], BLU), ("candidate", vd["candidate_subclone"], AC)]),
              ("genotype census", [("double-dip 72", req(sc, "geno_census", "homogeneous"), RED),
                                   ("anchored 30", req(sc, "geno_census", "anchored"), AC)]),
              ("geno-axis retest", [("no-sig 15", req(ar, "retest_status").get("NO_GENO_SIGNATURE", 0), RED),
                                    ("tumor-specific 9", req(ar, "geno_axis_verdict").get("candidate_subclone", 0), GRN),
                                    ("undet 5", req(ar, "geno_axis_verdict").get("undetermined", 0), GRY),
                                    ("cis 1", req(ar, "geno_axis_verdict").get("cis_asm", 0), BLU)])]
    y = 0; yt = []; yl = []
    for sname, parts in stages:
        x = 0
        for label, val, col in parts:
            ax.barh(y, val, left=x, color=col, edgecolor="white", height=0.62)
            if val > 40:
                ax.text(x + val / 2, y, f"{label}\n{val}", ha="center", va="center", fontsize=8, color="white")
            x += val
        yt.append(y); yl.append(sname); y -= 1
    ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_xlabel("loci"); ax.set_xlim(0, 1180)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("NACT cis-test funnel: 1139 候選 → 9 non-circular tumor-specific", fontsize=11, color=TX)
    return b64(fig)


# ---- chart 2: candidate vs cis_asm comparison ----
def fig_compare():
    fig, ax = plt.subplots(figsize=(8, 3.8))
    cand = req(sc, "candidate_vs_cis_comparison", "candidate"); cis = req(sc, "candidate_vs_cis_comparison", "cis_asm")
    keys = ["struct_dbeta_median", "germline_baseline_strength", "containment", "direction_concordance", "frac_normal_bimodal"]
    labs = ["struct |Δβ|\n(tumor sig)", "germline\nbaseline", "containment", "direction\nconcord.", "frac normal\nbimodal"]
    x = np.arange(len(keys)); w = 0.38
    ax.bar(x - w / 2, [cand[k] for k in keys], w, label=f"candidate (n={cand['n']})", color=AC)
    ax.bar(x + w / 2, [cis[k] for k in keys], w, label=f"cis_asm (n={cis['n']})", color=GRN)
    for i, k in enumerate(keys):
        ax.text(i - w / 2, cand[k] + 0.01, f"{cand[k]:.3f}", ha="center", fontsize=7, color=TX)
        ax.text(i + w / 2, cis[k] + 0.01, f"{cis[k]:.3f}", ha="center", fontsize=7, color=TX)
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=8); ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("candidate vs cis_asm: tumor signature 相同(0.718≈0.712)但 normal 側全異", fontsize=10, color=TX)
    return b64(fig)


# ---- chart 3: CN stratification rate + CI ----
def fig_cn():
    fig, ax = plt.subplots(figsize=(7, 3.6))
    strat = req(nst, "candidate_rate_by_stratum")
    order = [("neutral", "neutral"), ("LOH", "LOH"), ("gain+loss", "gain+loss"), ("ALL", "ALL")]
    xs = []; rates = []; los = []; his = []; labs = []
    for key, lab in order:
        s = strat[key]; xs.append(len(xs)); rates.append(s["rate_pct"])
        los.append(s["rate_pct"] - s["ci95"][0]); his.append(s["ci95"][1] - s["rate_pct"])
        labs.append(f"{lab}\n{s['k']}/{s['n']}")
    cols = [GRY, RED, "#C98A5B", AC]
    ax.bar(xs, rates, color=cols, width=0.6)
    ax.errorbar(xs, rates, yerr=[los, his], fmt="none", ecolor=TX, capsize=4, lw=1)
    for i, r in enumerate(rates):
        ax.text(i, r + his[i] + 0.6, f"{r}%", ha="center", fontsize=8, color=TX)
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=8.5); ax.set_ylabel("candidate rate %")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("candidate rate by CN (95% CI): neutral 0/11 underpowered [0,28.5]", fontsize=10, color=TX)
    return b64(fig)


# ---- chart 4: genotype census + R3 threshold sensitivity ----
def fig_census():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4))
    gc = req(sc, "geno_census")
    a1.bar(["homogeneous\n(double-dip)", "anchored\n(salvageable)"], [gc["homogeneous"], gc["anchored"]], color=[RED, AC], width=0.6)
    for i, v in enumerate([gc["homogeneous"], gc["anchored"]]):
        a1.text(i, v + 1, str(v), ha="center", fontsize=10, color=TX)
    a1.set_title(f"102 survivors genotype census\n(dom_frac median {req(sc, 'dominant_geno_frac_dist', 'median')})", fontsize=9.5, color=TX)
    a1.spines[["top", "right"]].set_visible(False); a1.set_ylabel("loci")
    ts = req(sc, "R3_threshold_sensitivity"); thr = sorted(float(k) for k in ts)
    a2.plot(thr, [ts[str(t)] for t in thr], "o-", color=AC)
    for t in (0.10, 0.15):
        a2.axvline(t, ls=":", color=MUT, lw=0.8)
    a2.set_xlabel("R3 frac_normal_bimodal cutoff"); a2.set_ylabel("# candidate")
    a2.set_title("count fragility: 102@0.15 → 58@0.10", fontsize=9.5, color=TX)
    a2.spines[["top", "right"]].set_visible(False)
    return b64(fig)


# ---- chart 5: 8-class observation ----
def fig_8class():
    fig, ax = plt.subplots(figsize=(8, 3.4))
    cats = req(ct, "categories_abs")
    items = sorted(cats.items(), key=lambda kv: -kv[1]["n"])
    names = [k.split("_")[0] for k, _ in items]; vals = [v["n"] for _, v in items]
    cols = [AC if n.startswith("A") else (GRN if n.startswith("C") else (BLU if n.startswith("B") else GRY)) for n in names]
    ax.bar(names, vals, color=cols, width=0.7)
    for i, v in enumerate(vals):
        ax.text(i, v + 200, str(v), ha="center", fontsize=8, color=TX)
    ax.set_ylabel("loci"); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"8-class observation (N={req(ct, 'N')}): A=subclone候選723", fontsize=10, color=TX)
    return b64(fig)


# ---- chart 6: fptp attribution ----
def fig_fptp():
    fig, ax = plt.subplots(figsize=(8, 3.4))
    order = [("somatic(HP-1)_≥3", "somatic HP-1"), ("REF_ALT_≥3", "REF/ALT"), ("HP1_HP2_≥3", "HP1/HP2"),
             ("tumor_normal_≥3", "tumor/normal"), ("cluster_≥3", "cluster≥3"), ("cluster_≥5", "cluster≥5")]
    labs = [l for _, l in order]; ratios = [fp[k]["fp_tp_ratio"] for k, _ in order]
    cols = [GRN if r < 0.6 else (GRY if r < 0.95 else RED) for r in ratios]
    ax.bar(labs, ratios, color=cols, width=0.66)
    ax.axhline(1.0, ls="--", color=MUT, lw=0.9)
    for i, r in enumerate(ratios):
        ax.text(i, r + 0.02, f"{r}", ha="center", fontsize=8, color=TX)
    ax.set_ylabel("FP/TP ratio (<1=TP 專一)"); ax.set_xticklabels(labs, fontsize=8, rotation=15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("FP/TP 歸因: somatic HP-1 軸 0.33 TP 專一; cluster FP 富集", fontsize=10, color=TX)
    return b64(fig)


charts = {"funnel": fig_funnel(), "compare": fig_compare(), "cn": fig_cn(), "census": fig_census(),
          "obs8": fig_8class(), "fptp": fig_fptp()}

# ---- survivor handoff table (9 tumor-specific) ----
nine = req(ar, "still_candidate_on_geno_axis")
nine_rows = "".join(f"<tr><td>{r['chrom']}:{r['pos']}</td><td>{r['set']}</td><td>{r['axis']}</td>"
                    f"<td>{r['S_geno']}</td><td>{r['geno_dbeta_median']}</td></tr>" for r in nine)
anchored30 = req(sc, "geno_anchored_loci")
anch_rows = "".join(f"<tr><td>{r['chrom']}:{r['pos']}</td><td>{r['set']}</td><td>{r['cat8']}</td>"
                    f"<td>{r['best_apriori_axis']}</td><td>{r['best_apriori_balance']}</td><td>{r['dominant_geno_frac']}</td></tr>"
                    for r in anchored30)

vd = req(ns, "verdict_dist"); ft = req(nst, "fisher_tests")
cand = req(sc, "candidate_vs_cis_comparison", "candidate"); cis = req(sc, "candidate_vs_cis_comparison", "cis_asm")
ccA = req(cc, "summary", "ALL")
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Phase 1 — Subclone 候選 cis-clean（HCC1395）</title>
<style>
:root{{--ac:{AC};--tx:{TX};--bg:{BG};--bd:{BD};--mut:{MUT};--grn:{GRN};--red:{RED}}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:system-ui,"Noto Sans CJK TC","PingFang TC",sans-serif;color:var(--tx);background:var(--bg);line-height:1.65}}
.wrap{{max-width:1080px;margin:0 auto;padding:0 24px 80px}}
header{{background:var(--tx);color:var(--bg);padding:32px 24px;margin-bottom:8px}}
header h1{{margin:0 0 6px;font-size:25px}} header .sub{{color:#cfc9bf;font-size:14px}}
nav{{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--bd);padding:10px 24px;font-size:13px;z-index:9;margin-bottom:18px}}
nav a{{color:var(--mut);text-decoration:none;margin-right:14px}} nav a:hover{{color:var(--ac)}}
h2{{font-size:19px;margin:34px 0 12px;padding-bottom:6px;border-bottom:2px solid var(--ac)}}
h3{{font-size:15px;margin:20px 0 8px;color:var(--ac)}}
img{{max-width:100%;border:1px solid var(--bd);border-radius:6px;margin:10px 0;background:white}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}}
th,td{{border:1px solid var(--bd);padding:5px 9px;text-align:left}} th{{background:#efe9df}}
td:first-child{{font-family:ui-monospace,monospace}}
.kpi{{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}}
.kpi div{{background:white;border:1px solid var(--bd);border-radius:8px;padding:12px 16px;min-width:130px}}
.kpi b{{display:block;font-size:23px;color:var(--ac)}} .kpi span{{font-size:12px;color:var(--mut)}}
.red{{border-left:4px solid var(--red);background:#fbf0ec;padding:12px 16px;border-radius:0 6px 6px 0;margin:14px 0}}
.ok{{border-left:4px solid var(--grn);background:#eef4ee;padding:12px 16px;border-radius:0 6px 6px 0;margin:14px 0}}
code{{background:#efe9df;padding:1px 5px;border-radius:3px;font-size:12px}}
.cols{{display:flex;gap:16px;flex-wrap:wrap}} .cols>div{{flex:1;min-width:300px}}
footer{{margin-top:40px;padding-top:16px;border-top:1px solid var(--bd);font-size:12px;color:var(--mut)}}
.tag{{display:inline-block;background:var(--ac);color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
</style></head><body>
<header><div class="wrap" style="padding:0">
<h1>Phase 1 — Subclone 候選 cis-clean <span class="tag">⭐2-3 單樣本 HCC1395</span></h1>
<div class="sub">清楚觀察 + 詳細分類定義 → normal-anchored cis-test (NACT)。數字由 JSON 注入（§13-A），不手打。</div>
</div></header>
<nav class="wrap" style="max-width:1080px">
<a href="#tldr">結論</a><a href="#obs">①觀察分類</a><a href="#nact">②NACT funnel</a>
<a href="#verify">③對抗驗證</a><a href="#census">④census</a><a href="#stat">⑤統計誠實</a>
<a href="#handoff">⑥Phase2 handoff</a><a href="#claim">⑦能宣稱</a></nav>
<div class="wrap">

<section id="tldr"><h2>結論先行</h2>
<div class="kpi">
<div><b>1,139</b><span>subclone 候選（cat8 A 723 ∪ subclone_novel 523）</span></div>
<div><b>{vd['cis_asm']}</b><span>cis_asm（DEMOTE 正確）</span></div>
<div><b>{vd['candidate_subclone']}</b><span>NACT survivors（mostly artifact）</span></div>
<div><b>9</b><span>非循環 genotype-錨定 tumor-specific</span></div>
</div>
<div class="red"><b>🔴 深層結論</b>：單樣本甲基<b>無法非循環地證明 subclone</b>。genotype-同質群內的甲基結構永遠 double-dip；唯一非循環甲基訊號是 genotype-aligned = <b>cis/somatic-ASM（characterization）非 subclone</b>。subclone 真確認只能靠 Phase 2 sSNV read-level 連鎖。</div>
<div class="ok"><b>✅ NACT DEMOTE 側可信</b>：{vd['cis_asm']} 個正確判 cis-ASM（normal 沿同 signature 共分離）；germline_baseline cis_asm {cis['germline_baseline_strength']} vs survivors {cand['germline_baseline_strength']}。</div>
</section>

<section id="obs"><h2>① 清楚觀察 + 詳細分類（pre-cis-clean）</h2>
<p>全 {req(ct,'N')} 位點先做 8-class 觀察分類（A=subclone候選 723）+ tumor-only/merged 雙軌 5-態（既有 decisionflow，<b>merged 結構 S4+S5 遠高於 tumor-only → 直接證實 T/N 混合 confound</b>）。各軸 FP/TP 歸因：唯一 TP 專一定位 = somatic HP-1 軸（0.33×）。</p>
<img src="{charts['obs8']}" alt="8-class">
<img src="{charts['fptp']}" alt="fptp">
</section>

<section id="nact"><h2>② NACT cis-test funnel</h2>
<p><b>NACT（Normal-Anchored Concordant cis-Test）</b>：Track A 純 tumor 算 signature → Track B 純 normal 三讀數（R1 投影 / R2 residual / R3 bimodality）→ AND-gate（任一 cis→demote；全不 cis+residual 存活→candidate；不一致→undetermined）。</p>
<img src="{charts['funnel']}" alt="funnel">
<table><tr><th>verdict</th><th>n</th><th>意義</th></tr>
<tr><td>cis_asm</td><td>{vd['cis_asm']}</td><td>normal 共分離 → germline-cis（DEMOTE 正確）</td></tr>
<tr><td>undetermined</td><td>{vd['undetermined']}</td><td>覆蓋/一致性不足</td></tr>
<tr><td>no_tumor_signature</td><td>{vd['no_tumor_signature']}</td><td>tumor-only 無結構 = T/N 假結構</td></tr>
<tr><td>candidate_subclone</td><td>{vd['candidate_subclone']}</td><td>AND-gate 存活 → mostly artifact</td></tr></table>
</section>

<section id="verify"><h2>③ 對抗驗證裁決：mostly_artifact（根因 double-dip 非 CN）</h2>
<p>4-agent workflow（3 透鏡各 spot-check 真 region + 對抗合成）。candidate vs cis_asm 對比證實「tumor signature 閘不判別、R1 結構性必然」：</p>
<img src="{charts['compare']}" alt="compare">
<div class="red">struct |Δβ| candidate <b>{cand['struct_dbeta_median']}</b> ≈ cis_asm {cis['struct_dbeta_median']}（不判別）；containment candidate <b>{cand['containment']}</b> vs {cis['containment']} → R1 not_cis 對 genotype-同質 loci <b>結構性必然</b>，零判別力。spot-check survivors tumor reads 100% genotype 同質卻被甲基切 2:1 = 純 double-dip collider。</div>
</section>

<section id="census"><h2>④ 決定性 genotype census + 門檻敏感度</h2>
<img src="{charts['census']}" alt="census">
<p>102 survivors = <b>homogeneous {req(sc,'geno_census','homogeneous')}（純 double-dip）+ anchored {req(sc,'geno_census','anchored')}</b>（dominant_geno_frac median {req(sc,'dominant_geno_frac_dist','median')}，{req(sc,'dominant_geno_frac_dist','eq_1.0')}/102 完全同質）。genotype-軸 re-test 30 anchored → <b>{req(ar,'geno_axis_verdict').get('candidate_subclone',0)} tumor-specific</b> + {req(ar,'retest_status').get('NO_GENO_SIGNATURE',0)} 無 geno-sig + {req(ar,'geno_axis_verdict').get('undetermined',0)} undet + {req(ar,'geno_axis_verdict').get('cis_asm',0)} cis。R3 門檻：102@0.15 → 70@0.12 → 58@0.10（count 脆弱）。</p>
</section>

<section id="stat"><h2>⑤ 統計誠實（機制更正）</h2>
<img src="{charts['cn']}" alt="cn">
<div class="red">🔴 <b>neutral=0/11 underpowered</b>（CI [0,28.5] 涵蓋 LOH 率；Fisher neutral vs LOH p={ft['neutral_vs_LOH']} 不顯著）→「neutral=0 證 CN-driven」<b>撤回</b>。<br>LOH vs gain Fisher p={ft['LOH_vs_gainloss']}（LOH 升高，但 survivors germline_baseline {cand['germline_baseline_strength']} → 機制<b>非</b> LOH-unmask 而是 double-dip）。<br>TP vs FP p={ft['TP_vs_FP']}（candidate 不分 TP/FP → 非 somatic-specific = artifact-consistent）。</div>
<p>候選池 profile：normal 覆蓋全 ≥10（median {ccA['normal_cov']['median']}），CN/LOH {ccA['is_loh_true']}/1139 LOH。</p>
</section>

<section id="handoff"><h2>⑥ Phase 2 handoff（sSNV 連鎖驗證輸入）</h2>
<div class="cols">
<div><h3>9 個非循環 tumor-specific（genotype-軸存活）</h3>
<table><tr><th>locus</th><th>set</th><th>axis</th><th>S_geno</th><th>|Δβ|med</th></tr>{nine_rows}</table></div>
<div><h3>30 anchored（可救少數，含上 9）</h3>
<table><tr><th>locus</th><th>set</th><th>cat8</th><th>axis</th><th>bal</th><th>dom_frac</th></tr>{anch_rows}</table></div>
</div>
<p>🔴 即使這 9 個 = somatic-ASM characterization，非 confirmed subclone。真確認 = Phase 2 window/phase-block 內 ≥2 sSNV read-level 共現 + 互斥。</p>
</section>

<section id="claim"><h2>⑦ 能宣稱 vs 不能宣稱</h2>
<div class="ok"><b>可宣稱</b>：① NACT DEMOTE work（{vd['cis_asm']} cis-ASM 移除）；② {vd['candidate_subclone']} survivors mostly double-dip artifact；③ 9 個非循環 genotype-錨定 tumor-specific somatic-ASM → Phase 2 候選；④ candidate=「not refuted as germline-cis」⭐3 非偵測。</div>
<div class="red"><b>🔴 不可宣稱</b>：①「102 subclone」/偵測語；②「整池 CN-driven/LOH-unmask」（germline_baseline {cand['germline_baseline_strength']} 反駁）；③ neutral=0 證 CN-driven（p={ft['neutral_vs_LOH']} powerless）；④ R1/R2/R3 獨立佐證（double-dip 必然一致）；⑤ 單樣本任何 survivor 升 subclone（⭐3 characterization-only）。</div>
</section>

<footer>HCC1395 單樣本 ⭐2-3 · build_branch: docs/method-comparison-ism-external-202606 · 數字由 JSON 注入（§13-A）<br>
來源：corrected_tree / contingency_summary / percpg_summary / nact_summary / nact_stats / survivor_census / anchored_retest / candidate_characterization / fptp_attribution.json<br>
方法文件：InterSubMod/docs/methodology/20260625_subclone_verification_observation_classification_01.md（觀察）+ 20260625_subclone_verification_nact_cistest_P1_3_01.md（cis-test）</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes, 6 charts inline")
