#!/usr/bin/env python3
"""
102 - Judgment & verification WORKSTATION (supersedes the locus-judgment display).

Covers ALL 30,350 ISM-analysed loci (every one has a figure after script 101). Lets the
user record include/exclude/undecided + a REASON (criteria / method / visual / other) on
EVERY card (PASS too), export/import the judgments for later correction, surfaces the
"possibly-missed" loci (over-strict structure but filtered), aggregates judgments into a
criteria-redesign signal, and shows live observation charts.

ALL numbers from real files: existence-scan summaries (classification + stats),
tier_assignment.json, validation3.json (independent validation), fig_index.json, funnel.

Output: display_v2/20260610_judgment_workstation_01.html
"""
import csv, json, datetime, os

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = f"{DV}/20260610_judgment_workstation_01.html"


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def cvmax(r):
    return max([num(r, k) or 0 for k in ("CramersV", "CramersV_HPFamily", "CramersV_HPFine")])


def cluster_real(r):
    return tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05 and not tb(r, "ClusterDispersionWarn")


def hp_consistent(r):
    return tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05 and not tb(r, "LabelHPDispersionWarn")


def main():
    tier_assign = json.load(open(f"{DV}/tier_assignment.json"))
    vloci = json.load(open(f"{DV}/validation3.json"))["loci"]
    fn = json.load(open(f"{DV}/funnel_numbers.json"))
    figidx = json.load(open(f"{DV}/fig_index.json")) if os.path.exists(f"{DV}/fig_index.json") else {}
    seqc2 = json.load(open(f"{DV}/seqc2_cn.json")) if os.path.exists(f"{DV}/seqc2_cn.json") else {}
    sqvs = json.load(open(f"{DV}/seqc2_vs_ism.json")) if os.path.exists(f"{DV}/seqc2_vs_ism.json") else {}

    chrs, chr_idx = [], {}
    rows = []
    counts = {0: 0, 1: 0, 2: 0}
    for cls_i, cls in enumerate(("tp", "fp")):
        for r in csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")):
            ch = r["Chr"]
            if ch not in chr_idx:
                chr_idx[ch] = len(chrs); chrs.append(ch)
            pos = int(r["Pos"]); key = f"{ch}_{pos}"
            cv = round(cvmax(r), 3); db = round(num(r, "HPMergedDelta") or 0.0, 3)
            reads = int(num(r, "NumReads") or 0)
            sig = tb(r, "Significant")
            overstrict = (reads >= 20 and tb(r, "PassedGating") and cluster_real(r) and hp_consistent(r) and cv < 0.1)
            cat = 0 if sig else (1 if overstrict else 2)
            counts[cat] += 1
            t = tier_assign.get(key, "")
            tval = {"A": 1, "B": 2}.get(t, 0)
            v = vloci.get(key, {})
            val = 1 if v.get("validated") else (0 if v.get("ok") else -1)
            mhn = int(min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0))
            covcat = ["Normal", "CNV_Gain", "Elevated", "Low", "High_Copy", "CNV_Loss"].index(r.get("Coverage_Category", "Normal")) \
                if r.get("Coverage_Category", "Normal") in ("Normal", "CNV_Gain", "Elevated", "Low", "High_Copy", "CNV_Loss") else 0
            lohsub = ["None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"].index(r.get("LOH_Subtype", "None")) \
                if r.get("LOH_Subtype", "None") in ("None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone") else 0
            qual = ["High", "Medium", "Low"].index(r.get("Quality_Tier", "High")) if r.get("Quality_Tier", "High") in ("High", "Medium", "Low") else 0
            sq = seqc2.get(key, [0, 2.0])
            rows.append([chr_idx[ch], pos, cls_i, cat, cv, db, reads, int(num(r, "NumCpGs") or 0),
                         1 if tb(r, "Potential_LOH") else 0, tval, val, mhn,
                         round(num(r, "ClusterPermanovaF") or 0, 1),
                         round(num(r, "LabelHPPermanovaF") or 0, 1),
                         1 if figidx.get(key, 0) else 0,
                         round(num(r, "GlobalP") or 1, 4),
                         round(num(r, "Coverage_Multiple") or 0, 2), covcat, lohsub, qual,
                         sq[0], sq[1]])

    nfig = sum(1 for x in rows if x[14])
    payload = json.dumps(rows, separators=(",", ":"))
    chrsj = json.dumps(chrs)
    funnel = json.dumps(dict(lps=fn["lps_tp_vcf_records"] + fn["lps_fp_vcf_records"],
                             ism=fn["ism_tp_analyzed"] + fn["ism_fp_analyzed"],
                             pass_=counts[0], missed=counts[1], filtered=counts[2], nfig=nfig))
    stamp = datetime.date.today().isoformat()

    H = TEMPLATE
    for tok, val_ in (("__DATA__", payload), ("__CHRS__", chrsj), ("__FUNNEL__", funnel),
                      ("__STAMP__", stamp), ("__NFIG__", str(nfig), ), ("__NTOTAL__", str(len(rows))),
                      ("__SQVS__", json.dumps(sqvs, ensure_ascii=False))):
        H = H.replace(tok, val_)
    with open(OUT, "w") as f:
        f.write(H)
    print(f"[102] wrote {OUT} ({os.path.getsize(OUT)//1024} KB)")
    print(f"     loci={len(rows)} with_fig={nfig} PASS={counts[0]} MISSED(可能漏掉)={counts[1]} FILTERED={counts[2]}")


TEMPLATE = r"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>位點判斷與驗證工作站 — HCC1395 ASM</title>
<style>
:root{--bg:#0f172a;--card:#1e293b;--mut:#94a3b8;--fg:#e2e8f0;--ac:#38bdf8;--pass:#22c55e;--miss:#f59e0b;--fil:#64748b;--rd:#ef4444}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.6 system-ui,-apple-system,'Segoe UI','Noto Sans CJK TC','Microsoft JhengHei',sans-serif}
header{padding:16px 24px;border-bottom:1px solid #334155;background:#0b1222;position:sticky;top:0;z-index:40}
h1{margin:0 0 3px;font-size:19px}h2{font-size:16px;margin:26px 0 8px;border-left:3px solid var(--ac);padding-left:9px}
.sub{color:var(--mut);font-size:12px}.wrap{max-width:1320px;margin:0 auto;padding:0 20px 90px}
.panel{background:var(--card);border:1px solid #334155;border-radius:10px;padding:14px 16px;margin:12px 0}
.row{display:flex;gap:16px;flex-wrap:wrap;align-items:flex-end}
.ctrl{display:flex;flex-direction:column;gap:3px}.ctrl label{font-size:11px;color:var(--mut)}
input[type=range]{width:150px}select,button,input[type=file]{background:#0b1222;color:var(--fg);border:1px solid #334155;border-radius:6px;padding:6px 9px;font-size:12.5px}
button{cursor:pointer}.v{font-weight:700;color:var(--ac)}
.stat4{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px}
.sb{background:#0b1222;border:1px solid #334155;border-radius:8px;padding:9px 12px}.sb .n{font-size:20px;font-weight:800}.sb .l{font-size:11px;color:var(--mut)}
.tabs{display:flex;gap:5px;flex-wrap:wrap;margin:8px 0}
.tab{cursor:pointer;padding:7px 13px;border-radius:7px;background:#16203a;border:1px solid #334155;font-size:12.5px}
.tab.on{background:var(--card);color:var(--ac);font-weight:700;border-color:var(--ac)}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(290px,1fr));gap:12px}
.card{background:var(--card);border:1px solid #334155;border-radius:9px;overflow:hidden}
.card.j-in{border-color:var(--pass)}.card.j-ex{border-color:var(--rd)}
.card .hd{padding:7px 10px;border-bottom:1px solid #334155;cursor:pointer}
.card .ttl{font-weight:700;font-size:12.5px}.card .meta{font-size:10.5px;color:var(--mut);margin-top:2px}
.figs{display:flex;gap:2px;background:#0b1222;cursor:pointer}.figs img{width:50%;display:block;background:#fff;min-height:40px}
.nofig{padding:14px;color:#64748b;font-size:11px;text-align:center}
.badge{display:inline-block;font-size:9.5px;padding:1px 6px;border-radius:9px;font-weight:700;margin-left:3px}
.b-pass{background:#14532d;color:#86efac}.b-miss{background:#78350f;color:#fcd34d}.b-fil{background:#1e293b;color:#cbd5e1}
.b-A{background:#052e16;color:#4ade80}.b-B{background:#0c2a4d;color:#60a5fa}.b-loh{background:#3b0764;color:#d8b4fe}
.b-vok{background:#052e16;color:#86efac}.b-vno{background:#3b1106;color:#fdba74}.b-fp{background:#450a0a;color:#fca5a5}
.jrow{display:flex;gap:4px;padding:6px 8px;border-top:1px solid #334155}
.jrow button{flex:1;font-size:10.5px;padding:4px 2px}
.jb-in.on{background:#14532d;color:#86efac;border-color:#22c55e}.jb-ex.on{background:#450a0a;color:#fca5a5;border-color:#ef4444}.jb-un.on{background:#1e293b;color:#cbd5e1}
.rsel{width:100%;padding:5px 8px;font-size:10.5px;border-top:1px solid #334155;border-left:0;border-right:0;border-bottom:0;border-radius:0}
.pager{display:flex;gap:8px;align-items:center;justify-content:center;margin:14px 0}
#modal{position:fixed;inset:0;background:rgba(0,0,0,.86);display:none;z-index:100;overflow:auto;padding:22px}
#modal .mc{max-width:1040px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:12px;padding:18px}
#modal img{width:100%;background:#fff;border-radius:6px}.mfigs{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin:12px 0}
.statg{display:grid;grid-template-columns:repeat(auto-fill,minmax(140px,1fr));gap:7px}
.st{background:#0b1222;border:1px solid #334155;border-radius:6px;padding:6px 9px}.st .l{font-size:10px;color:var(--mut)}.st .n{font-size:14px;font-weight:700}
.note{background:#1c1407;border:1px solid #78350f;border-radius:8px;padding:9px 13px;font-size:12px;color:#fcd34d;margin:10px 0}
.bar{height:15px;border-radius:3px;display:inline-block;vertical-align:middle}
.chart{display:flex;flex-direction:column;gap:6px}.chrow{display:flex;align-items:center;gap:8px;font-size:11.5px}
.chrow .lab{width:130px;color:var(--mut);text-align:right}.chrow .val{width:90px}
code{background:#0b1222;padding:1px 5px;border-radius:4px;font-size:11.5px;color:#fcd34d}
.kgreen{color:var(--pass)}.kred{color:var(--rd)}.kor{color:var(--miss)}
</style></head><body>
<header>
 <h1>位點判斷與驗證工作站 — HCC1395 ASM 篩選</h1>
 <div class="sub">全 __NTOTAL__ 個 ISM 分析位點 ｜ 每個都可記錄 應包含/確認排除/無判斷 + 原因 ｜ 可匯出存檔再校正 ｜ 生成 __STAMP__ ｜ 需與同目錄 figs/ 一起開</div>
</header>
<div class="wrap">

<h2>⓪ 圖例與數據說明（先看這裡）</h2>
<div class="panel">
 <div style="display:flex;gap:26px;flex-wrap:wrap">
  <div style="min-width:300px">
   <div class="sub" style="color:var(--ac);font-weight:700;margin-bottom:6px">兩張圖怎麼讀（HP / allele / 甲基 / 距離）</div>
   <svg viewBox="0 0 420 200" width="420" height="200" font-family="system-ui" style="background:#0b1222;border-radius:8px">
    <text x="8" y="16" fill="#cbd5e1" font-size="12">左圖：甲基 read×CpG（每列一條 read，依 HP 分群）</text>
    <rect x="14" y="26" width="12" height="54" fill="#1d4ed8"/><rect x="28" y="26" width="12" height="54" fill="#0ea5e9"/>
    <g font-size="11"><rect x="50" y="28" width="13" height="13" fill="#1d4ed8"/><text x="68" y="39" fill="#94a3b8">HP1（germline 單倍型1）</text>
    <rect x="50" y="45" width="13" height="13" fill="#16a34a"/><text x="68" y="56" fill="#94a3b8">HP1-1（HP1 的體細胞子單倍型）</text>
    <rect x="50" y="62" width="13" height="13" fill="#9333ea"/><text x="68" y="73" fill="#94a3b8">HP2　　<tspan fill="#ca8a04">▉</tspan> HP2-1</text></g>
    <text x="8" y="98" fill="#cbd5e1" font-size="11">Allele 側欄：<tspan fill="#dc2626">▉</tspan>ALT(變異)　<tspan fill="#0ea5e9">▉</tspan>REF(參考)</text>
    <text x="8" y="120" fill="#cbd5e1" font-size="11">甲基色：</text>
    <rect x="62" y="110" width="120" height="13" fill="url(#rb)"/><text x="186" y="120" fill="#94a3b8" font-size="11"><tspan fill="#dc2626">紅</tspan>=甲基 → <tspan fill="#2563eb">藍</tspan>=未甲基；灰=未覆蓋</text>
    <line x1="62" y1="128" x2="62" y2="140" stroke="#f59e0b" stroke-width="2" stroke-dasharray="3"/><text x="68" y="139" fill="#f59e0b" font-size="10">橘虛線=變異位置</text>
    <text x="8" y="162" fill="#cbd5e1" font-size="12">右圖：read×read 距離矩陣（依分群樹排序）</text>
    <rect x="14" y="170" width="120" height="13" fill="url(#mg)"/><text x="140" y="180" fill="#94a3b8" font-size="11">暗=近(相似) → 亮=遠；<b style="fill:#86efac">對角暗塊且對齊 HP 側欄 = 真分群</b></text>
    <defs><linearGradient id="rb"><stop offset="0" stop-color="#2563eb"/><stop offset="1" stop-color="#dc2626"/></linearGradient>
    <linearGradient id="mg"><stop offset="0" stop-color="#000004"/><stop offset="1" stop-color="#fcffa4"/></linearGradient></defs>
   </svg>
  </div>
  <div style="min-width:320px;flex:1">
   <div class="sub" style="color:var(--ac);font-weight:700;margin-bottom:6px">每個數據欄位是什麼（卡片/modal 會顯示）</div>
   <table style="font-size:11.5px"><tr><th>欄位</th><th>直覺</th><th>判讀</th></tr>
    <tr><td>CramersV</td><td>聚類群是否對應 HP 標籤（卡方）</td><td>0–1；≥0.1 才算；稀疏表會歸 0</td></tr>
    <tr><td>Δβ</td><td>HP 家族「間」減「內」的 read 距離分離量</td><td>大=兩家族甲基模式不同（距離量非平均差）</td></tr>
    <tr><td>PERMANOVA F</td><td>距離法測 HP 分離（對稀疏穩健）</td><td>≫1=分離；救回卡方歸 0 的真結構</td></tr>
    <tr><td>minHP</td><td>較小 HP 家族的 read 數</td><td>太小→卡方不可靠；越大越穩健</td></tr>
    <tr><td><b>SEQC2 CN（真值）</b></td><td>SEQC2 benchmark 拷貝數答案（gain/loss/loh/neutral + median CN）</td><td>權威 ground-truth；每位點 intersect <code>ngs_benchmark_cnvs_gain_loss_loh.bed</code></td></tr>
    <tr><td>ISM CN×（proxy）</td><td>ISM 自算 = NumReads / 53x(KDE 二倍體)</td><td><b>非 SEQC2</b>；read-深度估計，經 SEQC2 校準誤差 7.5%；受純度/GC 影響，可能與真值不一致</td></tr>
    <tr><td>獨立驗證</td><td>不用 TP/FP，換獨立子集能否複製</td><td>✓通過 / ⚠失敗(該肉眼) / 未測</td></tr>
   </table>
  </div>
 </div>
</div>

<h2>① 涵蓋與一致性（真實計數）</h2>
<div class="panel"><div class="stat4" id="cover"></div>
 <div class="sub" id="covnote" style="margin-top:8px"></div></div>

<h2>② 判斷紀錄 — 匯出 / 匯入 / 進度（PASS 也可記錄）</h2>
<div class="panel">
 <div class="row">
  <button onclick="exportJ()">⬇ 匯出判斷 (JSON)</button>
  <button onclick="exportCSV()">⬇ 匯出 CSV</button>
  <label class="ctrl"><span class="sub">⬆ 匯入存檔（再校正/重複確認）</span><input type="file" id="imp" accept=".json" onchange="importJ(event)"></label>
  <button onclick="if(confirm('清空所有判斷？'))clearJ()">✕ 清空</button>
 </div>
 <div class="stat4" id="jprog" style="margin-top:12px"></div>
</div>

<h2>③ 篩選標準健診（你的判斷 → 標準是否要重設計）</h2>
<div class="panel" id="redesign"></div>

<h2>④ 觀察圖表（即時，依當前資料）</h2>
<div class="panel"><div id="charts"></div></div>

<h2>④b SEQC2 CN 真值 × ISM 分類 關連分析（confound 檢查）</h2>
<div class="panel" id="sqvs"></div>

<h2>⑤ 即時門檻試算（全 __NTOTAL__）</h2>
<div class="panel">
 <div class="row">
  <div class="ctrl"><label>CramersV ≥ <span class="v" id="vcv">0.10</span></label><input type="range" id="cv" min="0" max="1" step="0.01" value="0.10"></div>
  <div class="ctrl"><label>|Δβ| ≥ <span class="v" id="vdb">0.20</span></label><input type="range" id="db" min="0" max="0.6" step="0.01" value="0.20"></div>
  <div class="ctrl"><label>reads ≥ <span class="v" id="vrd">20</span></label><input type="range" id="rd" min="0" max="80" step="1" value="20"></div>
  <div class="ctrl"><label>邏輯</label><select id="mode"><option value="cv">只 CramersV</option><option value="union">CramersV 或 Δβ</option><option value="and">兩者都要</option><option value="db">只 Δβ</option></select></div>
 </div>
 <div class="stat4" id="live" style="margin-top:12px"></div>
</div>

<h2>⑥ 位點圖庫 — 逐位點檢視與判斷（含「可能漏掉」）</h2>
<div class="tabs" id="cattabs"></div>
<div class="panel" style="border-radius:0 10px 10px 10px">
 <div class="row" style="font-size:12px;align-items:center">
  <span class="sub">排序</span><select id="sort"></select>
  <button id="sdir" title="正序/反序">▼ 降序</button>
  <select id="thf"><option value="off">門檻:不套用</option><option value="pass">只顯示「通過門檻」</option><option value="fail">只顯示「被篩掉」</option></select>
  <select id="cnf"><option value="all">SEQC2 CN:全部</option><option value="1">gain 增幅</option><option value="2">loss 缺失</option><option value="3">loh 雜合缺失</option><option value="0">neutral 正常</option></select>
  <label><input type="checkbox" id="onlyfig"> 只有圖的</label>
  <select id="tierf"><option value="all">Tier:全部</option><option value="A">★A</option><option value="B">B</option><option value="none">無</option></select>
  <select id="valf"><option value="all">驗證:全部</option><option value="1">✓通過</option><option value="0">⚠失敗</option><option value="-1">未測</option></select>
  <select id="jf"><option value="all">判斷:全部</option><option value="in">應包含</option><option value="ex">確認排除</option><option value="un">未判斷</option></select>
  <select id="clsf"><option value="all">TP+FP</option><option value="0">TP</option><option value="1">FP</option></select>
  <span id="gcount" class="sub"></span>
 </div>
 <div class="sub" id="thnote" style="margin-top:6px"></div>
 <div class="grid" id="grid"></div>
 <div class="pager" id="pager"></div>
</div>
</div>
<div id="modal" onclick="if(event.target.id==='modal')closeM()"><div class="mc" id="mc"></div></div>

<script>
const D=__DATA__, CHRS=__CHRS__, FN=__FUNNEL__, SQVS=__SQVS__;
function renderSqvs(){const el=document.getElementById('sqvs');if(!el||!SQVS.B_per_class_rates){el&&(el.innerHTML='(分析檔未生成)');return;}
 const R=SQVS.B_per_class_rates, order=['gain','loh','neutral','loss'];
 let h=`<div class="sub" style="margin-bottom:8px">關連強度 Cramér's V = <b>${SQVS.A_cramersV}</b>（弱-中度）。下表：每個 SEQC2 CN 類別中，ISM 各分類的比例與平均統計（全 ${SQVS.n.toLocaleString()} 位點）。</div>
 <table style="font-size:11.5px"><tr><th>SEQC2 CN (真值)</th><th>n</th><th>PASS%</th><th>可能漏掉%</th><th>Tier A%</th><th>ISM-LOH%</th><th>中位 reads</th><th>平均 CramersV</th><th>平均 Δβ</th><th>PASS 富集×</th></tr>`;
 for(const k of order){const r=R[k];const en=SQVS.D_pass_enrichment[k];
  const hi=k==='loh';
  h+=`<tr${hi?' style="color:#d8b4fe"':''}><td><b>${k}</b></td><td class="n">${r.n.toLocaleString()}</td><td class="n">${r.pass_pct}</td><td class="n">${r.missed_pct}</td><td class="n">${r.tierA_pct}</td><td class="n">${r.ismLOH_pct}</td><td class="n">${r.median_reads}</td><td class="n">${r.mean_cv}</td><td class="n">${r.mean_db}</td><td class="n">${en}×</td></tr>`;}
 h+=`</table>
 <div class="note" style="margin-top:10px"><b>結論（confound 檢查）：</b>ASM 結構在 <b>SEQC2 LOH 區最低</b>（PASS 2.0% / Tier A 0.0% / 平均 CramersV 0.018 / Δβ 0.006，富集 0.46×），在 neutral/gain 最高（富集 1.23×）。
 <b>方向符合生物學</b>：LOH＝一條單倍型丟失 → 只剩單一 haplotype → 無法有 HP1-vs-HP2 甲基差異 → ASM 結構本就該少。<b>這不是 CN 製造假 ASM，而是 CN 狀態決定 ASM 是否「可能存在」（需兩條單倍型）。</b>
 過嚴格(稀疏歸零)位點 <b>89% 在 gain</b>（reads 多但 HP-fine 子群不平衡）、僅 4% 在 LOH（基線 29%）→ 過嚴格非 LOH 驅動。獨立驗證通過率各 CN 類別大致均勻（71-92%）→ 驗證過的 ASM 非 CN artifact。
 ⚠ 單樣本 HCC1395、關連屬描述性（L3 機制推論）；對齊既有 CN-confound pilot（甲基分群差異主要 mutation-linked，confound 通道大致排除）。</div>`;
 el.innerHTML=h;}
// column index: 0 ck,1 pos,2 cls,3 cat,4 cv,5 db,6 rd,7 cg,8 loh,9 tier,10 val,11 mhn,12 pf,13 hpf,14 fig,15 gp
const C={ck:0,ps:1,cl:2,cat:3,cv:4,db:5,rd:6,cg:7,loh:8,tier:9,val:10,mhn:11,pf:12,hpf:13,fig:14,gp:15,cnv:16,covcat:17,lohsub:18,qual:19,sqc:20,sqcn:21};
const $=id=>document.getElementById(id);
const key=r=>CHRS[r[C.ck]]+'_'+r[C.ps];
const CATN=['PASS 通過','可能漏掉(有結構被篩)','確認篩掉(無結構)'];
const COVCAT=['Normal 正常','CNV_Gain 增幅','Elevated 偏高','Low 偏低','High_Copy 高拷貝','CNV_Loss 缺失'];
const LOHSUB=['無','LOH_Noise 雜訊','LOH_Weak 弱','LOH_Strong 強','LOH_Subclone 亞克隆'];
const QUAL=['High','Medium','Low'];
const SQCLS=['neutral 正常','gain 增幅','loss 缺失','loh 雜合缺失'];  // SEQC2 ground-truth class

// ---- judgments (localStorage) ----
const JK='ism_judge_ws_v1';
let J=JSON.parse(localStorage.getItem(JK)||'{}');  // key -> {c:'in'/'ex'/'un', r:'criteria'/'method'/'visual'/'other'/'', t:iso}
function saveJ(){localStorage.setItem(JK,JSON.stringify(J));}
function setJ(k,c){const cur=J[k]||{};cur.c=(cur.c===c?'un':c);cur.t=new Date().toISOString().slice(0,19);J[k]=cur;saveJ();draw();prog();redesign();}
function setR(k,r){const cur=J[k]||{c:'un'};cur.r=r;cur.t=new Date().toISOString().slice(0,19);J[k]=cur;saveJ();redesign();}
window.setJ=setJ;window.setR=setR;

// ---- coverage ----
function coverage(){
 $('cover').innerHTML=[
  ['caller→longphase-S 輸出',FN.lps.toLocaleString(),''],
  ['ISM 實際分析(本工作站)',FN.ism.toLocaleString(),`丟棄 ${FN.lps-FN.ism}（reads太少）`],
  ['有對應圖片',FN.nfig.toLocaleString(),`${(100*FN.nfig/FN.ism).toFixed(1)}%`],
  ['PASS / 可能漏掉 / 確認篩掉',`${FN.pass_} / ${FN.missed} / ${FN.filtered}`,'三分類']
 ].map(s=>`<div class="sb"><div class="n">${s[1]}</div><div class="l">${s[0]}</div>${s[2]?'<div class="l kor">'+s[2]+'</div>':''}</div>`).join('');
 $('covnote').innerHTML=`ISM 比 caller 少 ${FN.lps-FN.ism} 個（read/CpG 太少無法分析）。「可能漏掉」= 有乾淨 HP 結構(PERMANOVA)卻被 gate 篩掉，是你最該逐一檢視是否該納回的一群。`;
}

// ---- judgment progress ----
function prog(){
 const v=Object.values(J);
 const ji=v.filter(x=>x.c==='in').length, je=v.filter(x=>x.c==='ex').length, jn=D.length-ji-je;
 const withR=v.filter(x=>x.c!=='un'&&x.c&&x.r).length, need=v.filter(x=>x.c!=='un'&&x.c&&!x.r).length;
 $('jprog').innerHTML=[
  ['應包含',ji,'kgreen'],['確認排除',je,'kred'],['未判斷',jn,'kmut'],['已標原因 / 待標',`${withR} / ${need}`,'']
 ].map(s=>`<div class="sb"><div class="n ${s[2]||''}">${typeof s[1]==='number'?s[1].toLocaleString():s[1]}</div><div class="l">${s[0]}</div></div>`).join('');
}

// ---- criteria redesign signal ----
function redesign(){
 let passEx=0,missIn=0,filIn=0; const byR={criteria:0,method:0,visual:0,other:0};
 for(const r of D){const j=J[key(r)];if(!j||j.c==='un'||!j.c)continue;
   if(r[C.cat]===0&&j.c==='ex')passEx++;          // PASS 卻判排除 = 誤收(標準太寬)
   if(r[C.cat]===1&&j.c==='in')missIn++;          // 可能漏掉 判包含 = 漏收(標準太嚴)
   if(r[C.cat]===2&&j.c==='in')filIn++;           // 確認篩掉 判包含 = 更嚴重漏收
   if(j.r&&byR[j.r]!==undefined)byR[j.r]++;
 }
 const tot=passEx+missIn+filIn;
 let verdict='judgments 還太少，先逐位點標註再評估。';
 if(tot>=8){
   const parts=[];
   if(missIn+filIn>passEx*1.5) parts.push('<b class="kor">標準偏嚴（漏收多）→ 考慮放寬（如改用 PERMANOVA 納回、降 CramersV 可靠性門檻）</b>');
   else if(passEx>(missIn+filIn)*1.5) parts.push('<b class="kred">標準偏寬（誤收多）→ 考慮收緊</b>');
   else parts.push('誤收/漏收相當 → 標準方向大致合理，差異或在個別位點。');
   if(byR.method>=Math.max(byR.criteria,byR.visual)) parts.push('原因多為「方法設計」→ 可能要改 ISM 統計（如 gate 改 PERMANOVA 主軸）。');
   else if(byR.criteria>=byR.visual) parts.push('原因多為「標準不好」→ 調門檻即可。');
   else parts.push('原因多為「沒看清楚」→ 先複查視覺再下結論。');
   verdict=parts.join(' ');
 }
 $('redesign').innerHTML=`
  <div class="stat4">
   <div class="sb"><div class="n kred">${passEx}</div><div class="l">PASS 判為「應排除」<br>= 誤收（標準太寬）</div></div>
   <div class="sb"><div class="n kor">${missIn}</div><div class="l">可能漏掉 判為「應包含」<br>= 漏收（標準太嚴）</div></div>
   <div class="sb"><div class="n kor">${filIn}</div><div class="l">確認篩掉 判為「應包含」<br>= 更嚴重漏收</div></div>
   <div class="sb"><div class="n">${byR.criteria}/${byR.method}/${byR.visual}</div><div class="l">原因：標準不好 / 方法設計 / 沒看清楚</div></div>
  </div>
  <div class="note" style="margin-top:10px"><b>建議：</b>${verdict}</div>`;
}

// ---- charts (live SVG from current D) ----
function bars(title,items){ // items=[[label,count,color]]
 const max=Math.max(1,...items.map(i=>i[1]));
 return `<div style="min-width:260px"><div class="sub" style="margin-bottom:5px">${title}</div><div class="chart">`+
  items.map(i=>`<div class="chrow"><div class="lab">${i[0]}</div><div class="bar" style="width:${Math.max(2,180*i[1]/max)}px;background:${i[2]}"></div><div class="val">${i[1].toLocaleString()}</div></div>`).join('')+`</div></div>`;
}
function charts(){
 const cat=[0,0,0],tier=[0,0,0],val=[0,0,0],sqc=[0,0,0,0],covc=[0,0,0,0,0,0];
 let cvbin=[0,0,0,0,0];
 for(const r of D){cat[r[C.cat]]++;tier[r[C.tier]]++;val[r[C.val]+1]++;sqc[r[C.sqc]]++;covc[r[C.covcat]]++;
   cvbin[Math.min(4,Math.floor(r[C.cv]*5))]++;}
 $('charts').innerHTML=`<div style="display:flex;gap:30px;flex-wrap:wrap">
  ${bars('分類組成',[['PASS',cat[0],'#22c55e'],['可能漏掉',cat[1],'#f59e0b'],['確認篩掉',cat[2],'#64748b']])}
  ${bars('Tier 分層',[['★Tier A',tier[1],'#4ade80'],['Tier B',tier[2],'#60a5fa'],['無',tier[0],'#475569']])}
  ${bars('獨立驗證(已測者)',[['✓通過',val[2],'#16a34a'],['⚠失敗',val[1],'#f59e0b'],['未測',val[0],'#334155']])}
  ${bars('CramersV 分布',cvbin.map((c,i)=>['['+(i/5).toFixed(1)+'-'+((i+1)/5).toFixed(1)+')',c,'#38bdf8']))}
  ${bars('★ SEQC2 CN 類別 (真值)',[['gain 增幅',sqc[1],'#fdba74'],['loh 雜合缺失',sqc[3],'#d8b4fe'],['neutral 正常',sqc[0],'#64748b'],['loss 缺失',sqc[2],'#7dd3fc']])}
  ${bars('ISM CNV 類別 (read-深度 proxy)',[['Normal',covc[0],'#64748b'],['Gain',covc[1],'#d97706'],['Elevated',covc[2],'#b45309'],['High_Copy',covc[4],'#92400e'],['Low',covc[3],'#0ea5e9'],['Loss',covc[5],'#0369a1']])}
 </div><div class="sub" style="margin-top:8px">圖表即時由當前 __NTOTAL__ 位點計算。<b>SEQC2 CN = 真值</b>（ngs_benchmark_cnvs_gain_loss_loh.bed）；<b>ISM CN× = read-深度 proxy</b>（NumReads/53x KDE 二倍體，經 SEQC2 校準誤差 7.5%）。兩者可能不一致。</div>`;
}

// ---- live threshold ----
function live(){
 const cv=+$('cv').value,db=+$('db').value,rd=+$('rd').value,m=$('mode').value;
 $('vcv').textContent=cv.toFixed(2);$('vdb').textContent=db.toFixed(2);$('vrd').textContent=rd;
 let ptp=0,pfp=0,ttp=0,tfp=0;
 for(const r of D){const tp=r[C.cl]===0;if(tp)ttp++;else tfp++;if(r[C.rd]<rd)continue;
   let p=m==='cv'?r[C.cv]>=cv:m==='db'?Math.abs(r[C.db])>=db:m==='and'?(r[C.cv]>=cv&&Math.abs(r[C.db])>=db):(r[C.cv]>=cv||Math.abs(r[C.db])>=db);
   if(p){if(tp)ptp++;else pfp++;}}
 const ratio=pfp>0?(ptp/pfp).toFixed(1):'∞';
 $('live').innerHTML=[['通過 TP',ptp.toLocaleString(),'kgreen'],['通過 FP',pfp.toLocaleString(),'kred'],
  ['TP:FP',ratio+':1',''],['TP 富集 vs 47:1',((ptp/ttp)/((pfp/tfp)||1e-9)/1).toFixed(2)+'×','']]
  .map(s=>`<div class="sb"><div class="n ${s[2]||''}">${s[1]}</div><div class="l">${s[0]}</div></div>`).join('');
}
['cv','db','rd','mode'].forEach(id=>$(id).addEventListener('input',()=>{live();if($('thf').value!=='off'){page=0;draw();}}));

// ---- gallery ----
const SORTKEY={'reads':r=>r[C.rd],'CramersV':r=>r[C.cv],'|Δβ|':r=>Math.abs(r[C.db]),'PERMANOVA F':r=>r[C.pf],
 'HP-PERMANOVA F':r=>r[C.hpf],'minHP':r=>r[C.mhn],'SEQC2 median CN':r=>r[C.sqcn],'ISM CN×(proxy)':r=>r[C.cnv],'global_p':r=>r[C.gp],'位置':r=>r[C.ck]*1e10+r[C.ps]};
$('sort').innerHTML=Object.keys(SORTKEY).map(k=>`<option>${k}</option>`).join('');
let curCat=1,page=0,PER=24,sdir=-1;  // sdir=-1 降序 / 1 升序
$('sdir').onclick=()=>{sdir*=-1;$('sdir').textContent=sdir===-1?'▼ 降序':'▲ 升序';page=0;draw();};
$('cattabs').innerHTML=CATN.map((n,i)=>`<div class="tab${i===1?' on':''}" data-c="${i}">${n} <span id="cc${i}"></span></div>`).join('');
document.querySelectorAll('.tab').forEach(t=>t.onclick=()=>{curCat=+t.dataset.c;page=0;document.querySelectorAll('.tab').forEach(x=>x.classList.toggle('on',x===t));draw();});
['sort','onlyfig','tierf','valf','jf','clsf','thf','cnf'].forEach(id=>$(id).addEventListener('input',()=>{page=0;draw();}));

function passThreshold(r){ // 依 §⑤ 滑桿判斷此位點是否「通過門檻」
 const cv=+$('cv').value,db=+$('db').value,rd=+$('rd').value,m=$('mode').value;
 if(r[C.rd]<rd)return false;
 return m==='cv'?r[C.cv]>=cv:m==='db'?Math.abs(r[C.db])>=db:m==='and'?(r[C.cv]>=cv&&Math.abs(r[C.db])>=db):(r[C.cv]>=cv||Math.abs(r[C.db])>=db);
}
function filtered(){
 let a=D.filter(r=>r[C.cat]===curCat);
 if($('onlyfig').checked)a=a.filter(r=>r[C.fig]);
 const th=$('thf').value;
 if(th==='pass')a=a.filter(passThreshold); else if(th==='fail')a=a.filter(r=>!passThreshold(r));
 const cn=$('cnf').value;
 if(cn!=='all')a=a.filter(r=>r[C.sqc]===+cn);  // SEQC2 ground-truth class
 const tf=$('tierf').value;if(tf!=='all')a=a.filter(r=>({A:1,B:2,none:0})[tf]===r[C.tier]);
 const vf=$('valf').value;if(vf!=='all')a=a.filter(r=>String(r[C.val])===vf);
 const jf=$('jf').value;if(jf!=='all')a=a.filter(r=>{const c=(J[key(r)]||{}).c||'un';return jf==='un'?(c==='un'):c===jf;});
 const cf=$('clsf').value;if(cf!=='all')a=a.filter(r=>String(r[C.cl])===cf);
 const k=SORTKEY[$('sort').value]||SORTKEY['reads'];
 a.sort((x,y)=>sdir*(k(x)-k(y)));
 return a;
}
function seqc2Badge(r){const c=r[C.sqc],cn=r[C.sqcn];
 const m={1:['#3f1d0a','#fdba74','gain'],2:['#0a2540','#7dd3fc','loss'],3:['#3b0764','#d8b4fe','LOH'],0:['#1e293b','#94a3b8','neutral']}[c];
 return `<span class="badge" style="background:${m[0]};color:${m[1]}" title="SEQC2 真值">SEQC2:${m[2]}${c===1||c===2?' CN'+cn:''}</span>`;}
function cnvBadge(r){const cn=r[C.cnv];  // ISM read-depth proxy (NumReads/53x)
 const col=cn>1.3?['#2a1505','#d97706']:cn>0&&cn<0.7?['#06182e','#0ea5e9']:['#171f2e','#64748b'];
 return `<span class="badge" style="background:${col[0]};color:${col[1]}" title="ISM read-深度 proxy = NumReads/53x">ISM CN×${cn.toFixed(1)}</span>`;}
function badges(r){let b='';
 b+=['<span class="badge b-pass">PASS</span>','<span class="badge b-miss">可能漏掉</span>','<span class="badge b-fil">篩掉</span>'][r[C.cat]];
 if(r[C.tier]===1)b+='<span class="badge b-A">★A</span>';else if(r[C.tier]===2)b+='<span class="badge b-B">B</span>';
 if(r[C.val]===1)b+='<span class="badge b-vok">✓驗</span>';else if(r[C.val]===0)b+='<span class="badge b-vno">⚠</span>';
 if(r[C.cl]===1)b+='<span class="badge b-fp">FP</span>';
 b+=seqc2Badge(r)+cnvBadge(r);return b;}
function jbtn(k,c,lab,cls){const on=((J[k]||{}).c===c)?' on':'';return `<button class="jb-${cls}${on}" onclick="event.stopPropagation();setJ('${k}','${c}')">${lab}</button>`;}
function card(r){
 const k=key(r),ch=CHRS[r[C.ck]],j=J[k]||{};
 const fig=r[C.fig]?`<div class="figs" onclick="openM('${k}')"><img loading="lazy" src="figs/${k}_meth.png"><img loading="lazy" src="figs/${k}_dist.png"></div>`
   :`<div class="nofig">無圖（reads/CpG 太少）</div>`;
 const reasons=['','criteria','method','visual','other'],rlab={'':'— 原因（誤判時標）—',criteria:'標準不好',method:'方法設計',visual:'沒看清楚',other:'其他'};
 const rsel=`<select class="rsel" onchange="event.stopPropagation();setR('${k}',this.value)">`+reasons.map(x=>`<option value="${x}"${(j.r||'')===x?' selected':''}>${rlab[x]}</option>`).join('')+`</select>`;
 const cardcls=j.c==='in'?' j-in':j.c==='ex'?' j-ex':'';
 return `<div class="card${cardcls}">
  <div class="hd" onclick="openM('${k}')"><div class="ttl">${ch}:${r[C.ps]} ${badges(r)}</div>
   <div class="meta">CV=${r[C.cv].toFixed(2)} ｜ Δβ=${r[C.db]>=0?'+':''}${r[C.db].toFixed(2)} ｜ reads=${r[C.rd]} ｜ permF=${r[C.pf]} ｜ minHP=${r[C.mhn]}</div></div>
  ${fig}
  <div class="jrow">${jbtn(k,'in','應包含','in')}${jbtn(k,'ex','確認排除','ex')}${jbtn(k,'un','無判斷','un')}</div>
  ${rsel}</div>`;
}
function draw(){
 const a=filtered();
 for(let i=0;i<3;i++)$('cc'+i).textContent='('+D.filter(r=>r[C.cat]===i).length.toLocaleString()+')';
 $('gcount').textContent=a.length.toLocaleString()+' 個';
 if($('thf').value!=='off'){const cat=D.filter(r=>r[C.cat]===curCat);const p=cat.filter(passThreshold).length;
   $('thnote').innerHTML=`此分類「${CATN[curCat]}」共 ${cat.length.toLocaleString()} 個，在目前門檻(CV≥${(+$('cv').value).toFixed(2)}、|Δβ|≥${(+$('db').value).toFixed(2)}、reads≥${$('rd').value}、${$('mode').value})下：<b class="kgreen">通過 ${p.toLocaleString()}</b> ／ <b class="kor">被篩掉 ${(cat.length-p).toLocaleString()}</b>。拉動 §⑤ 滑桿即時改變。`;}
 else $('thnote').textContent='';
 const pages=Math.max(1,Math.ceil(a.length/PER));page=Math.min(page,pages-1);
 $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
 $('pager').innerHTML=`<button ${page<=0?'disabled':''} onclick="page--;draw()">‹</button><span>第 ${page+1}/${pages} 頁</span><button ${page>=pages-1?'disabled':''} onclick="page++;draw()">›</button>`;
}

// ---- modal ----
window.openM=k=>{const r=D.find(x=>key(x)===k);if(!r)return;
 const ch=CHRS[r[C.ck]];
 const st=[['分類',CATN[r[C.cat]]],['CramersV (gate)',r[C.cv].toFixed(3)],['Δβ (HPMerged)',(r[C.db]>=0?'+':'')+r[C.db].toFixed(3)],
  ['NumReads',r[C.rd]],['NumCpG',r[C.cg]],['global_p',r[C.gp]],['Cluster PERMANOVA F',r[C.pf]],['HP PERMANOVA F',r[C.hpf]],
  ['Tier',['—','A','B'][r[C.tier]]],['獨立驗證',r[C.val]===1?'✓通過':r[C.val]===0?'⚠失敗':'未測'],['min HP group',r[C.mhn]],
  ['SEQC2 CN 類別 (真值)',SQCLS[r[C.sqc]]],['SEQC2 median CN',r[C.sqcn]],
  ['ISM CN× (read-深度 proxy)',r[C.cnv].toFixed(2)+' (=reads/53x)'],['ISM CNV 類別 (proxy)',COVCAT[r[C.covcat]]],
  ['ISM LOH (proxy)',r[C.loh]?LOHSUB[r[C.lohsub]]||'是':'否'],['品質',QUAL[r[C.qual]]]];
 const fig=r[C.fig]?`<div class="mfigs"><div><div class="sub" style="text-align:center">甲基 read×CpG</div><img src="figs/${k}_meth.png"></div><div><div class="sub" style="text-align:center">read×read 距離（樹狀+HP側欄）</div><img src="figs/${k}_dist.png"></div></div>`:'<div class="nofig">無圖</div>';
 $('mc').innerHTML=`<div style="display:flex;justify-content:space-between"><h2 style="border:none;margin:0">${ch}:${r[C.ps]} ${badges(r)}</h2><button onclick="closeM()">關閉 ✕</button></div>
  ${fig}<div class="statg">${st.map(s=>`<div class="st"><div class="l">${s[0]}</div><div class="n">${s[1]}</div></div>`).join('')}</div>
  <div class="jrow" style="margin-top:10px">${jbtn(k,'in','應包含','in')}${jbtn(k,'ex','確認排除','ex')}${jbtn(k,'un','無判斷','un')}</div>
  <div class="sub" style="margin-top:8px">判讀：距離圖對角塊乾淨且對齊 HP 側欄 = 真分群（高信心）；無塊或不對齊 = 弱/假。可能漏掉(橘)請特別看是否該納回。</div>`;
 $('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';
document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});

// ---- export / import ----
function rows4export(){return D.filter(r=>{const c=(J[key(r)]||{}).c;return c&&c!=='un';}).map(r=>{const k=key(r),j=J[k];
 return {locus:k,chr:CHRS[r[C.ck]],pos:r[C.ps],cls:r[C.cl]===0?'TP':'FP',category:CATN[r[C.cat]],cv:r[C.cv],dbeta:r[C.db],tier:['','A','B'][r[C.tier]],validated:r[C.val],choice:j.c,reason:j.r||'',ts:j.t||''};});}
window.exportJ=()=>{const o={generated:new Date().toISOString(),n:rows4export().length,judgments:rows4export()};
 dl('ism_judgments_'+new Date().toISOString().slice(0,10)+'.json',JSON.stringify(o,null,1),'application/json');};
window.exportCSV=()=>{const r=rows4export();const h=['locus','chr','pos','cls','category','cv','dbeta','tier','validated','choice','reason','ts'];
 const csv=[h.join(',')].concat(r.map(x=>h.map(k=>x[k]).join(','))).join('\n');dl('ism_judgments_'+new Date().toISOString().slice(0,10)+'.csv',csv,'text/csv');};
function dl(name,txt,mime){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([txt],{type:mime}));a.download=name;a.click();}
window.importJ=e=>{const f=e.target.files[0];if(!f)return;const rd=new FileReader();
 rd.onload=()=>{try{const o=JSON.parse(rd.result);const arr=o.judgments||o;let n=0;
  for(const x of arr){if(x.locus&&x.choice){J[x.locus]={c:x.choice,r:x.reason||'',t:x.ts||''};n++;}}
  saveJ();draw();prog();redesign();alert('匯入 '+n+' 筆判斷');}catch(err){alert('匯入失敗：'+err);}};rd.readAsText(f);};
window.clearJ=()=>{J={};saveJ();draw();prog();redesign();};

coverage();prog();redesign();charts();renderSqvs();live();draw();
</script></body></html>"""


if __name__ == "__main__":
    main()
