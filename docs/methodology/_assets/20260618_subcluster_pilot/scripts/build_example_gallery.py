#!/usr/bin/env python3
"""[分類範例圖庫 standalone HTML] 7 類重點分類各挑乾淨範例位點, 把 tumor phylo + T/N 兩圖 base64 內嵌
(自包含, 不依賴 symlink, 任意處可顯示)。stats 由 records_v5 注入(§13 不手打), 視覺特徵=策展說明文字。
輸出 20260624_classification_examples_gallery.standalone.html。"""
import json, os, base64

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUTD = f"{A}/obs_ws/cpp_wg"
FIGS = f"{A}/figs_cpp_wg_full"; FIGS_TN = f"{A}/figs_cpp_wg_full_tn"
rec = {f"{r['chrom']}:{int(r['pos'])}": r for r in json.load(open(f"{A}/phylo_cpp_wg_full_records_v5.json"))}

# 策展 7 例: (key(chrom:windowstart), 標題, 視覺特徵說明)
EX = [
    ("chr1:240966025", "① C-S1 cis-ASM 對齊 (1:1)",
     "甲基熱圖分乾淨群，<b>群邊界與右側 HP｜ALT 側欄完全對齊</b> → 結構由既有單倍型/等位解釋 = cis-ASM（非 subclone）。距離矩陣對角塊乾淨。"),
    ("chr17:32608627", "② A SUBCLONE 候選 ⭐（單標籤 somatic + 多結構）",
     "甲基熱圖<b>明顯分群（大 Δβ）</b>，<b>但 HP｜ALT 側欄全同色（單一克隆標籤）</b> → 結構<b>不</b>由標籤解釋 = 同一 somatic 克隆內再分群 = subclone 候選。"),
    ("chr9:26269295", "③ ①many:1 subclone-like（S2 未對齊・結構>標籤）",
     "切多群（coarse=4）<b>但 V 低、不對齊</b> → 無監督結構比 a-priori 標籤更細 = subclone 訊號（非 cis-ASM）。"),
    ("chr16:22539343", "④ ②1:many 跨標籤（無 ASM / trans）",
     "高覆蓋，<b>同一甲基塊內 HP 側欄混色</b>（HP1+HP2 共享同甲基）→ 甲基不分單倍型 = 無等位特異甲基。"),
    ("chr8:133069620", "⑤ B1 LOH 無法區分（單標籤）",
     "單一均勻甲基塊、<b>側欄單色（LOH 丟失一單倍型）</b> → 沒有可比軸 = unevaluable（測不了，非「確認單群」）。"),
    ("chr14:19781114", "⑥ D 可測無結構（有軸但單群）",
     "<b>側欄有 ≥2 標籤（可測）</b>，但甲基<b>單一均勻塊（無分群）</b> → 真單一甲基群體（測了沒差異）。"),
    ("chr1:58056920", "⑦ LOH + somatic 軸 + 結構（癌突變對齊好例）",
     "LOH 區，<b>HP1(germline) vs HP1-1(somatic/ALT) 分群且對齊</b> → 留存單倍型上 somatic 變異連鎖甲基。"),
]


def uri(path):
    if not os.path.exists(path):
        return None
    return "data:image/png;base64," + base64.b64encode(open(path, "rb").read()).decode()


def badges(r):
    b = [("TP" if r["set"] == "TP" else "FP", "bt" if r["set"] == "TP" else "bf"),
         (r.get("cat8", "?"), "bc"), (str(r.get("ctype", ""))[:7], "bx"),
         (r["cn_state"] + ("" if r["cn_state"] == "loh" else f" CN{r['cn_value']}" if isinstance(r.get("cn_value"), (int, float)) else ""), "bn"),
         (f"coarse {r['coarse_ng']}", "bg"), (f"V {r.get('V_hp',0)}/{r.get('V_allele',0)}", "bg"),
         (f"Δβ群 {r.get('m_dbeta_group')}" if r.get("m_dbeta_group") is not None else "Δβ群 —", "bm"),
         (f"n {r['n']} (T{r.get('n_tumor','?')}/N{r.get('n_normal','?')})", "bg")]
    return "".join(f'<span class="bdg {c}">{v}</span>' for v, c in b)


cards = []
for k, title, sig in EX:
    r = rec.get(k)
    if not r:
        continue
    p = int(r["pos"]); base = f"cpp_{r['chrom']}_{p}_{p+10000}"
    ph = uri(f"{FIGS}/{base}.png"); tn = uri(f"{FIGS_TN}/{base}_tn.png")
    sk = f"{r['chrom']}:{p+5000}"
    imgs = ""
    if ph: imgs += f'<figure><figcaption>① tumor phylo（樹·甲基·距離, 側欄 HP｜ALT）</figcaption><img src="{ph}"></figure>'
    if tn: imgs += f'<figure><figcaption>② tumor-vs-normal 甲基對照（上 tumor／下 normal, 側欄 T/N｜HP）</figcaption><img src="{tn}"></figure>'
    cards.append(f'''<section class="ex"><h2>{title}</h2>
<div class="meta">搜尋鍵 <code>{sk}</code>（貼進儀表板 🔍 可跳到此位點）</div>
<div class="bdgs">{badges(r)}</div>
<div class="sig">👁 <b>圖中該看什麼：</b>{sig}</div>
<div class="imgs">{imgs}</div></section>''')

HTML = f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>subclone 分類範例圖庫</title><style>
:root{{--bg:#0a0e16;--card:#111826;--fg:#e7edf5;--mut:#8b9bb4;--line:#2a3344;--ac:#D97757}}
*{{box-sizing:border-box}}body{{font-family:system-ui,"Noto Sans CJK TC",sans-serif;background:var(--bg);color:var(--fg);margin:0;line-height:1.6}}
.wrap{{max-width:1080px;margin:0 auto;padding:16px}}
h1{{font-size:20px;margin:4px 0}}.sub{{color:var(--mut);font-size:12.5px}}
.note{{background:#2a2410;border:1px solid #d97706;border-radius:8px;padding:9px 13px;font-size:12px;margin:10px 0}}.note b{{color:#fbbf24}}
.legend{{background:#0b1222;border:1px solid var(--line);border-radius:8px;padding:9px 13px;font-size:12px;color:var(--mut);margin:8px 0}}
.ex{{background:var(--card);border:1px solid var(--line);border-radius:11px;padding:13px 15px;margin:14px 0}}
.ex h2{{font-size:15.5px;margin:0 0 6px}}.meta{{font-size:12px;color:var(--mut);margin-bottom:6px}}code{{background:#0b1222;border:1px solid var(--line);border-radius:4px;padding:1px 6px;color:#7fb4e8}}
.bdgs{{display:flex;flex-wrap:wrap;gap:5px;margin:6px 0}}.bdg{{font-size:10.5px;padding:2px 7px;border-radius:9px;border:1px solid var(--line)}}
.bt{{background:#10243a;color:#7cc4ff}}.bf{{background:#3a1a1a;color:#ffa3a3}}.bc{{background:#3a2e0a;color:#fcd34d;border-color:#eab308;font-weight:700}}.bx{{background:#241038;color:#d6b3ff;border-color:#a855f7}}.bn{{background:#1a2130;color:var(--mut)}}.bg{{background:#0f1a2a;color:#9fc6ff}}.bm{{background:#06262e;color:#5fd6c8;border-color:#14b8a6}}
.sig{{background:#0d1422;border-left:3px solid var(--ac);border-radius:0 6px 6px 0;padding:7px 12px;font-size:12.5px;margin:8px 0}}.sig b{{color:#ffd9a0}}
.imgs{{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:8px}}@media(max-width:760px){{.imgs{{grid-template-columns:1fr}}}}
figure{{margin:0}}figcaption{{font-size:11px;color:var(--mut);margin-bottom:3px}}img{{width:100%;border-radius:6px;background:#fff;display:block}}
.cmp{{background:#0b1222;border:1px solid var(--line);border-radius:8px;padding:10px 14px;font-size:12.5px;margin:14px 0}}.cmp b{{color:#ffd9a0}}
footer{{color:var(--mut);font-size:11px;margin-top:18px;border-top:1px solid var(--line);padding-top:10px}}
</style></head><body><div class="wrap">
<h1>Subclone 分類重點範例圖庫 — HCC1395 single-sample</h1>
<div class="sub">每類挑一個乾淨範例，內嵌 tumor phylo + tumor-vs-normal 兩張圖，附「圖中該看什麼」。⭐2-3 觀察/刻畫層，非 subclone caller。</div>
<div class="note">🔴 <b>紅線</b>：① 對齊 = cis-ASM 跡象<b>非 subclone</b>；② confident-multi 在 FP 更多 = <b>反判別</b>；③ 「subclone 候選」是<b>候選非確認</b>（單樣本無單細胞真值，確認需 multi-region/single-cell）。</div>
<div class="legend">📖 <b>怎麼讀</b>：① tumor phylo = 左 UPGMA 樹（群色）｜中 甲基 read×CpG（RdBu_r：紅=甲基/藍=未甲基/灰=NaN）｜右 read×read 距離（暗=近），側欄 HP（藍HP1/紫HP2）｜ALT（紅ALT/琥珀REF）。② T/N 對照 = 上 tumor／下 normal 甲基熱圖，側欄 T/N（橘tumor/綠normal）｜HP。<b>關鍵：看甲基分群塊邊界與側欄標籤對不對齊</b>。</div>
{"".join(cards)}
<div class="cmp">🔬 <b>對照觀察（最有教學性）</b>：<br>
• <b>①(對齊) vs ②(單標籤多結構)</b>：①甲基塊<b>跟側欄對齊</b>(cis-ASM)；②甲基塊<b>跟側欄不對齊、側欄單色</b>(克隆內 subclone)。<br>
• <b>③(結構>標籤) vs ④(結構<標籤)</b>：③切多群但側欄對不上；④側欄混進同一塊。<br>
• <b>⑤(無法區分) vs ⑥(真單群)</b>：⑤側欄<b>單色</b>(沒得比)；⑥側欄<b>多色但甲基不分</b>(比了沒差異)。</div>
<footer>來源 phylo_cpp_wg_full_records_v5.json（stats 注入）+ figs_cpp_wg_full/＋_tn/（base64 內嵌, 自包含）。分類定義+各類數量比例: InterSubMod/docs/methodology/20260624_subclone_classification_decision_tree_corrected_01.md。commit d360316 · feat/summary-nreadsvalid · 單樣本 HCC1395 ⭐2-3。</footer>
</div></body></html>'''

outp = f"{OUTD}/20260624_classification_examples_gallery.standalone.html"
open(outp, "w").write(HTML)
print(f"WROTE {outp} ({len(HTML)//1024} KB) | {len(cards)}/7 examples embedded")
