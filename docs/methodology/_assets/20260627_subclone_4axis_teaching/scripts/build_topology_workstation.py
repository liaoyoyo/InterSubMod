#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互動式拓樸工作站 HTML v2：
  - 每區 genotype 向量(S1/S2..標籤) + 克隆樹 SVG(節點=S-mut-set·reads·%·HP) + 分支(V/H)
  - 篩選: chr / genome_ctx(telo/centro/arm) / 拓樸型 / determinacy / 含FP / 最少群數 / 搜尋
  - 排序: 座標 / 複雜度(n_sSNV) / 群數 / region
  - region badge: genome_ctx + TP/FP 組成
  - chr17 完整 worked panel(S/r/m 一致標籤 + 樹 + 16 sig CpG)
§13-A: 數字由 topology_per_region.json 注入。
用法: python3 build_topology_workstation.py
"""
import json, os, glob
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
MULTI_OUT = os.environ.get("SM_OUT") or os.path.normpath(os.path.join(HERE, "..", "..", "20260629_multisample_topology_workstation.standalone.html"))
# OUT(舊單樣本 20260628)已 deprecated:build 只產多樣本 MULTI_OUT(=主結果);舊單樣本檔不再寫出

# 多樣本(2026-06-29):SM_SAMPLES="name:dir,name:dir" → 多分頁;預設納入已完成樣本(HCC1395 凍結 + multisample_subclone 下有 topology 的)。
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
def _sample_dirs():
    env = os.environ.get("SM_SAMPLES", "")
    if env:
        return [tuple(x.split(":", 1)) for x in env.split(",") if ":" in x]
    pairs = [("HCC1395", DATA)]  # 凍結主樣本
    if os.path.isdir(MSROOT):
        for s in sorted(os.listdir(MSROOT)):
            tp = os.path.join(MSROOT, s, "topology_per_region.json")
            if os.path.exists(tp):
                pairs.append((s, os.path.join(MSROOT, s)))
    return pairs

def _load_sample(dr):
    d = json.load(open(os.path.join(dr, "topology_per_region.json"), encoding="utf-8"))
    rec = {"stats": d["stats"], "detail": d["detail"], "chr17": d.get("chr17_worked"),
           "chroms": d.get("chroms", []), "provenance": d.get("provenance", {})}
    csp = os.path.join(dr, "candidate_scoring.json")
    if os.path.exists(csp):
        cs = json.load(open(csp, encoding="utf-8"))
        rec["scoring"] = {"summary": {k: cs.get(k) for k in ("n_total","n_need_confirm","score_formula","situation_dist","resolution_dist","score_buckets","needs_methyl_n")},
                          "queue": cs.get("queue", [])}
    else:
        rec["scoring"] = {"summary": {}, "queue": []}
    gp = os.path.join(dr, "region_gene_annotation.json")
    rec["gene"] = json.load(open(gp, encoding="utf-8")).get("regions", {}) if os.path.exists(gp) else {}
    accp = os.path.join(dr, "single_snv_accounting.json")
    rec["accounting"] = json.load(open(accp, encoding="utf-8")) if os.path.exists(accp) else None
    ctp = os.path.join(dr, "candidate_trees.json")  # R6/Part B: enumerate_candidate_trees 誠實版
    rec["candtrees"] = {x["region"]: x for x in json.load(open(ctp, encoding="utf-8")).get("candidate_trees", [])} if os.path.exists(ctp) else {}
    rdp = os.path.join(dr, "rd_perregion.json")  # read-driven per-region 交叉確認(22-way 平行遍歷 read)
    rec["rd"] = json.load(open(rdp, encoding="utf-8")) if os.path.exists(rdp) else {}
    idp = os.path.join(dr, "ideogram_data.json")  # per-sample HG38 ideogram(census+topology 衍生,免BAM)
    rec["ideogram"] = json.load(open(idp, encoding="utf-8")) if os.path.exists(idp) else None
    # LOH 底帶(per-sample):multisample 樣本用 {dr}/{S}_clp_cn.bed 篩 loh;HCC1395 在 SAMPLES 後 patch 用 longphase-TO(SEQC2-validated)
    rec["loh"], rec["loh_source"] = [], None
    _beds = glob.glob(os.path.join(dr, "*_clp_cn.bed"))
    if _beds:
        rec["loh_source"] = "clp_cn (CLP CN caller)"
        seen = set()
        for _ln in open(_beds[0], encoding="utf-8"):
            _p = _ln.split()
            if len(_p) >= 4 and _p[3].lower() == "loh":
                try:
                    iv = (_p[0], int(_p[1]), int(_p[2]))
                    if iv not in seen:
                        seen.add(iv); rec["loh"].append([iv[0], iv[1], iv[2]])
                except ValueError:
                    pass
    return rec

SAMPLES = {name: _load_sample(dr) for name, dr in _sample_dirs()}
SAMPLE_NAMES = list(SAMPLES.keys())
# HCC1395 凍結樣本 ideogram 在 20260630_perregion_workstation/data(非其 frozen dir)→ 補上(per-sample 接線)
_ideo = os.path.normpath(os.path.join(HERE, "..", "..", "20260630_perregion_workstation", "data", "ideogram_data.json"))
if "HCC1395" in SAMPLES and SAMPLES["HCC1395"].get("ideogram") is None and os.path.exists(_ideo):
    SAMPLES["HCC1395"]["ideogram"] = json.load(open(_ideo, encoding="utf-8"))
# 結構解析 + 所有子read(2026-07-04 gained-pair 定序 + block 邊 + provenance;全 3885 區)
_res = os.path.join(DATA, "subreads_all.json")
_res0 = os.path.join(DATA, "resolution_subreads.json")  # fallback(僅 42 多突變區)
if os.path.exists(_res):
    RESOLUTION_JSON = json.dumps(json.load(open(_res, encoding="utf-8")), ensure_ascii=False)
elif os.path.exists(_res0):
    RESOLUTION_JSON = json.dumps(json.load(open(_res0, encoding="utf-8")).get("regions", {}), ensure_ascii=False)
else:
    RESOLUTION_JSON = "{}"
# >8 sSNV 全 pairwise 樹(2026-07-04 去 8-cap;CN-clean 建樹·CN-gain multiplicity-caveat)
_gt8 = os.path.join(DATA, "gt8_trees.json")
GT8_JSON = json.dumps(json.load(open(_gt8, encoding="utf-8")).get("trees", {}), ensure_ascii=False) if os.path.exists(_gt8) else "{}"
# LOH 全 7 樣本用 longphase-TO tumor_phased_LOH.bed(SEQC2-validated;取代 clp_cn)。DORADO=同細胞株用 HCC1395。
_LPT = "/big7_disk/liaoyoyo2001/longphase-to-mod/output"
LOH_BEDS = {
    "HCC1395":        _LPT + "/baseline/tumor_phased_LOH.bed",
    "HCC1395_DORADO": _LPT + "/baseline/tumor_phased_LOH.bed",  # 同細胞株(Dorado 重 basecall,同 DNA→同 LOH)
    "COLO829":        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437":          _LPT + "/v6_5sample_extension/H1437/tumor_phased_LOH.bed",
    "H2009":          _LPT + "/v6_5sample_extension/H2009/tumor_phased_LOH.bed",
    "HCC1937":        _LPT + "/v6_5sample_extension/HCC1937/tumor_phased_LOH.bed",
    "HCC1954":        _LPT + "/v6_5sample_extension/HCC1954/tumor_phased_LOH.bed",
}
def _read_loh_bed(path):
    loh, seen = [], set()
    for ln in open(path, encoding="utf-8"):
        p = ln.split()
        if len(p) >= 3:
            try:
                iv = (p[0], int(p[1]), int(p[2]))
                if iv not in seen:
                    seen.add(iv); loh.append([iv[0], iv[1], iv[2]])
            except ValueError:
                pass
    return loh
for _nm, _rec in SAMPLES.items():
    _lp = LOH_BEDS.get(_nm)
    if _lp and os.path.exists(_lp):
        _rec["loh"] = _read_loh_bed(_lp)
        _rec["loh_source"] = "longphase-TO (SEQC2-validated)" + ("·同細胞株借 HCC1395" if _nm == "HCC1395_DORADO" else "")
    # 否則保留 _load_sample 載入的 clp_cn 作為 fallback
SAMPLES_JSON = json.dumps(SAMPLES, ensure_ascii=False)
# R3+R7: chr17 read×read 距離矩陣 + 分群×lineage 交叉表(HCC1395 worked;固定教學甲基展板用)
_c17t = os.path.join(DATA, "chr17_tree_data.json")
CHR17TREE_JSON = json.dumps(json.load(open(_c17t, encoding="utf-8")), ensure_ascii=False) if os.path.exists(_c17t) else "null"
# (死碼移除 2026-07-01 審查 p2:FIRST/d/acc/B/DJ/UNIVERSE_BANNER 從不被 HTML 引用;宇宙帳本由 JS renderUniverse 逐樣本 runtime 渲染。rebuild byte-identical 已證移除無影響)
# HG38 ideogram 改為 per-sample 嵌入 __SAMPLES__[s].ideogram(見上 SAMPLES 補丁);舊全域 __IDEOGRAM__ 已移除(merge 2026-07-01)

GLOSSARY = [
 ("sSNV / S1·S2·S3", "體細胞單核苷酸變異；S1..Sk = 區內依座標排序的 sSNV（基因型向量第 i 位 = Si）。", "癌細胞才有、正常細胞沒有的點突變。一個區域有 k 個就標 S1..Sk，順序按基因座位置。"),
 ("read 群 / r / population", "同一基因型向量(如 RAR)的 read 集合 = 一種『細胞狀態』。", "一條 read=一條 DNA 分子=一個細胞的一條染色體。攜帶相同突變組合的 read 歸為同一群（population），代表一個 lineage 節點。"),
 ("m / 甲基位點", "顯著差異的 CpG 甲基化位點（chr17: m1..m16）。", "DNA 甲基化標記；m 是 L1 vs L2 lineage 間 |Δβ| 大的 CpG。⚠ 實測對齊 cis-genotype 軸，非獨立 lineage。"),
 ("HP / H1·H2·H3", "germline 單倍型（哪條親代染色體）。H1/H2=兩條;H3?=未定相(somatic-ALT read)。", "由 longphase-S haplotag 決定。正常人 2 條 haplotype 根 H1、H2。HP tag『1-1』→H1、『2-1』→H2、『3』→H3?。"),
 ("HP{h}-path", "lineage 標籤 = HP{根}-{分支1}-{分支2}…（Dewey 路徑）。", "如 H1-1=H1 上第一個 somatic 事件;H1-1-1=其後代;H1-2=姊妹分支。分支編號=VAF 遞減。"),
 ("vertical 直系 / nested", "祖先→後代（一個細胞先後累積兩突變）。", "2×2 有 AA 格(兩突變同 read) + 一側零格 → 巢狀。樹上往下一層。"),
 ("horizontal 姊妹 / sibling", "兩突變從不共現但同 haplotype → 不同 subclone 平行分支。", "2×2 的 AA=0(兩 ALT 從不同 read) + same-HP。樹上同層分叉。"),
 ("co_linked", "兩突變完美共現(只見 AA) → 同一 lineage 事件(同節點)。", "RA=AR=0、AA≥2。兩突變永遠同進退，無法內部排序。"),
 ("mutual_excl 互斥", "兩 ALT 從不共現(AA=0)。", "diff-HP→allelic(兩條染色體各自突變,非 subclone);same-HP→sibling subclone。"),
 ("independent / four-gamete 違反", "RR/RA/AR/AA 四格全有 → 不相容單一樹。", "違反無限位點假設(回復/重複突變/CNV multiplicity/偽影)。"),
 ("perfect-phylogeny 完美系統發生樹", "二元字元(REF/ALT)的樹:每位點只突變一次。", "古典定理:二元字元『每對相容』即足以保證整棵樹存在 → pairwise 拼接合法，不需單分子整跨。"),
 ("2×2 共現 (RR/RA/AR/AA)", "對每對 sSNV 數共讀 read 在兩位點的 REF/ALT 組合。", "RR=都不帶、AA=都帶、RA/AR=只一帶。哪格為零決定關係(co_linked/nested/互斥)。"),
 ("ε=2% 噪聲底線", "cell 為真 ⟺ count > coread×2%（ONT 錯誤率）。", "保留最低 1 條(低 coread);高 coread 單讀(1≤coread×2%)判噪聲。經 FP 裁判+結構穩定+塌陷集中三路定案。"),
 ("coread 共讀", "同時覆蓋兩個位點的 read 數。", "≥6 才算 powered。決定一個零格是否可信。"),
 ("VAF / CCF", "VAF=變異等位基因頻率;CCF=癌細胞比例(clonal prevalence)。", "高 VAF=clonal(早/大);低=subclonal(晚/小)。只在 CN-clean 可信。用於分支編號 + 單位點刻畫。"),
 ("determinacy", "拓樸能否唯一辨識:A_determined(單分子向量) / A_ambiguous(缺中間群順序未定) / B_pairwise(拼接非單分子) / C_underdetermined(單ALT群·缺連接read·只觀測到一個 ALT 群無第二群可定序;非「多棵樹相容」之意) / incompatible(成環)。", "🔴『樹存在』(絕大多數位點相容)≠『能辨識是哪一棵』(僅約半數 A_determined;其餘 B拼接/C欠定/缺中間群/成環)。實際比例見上方 determinacy 圓餅(隨樣本變)。"),
 ("situation tier", "A 單分子整跨 / B 可整跨pairwise / C 必鏈接(span>read)。", "≥3 位點先分 situation 再處理:有沒有一條 read 穿過全部決定證據強度。"),
 ("genome_ctx", "telomere(端粒,≤3Mb 端)/ centromere(著絲點±3Mb)/ arm(染色體臂)。", "hg38 染色體長度+centromere 近似。centromere/telomere 區偽影風險高。"),
 ("TP / FP", "真/假陽性標籤(來源隨樣本)。🔴 HCC1395=SEQC2 truth set;其餘樣本=per_sSNV_census 衍生(無外部 truth set、標籤弱、勿當判別力佐證,見宇宙帳本 caveat)。", "🔴 只用於觀察評估,絕不進前處理/定義(build 用 TP∪FP union + normal 比對定 somatic)。"),
 ("Tier-R / Tier-PS", "Tier-R=same-read(≤50kb,同分子);Tier-PS=same phase-set(>50kb,統計相位,未做)。", "克隆連鎖只認 Tier-R;isolated 區可能有 Tier-PS partner 待救。"),
 ("cluster-count (c, 上界 ≤ k)", "c = 含≥1 ALT 的 distinct genotype 向量數(germline 不計);perfect-phylogeny 下 c ≤ k(非 2^k;含 germline 的 population 數才 ≤ k+1)。", "實測 c 多為 1-2 → 拓樸搜尋空間極小。先定 c 再縮限拓樸。"),
 ("ambiguous-parentage 缺中間群", "節點突變集跳>1(中間 population 沒觀察到)→累積順序未定。", "76 區。如 {0,3,4} 缺 {0,3} 等中間群 → 0,3,4 哪個先未定。"),
 ("linked / underpowered / isolated", "全 sSNV 三桶:可建樹 / 有 partner 無共讀(可救) / 無 partner(Tier-R 樹外)。", "三桶比例隨樣本變(見上方宇宙帳本,隨分頁更新)。單位點非全無法處理:underpowered 有 CCF、isolated 有 caller VAF+可能 Tier-PS。"),
 ("cis-ASM / double-dip", "甲基隨突變的 cis 局部效應 / 用同量定群又驗群的循環。", "chr17 證甲基分群對齊突變 genotype 軸(cis)非獨立 lineage → 甲基不能當獨立驗證器。06-28 normal cis-control 已測:CROSS-HP 35.4% 可控、SAME-HP 多數區 normal 無對應 within-HP 軸=結構性無法 control(需 single-cell)。"),
 ("bounded-auxiliary 甲基定位", "甲基=corroborate 非 detect 的有界輔助(Tier-3 機率層)。", "排序:genetic 共現>HP 定根>甲基。🔴 06-28 cis-control pilot 已測定案:cis-control 只對 CROSS-HP 區有效(約三成)、SAME-HP 多數區結構性無解→甲基不能升 resolver;最終角色=cluster-count sanity + 43 區 CROSS-HP 弱排序 PoC。"),
 ("⭐3 / single-pipeline", "單樣本 HCC1395 單一 pipeline 的證據上限。", "升 ⭐4 需 ≥5/7 樣本+COLO829+single-cell 正交確認。"),
]
GLOSSARY_CHAPTERS = [("① 基本標籤", [0,1,2,3,4]), ("② SNV 關係 / 拓樸型", [5,6,7,8,9,10]), ("③ 共現與證據量", [11,12,13,14]), ("④ 可辨識性 determinacy", [15,16,20,21]), ("⑤ 基因體與真值", [17,18,19]), ("⑥ 全 sSNV 帳本與甲基", [22,23,24,25])]
def _gterm(i):
    t, s, dd = GLOSSARY[i]
    return f'<details style="border:1px solid #f1f3f5;border-radius:6px;padding:5px 9px;font-size:12px"><summary style="cursor:pointer"><b>{t}</b></summary><div style="margin-top:4px;color:#343a40">{s}</div><div style="margin-top:3px;color:#868e96;font-size:11px">{dd}</div></details>'
GLOSSARY_HTML = ('<details style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 14px;margin:10px 0"><summary style="cursor:pointer;font-weight:600;color:#1971c2">📖 名詞與概念解釋（分 6 章節；點章節→點詞展開）</summary>'
 + "".join('<details style="border:1px solid #e9ecef;border-radius:6px;padding:6px 10px;margin-top:6px;background:#f8f9fa"><summary style="cursor:pointer;font-weight:600;color:#495057">' + ch + f'（{len(idxs)} 詞）</summary><div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(290px,1fr));gap:6px;margin-top:8px">' + "".join(_gterm(i) for i in idxs) + '</div></details>' for ch, idxs in GLOSSARY_CHAPTERS)
 + '</details>')

CSS = """
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa}
.wrap{max-width:1320px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:.2em 0}h3{margin:.3em 0}.sub{color:#868e96;font-size:12.5px}
.stats{display:flex;gap:12px;flex-wrap:wrap;margin:10px 0}.scard{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:9px 11px;min-width:200px}
.scard h4{margin:0 0 5px;font-size:11.5px;color:#495057}.bar{display:flex;align-items:center;gap:5px;font-size:10.5px;margin:2px 0}.bar i{height:10px;background:#1c7ed6;border-radius:2px;display:inline-block}
.scard{cursor:pointer;transition:box-shadow .12s}.scard:hover{box-shadow:0 2px 12px rgba(0,0,0,.12)}.scard h4 .more{float:right;font-size:9px;color:#adb5bd;font-weight:400}
.smbg{position:fixed;inset:0;background:rgba(0,0,0,.45);z-index:50;display:none;align-items:center;justify-content:center;padding:20px}
.smbox{background:#fff;border-radius:12px;max-width:780px;width:100%;max-height:86vh;overflow:auto;padding:18px 22px;position:relative}
.smbox .x{position:absolute;top:10px;right:15px;cursor:pointer;font-size:22px;color:#868e96;background:none;border:none;line-height:1}
details.c17{background:#fff;border:1px solid #ffd8a8;border-radius:8px;padding:10px 14px;margin:10px 0}details.c17 summary{cursor:pointer;font-weight:600;color:#d9480f}
.ctrl{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:10px 0;font-size:12.5px;background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:9px}
.ctrl select,.ctrl input{padding:3px 7px;border:1px solid #ced4da;border-radius:5px;font-size:12.5px}
.main{display:grid;grid-template-columns:400px 1fr;gap:12px}@media(max-width:860px){.main{grid-template-columns:1fr}}
.list{background:#fff;border:1px solid #dee2e6;border-radius:8px;max-height:76vh;overflow:auto}
.row{padding:6px 10px;border-bottom:1px solid #f1f3f5;cursor:pointer;font-size:12px}.row:hover{background:#e7f5ff}.row.sel{background:#d0ebff}.row b{color:#1c7ed6}
.tag{font-size:9.5px;padding:1px 6px;border-radius:9px;margin-left:3px}
.t_linear{background:#d3f9d8;color:#2b8a3e}.t_branched{background:#e5dbff;color:#5f3dc4}.t_star{background:#fff3bf;color:#b08900}.t_single{background:#f1f3f5;color:#868e96}
.ctx_telomere{background:#d0ebff;color:#1971c2}.ctx_centromere{background:#ffe3e3;color:#c92a2a}.ctx_arm{background:#f1f3f5;color:#868e96}
.facets{display:flex;gap:10px;flex-wrap:wrap;margin:8px 0}
.fcard{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:8px 11px;flex:1;min-width:230px}
.fcard .fh{display:flex;align-items:center;gap:6px;font-size:11.5px;font-weight:700;color:#495057;margin-bottom:7px;padding-bottom:5px;border-bottom:1px solid #f1f3f5}
.fcard .fh .fhd{font-weight:400;color:#adb5bd;font-size:10.5px;margin-left:auto}
.chips{display:flex;flex-wrap:wrap;gap:5px}
.chip{display:inline-flex;align-items:center;gap:5px;padding:3px 9px;border-radius:14px;border:1px solid transparent;cursor:pointer;font-size:11px;user-select:none;white-space:nowrap}
.chip:hover{filter:brightness(.95)}.chip input{margin:0;cursor:pointer}
.chip.on{box-shadow:inset 0 0 0 2px currentColor;font-weight:600}
.chip .cnt{font-size:9.5px;background:rgba(0,0,0,.09);border-radius:8px;padding:0 5px;font-weight:700}
.det_A{background:#d3f9d8;color:#2b8a3e}.det_amb{background:#fff3bf;color:#b08900}.det_B{background:#d0ebff;color:#1971c2}.det_C{background:#ffe3e3;color:#e8590c}.det_incompat{background:#ffd6e7;color:#c2255c}.det_other{background:#f1f3f5;color:#868e96}
.detail{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:15px;min-height:76vh;overflow-x:auto}
.detail table{max-width:100%}
.zone{margin:20px 0 6px;padding:8px 14px;border-radius:7px;font-size:14px;font-weight:700;background:#e7f5ff;color:#1971c2;border-left:5px solid #1971c2}
.detail h3{position:sticky;top:0;background:#fff;margin:-15px -15px 8px;padding:12px 15px 8px;z-index:5;border-bottom:1px solid #f1f3f5;border-radius:8px 8px 0 0}
.kv{display:flex;gap:10px;flex-wrap:wrap;margin:6px 0}.kv .b{background:#f1f3f5;border-radius:6px;padding:5px 10px;font-size:11.5px}
table{border-collapse:collapse;font-size:11.5px;margin:7px 0}th,td{border:1px solid #dee2e6;padding:3px 8px}th{background:#f1f3f5}
.note{font-size:11.5px;color:#868e96}.mono{font-family:ui-monospace,Menlo,monospace}
"""

JS = r"""
function bootWS(){
const D=window.__DATA__;
const TT={'linear(全直系)':'t_linear','branched(直系+姊妹)':'t_branched','star(全姊妹)':'t_star','single':'t_single','germline_only':'t_single'};
const el=id=>document.getElementById(id);
const PAL=['#1c7ed6','#37b24d','#f59f00','#ae3ec9','#e8590c','#1098ad','#f03e3e','#adb5bd','#7048e8'];
function bars(o,m){let vs=Object.values(o),tot=vs.reduce((a,b)=>a+b,0)||1,mx=Math.max(...vs,1);return Object.entries(o).sort((a,b)=>b[1]-a[1]).slice(0,m||9).map(([k,v],idx)=>`<div class="bar"><i style="width:${Math.max(3,78*v/mx)}px;background:${PAL[idx%PAL.length]}"></i>${k}: <b>${v}</b> (${(100*v/tot).toFixed(1)}%)</div>`).join('')+`<div class="bar" style="color:#868e96">— 合計 ${tot} (100%)</div>`}
function pie(o){let es=Object.entries(o).sort((a,b)=>b[1]-a[1]),tot=es.reduce((a,b)=>a+b[1],0)||1,cx=30,cy=30,r=26,ri=14,a0=-Math.PI/2,p='';
 if(es.filter(e=>e[1]>0).length<=1){return '<svg width="60" height="60" viewBox="0 0 60 60" style="flex:0 0 auto"><circle cx="30" cy="30" r="26" fill="'+PAL[0]+'"/><circle cx="30" cy="30" r="14" fill="#fff"/></svg>'} // 100% 單類別→整環(避退化弧消失)
 es.forEach(([k,v],i)=>{let a1=a0+2*Math.PI*v/tot,lg=(a1-a0)>Math.PI?1:0,x0=cx+r*Math.cos(a0),y0=cy+r*Math.sin(a0),x1=cx+r*Math.cos(a1),y1=cy+r*Math.sin(a1),xi1=cx+ri*Math.cos(a1),yi1=cy+ri*Math.sin(a1),xi0=cx+ri*Math.cos(a0),yi0=cy+ri*Math.sin(a0);p+=`<path d="M${x0.toFixed(1)} ${y0.toFixed(1)} A${r} ${r} 0 ${lg} 1 ${x1.toFixed(1)} ${y1.toFixed(1)} L${xi1.toFixed(1)} ${yi1.toFixed(1)} A${ri} ${ri} 0 ${lg} 0 ${xi0.toFixed(1)} ${yi0.toFixed(1)} Z" fill="${PAL[i%PAL.length]}"><title>${k}: ${v} (${(100*v/tot).toFixed(1)}%)</title></path>`;a0=a1});return `<svg width="60" height="60" viewBox="0 0 60 60" style="flex:0 0 auto">${p}</svg>`}
function pieBars(o,m){return `<div style="display:flex;gap:9px;align-items:flex-start">${pie(o)}<div style="flex:1;min-width:0">${bars(o,m)}</div></div>`}
el('s_topo').innerHTML=pieBars(D.stats.topology_type);el('s_clust').innerHTML=pieBars(Object.fromEntries(Object.entries(D.stats.n_clusters).map(([k,v])=>['c='+k,v])));
el('s_det').innerHTML=pieBars(D.stats.determinacy);el('s_root').innerHTML=pieBars(D.stats.n_roots);
// R4: 建樹位點分布(第5卡) — 區級從 detail 算 + sSNV 位點宇宙從 accounting
(function(){let dt=D.detail,tot=dt.length;
 let buckets=[['2',dt.filter(r=>r.n_sSNV==2).length],['3',dt.filter(r=>r.n_sSNV==3).length],['4',dt.filter(r=>r.n_sSNV==4).length],['5',dt.filter(r=>r.n_sSNV==5).length],['6',dt.filter(r=>r.n_sSNV==6).length],['7',dt.filter(r=>r.n_sSNV==7).length],['8',dt.filter(r=>r.n_sSNV==8).length],['>8',dt.filter(r=>r.n_sSNV>8).length]];
 let mx=Math.max(...buckets.map(b=>b[1]),1);
 let html='<div class="note" style="margin-bottom:3px">每區 sSNV 數（縱）→ 區域個數（橫）；共 <b>'+tot+'</b> 可建樹區（n_sSNV≥2）。多數為 2-sSNV 對。</div>';
 html+=buckets.map(([k,v])=>{let pct=(100*v/tot).toFixed(1);let c=k=='>8'?'#fa5252':(+k>=4?'#37b24d':'#1c7ed6');return `<div class="bar" title="${v} 區 = 可建樹區 ${pct}%"><span style="display:inline-block;width:30px;text-align:right;font-weight:600">${k}</span> <i style="width:${Math.max(2,150*v/mx)}px;background:${c}"></i> ${v} <span style="color:#868e96">(${pct}%)</span></div>`;}).join('');
 let a=D.accounting,snv_in_trees=dt.reduce((s,r)=>s+r.n_sSNV,0);
 html+='<div class="bar" style="color:#495057;margin-top:5px;border-top:1px solid #f1f3f5;padding-top:4px">相較總量：可建樹區涵蓋 <b>'+snv_in_trees.toLocaleString()+'</b> sSNV-in-region';
 if(a){html+='　/　全宇宙 <b>'+(a.universe_total||0).toLocaleString()+'</b>（linked <b>'+a.buckets.linked.pct+'%</b> 可建樹·單位點 <b>'+a.single_pct+'%</b>·其餘 isolated/underpowered）';}
 html+='</div><div class="note" style="margin-top:2px">藍=2-3 sSNV·綠=≥4(較豐富樹)·紅=>8 截斷</div>';
 el('s_nsnv').innerHTML=html;})();
// situation 分類 — 🔴 完全複製 candidate_scoring.py situation()(L42-48)保證計數與評分佇列一致;Level1+Level2 共用
function regSit(r){ if(r.determinacy=='incompatible')return['衝突(成環)','#f03e3e'];
 if((r.ambig_nodes||0)>0)return['順序 2-3 順位待定','#ffd43b'];
 if(r.determinacy=='C_underdetermined')return['多樹相容(欠定)','#adb5bd'];
 var hp=r.haplotypes; if(hp!='H1'&&hp!='H2'&&hp!='?')return['跨HP(兩棵樹)','#f59f00'];
 if(r.determinacy=='B_pairwise_structure')return['pairwise 拼接','#74c0fc'];
 return['已確定','#37b24d']; }
// 逐區 LOH 判定:longphase-TO LOH.bed(D.loh)overlap 或 cn=loh(補 6 樣本 cn=unknown 看不到 bed LOH 的洞)
let LOHIV={};(D.loh||[]).forEach(function(iv){(LOHIV[iv[0]]=LOHIV[iv[0]]||[]).push(iv)});
function regLOH(r){if(r.cn=='loh')return true;let L=LOHIV[r.chrom];if(!L)return false;let a=r.start,b=r.start+(r.span||0);return L.some(function(iv){return iv[1]<b&&iv[2]>a})}
// topology_type 顯示正名(③):row-laminar 建樹結構上只能「多根 lineage」(各自從 germline 分出),畫不出 subclone-of-subclone
const TTLABEL={'branched(直系+姊妹)':'多根 lineage','linear(全直系)':'linear 直系','star(全姊妹)':'star 多根','single':'single','germline_only':'germline'};
// 癌基因命中(#4):該區 gene 註釋含 COSMIC cancer_genes
function geneHit(region){var g=(D.gene||{})[region];return !!(g&&g.cancer_genes&&Object.keys(g.cancer_genes).length)}
// HP 單倍型家族(ideogram HP 模式):此區 somatic 落在哪條 germline HP
function hpFamily(r){var h=r.haplotypes||'';var h1=h.indexOf('H1')>=0,h2=h.indexOf('H2')>=0;
 if(h1&&h2)return['H1H2(跨HP)','#f59f00']; if(h1)return['H1','#1c7ed6']; if(h2)return['H2','#7048e8']; return['未定相(H3?/?)','#adb5bd'];}
// situation 顯示正名(#8):regSit/JSON key 不動(計數與評分一致),只改給人看的字
var SIT_LABEL={'多樹相容(欠定)':'單ALT群·缺連接read(欠定)','跨HP(兩棵樹)':'跨HP·兩棵樹(allelic)'};
function sitDisp(n){return SIT_LABEL[n]||n;}
// HG38 ideogram: 結果(situation)+LOH 整合單圖(預設) / 樹形(次) + LOH 半透明底帶
window.__ideoMode=window.__ideoMode||'situation';
window.toggleIdeo=function(m){window.__ideoMode=m;renderIdeo();};
function renderIdeo(){var ID=D.ideogram;var host=el('ideogram');if(!host)return;if(!ID){host.innerHTML='';return;}
 var chroms=Object.keys(ID.per_chrom);if(!chroms.length){host.innerHTML='';return;}
 var mode=window.__ideoMode||'situation';
 var maxlen=Math.max.apply(null,chroms.map(function(c){return ID.per_chrom[c].len;}));
 var PXW=860,shapeCol={F:'#2f9e44',S:'#1c7ed6',I:'#e03131',N:'#adb5bd'},t=ID.totals||{};
 var X0=64,X=function(p){return X0+PXW*p/maxlen;};
 var rowOf={};chroms.forEach(function(c,i){rowOf[c]=i;});
 var lohBy={},lohN=(D.loh||[]).length;(D.loh||[]).forEach(function(iv){(lohBy[iv[0]]=lohBy[iv[0]]||[]).push(iv);});
 var sitByChrom={},sitCnt={},snvByChrom={};
 if(mode=='situation'){(D.detail||[]).forEach(function(r,idx){var c=r.chrom;if(rowOf[c]==null)return;var sv=regSit(r);if(!sv)return;
   var st=(r.start!=null?r.start:parseInt(((r.region||'').split(':')[1]||'0').split('-')[0])||0);
   (sitByChrom[c]=sitByChrom[c]||[]).push([st,sv[0],sv[1],idx,r.region,r.cn]);sitCnt[sv[0]]=(sitCnt[sv[0]]||0)+1;});}
 if(mode=='snv'){(D.detail||[]).forEach(function(r,idx){var c=r.chrom;if(rowOf[c]==null)return;var sv=regSit(r);if(!sv)return;
   var ps=Object.keys(r.node_hp||{}).map(function(k){return parseInt(k.split(':')[1])||0}).filter(function(p){return p>0}).sort(function(a,b){return a-b});
   if(!ps.length)return;(snvByChrom[c]=snvByChrom[c]||[]).push({ps:ps,col:sv[1],idx:idx,region:r.region,sit:sv[0]});sitCnt[sv[0]]=(sitCnt[sv[0]]||0)+ps.length;});}
 var hpByChrom={},hpCnt={};
 if(mode=='hp'){(D.detail||[]).forEach(function(r,idx){var c=r.chrom;if(rowOf[c]==null)return;var hf=hpFamily(r);
   (hpByChrom[c]=hpByChrom[c]||[]).push([(r.start||0),hf[0],hf[1],idx,r.region,r.haplotypes]);hpCnt[hf[0]]=(hpCnt[hf[0]]||0)+1;});}
 // CN 帶(gain/loss;loh 已由 LOH 背景紫帶顯示):HCC1395 cn 欄完整,6 樣本 cn=unknown→空。非樹形模式才畫
 var cnByChrom={},cnN=0;if(mode!='shape')(D.detail||[]).forEach(function(r){if(r.cn!='gain'&&r.cn!='loss')return;var c=r.chrom;if(rowOf[c]==null)return;(cnByChrom[c]=cnByChrom[c]||[]).push([(r.start||0),(r.start||0)+(r.span||0),r.cn]);cnN++;});
 var bt=function(m,lab){return '<button onclick="toggleIdeo(\''+m+'\')" style="font-size:11px;padding:2px 10px;border:1px solid #1971c2;cursor:pointer;background:'+(mode==m?'#1971c2':'#fff')+';color:'+(mode==m?'#fff':'#1971c2')+'">'+lab+'</button>';};
 var toggle='<span style="margin-left:8px;display:inline-flex">'+bt('situation','結果+LOH')+bt('snv','結果/每sSNV')+bt('hp','HP單倍型')+bt('shape','樹形')+'</span>';
 var lohLeg=lohN?'　<span style="background:#d0bfff;padding:0 6px;border-radius:2px">▬ LOH '+lohN+' 段</span>('+(D.loh_source||'?')+')':'　<span class="note">LOH:無資料</span>';
 var cnLeg=cnN?'　<span style="color:#ff8787">▮gain</span> <span style="color:#4dabf7">▮loss</span>('+cnN+' 區·HCC1395 cn欄;6樣本 unknown 僅 LOH)':'';
 var sitLegItems='<span style="color:#37b24d">▏已確定 '+(sitCnt['已確定']||0)+'</span> <span style="color:#74c0fc">▏pairwise '+(sitCnt['pairwise 拼接']||0)+'</span> <span style="color:#adb5bd">▏單群欠定 '+(sitCnt['多樹相容(欠定)']||0)+'</span> <span style="color:#f59f00">▏跨HP '+(sitCnt['跨HP(兩棵樹)']||0)+'</span> <span style="color:#ffd43b">▏順序待定 '+(sitCnt['順序 2-3 順位待定']||0)+'</span> <span style="color:#f03e3e">▏衝突 '+(sitCnt['衝突(成環)']||0)+'</span>';
 var hpLegItems='<span style="color:#1c7ed6">▏H1 '+(hpCnt['H1']||0)+'</span> <span style="color:#7048e8">▏H2 '+(hpCnt['H2']||0)+'</span> <span style="color:#f59f00">▏H1H2跨HP '+(hpCnt['H1H2(跨HP)']||0)+'</span> <span style="color:#adb5bd">▏未定相 '+(hpCnt['未定相(H3?/?)']||0)+'</span>';
 var legend = mode=='shape'
   ? '樹位點 <span style="color:#2f9e44">▏full_tree</span> <span style="color:#1c7ed6">▏結構(linear/sibling/co_linked)</span> <span style="color:#e03131">▏成環</span>　密度軌 <span style="color:#fa8c16">▮underpowered '+(t.underpowered||0)+'</span> <span style="color:#868e96">▮isolated '+(t.isolated||0)+'</span>'+lohLeg
   : mode=='snv'
   ? '<b>每個已定相 sSNV 一點</b>·色=該區樹結果;<b>同一區的 sSNV 連一條底線=同一棵樹</b>('+sitLegItems+')'+lohLeg+cnLeg+'　<span class="note">⚠ genome 尺度多數區極小(sSNV 幾乎重疊)→看某區內請點該區用下方 locus track</span>'
   : mode=='hp'
   ? '<b>每區依 germline HP 家族上色</b>（此區 somatic 落在哪條親代染色體;點 tick→跳 detail）：'+hpLegItems+lohLeg+cnLeg
   : '每區依結果上色(點 tick→跳 detail)：'+sitLegItems+lohLeg+cnLeg;
 var s='<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 14px;margin:10px 0;font-size:12px">';
 s+='<b>🗺️ 分布圖</b>'+toggle+'<br><span style="line-height:1.95">'+legend+'</span><span class="note">(染色體比例=GRCh38 真實長度;LOH=半透明紫底帶;hover 看數)</span>';
 s+='<svg viewBox="0 0 '+(PXW+90)+' '+(chroms.length*24+16)+'" width="100%" style="margin-top:6px;font-family:ui-monospace,monospace">';
 chroms.forEach(function(c,i){var pc=ID.per_chrom[c],y=14+i*24,w=PXW*pc.len/maxlen,bin=pc.bin;
   s+='<text x="'+(X0-6)+'" y="'+(y+4)+'" text-anchor="end" font-size="10" font-weight="600">'+c+'</text>';
   s+='<rect x="'+X0+'" y="'+(y-5)+'" width="'+w+'" height="10" rx="2" fill="#f8f9fa" stroke="#dee2e6"/>';
   (lohBy[c]||[]).forEach(function(iv){var x1=X(iv[1]),x2=X(iv[2]);if(x2-x1<0.4)x2=x1+0.4;s+='<rect x="'+x1.toFixed(1)+'" y="'+(y-5)+'" width="'+(x2-x1).toFixed(1)+'" height="10" fill="#9775fa" opacity="0.35"><title>'+c+' LOH '+iv[1]+'-'+iv[2]+'</title></rect>';});
   if(mode!='shape')(cnByChrom[c]||[]).forEach(function(cv){var x1=X(cv[0]),x2=X(cv[1]);if(x2-x1<0.4)x2=x1+0.4;s+='<rect x="'+x1.toFixed(1)+'" y="'+(y+6)+'" width="'+(x2-x1).toFixed(1)+'" height="3" fill="'+(cv[2]=='gain'?'#ff8787':'#4dabf7')+'" opacity="0.85"><title>'+c+' CN '+cv[2]+' '+cv[0]+'-'+cv[1]+'</title></rect>';});
   if(mode=='shape'){var bw=Math.max(1,PXW*bin/maxlen);
     (pc.und_bins||[]).forEach(function(v,bi){if(v>0)s+='<rect x="'+X(bi*bin)+'" y="'+(y+6)+'" width="'+bw+'" height="'+Math.min(6,1+v/2)+'" fill="#fa8c16" opacity="0.75"><title>'+c+' ~'+(bi*2)+'Mb underpowered '+v+'</title></rect>';});
     (pc.iso_bins||[]).forEach(function(v,bi){if(v>0)s+='<rect x="'+X(bi*bin)+'" y="'+(y+13)+'" width="'+bw+'" height="'+Math.min(5,1+v/4)+'" fill="#868e96" opacity="0.5"><title>'+c+' ~'+(bi*2)+'Mb isolated '+v+'</title></rect>';});
     (pc.trees||[]).forEach(function(tr){var x=X(tr[0]);s+='<line x1="'+x+'" y1="'+(y-8)+'" x2="'+x+'" y2="'+(y-1)+'" stroke="'+(shapeCol[tr[1]]||'#adb5bd')+'" stroke-width="1.1"><title>'+c+':'+tr[0]+' '+tr[1]+'</title></line>';});}
   else if(mode=='situation'){(sitByChrom[c]||[]).sort(function(a,b){var o={'已確定':0,'pairwise 拼接':1,'多樹相容(欠定)':2,'順序 2-3 順位待定':3,'跨HP(兩棵樹)':4,'衝突(成環)':5};return (o[a[1]]||0)-(o[b[1]]||0);}).forEach(function(g){var x=X(g[0]);s+='<line x1="'+x+'" y1="'+(y-9)+'" x2="'+x+'" y2="'+(y-1)+'" stroke="'+g[2]+'" stroke-width="1.5" style="cursor:pointer" data-i="'+g[3]+'"><title>'+g[4]+' · '+sitDisp(g[1])+(g[5]&&g[5]!='unknown'?' · CN '+g[5]:'')+'</title></line>';});}
   else if(mode=='snv'){(snvByChrom[c]||[]).forEach(function(g){var x0=X(g.ps[0]),x1=X(g.ps[g.ps.length-1]);
     if(g.ps.length>=2){if(x1-x0<0.6)x1=x0+0.6;s+='<line x1="'+x0.toFixed(1)+'" y1="'+(y-2)+'" x2="'+x1.toFixed(1)+'" y2="'+(y-2)+'" stroke="'+g.col+'" stroke-width="1" opacity="0.55"/>';}
     g.ps.forEach(function(p){var x=X(p);s+='<circle cx="'+x.toFixed(1)+'" cy="'+(y-6)+'" r="1.5" fill="'+g.col+'" style="cursor:pointer" data-i="'+g.idx+'"><title>'+g.region+' · '+sitDisp(g.sit)+' · '+g.ps.length+' sSNV 同樹 · pos '+p+'</title></circle>';});});}
   else if(mode=='hp'){(hpByChrom[c]||[]).forEach(function(g){var x=X(g[0]);s+='<line x1="'+x+'" y1="'+(y-9)+'" x2="'+x+'" y2="'+(y-1)+'" stroke="'+g[2]+'" stroke-width="1.5" style="cursor:pointer" data-i="'+g[3]+'"><title>'+g[4]+' · HP '+g[5]+'</title></line>';});}
 });
 s+='</svg></div>';host.innerHTML=s;
 if(mode=='situation'||mode=='snv'||mode=='hp')host.onclick=function(e){var el2=e.target.closest('[data-i]');if(!el2)return;show(+el2.dataset.i,null);var d=el('detail');if(d)d.scrollIntoView({behavior:'smooth',block:'center'});};
 else host.onclick=null;
}
renderIdeo();
// #10 每樣本記分卡:結果分類(同 regSit→與評分一致) + 截斷 + LOH + 癌基因命中
(function(){var sc=el('scorecard');if(!sc)return;var dt=D.detail||[];var tot=dt.length||1;
 var cnt={};dt.forEach(function(r){var sv=regSit(r);if(sv)cnt[sv[0]]=(cnt[sv[0]]||0)+1;});
 var trunc=dt.filter(function(r){return r.truncated}).length;
 var lohIv={};(D.loh||[]).forEach(function(iv){(lohIv[iv[0]]=lohIv[iv[0]]||[]).push(iv);});
 var lohN=dt.filter(function(r){if(r.cn=='loh')return true;var a=r.start,b=r.start+(r.span||0),L=lohIv[r.chrom];return !!(L&&L.some(function(iv){return iv[1]<b&&iv[2]>a}))}).length;
 var geneN=Object.keys(D.gene||{}).filter(function(rg){var g=(D.gene||{})[rg];return g&&g.cancer_genes&&Object.keys(g.cancer_genes).length}).length;
 var card=function(lab,n,col){return '<div style="flex:1 1 86px;text-align:center;background:#fff;border:1px solid #dee2e6;border-top:3px solid '+col+';border-radius:6px;padding:6px 4px"><div style="font-size:18px;font-weight:700;color:'+col+'">'+n.toLocaleString()+'</div><div style="font-size:10px;color:#666">'+lab+'</div><div style="font-size:9px;color:#aaa">'+(100*n/tot).toFixed(0)+'%</div></div>'};
 sc.innerHTML='<div style="display:flex;gap:6px;flex-wrap:wrap;margin:8px 0">'
  +card('已確定',cnt['已確定']||0,'#37b24d')+card('pairwise',cnt['pairwise 拼接']||0,'#74c0fc')+card('單群欠定',cnt['多樹相容(欠定)']||0,'#adb5bd')+card('跨HP',cnt['跨HP(兩棵樹)']||0,'#f59f00')+card('順序待定',cnt['順序 2-3 順位待定']||0,'#ffd43b')+card('衝突',cnt['衝突(成環)']||0,'#f03e3e')+card('截斷>8',trunc,'#e8590c')+card('LOH 區',lohN,'#9775fa')+card('癌基因命中',geneN,'#e64980')
  +'</div><div class="note" style="margin:-2px 0 6px">📋 <b>每樣本記分卡</b>（共 '+tot+' 區）：結果分類數(與評分佇列同 regSit) + 截斷 + LOH('+(D.loh_source||'無資料')+') + 癌基因命中。</div>';})();
// 三軸觀察(c×拓撲×樣本):characterization·partly-artifact·L3 — 依 20260702 分析報告 + 對抗驗證
(function(){var host=el('threeaxis');if(!host)return;var dt=D.detail||[];
 var TTC={single:'#adb5bd',linear:'#1c7ed6',branched:'#f59f00',incompatible:'#f03e3e'},TTL={single:'single',linear:'linear 直系',branched:'多根 lineage(非嵌套)',incompatible:'成環'};
 var CONF={full_tree:1,linear_nested:1,sibling_only:1,co_linked_lineage:1},order=['single','linear','branched','incompatible'],CB=['1','2','3','4+'];
 var joint={};CB.forEach(function(c){joint[c]={};});
 dt.forEach(function(r){var c=r.n_clusters<=3?String(r.n_clusters):'4+';var t=(r.topology_type||'').split('(')[0];joint[c][t]=(joint[c][t]||0)+1;});
 var maxc=Math.max.apply(null,CB.map(function(c){return Object.values(joint[c]).reduce(function(a,b){return a+b},0)}).concat([1]));
 var bars=CB.map(function(c){var tot=Object.values(joint[c]).reduce(function(a,b){return a+b},0);if(!tot)return '';
   var segs=order.filter(function(t){return joint[c][t]}).map(function(t){var v=joint[c][t];return '<div title="c='+c+' '+TTL[t]+' '+v+'('+Math.round(100*v/tot)+'%)" style="height:'+Math.max(1,140*v/maxc)+'px;background:'+TTC[t]+'"></div>';}).join('');
   return '<div style="text-align:center"><div style="display:flex;flex-direction:column-reverse;width:44px;height:145px;border-bottom:1px solid #ccc">'+segs+'</div><div style="font-size:10px;margin-top:2px">c='+c+'</div><div style="font-size:9px;color:#868e96">'+tot+'</div></div>';}).join('');
 // canonical shape 一致性(此樣本;樹同構標準式) — 答「是否所有樹一致/有效」
 var canonOf=function(edges){var ch={};edges.forEach(function(e){(ch[e[0]]=ch[e[0]]||[]).push(e[1])});var seen={};var rec=function(n){if(seen[n])return '(X)';seen[n]=1;return '('+(ch[n]||[]).map(rec).sort().join('')+')'};return rec('ROOT')};
 var shapeSet={},ntree=0,ninc=0;
 dt.forEach(function(r){if(r.determinacy=='incompatible')ninc++;var e=r.edges||[];if(!e.length)return;ntree++;shapeSet[canonOf(e)]=1;});
 var incpct=(100*ninc/(dt.length||1)).toFixed(1);
 var consist='<div style="background:#fff9db;border:1px solid #ffe08a;border-radius:5px;padding:4px 9px;margin-bottom:6px;font-size:10.5px;color:#7a5c00">📐 <b>枚舉完整</b> ✅：此樣本 <b>'+Object.keys(shapeSet).length+'</b> 種 canonical 樹形（全 7 樣本共 <b>11 種</b>·0 未分類·樹同構式）。 <b>有效性</b>：incompatible（四配子/perfect-phylogeny 違反）<b>'+incpct+'%</b>（'+ninc+'/'+dt.length+'）— 跨樣本 0.4%–19.1%，<b>非「100% 有效」</b>（R1 審計:原「成環 0」是 tautology）。樹空間小=經驗淺薄（深度多 1-2）,非壓縮 n^n。</div>';
 var cross=[];Object.keys(window.__SAMPLES__||{}).forEach(function(s){var sd=window.__SAMPLES__[s];if(!sd||!sd.detail)return;
   var c2=sd.detail.filter(function(r){return r.n_clusters==2});var n=c2.length||1;
   var br=c2.filter(function(r){return (r.topology_type||'').indexOf('branched')>=0}).length,conf=c2.filter(function(r){return CONF[r.tree_shape]}).length;
   cross.push({s:s,n:c2.length,br:(100*br/n).toFixed(1),conf:(100*conf/n).toFixed(1)});});
 cross.sort(function(a,b){return b.br-a.br});
 var FLAG={COLO829:'⚠低coread artifact',HCC1954:'⚠undetermined',H2009:'⚠資料品質'},active=(document.querySelector('.stab.active')||{}).textContent;
 var crossRows=cross.map(function(x){var cur=x.s===active;
   return '<div style="display:flex;align-items:center;gap:5px;margin:2px 0;font-size:11px'+(cur?';font-weight:700;background:#fff9db':'')+'"><span style="width:150px;flex:0 0 auto">'+x.s+(FLAG[x.s]?' <span style="color:#c2255c;font-size:9px">'+FLAG[x.s]+'</span>':'')+'</span><span style="width:34px;text-align:right;color:#f08c00;flex:0 0 auto">'+x.br+'%</span><div style="flex:0 0 100px;background:#f1f3f5;border-radius:3px"><div style="width:'+x.br+'%;height:9px;background:#f59f00;border-radius:3px"></div></div><span style="width:34px;text-align:right;color:#2b8a3e;flex:0 0 auto">'+x.conf+'%</span><div style="flex:0 0 100px;background:#f1f3f5;border-radius:3px"><div style="width:'+x.conf+'%;height:9px;background:#37b24d;border-radius:3px"></div></div><span class="note" style="flex:0 0 auto">n='+x.n+'</span></div>';}).join('');
 host.innerHTML='<div style="background:#fffbf5;border:1px solid #ffe0a3;border-radius:5px;padding:4px 9px;margin-bottom:6px;font-size:10.5px;color:#a37200">🏷️ characterization · single-bulk · 6樣本標籤弱(census) · <b>無 CN 控制</b> · pseudoreplication(真 n=7) · <b>L3/⭐3</b></div>'+consist
  +'<div style="display:flex;gap:22px;flex-wrap:wrap"><div><b style="font-size:12px">本樣本 c × 拓撲(堆疊)</b><div class="note"><span style="color:#adb5bd">■single</span> <span style="color:#1c7ed6">■linear直系</span> <span style="color:#f59f00">■多根lineage</span> <span style="color:#f03e3e">■成環</span></div><div style="display:flex;gap:8px;align-items:flex-end;margin-top:4px">'+bars+'</div></div>'
  +'<div style="flex:1;min-width:360px"><b style="font-size:12px">c=2：<span style="color:#f08c00">非嵌套(branched·幾何上界)%</span> vs <span style="color:#2b8a3e">read-confirmed%</span> — 7 樣本</b>'+crossRows+'</div></div>'
  +'<div class="note" style="margin-top:6px;line-height:1.7">🔴 <b>branched%=幾何上界(2 條非嵌套 ALT 向量)非 read 驗證</b>→看右邊 confirmed%(COLO829 87%→僅 39% 驗證)。跨樣本 CramérV=<b>0.227</b>(small-med;<b>不報 p 值</b>:pseudoreplication·真 n=7)。<b>CN 決定性混淆</b>:HCC1395(唯一有 CN)c=2 branched% 穩健對比 gain <b>67.6%</b> vs LOH <b>27.5%</b>(差 40 點;neutral n=65 小樣本參考),6/7 樣本 cn=unknown 無法 adjust。<b>兩離群成因不同</b>:COLO829=低coread artifact、HCC1954=undetermined。已排除:sSNV 密度、basecaller。此為描述性觀察·<b>非</b>平行 subclone 或生物差異證據。詳見 InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md</div>';})();
// R4: 統計卡 popup(放大 pie+全 bin+名詞字典),onclick 在 bootWS 閉包內→永遠對應當前 D(per-sample)
const STAT_DICT={topology_type:{title:'拓樸型態',desc:'每區 read 群在系統發生樹上的形狀(只計 n_sSNV≥2 區)。🔴 <b>表達力上限</b>:此建樹用 row-laminar(每 read 的 ALT 集合兩兩巢狀/互斥),<b>結構上畫不出「subclone 的 subclone」(深分支)</b> — 實測 1112/1113 個「branched」只在 germline 根分支。故 branched 顯示為「<b>多根 lineage</b>」(各自從 germline 分出),非真正的多層亞克隆樹。要建深分支需單分子整跨全部位點或多樣本/單細胞。',items:{'single':'單群:reads 全塌成一個基因型,無分支','linear(全直系)':'全直系鏈 germline→A→AB→…','branched(直系+姊妹)':'<b>多根 lineage</b>:多條從 germline 各自分出(⚠ 非 subclone-of-subclone·root-laminar 畫不出深分支)','star(全姊妹)':'全姊妹:多條從 germline 各自分出','germline_only':'只有 germline'}},n_clusters:{title:'群數 c（=觀測到的不同 ALT 細胞群數）',desc:'c = 該區內「含 ≥1 個 somatic 突變的不同 genotype 向量」數 = 觀測到的不同 ALT 細胞群 / lineage 狀態（germline 全-R 向量不計）。c=1 即只觀測到單一 ALT 群（無第二群可比較定序→易落 single/欠定）；c=2 是最乾淨的兩 lineage 分離；c 越大分支越複雜。perfect-phylogeny 下 c ≤ k（觀察到 germline root 時，root 不計；k=該區 sSNV 數），實測 c 多為 1-2，故拓樸搜尋空間極小。',items:{}},determinacy:{title:'determinacy 可辨識性',desc:'樹「存在」≠「能辨識是哪棵」。',items:{'A_determined(單分子向量)':'單分子向量唯一可辨識','A_ambiguous_order(缺中間群)':'缺中間群→累積順序未定','B_pairwise_structure':'pairwise 拼接,非單分子整跨','C_underdetermined':'單ALT群·缺連接read(只觀測到單一 ALT 群、無第二群可定序→需深覆蓋;非「多棵樹相容」之意)','incompatible':'四配子違反→成環,無法成單一樹','other':'單群無分支'}},n_roots:{title:'HP 根數',desc:'somatic 事件散在幾條 germline 單倍型。≥2 = 跨 HP(allelic,非 subclone)。',items:{}}};
window.openStatModal=function(which){let o=D.stats[which];if(which==='n_clusters')o=Object.fromEntries(Object.entries(o).map(([k,v])=>['c='+k,v]));let d=STAT_DICT[which]||{title:which,desc:'',items:{}};let tot=Object.values(o).reduce((a,b)=>a+b,0)||1;let dict=Object.entries(d.items||{}).filter(([k])=>o[k]!=null).map(([k,v])=>`<div style="font-size:11.5px;margin:3px 0"><b class="mono">${k}</b> — ${v}</div>`).join('');el('statmodal_body').innerHTML=`<h3 style="margin-top:0">${d.title}（${(window.__DATA__&&document.querySelector('.stab.active')?document.querySelector('.stab.active').dataset.s:'')}；合計 ${tot}）</h3><div class="note" style="margin-bottom:10px">${d.desc}</div><div style="display:flex;gap:24px;align-items:flex-start;flex-wrap:wrap"><div style="transform:scale(1.7);transform-origin:top left;margin:18px 40px 40px 8px">${pie(o)}</div><div style="flex:1;min-width:280px">${bars(o,99)}</div></div>${dict?`<div style="margin-top:12px;border-top:1px solid #f1f3f5;padding-top:9px"><b style="font-size:12.5px">類別說明</b><div style="margin-top:5px">${dict}</div></div>`:''}`;el('statmodal').style.display='flex'};
window.closeStatModal=function(){el('statmodal').style.display='none'};
[['s_topo','topology_type'],['s_clust','n_clusters'],['s_det','determinacy'],['s_root','n_roots']].forEach(([id,key])=>{let sc=el(id).closest('.scard');if(sc)sc.onclick=()=>openStatModal(key)});
(function(){let a=D.accounting,u=el('universe');if(!u)return;if(!a){u.innerHTML='';return;}let B=a.buckets||{},gg=(o,k)=>(o&&o[k]!=null?o[k]:0),nf=x=>(x||0).toLocaleString();
 u.innerHTML=`<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:11px 14px;margin:10px 0;font-size:12.5px"><b>全 sSNV 宇宙帳本（${nf(a.universe_total)} = TP ${nf(a.tp)}（${(100*a.tp/(a.universe_total||1)).toFixed(1)}%）+ FP ${nf(a.fp)}（${(100*a.fp/(a.universe_total||1)).toFixed(1)}%）；隨樣本變）</b><br><div style="display:flex;gap:18px;flex-wrap:wrap;margin-top:6px"><span>🟢 <b>linked ${nf(gg(B.linked,'n'))}（${gg(B.linked,'pct')}%）</b> 可建樹</span><span>🟡 <b>underpowered ${nf(gg(B.underpowered,'n'))}（${gg(B.underpowered,'pct')}%）</b> 有 partner 無共讀→加深覆蓋可救</span><span>⚪ <b>isolated ${nf(gg(B.isolated,'n'))}（${gg(B.isolated,'pct')}%）</b> read-span 內無 partner</span><span class="note">（三桶加總＝宇宙 ✓）</span></div>${(a.derived_from_census||a.fp_truth_sparse||a.cn_all_neutral)?`<div class="note" style="margin-top:6px;color:#a37200">⚠ ${a.derived_from_census?'帳本由 per_sSNV_census 衍生':''}${a.fp_truth_sparse?'；FP truth 稀疏（無外部 truth set→TP/FP 標籤弱，勿當判別力佐證）':''}${a.cn_all_neutral?'；CN 未併 census（by_cn 全 neutral）':''}</div>`:''}</div>`;})();
// 克隆樹判讀圖鑑(R3+R7) — 固定取 HCC1395,渲染一次,與分頁無關
function renderTeaching(){
 let TS=window.__SAMPLES__||{}, H=TS['HCC1395']||TS[Object.keys(TS).find(k=>TS[k]&&TS[k].chr17)];
 let host=el('teachbody'); if(!host)return;
 if(!H||!H.chr17){host.innerHTML='<div class="note">教學資料(HCC1395 chr17)不可用</div>';return;}
 let c=H.chr17, byreg={}; (H.detail||[]).forEach(r=>byreg[r.region]=r);
 let germ=(c.populations.find(p=>p.muts=='germline')||{}).reads||0, germVec=(c.populations.find(p=>p.muts=='germline')||{}).vec;
 let treeVecs=new Set((c.edges||[]).reduce((a,e)=>a.concat(e),[]));
 let st=c.snvs.map(s=>`<tr><td class="mono"><b>${s.S}</b></td><td class="mono">${s.pos}</td><td>${s.change}</td><td>${s.role}</td><td>VAF ${s.vaf}</td><td>${s.hp}</td><td>${s.src}</td><td>${s.somatic_confirmed?'✓':'✗ normal有ALT'}</td></tr>`).join('');
 let pp=c.populations.map(p=>{let onT=treeVecs.has(p.vec)||p.vec==germVec;return `<tr style="${onT?'':'background:#fff0f6'}"><td class="mono">${p.vec}</td><td>${p.muts}</td><td>${p.reads}</td><td>${p.pct}%</td><td class="note">${onT?'':'⚠ 噪聲丟棄(無上樹)'}</td></tr>`}).join('');
 let mm=c.sig_cpg.slice(0,8).map(x=>`<tr><td class="mono">${x.m}</td><td>${x.cpg}</td><td>${x.L1}</td><td>${x.L2}</td><td>${x.dbeta}</td></tr>`).join('');
 let hero=`<div class="kv"><div class="b">locus ${c.locus}</div><div class="b">ctx ${c.genome_ctx}</div><div class="b">拓樸 ${c.topology_type}</div><div class="b">噪聲丟 ${c.dropped_noise} reads</div></div>
  <div style="background:#e7f5ff;border:1px solid #a5d8ff;border-radius:6px;padding:7px 11px;font-size:11.5px;line-height:1.85;margin:6px 0">① <b>germline 根</b>(灰框·全R向量·無somatic·起點) ② <b>直系鏈</b> germline→RRAR(+S3 α祖先)→RAAA(+S2+S4 β後代),S2/S4 <b>co_linked</b> 同步加(順序未定) ③ <b>姊妹分支</b> germline→ARRR(+S1):S1=48357368 是 <b>ClairS FP</b>;因與 α/β 鏈<b>無共現 read</b>→拓樸上自然落為獨立 sibling(演算法<b>純看共現、未用 FP 標籤</b>判別) ④ <b>+S</b> 對到下方 S 表 ⑤ <b>甲基 cis-ASM</b>:m 對齊 genotype 軸(L1=RRAR vs L2=RAAA)非獨立 lineage(見下方展板)</div>
  <b>克隆樹</b>${tree(c.edges,Object.fromEntries(c.populations.map(p=>[p.vec,p.reads])),c.populations.length,'H1',germ,0)}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:6px"><div><b>S 位點(somatic sSNV)</b><table><tr><th>S</th><th>pos</th><th>變異</th><th>角色</th><th>VAF</th><th>HP</th><th>TP/FP</th><th>som</th></tr>${st}</table></div><div><b>細胞群(粉底=噪聲未上樹)</b><table><tr><th>向量</th><th>突變</th><th>reads</th><th>%</th><th></th></tr>${pp}</table></div></div>
  <b>甲基差異位點 m(⚠ cis-ASM:對齊 genotype 軸非獨立 lineage)</b><table><tr><th>m</th><th>CpG</th><th>L1 β</th><th>L2 β</th><th>Δβ</th></tr>${mm}</table>`;
 let card=(t,note,extra)=>`<div style="border:1px solid #dee2e6;border-radius:8px;padding:9px 11px;background:#fff"><div style="font-weight:700;font-size:12.5px;color:#1971c2">${t}</div><div class="note" style="margin:3px 0 6px">${note}</div>${extra}</div>`;
 let TR=r=>tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes);
 let cards='';
 let rA=byreg['chr1:24630300-24635403']; if(rA)cards+=card('卡A · 姊妹分支 sibling','AA=0(兩突變從不共現)+ 同 HP → 平行 subclone、樹上同層分叉(橙框)。',TR(rA));
 let rB=byreg['chr1:5816053-5816054']; if(rB)cards+=card('卡B · co_linked 完美共現','只見 AA、RA=AR=0 → 兩突變永遠同進退、綁同一事件、無法內部排序。',TR(rB));
 let rC=byreg['chr1:12724814-12732374']; if(rC){let cd=enumCandidates(rC);let cl=cd?'<div style="margin-top:6px;font-size:11px;background:#f8f0fc;border-radius:5px;padding:5px 8px">候選累積序(共 '+cd.trueCount+' 種·等機率·中間群未觀察):<br>'+cd.cands.map((cc,ix)=>(ix+1)+'. '+cc.map(a=>((rC.node_paths||{})[a.node]||a.node)+': '+a.order.join('→')).join('；')).join('<br>')+'</div>':'';cards+=card('卡C · 缺中間群 ambiguous','節點一次獲≥2 變異、中間群未觀察 → 累積順序未定(黃框)。可能順序等機率列出:',TR(rC)+cl);}
 let rD=byreg['chr2:97211773-97229016']; if(rD){let cf=fourGamete(rD.populations);cards+=card('卡D · 四配子衝突 four-gamete','RR/RA/AR/AA 四格全有 → 不相容單一樹。錨 RR=germline → AA=雙突變最遠。',cf.length?'<table><tr><th>對</th><th>RR</th><th>RA</th><th>AR</th><th>AA</th></tr>'+cf.map(x=>'<tr><td class="mono"><b>'+x.pair+'</b></td><td>'+x.g.RR+'</td><td>'+x.g.RA+'</td><td>'+x.g.AR+'</td><td>'+x.g.AA+'</td></tr>').join('')+'</table>':'<div class="note">(此區四配子未觸發)</div>');}
 let rE=byreg['chr7:100784203-100798014']; if(rE)cards+=card('卡E · CROSS-HP allelic(對比卡A)','n_roots≥2:突變散在 H1/H2 兩單倍型=各自染色體突變(allelic)非 subclone → 畫兩棵分開 HP 樹。',posTree(rE)||TR(rE));
 let mt=window.__CHR17TREE__, exhibit='';
 if(mt&&mt.cross_clu2_x_lineage){let cx=mt.cross_clu2_x_lineage,lins=['L0','L1','L2','other'];
  let xrows=Object.keys(cx).sort().map(cl=>'<tr><td><b>甲基群 '+cl+'</b></td>'+lins.map(L=>'<td style="'+((cx[cl][L]||0)>=10?'background:#ffe3e3;font-weight:700':'')+'">'+(cx[cl][L]||0)+'</td>').join('')+'</tr>').join('');
  let ax=mt.axis_sig_cpg_count||{},axmax=Math.max.apply(null,Object.values(ax).concat([1]));
  let axbars=Object.entries(ax).map(([k,v])=>'<div class="bar"><i style="width:'+Math.max(3,90*v/axmax)+'px;background:'+(k.indexOf('lineage')>=0?'#e8590c':(k.indexOf('HP')>=0?'#868e96':'#1c7ed6'))+'"></i>'+k+': <b>'+v+'</b></div>').join('');
  exhibit='<div style="background:#fff5f5;border:1px solid #ffc9c9;border-radius:8px;padding:10px 13px;margin-top:12px"><b>🧬 甲基 read 距離+分群「能/不能做什麼」(chr17 worked,'+mt.n_reads+' reads × '+mt.n_cpg+' CpG)</b><div class="note" style="margin:4px 0">把 read 用甲基距離分 k=2 群,再對遺傳 lineage 交叉:</div><table><tr><th>甲基分群＼遺傳lineage</th>'+lins.map(L=>'<th>'+L+'</th>').join('')+'</tr>'+xrows+'</table><div style="color:#c92a2a;font-weight:600;margin:6px 0;font-size:11.5px">🔴 甲基群「1」同時含 L1('+((cx['1']&&cx['1'].L1!=null)?cx['1'].L1:'?')+') 與 L2('+((cx['1']&&cx['1'].L2!=null)?cx['1'].L2:'?')+') → <b>甲基分群 ≠ 遺傳 lineage</b>(cis-ASM double-dip 本體);甲基不能 recover subclone。</div><div class="note">各軸顯著 CpG 數(甲基「對齊」哪個軸):</div>'+axbars+'<div class="note" style="margin-top:5px">→ α(genotype)軸最強、lineage 軸弱、HP 軸 0 ⇒ 甲基是 <b>ASM 存在性偵測器</b>非 lineage 排序器。基因組級 PERMANOVA 幾乎無 testable 區、recover≈0(統計上無力 recover subclone)。<b>可用=負篩/L3弱旗標/教學;不可用=分群器/排序器/定群器</b>(06-28 cis-control 裁決,見 20260628_cis_control_scope_pilot_verdict)。</div></div>';}
 host.innerHTML='<div class="note" style="margin-bottom:8px">節點=一種細胞狀態(read 群);往下=直系、同層分叉=姊妹;邊上 +S=新增 somatic 變異;根=germline。<b>以下範例固定取自 HCC1395,與上方分頁樣本無關。</b></div><h4 style="margin:6px 0">① HERO:chr17:48357368-48365161 逐元素判讀（綁 production 真實區·4 sSNV）</h4>'+hero+'<h4 style="margin:16px 0 6px">② 五種 SNV 關係配套真例(各一真實區)</h4><div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:8px">'+cards+'</div>'+exhibit+'<h4 style="margin:16px 0 6px">③ 標籤圖例</h4><div class="note" style="line-height:1.9"><span style="color:#1565c0;font-weight:700">●直系 vertical</span> 後代多帶1變異 ｜ <span style="color:#d9480f;font-weight:700">●姊妹 sibling</span> 同父平行 subclone ｜ <span style="color:#5f3dc4;font-weight:700">co_linked</span> 完美共現綁同事件 ｜ <span style="color:#e8590c;font-weight:700">⚠缺中間群</span> 順序未定 ｜ <span style="color:#c2255c;font-weight:700">四配子違反</span> 成環無法成樹 ｜ <span style="color:#2b8a3e;font-weight:700">+S</span> 新增變異';
}
// R6 草圖: 無法定位群偵測總覽(per-sample;「群存在但定不出位置/順序」)
(function(){let u=el('unlocatable');if(!u)return;let sm=(D.scoring&&D.scoring.summary)||{},sd=sm.situation_dist;if(!sd){u.innerHTML='';return;}
 let tot=Object.values(sd).reduce((a,b)=>a+b,0)||1;
 let order=[['已確定','#37b24d'],['pairwise 拼接','#74c0fc'],['多樹相容(欠定)','#adb5bd'],['跨HP(兩棵樹)','#f59f00'],['順序 2-3 順位待定','#ffd43b'],['衝突(成環)','#f03e3e']];
 let seg=order.filter(([k])=>sd[k]).map(([k,c])=>`<div title="${k}: ${sd[k]}（${(100*sd[k]/tot).toFixed(1)}%）" style="width:${100*sd[k]/tot}%;background:${c}"></div>`).join('');
 let unloc=[['衝突(成環)','#f03e3e','四配子違反 → 無單一樹(多為 CN-gain 假環)'],['跨HP(兩棵樹)','#f59f00','突變散在 H1/H2 → 非單一譜系、是兩棵樹(allelic)'],['多樹相容(欠定)','#adb5bd','只觀測到單一 ALT 群(無第二群可定序;非「多棵樹相容」之意)→需深覆蓋'],['順序 2-3 順位待定','#ffd43b','缺中間群 → 累積順序未定']];
 let rows=unloc.filter(([k])=>sd[k]!=null).map(([k,c,d])=>`<div onclick="filterSituation('${k}')" style="display:flex;gap:8px;align-items:baseline;margin:3px 0;font-size:11.5px;cursor:pointer;border-radius:4px;padding:2px 6px" onmouseover="this.style.background='#fff0f6'" onmouseout="this.style.background=''"><span style="display:inline-block;width:11px;height:11px;border-radius:3px;background:${c};flex:0 0 auto"></span><b style="width:150px;flex:0 0 auto">${sitDisp(k)}</b><b style="color:${c};width:46px;flex:0 0 auto">${sd[k]}</b><span class="note">${d}</span><span class="note" style="margin-left:auto;color:#1971c2;flex:0 0 auto;white-space:nowrap">▸ 篩選佇列</span></div>`).join('');
 window.filterSituation=function(sit){let s=el('q_sit');if(!s)return;if([...s.options].some(o=>o.value===sit)){s.value=sit;if(typeof renderQ==='function')renderQ();}let qq=el('queue');if(qq)qq.scrollIntoView({behavior:'smooth',block:'center'});};
 let q=(D.scoring.queue||[]),withp=q.filter(x=>x.parsimony_first_rank_prob!=null),hi=withp.filter(x=>x.parsimony_first_rank_prob>=0.7).length;
 let nunloc=(sd['衝突(成環)']||0)+(sd['跨HP(兩棵樹)']||0)+(sd['多樹相容(欠定)']||0)+(sd['順序 2-3 順位待定']||0);
 u.innerHTML=`<div class="zone" style="background:#fff0f6;color:#c2255c;border-left-color:#f06595">🔍 無法定位群偵測（草圖）—「群存在但定不出位置/順序」(此樣本 ${nunloc} 區)</div>
  <div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 13px;font-size:12.5px">
  <div style="display:flex;height:16px;border-radius:5px;overflow:hidden;margin-bottom:7px">${seg}</div>
  <div class="note" style="margin-bottom:6px">合計 ${tot} 區｜🟢已確定 ${sd['已確定']||0}｜🔵pairwise 拼接 ${sd['pairwise 拼接']||0}(可建樹非單分子)｜<b style="color:#c2255c">↓ 4 種無法定位</b></div>${rows}
  <div style="margin-top:8px;border-top:1px solid #f1f3f5;padding-top:7px;font-size:11.5px"><b>機率(誠實兩軌)</b>：遺傳 parsimony 高信度(≥0.7)= <b>${withp.length?hi:'未回填'}</b>${withp.length?` 區 / ${withp.length} 個有 parsimony 值(順序待定類·非全 ${nunloc} 無法定位區) → ${hi?'':'連遺傳都信心不足'}`:'（此樣本待上游回填）'}　｜　🧬 甲基：SAME-HP <b>不給機率</b>(cis-ASM double-dip)、CROSS-HP 弱(~35%)、乾淨可用≈0（06-28）</div>
  <div class="note" style="margin-top:5px;color:#a37200">⚠ 「不同的樹」全列舉(替代整樹)= 上游 <b>enumerate_candidate_trees</b>(已實作;B1 缺中間群區 detail 有候選樹 carousel·等機率);成環區 detail 另有「可能固定樹」carousel(假環示意)。此面板誠實顯示「什麼定不出來」,不假裝在列舉答案。</div></div>`;})();
if(!window.__teachDone){renderTeaching();window.__teachDone=true;}
// S-label 化基因型向量
function sLabels(g){let s=[...g].map((c,i)=>c=='A'?('S'+(i+1)):null).filter(Boolean);return s.length?s.join('+'):'germline'}
function gainedS(parent,child){let g=[];for(let i=0;i<child.length;i++){if(child[i]=='A'&&(!parent||parent[i]!='A'))g.push('S'+(i+1))}return g}
// 四配子檢定:對每對 sSNV(i,j)數 read 在兩位點的 RR/RA/AR/AA。四格全有=incompatible(無限位點違反)。RR=germline(normal錨REF)、AA=雙突變(最遠)。
function fourGamete(pops){let vs=Object.keys(pops||{}).filter(v=>v);if(!vs.length)return [];let L=vs[0].length,out=[];for(let i=0;i<L;i++)for(let j=i+1;j<L;j++){let g={RR:0,RA:0,AR:0,AA:0};vs.forEach(v=>{let k=v[i]+v[j];if(g[k]!=null)g[k]+=pops[v]});if(g.RR&&g.RA&&g.AR&&g.AA)out.push({pair:'S'+(i+1)+'–S'+(j+1),i:i+1,j:j+1,g:g})}return out}
// 候選累積順序列舉:缺中間群節點(獲得≥2突變)→列舉可能累積序(等機率,中間群未觀察→無讀數可分)。前端可做;ranked 替代整樹需上游 enumerate_candidate_trees。
function fact(n){let r=1;for(let i=2;i<=n;i++)r*=i;return r}
function enumCandidates(r){let par={};(r.edges||[]).forEach(([p,c])=>{if(p!='ROOT')par[c]=p});let amb=[];Object.keys(r.populations).forEach(n=>{let g=gainedS(par[n]||null,n);if(g.length>=2)amb.push({node:n,gained:g})});if(!amb.length)return null;function perms(a){if(a.length<=1)return [a];let o=[];a.forEach((x,i)=>{perms(a.slice(0,i).concat(a.slice(i+1))).forEach(p=>o.push([x].concat(p)))});return o}let trueCount=amb.reduce((acc,a)=>acc*fact(a.gained.length),1);let perNode=amb.map(a=>({node:a.node,orders:a.gained.length<=4?perms(a.gained):[a.gained]}));let cands=[[]];perNode.forEach(pn=>{let nx=[];cands.forEach(c=>pn.orders.forEach(o=>{if(nx.length<24)nx.push(c.concat([{node:pn.node,order:o}]))}));cands=nx});return {cands:cands,trueCount:trueCount,bigNode:amb.some(a=>a.gained.length>4)}}
function candCard(cs,idx,r){let c=cs[idx],np=r.node_paths||{};let body=c.map(a=>`<div style="margin:3px 0"><b class="mono" style="color:#1971c2">${np[a.node]||a.node}</b>（${a.node}）：${a.order.map((s,i)=>`<span style="background:#ebfbee;border:1px solid #2f9e44;border-radius:9px;padding:1px 7px;color:#2b8a3e;font-weight:700;margin:0 1px">${i+1}.${s}</span>`).join('→')}</div>`).join('');return `<div style="display:flex;align-items:center;gap:10px"><button onclick="candNav(-1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px">◀</button><div style="flex:1">${body}<div class="note" style="margin-top:3px">第 ${idx+1} / ${cs.length} 個候選</div></div><button onclick="candNav(1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px">▶</button></div>`}
window.candNav=function(d){let s=window.__cand;if(!s)return;s.idx=(s.idx+d+s.cands.length)%s.cands.length;let box=document.getElementById('candbox');if(box)box.innerHTML=candCard(s.cands,s.idx,s.r)}
// Part B: 上游 enumerate 的完整候選樹 carousel(虛擬中間節點;equiprobable 誠實標)
function candTreeCard(ctr,idx,r){let c=ctr.candidate_set[idx];let th=tree(c.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes);let single=ctr.candidate_set.length<=1;let bs='font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px;margin-top:46px';let nav1=single?'':`<button onclick="ctNav(-1)" style="${bs}">◀</button>`,nav2=single?'':`<button onclick="ctNav(1)" style="${bs}">▶</button>`;let lab=single?`<b style="color:#2b8a3e">✅ 唯一確認樹</b>（resolution=${ctr.resolution||'—'}·非多候選）`:`候選樹 <b>${idx+1} / ${ctr.candidate_set.length}</b>　softmax ${(c.softmax_prob*100).toFixed(0)}%${ctr.equiprobable?'　⚖ <b>等機率</b>(無 read 證據可排)':''}`;return `<div style="display:flex;align-items:flex-start;gap:8px">${nav1}<div style="flex:1;min-width:0">${th}<div class="note" style="margin-top:4px">${lab}　虛擬中間節點(0 reads·未觀察): ${c.virtual_nodes.length?c.virtual_nodes.join('、'):'無'}</div></div>${nav2}</div>`}
window.ctNav=function(d){let s=window.__ctr;if(!s)return;let n=s.cs.candidate_set.length;s.idx=(s.idx+d+n)%n;let box=document.getElementById('cttreebox');if(box)box.innerHTML=candTreeCard(s.cs,s.idx,s.r)}
function tree(edges,popcount,nc,hp,germR,np,ambig){np=np||{};ambig=ambig||0;
 let germKey=Object.keys(popcount).find(k=>k&&/^R+$/.test(k)); // germline=全-R向量;germline_reads欄位不可靠(很多區為0),改用 populations 全-R count 與下方表格一致
 let germN=germKey!=null?popcount[germKey]:(germR||0);
 let ch={},par={},all=new Set();
 (edges||[]).forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT'){all.add(p);par[c]=p}});
 if(!all.size)return `<div class="note" style="background:#fff9db;border:1px solid #ffe0a3;border-radius:5px;padding:6px 9px">⚠ <b>此區無 genotype-向量克隆樹</b>：觀測 read 基因型全塌成單一群（germline ${germN} reads）、向量層無分支。常見原因:① 只觀測到單一 ALT 群(缺第二群可分支) ② <b>genotype 向量截斷</b>(僅取前 8 位點→多位點塌成同一截斷向量,如全 A)。→ 改看上方 <b>locus 排列</b> 與 <b>全位點分群樹</b>(位點層·突破 8 上限)。</div>`;
 let sib={};Object.keys(ch).forEach(p=>{(ch[p]||[]).forEach(c=>{sib[c]=(ch[p].length>=2)})});
 let depth={},pos={},leaf=0,seen={};
 function lay(n,dp){if(seen[n])return pos[n]!=null?pos[n]:0;seen[n]=1;depth[n]=dp;let k=(ch[n]||[]).filter(x=>!seen[x]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]} // seen 防護:成環邊跳過,避免無限遞迴 stack overflow
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));[...all].forEach(n=>{if(pos[n]==null)lay(n,1)}); // 孤兒節點補位(與 posSVG 一致硬化,防 NaN)
 let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]||0),1);
 let NW=210,NH=92,GX=52,GY=140;
 let W=Math.max(450,(leaf||1)*(NW+GX)),H=108+md*GY;
 let X=p=>34+p*(NW+GX)+NW/2,Y=dp=>56+dp*GY;
 let totR=Object.values(popcount).reduce((a,b)=>a+b,0)||1; // popcount 已含 germline 向量,勿再加 germR(否則 germline 雙算→比例與下方表格對不上)
 let relSet=new Set();
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" style="font-family:ui-monospace,Menlo,monospace">`;
 roots.forEach(r=>{let g=gainedS(null,r),mx=(X(gx)+X(pos[r]))/2,my=(Y(0)+Y(1))/2,w=Math.max(44,g.join('+').length*7+12);
  s+=`<line x1="${X(gx)}" y1="${Y(0)+NH/2}" x2="${X(pos[r])}" y2="${Y(1)-NH/2}" stroke="${g.length>=2?'#f08c00':'#ced4da'}" stroke-width="${g.length>=2?(2.2+1.5*g.length).toFixed(1):'1.6'}"/>`;
  if(g.length)s+=`<rect x="${mx-w/2}" y="${my-10}" width="${w}" height="19" rx="9" fill="#ebfbee" stroke="#2f9e44"/><text x="${mx}" y="${my+4}" text-anchor="middle" font-size="11" fill="#2b8a3e" font-weight="700">+${g.join('+')}</text>`;}); // germline→第一代 edge 也標 +S(原缺)
 (edges||[]).forEach(([p,c])=>{if(p!='ROOT'){
  let g=gainedS(p,c),mx=(X(pos[p])+X(pos[c]))/2,my=(Y(depth[p])+Y(depth[c]))/2;
  s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+NH/2}" x2="${X(pos[c])}" y2="${Y(depth[c])-NH/2}" stroke="${g.length>=2?'#f08c00':'#ced4da'}" stroke-width="${g.length>=2?(2.2+1.5*g.length).toFixed(1):'1.6'}"/>`;
  let w=Math.max(44,g.join('+').length*7+12);if(g.length)s+=`<rect x="${mx-w/2}" y="${my-10}" width="${w}" height="19" rx="9" fill="#ebfbee" stroke="#2f9e44"/><text x="${mx}" y="${my+4}" text-anchor="middle" font-size="11" fill="#2b8a3e" font-weight="700">+${g.join('+')}</text>`;
 }});
 let gp=(100*germN/totR).toFixed(0);
 let gtop=Y(0)-NH/2;
 s+=`<rect x="${X(gx)-NW/2}" y="${gtop}" width="${NW}" height="${NH}" rx="12" fill="#f1f3f5" stroke="#868e96" stroke-width="2.5"/>`;
 s+=`<path d="M${X(gx)-NW/2} ${gtop+24} v-12 a12 12 0 0 1 12 -12 h${NW-24} a12 12 0 0 1 12 12 v12 Z" fill="#868e96"/>`;
 s+=`<text x="${X(gx)}" y="${gtop+17}" text-anchor="middle" font-size="14" font-weight="800" fill="#fff">⌂ germline（${hp||'根'}）</text>`;
 s+=`<text x="${X(gx)}" y="${gtop+46}" text-anchor="middle" font-size="10.5" fill="#868e96">無 somatic 變異 · 起點</text>`;
 s+=`<line x1="${X(gx)-NW/2+12}" y1="${gtop+78}" x2="${X(gx)+NW/2-12}" y2="${gtop+78}" stroke="#868e96" stroke-opacity="0.28"/>`;
 s+=`<text x="${X(gx)}" y="${gtop+89}" text-anchor="middle" font-size="12.5" font-weight="700" fill="#212529">${germN} reads · ${gp}%</text>`;
 nodes.forEach(n=>{
  let cnt=popcount[n]||0,pct=(100*cnt/totR).toFixed(0),lab=np[n]||'—';
  let g=gainedS(par[n],n);
  let isSib=sib[n], isMulti=g.length>=2, isAmbig=isMulti&&ambig>0, isCo=isMulti&&!isAmbig;
  if(isSib)relSet.add('姊妹(sibling)'); else if(g.length)relSet.add('直系(vertical)');
  if(isAmbig)relSet.add('⚠缺中間群(順序未定)'); else if(isCo)relSet.add('co_linked(完美共現)');
  let multiTag=isAmbig?' · ⚠缺中間群(順序未定)':(isCo?' · co_linked(完美共現)':'');
  let rel=g.length?((isSib?'姊妹分支(sibling)':'直系(vertical)')+multiTag):'（與父同型）';
  let gtxt=g.length?('獲得 +'+g.join('+')):'';
  let x=X(pos[n]),y=Y(depth[n]);
  let fill=isSib?'#fff4e6':'#e7f5ff', stroke=isSib?'#e8590c':'#1c7ed6';
  if(isAmbig){fill='#fff9db';stroke='#f08c00';}
  let top=y-NH/2,relCol=isSib?'#d9480f':(isAmbig?'#e8590c':'#2b8a3e');
  s+=`<rect x="${x-NW/2}" y="${top}" width="${NW}" height="${NH}" rx="12" fill="${fill}" stroke="${stroke}" stroke-width="2.5"/>`;
  // lineage 標籤頭條(圓角上緣填色)
  s+=`<path d="M${x-NW/2} ${top+24} v-12 a12 12 0 0 1 12 -12 h${NW-24} a12 12 0 0 1 12 12 v12 Z" fill="${stroke}"/>`;
  s+=`<text x="${x}" y="${top+17}" text-anchor="middle" font-size="14" font-weight="800" fill="#fff">${lab}</text>`;
  s+=`<text x="${x}" y="${top+40}" text-anchor="middle" font-size="10" font-weight="700" fill="${relCol}">${rel}</text>`;
  if(gtxt)s+=`<rect x="${x-Math.max(30,gtxt.length*5.4)}" y="${top+46}" width="${2*Math.max(30,gtxt.length*5.4)}" height="16" rx="8" fill="#ebfbee" stroke="#40c057"/><text x="${x}" y="${top+58}" text-anchor="middle" font-size="10.5" font-weight="700" fill="#2b8a3e">${gtxt}</text>`;
  s+=`<text x="${x}" y="${top+74}" text-anchor="middle" font-size="9" fill="#868e96" font-family="ui-monospace,monospace">${n}（=${sLabels(n)}）</text>`;
  s+=`<line x1="${x-NW/2+12}" y1="${top+78}" x2="${x+NW/2-12}" y2="${top+78}" stroke="${stroke}" stroke-opacity="0.28"/>`;
  s+=`<text x="${x}" y="${top+89}" text-anchor="middle" font-size="12.5" font-weight="700" fill="#212529">${cnt} reads · ${pct}%</text>`;
 });
 s+='</svg>';
 let leg=`<details style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;padding:5px 11px;margin-top:5px;font-size:11px"><summary style="cursor:pointer;font-weight:600">🔖 SNV 關係圖例（點開）${relSet.size?'：此樹有 '+[...relSet].join('、'):''}</summary><div style="line-height:1.7;margin-top:5px">
 <span style="color:#1565c0;font-weight:700">●直系 vertical</span> 往下一層、後代多帶 1 變異(+S，藍框)<br>
 <span style="color:#d9480f;font-weight:700">●姊妹分支 sibling</span> 同一父不同分支、平行 subclone(橙框)<br>
 <span style="color:#5f3dc4;font-weight:700">co_linked(完美共現)</span> 一節點獲≥2 變異且區無 ambiguous → 兩變異綁同一事件<br>
 <span style="color:#e8590c;font-weight:700">⚠缺中間群(順序未定)</span> 一節點獲≥2 變異但區有 ambiguous → 跳>1突變、未觀察到中間群、累積順序未定（黃框；<b>非</b>確定 co_linked）<br>
 <span style="color:#2b8a3e;font-weight:700">+S</span> 該分支新增 somatic 變異<br>
 <span style="color:#f08c00;font-weight:700">━ 粗橙線</span> 一次獲 ≥2 突變（缺中間群／co_linked；線越粗代表該 edge 增加越多突變）</div></details>`;
 return s+leg;
}
// 2-root: 位置樹按 HP 分兩棵
function hasCycleEdges(edges){let ch={};edges.forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c)});let col={},cyc=false;function dfs(u){col[u]=1;for(let v of (ch[u]||[])){if(col[v]==1){cyc=true;return}if(!col[v])dfs(v);if(cyc)return}col[u]=2}for(let n of Object.keys(ch)){if(!col[n])dfs(n);if(cyc)break}return cyc}
function posTable(ns,vaf,hp,nested){let parent={};nested.forEach(([a,b])=>{(parent[b]=parent[b]||[]).push(a.split(':')[1])});let rows=[...ns].sort().map(n=>{let v=vaf[n],par=(parent[n]||[]).join('、');return `<tr><td class="mono">${n.split(':')[1]}</td><td>${v!=null?v:'?'}</td><td class="note">${par?'巢狀於 '+par:'—(根/無上游)'}</td></tr>`}).join('');return `<table style="font-size:10.5px"><tr><th>${hp} 位點</th><th>VAF</th><th>巢狀於</th></tr>${rows}</table>`}
function posBarStrip(ns,vaf){let arr=[...ns].sort((a,b)=>(+a.split(':')[1])-(+b.split(':')[1]));let mx=Math.max(...arr.map(n=>vaf[n]||0),0.01);return '<div style="display:flex;gap:3px;overflow-x:auto;padding:4px 0;max-width:100%">'+arr.map(n=>{let v=vaf[n]||0,h=Math.max(4,44*v/mx),p=(n.split(':')[1]||n);return '<div style="flex:0 0 auto;width:32px;text-align:center" title="'+p+' VAF '+v+'"><div style="height:48px;display:flex;align-items:flex-end;justify-content:center"><div style="width:16px;height:'+h+'px;background:#e8590c;border-radius:2px 2px 0 0"></div></div><div style="font-size:7.5px;color:#868e96;white-space:nowrap;overflow:hidden">'+p.slice(-5)+'</div><div style="font-size:8px;color:#495057">'+v+'</div></div>'}).join('')+'</div>'}
function posTree(r){
 let byhp={};Object.entries(r.node_hp||{}).forEach(([p,h])=>{if(h=='H1'||h=='H2')(byhp[h]=byhp[h]||[]).push(p)});
 let dropH3=Object.values(r.node_hp||{}).filter(h=>h!='H1'&&h!='H2').length;
 let hps=Object.keys(byhp).sort();if(hps.length<2)return '';
 let gk=Object.keys(r.populations||{}).find(k=>k&&/^R+$/.test(k)),germN=gk?r.populations[gk]:0,totR=Object.values(r.populations||{}).reduce((a,b)=>a+b,0)||1; // #5 germline/全區 read 數
 let head='<div class="note" style="width:100%;margin-bottom:2px">節點顯示 <b>VAF</b>(該位點變異分率;位點層無 per-位點 read 數)；germline <b>'+germN+' reads</b>·全區 <b>'+totR+' reads</b>；各 lineage 絕對 read 數見下方「細胞群」表。</div>';
 return '<div style="display:flex;gap:20px;flex-wrap:wrap">'+head+hps.map(h=>{
  let ns=new Set(byhp[h]);let ned=(r.pos_nested||[]).filter(e=>ns.has(e[0])&&ns.has(e[1]));
  let body;
  if(ns.size>24){body='<div class="note" style="color:#a37200;margin:3px 0">'+h+' 共 '+ns.size+' 位點(過多→改左右滑桿看各位點 VAF·基因組順序):</div>'+posBarStrip(ns,r.pos_vaf||{});}
  else if(hasCycleEdges(ned)){body='<div class="note" style="color:#c2255c;margin:3px 0">⚠ 此 '+h+' 位點 pairwise nested <b>成環/互指(incompatible)</b>→ 無法成單一樹,改列位點+VAF 表:</div>'+posTable(ns,r.pos_vaf||{},h,ned);}
  else{let hasp=new Set(ned.map(e=>e[1]));let edges=ned.map(e=>[e[0],e[1]]);[...ns].forEach(n=>{if(!hasp.has(n))edges.unshift(['ROOT',n])});body=posSVG(edges,ns,r.pos_vaf||{},h,null,germN);}
  return '<div style="min-width:0;max-width:100%"><b style="color:#d9480f">'+h+' 樹（'+ns.size+' 位點）</b>'+body+'</div>';
 }).join('')+(dropH3?'<div class="note" style="color:#a37200;width:100%;margin-top:3px">⚠ 另有 '+dropH3+' 個未定相(H3?/somatic-ALT)位置未納入上面 HP 樹</div>':'')+'</div>';
}
// 位置層克隆樹(豐富樣式)。hpLabel=整棵 HP;hpMap=node→HP 多HP上色;germN=germline read 數(#5)。過寬→水平捲動(#4)
function posSVG(edges,ns,vaf,hpLabel,hpMap,germN){
 let HPC={H1:'#1c7ed6',H2:'#7048e8'},rootCol=HPC[hpLabel]||'#495057';
 let ch={},all=new Set();edges.forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT')all.add(p)});
 if(!all.size)return '<div class="note">無結構</div>';
 let depth={},pos={},leaf=0,seen={};function lay(n,dp){if(seen[n])return pos[n]!=null?pos[n]:0;seen[n]=1;depth[n]=dp;let k=(ch[n]||[]).filter(x=>!seen[x]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]} // seen 防護:pos_nested 成環(incompatible)時跳過,避免 stack overflow
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));[...all].forEach(n=>{if(pos[n]==null)lay(n,1)});let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]),1),NW=124,NH=52,GX=20,GY=84;
 let W=Math.max(220,(leaf||1)*(NW+GX)),H=80+md*GY,X=p=>22+p*(NW+GX)+NW/2,Y=dp=>34+dp*GY;
 let s=`<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" style="font-family:ui-monospace,monospace">`;
 s+=`<rect x="${X(gx)-54}" y="${Y(0)-20}" width="108" height="42" rx="9" fill="#f1f3f5" stroke="${rootCol}" stroke-width="2"/><text x="${X(gx)}" y="${Y(0)-5}" text-anchor="middle" font-size="11" font-weight="800" fill="${rootCol}">${hpLabel||'germ'} germline</text><text x="${X(gx)}" y="${Y(0)+7}" text-anchor="middle" font-size="8" fill="#868e96">起點·無 somatic</text>${germN!=null?`<text x="${X(gx)}" y="${Y(0)+18}" text-anchor="middle" font-size="9.5" font-weight="700" fill="#212529">${germN} reads</text>`:''}`;
 roots.forEach(rt=>s+=`<line x1="${X(gx)}" y1="${Y(0)+22}" x2="${X(pos[rt])}" y2="${Y(1)-NH/2}" stroke="#ced4da" stroke-width="1.6"/>`);
 edges.forEach(([p,c])=>{if(p!='ROOT')s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+NH/2}" x2="${X(pos[c])}" y2="${Y(depth[c])-NH/2}" stroke="#ced4da" stroke-width="1.6"/>`});
 nodes.forEach(n=>{let v=vaf[n],hp=(hpMap&&hpMap[n])||hpLabel,col=HPC[hp]||'#e8590c',fill=hp=='H1'?'#dbe9f7':hp=='H2'?'#ece0f5':'#fff4e6',x=X(pos[n]),y=Y(depth[n]);
  s+=`<rect x="${x-NW/2}" y="${y-NH/2}" width="${NW}" height="${NH}" rx="9" fill="${fill}" stroke="${col}" stroke-width="2"/>`;
  s+=`<rect x="${x-NW/2}" y="${y-NH/2}" width="${NW}" height="17" rx="9" fill="${col}"/><rect x="${x-NW/2}" y="${y-NH/2+8}" width="${NW}" height="9" fill="${col}"/>`;
  s+=`<text x="${x}" y="${y-NH/2+13}" text-anchor="middle" font-size="10" font-weight="800" fill="#fff">${hpMap?(hp||'HP?')+'·':''}${n.split(':')[1]}</text>`;
  s+=`<text x="${x}" y="${y+8}" text-anchor="middle" font-size="10.5" font-weight="700" fill="#212529">VAF ${v!=null?v:'?'}</text>`;
  s+=`<text x="${x}" y="${y+19}" text-anchor="middle" font-size="7.5" fill="#aaa">${n.split(':')[0]}</text>`;});
 return '<div style="overflow-x:auto;max-width:100%;padding-bottom:4px">'+s+'</svg></div>';
}
// #5 截斷區全位點分群樹(位點層 pos_nested,突破 8 向量上限;只已定相位點) + #4 成環候選固定樹(移除最小衝突邊)
function allPosTree(r){
 let nh=r.node_hp||{},ns=Object.keys(nh);if(ns.length<2)return '';
 let nsSet=new Set(ns),ned=(r.pos_nested||[]).filter(e=>nsSet.has(e[0])&&nsSet.has(e[1]));
 let gk=Object.keys(r.populations||{}).find(k=>k&&/^R+$/.test(k)),germN=gk?r.populations[gk]:0;
 let body,note='';
 if(ns.length>40){body=posBarStrip(nsSet,r.pos_vaf||{});note='位點>40→改左右滑桿 VAF bar(基因組順序)';}
 else if(hasCycleEdges(ned)){body=posTable(nsSet,r.pos_vaf||{},'全',ned);note='位點層亦成環→位點+VAF 表(可能樹見下方 carousel)';}
 else{let hasp=new Set(ned.map(e=>e[1])),edges=ned.map(e=>[e[0],e[1]]);[...nsSet].forEach(n=>{if(!hasp.has(n))edges.unshift(['ROOT',n])});body=posSVG(edges,nsSet,r.pos_vaf||{},null,nh,germN)+'<div class="note">↔ 樹過寬可左右捲動</div>';}
 let dropHP=Object.values(nh).filter(h=>h!='H1'&&h!='H2').length;
 return `<details open style="background:#fff;border:1px solid #ffd8a8;border-radius:6px;padding:8px;margin:6px 0"><summary style="cursor:pointer;font-weight:700;color:#e8590c">🧬 全位點分群樹（截斷區·突破前8上限·${ns.length} 已定相位點${note?'·'+note:''}）</summary><div class="note" style="margin:4px 0">genotype 向量截斷在 8 位點;此樹改用<b>位點層 pos_nested 巢狀</b>畫出全 <b>${ns.length}</b> 個<b>已定相</b> sSNV(色=HP·藍H1/紫H2/灰H3?)。共 ${r.n_sSNV} sSNV,未定相 <b>${r.n_sSNV-ns.length}</b> 個無座標不納入(caveat)。</div>${body}${dropHP?'<div class="note" style="color:#a37200">⚠ 含 '+dropHP+' 個未定相(H3?)位點(灰)</div>':''}</details>`;
}
function cycleCandidates(r){
 let nh=r.node_hp||{},nsSet=new Set(Object.keys(nh));
 let ned=(r.pos_nested||[]).filter(e=>nsSet.has(e[0])&&nsSet.has(e[1]));
 if(ned.length<2||!hasCycleEdges(ned))return null;
 let cands=[],seen={};
 function addCand(sub,label){let hasp=new Set(sub.map(e=>e[1])),edges=sub.map(e=>[e[0],e[1]]);[...nsSet].forEach(n=>{if(!hasp.has(n))edges.unshift(['ROOT',n])});
  let key=JSON.stringify(sub.slice().map(e=>e.join('>')).sort());if(!seen[key]&&cands.length<6){seen[key]=1;cands.push({edges:edges,label:label});}}
 // 1) 移除單一最小衝突邊 → 乾淨候選
 for(let i=0;i<ned.length&&cands.length<6;i++){let sub=ned.filter((e,j)=>j!=i);if(!hasCycleEdges(sub))addCand(sub,'移除 '+ned[i][0].split(':')[1]+'↔'+ned[i][1].split(':')[1]);}
 // 2) fallback(多環區):貪婪生成樹(輪轉起始序)→ 多環也能給固定樹示意
 if(!cands.length){for(let s=0;s<Math.min(ned.length,4);s++){let order=ned.slice(s).concat(ned.slice(0,s)),keep=[];order.forEach(e=>{if(!hasCycleEdges(keep.concat([e])))keep.push(e)});if(keep.length)addCand(keep,'貪婪生成樹 #'+(s+1)+'(多環→需移除多邊)');}}
 return cands.length?cands:null;
}
function cycleCard(cs,idx,r){let c=cs[idx],nh=r.node_hp||{};let th=posSVG(c.edges,new Set(Object.keys(nh)),r.pos_vaf||{},null,nh);
 return `<div style="display:flex;align-items:flex-start;gap:8px"><button onclick="cycNav(-1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px;margin-top:40px">◀</button><div style="flex:1;min-width:0">${th}<div class="note" style="margin-top:4px">候選固定樹 <b>${idx+1} / ${cs.length}</b>（${c.label}）</div></div><button onclick="cycNav(1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px;margin-top:40px">▶</button></div>`}
window.cycNav=function(d){let s=window.__cyc;if(!s)return;s.idx=(s.idx+d+s.cs.length)%s.cs.length;let box=document.getElementById('cycbox');if(box)box.innerHTML=cycleCard(s.cs,s.idx,s.r)};
// #7 甲基裁決(cis-control)收斂成單一函式(原散在 undefined/four-gamete 兩處)
function methylVerdict(r){let cross=(r.n_roots>=2)||((r.haplotypes||'').indexOf('H1')>=0&&(r.haplotypes||'').indexOf('H2')>=0);
 return `<div style="margin:6px 0;padding:7px 10px;border-radius:6px;font-size:11.5px;background:${cross?'#e7f5ff':'#fff5f5'};border:1px solid ${cross?'#74c0fc':'#ffc9c9'}">🧬 <b>甲基能否裁決此區的樹/順序：</b>${cross?'<b style="color:#1971c2">CROSS-HP(跨 H1/H2)</b> → 衝突/分支屬 <b>allelic</b>(兩突變在不同染色體)、本就無 subclone 順序可排;甲基只是 germline-ASM、<b>不裁決</b>。':'<b style="color:#c92a2a">SAME-HP(單一 germline HP)</b> → 甲基隨 genotype cis 共變(cis-ASM double-dip)、normal 無對應 within-HP 軸→<b>結構性無法裁決</b>。'} 🔴 <b>06-28 cis-control 定案</b>:甲基排序 ρ≈0.18 p≈0.06-0.08 <b>不顯著</b> → <b>不能裁決哪棵樹/哪個順序更可能</b>(bounded-auxiliary);只能靠讀數/VAF 弱先驗 + single-cell/multi-region。</div>`;
}
// Level2: 逐區 locus lollipop — 線=該區基因組範圍;每 sSNV 一根 lollipop(位置·HP色·VAF高);底色=結果
function locusTrack(r){
 var nh=r.node_hp||{},vaf=r.pos_vaf||{};
 // 修 bug:改用 pos_vaf ∪ node_hp 的位點(未定相區 node_hp 空但 pos_vaf 有位點→原本整條不顯示);未定相 hp 標 H?(灰)
 var allk={};Object.keys(vaf).forEach(function(k){allk[k]=1});Object.keys(nh).forEach(function(k){allk[k]=1});
 var keys=Object.keys(allk);
 if(!keys.length)return '';
 var pos=keys.map(function(k){return {k:k,p:parseInt(k.split(':')[1])||0,hp:nh[k]||'H?',v:vaf[k]};}).filter(function(o){return o.p>0}).sort(function(a,b){return a.p-b.p});
 if(!pos.length)return '';
 var rs=(r.start!=null?r.start:pos[0].p),re=(r.start!=null?(r.start+(r.span||0)):pos[pos.length-1].p);
 var lo=Math.min(rs,pos[0].p),hi=Math.max(re,pos[pos.length-1].p),rng=Math.max(1,hi-lo),pad=rng*0.04;
 var W=760,AX=46,X0=AX,XW=W-AX-16,X=function(p){return X0+XW*(p-(lo-pad))/(rng+2*pad)};
 var hpc={H1:'#1c7ed6',H2:'#7048e8'};
 var sv=regSit(r),bg=sv?sv[1]:'#adb5bd',sname=sv?sv[0]:'單群/其他';
 var maxv=Math.max.apply(null,pos.map(function(o){return o.v||0}).concat([0.01]));
 var BASE=94,HMAX=60;
 var s='<svg viewBox="0 0 '+W+' 150" width="100%" style="font-family:ui-monospace,monospace;max-width:'+W+'px">';
 // #7 Y 軸(VAF 高度) + 中線
 s+='<text x="'+(AX-8)+'" y="'+(BASE-HMAX-8)+'" text-anchor="middle" font-size="9" font-weight="700" fill="#495057">VAF</text>';
 s+='<line x1="'+(AX-8)+'" y1="'+(BASE-HMAX)+'" x2="'+(AX-8)+'" y2="'+BASE+'" stroke="#adb5bd"/>';
 s+='<text x="'+(AX-11)+'" y="'+(BASE+3)+'" text-anchor="end" font-size="8" fill="#868e96">0</text>';
 s+='<text x="'+(AX-11)+'" y="'+(BASE-HMAX/2+3)+'" text-anchor="end" font-size="8" fill="#868e96">'+(maxv/2).toFixed(2)+'</text>';
 s+='<text x="'+(AX-11)+'" y="'+(BASE-HMAX+6)+'" text-anchor="end" font-size="8" fill="#868e96">'+maxv.toFixed(2)+'</text>';
 s+='<line x1="'+(AX-6)+'" y1="'+(BASE-HMAX/2)+'" x2="'+X(hi)+'" y2="'+(BASE-HMAX/2)+'" stroke="#e9ecef" stroke-dasharray="2 3"/>';
 s+='<rect x="'+X(lo)+'" y="'+(BASE-2)+'" width="'+(X(hi)-X(lo))+'" height="5" rx="2" fill="'+bg+'" opacity="0.30"/>';
 s+='<line x1="'+X(lo)+'" y1="'+BASE+'" x2="'+X(hi)+'" y2="'+BASE+'" stroke="'+bg+'" stroke-width="2"/>';
 s+='<text x="'+X(rs)+'" y="'+(BASE+15)+'" font-size="9" fill="#495057">'+rs.toLocaleString()+'</text>';
 s+='<text x="'+X(re)+'" y="'+(BASE+15)+'" font-size="9" fill="#495057" text-anchor="end">'+re.toLocaleString()+'</text>';
 pos.forEach(function(o){var x=X(o.p),h=Math.max(5,HMAX*(o.v||0)/maxv),col=hpc[o.hp]||'#adb5bd';
   s+='<line x1="'+x+'" y1="'+BASE+'" x2="'+x+'" y2="'+(BASE-h)+'" stroke="'+col+'" stroke-width="1.4"/>';
   s+='<circle cx="'+x+'" cy="'+(BASE-h)+'" r="4" fill="'+col+'"><title>'+o.k+' · '+(o.hp||'HP?')+' · VAF '+(o.v!=null?o.v:'?')+'</title></circle>';});
 // #7 距離比例尺(bp/kb)
 var niceLen=Math.pow(10,Math.max(0,Math.floor(Math.log10(rng/3))));var sbPx=XW*niceLen/(rng+2*pad);var sbLab=niceLen>=1000?(niceLen/1000)+' kb':niceLen+' bp';
 s+='<line x1="'+AX+'" y1="'+(BASE+30)+'" x2="'+(AX+sbPx).toFixed(1)+'" y2="'+(BASE+30)+'" stroke="#495057" stroke-width="2"/><line x1="'+AX+'" y1="'+(BASE+26)+'" x2="'+AX+'" y2="'+(BASE+34)+'" stroke="#495057"/><line x1="'+(AX+sbPx).toFixed(1)+'" y1="'+(BASE+26)+'" x2="'+(AX+sbPx).toFixed(1)+'" y2="'+(BASE+34)+'" stroke="#495057"/><text x="'+(AX+sbPx+7).toFixed(1)+'" y="'+(BASE+34)+'" font-size="9" fill="#495057">↔ '+sbLab+' 比例尺</text>';
 s+='</svg>';
 var dropHP=pos.filter(function(o){return o.hp!='H1'&&o.hp!='H2'}).length;
 return '<div style="background:#fff;border:1px solid #dee2e6;border-radius:6px;padding:8px 10px;margin:6px 0"><b>🧬 locus 突變排列（每 sSNV 一根 lollipop·桿高=VAF·色=HP）</b> <span class="note">底色帶=該區結果「<b style="color:'+bg+'">'+sname+'</b>」·範圍 '+rs.toLocaleString()+'–'+re.toLocaleString()+'（'+(((re-rs)/1000).toFixed(1))+'kb）·'+pos.length+' sSNV'+(dropHP?'（其中 '+dropHP+' 未定相·灰）':'')+(r.truncated?'（共 '+r.n_sSNV+'，截斷）':'')+'</span> <span style="margin-left:6px"><span style="color:#1c7ed6">●H1</span> <span style="color:#7048e8">●H2</span>'+(dropHP?' <span style="color:#adb5bd">●HP?</span>':'')+'</span>'+s+'<div class="note">對照右側克隆樹：同一批 sSNV 在基因組上的實際位置與 HP 分群（跨 HP＝兩條染色體各自突變＝allelic）。</div></div>';
}
// filters + sort
let det=D.detail;
let chrsel=el('f_chr');chrsel.innerHTML='<option value="">全</option>';D.chroms.forEach(c=>{let o=document.createElement('option');o.value=c;o.textContent=c;chrsel.appendChild(o)});
const FACET_CFG={topology_type:{order:['linear(全直系)','branched(直系+姊妹)','star(全姊妹)','single','germline_only'],label:{'linear(全直系)':'linear 全直系','branched(直系+姊妹)':'多根 lineage','star(全姊妹)':'star 全姊妹','single':'single 單群','germline_only':'germline only'},cls:{'linear(全直系)':'t_linear','branched(直系+姊妹)':'t_branched','star(全姊妹)':'t_star','single':'t_single','germline_only':'t_single'}},determinacy:{order:['A_determined(單分子向量)','A_ambiguous_order(缺中間群)','B_pairwise_structure','C_underdetermined','incompatible','other'],label:{'A_determined(單分子向量)':'A 唯一·單分子','A_ambiguous_order(缺中間群)':'A 順序未定','B_pairwise_structure':'B 拼接','C_underdetermined':'C 單ALT群欠定','incompatible':'✗ 成環','other':'— 單群無分支'},cls:{'A_determined(單分子向量)':'det_A','A_ambiguous_order(缺中間群)':'det_amb','B_pairwise_structure':'det_B','C_underdetermined':'det_C','incompatible':'det_incompat','other':'det_other'}},genome_ctx:{order:['arm','telomere','centromere'],label:{'arm':'arm 臂','telomere':'telomere 端粒','centromere':'centromere 著絲點'},cls:{'arm':'ctx_arm','telomere':'ctx_telomere','centromere':'ctx_centromere'}}};
['topology_type','determinacy','genome_ctx'].forEach(k=>{let c=el('cb_'+k);c.innerHTML='';let cfg=FACET_CFG[k],cnt={};det.forEach(r=>{cnt[r[k]]=(cnt[r[k]]||0)+1});let vals=cfg.order.filter(v=>cnt[v]!=null).concat(Object.keys(cnt).filter(v=>!cfg.order.includes(v)).sort());vals.forEach(o=>{let lab=document.createElement('label');lab.className='chip '+(cfg.cls[o]||'det_other');lab.title=o;lab.innerHTML='<input type="checkbox" value="'+o+'">'+(cfg.label[o]||o)+'<span class="cnt">'+cnt[o]+'</span>';c.appendChild(lab)});c.onchange=()=>{c.querySelectorAll('.chip').forEach(ch=>ch.classList.toggle('on',ch.querySelector('input').checked));render()}});
function cset(k){return new Set([...el('cb_'+k).querySelectorAll('input:checked')].map(x=>x.value))}
const SORT={coord:(a,b)=>a.chrom.localeCompare(b.chrom,undefined,{numeric:true})||a.start-b.start,nsnv:(a,b)=>b.n_sSNV-a.n_sSNV,nclust:(a,b)=>b.n_clusters-a.n_clusters||b.n_sSNV-a.n_sSNV,glen:(a,b)=>((Object.keys(b.populations||{})[0]||'').length)-((Object.keys(a.populations||{})[0]||'').length)||b.n_sSNV-a.n_sSNV,region:(a,b)=>a.region.localeCompare(b.region)};
function render(){
 let ch=el('f_chr').value,mc=+el('f_minc').value,tf=el('f_tpfp').value,loh=el('f_loh').checked,undef=el('f_undef').checked,gene=el('f_gene')?el('f_gene').checked:false,q=el('f_q').value.trim(),so=el('f_sort').value;
 let tt=cset('topology_type'),dd=cset('determinacy'),gc=cset('genome_ctx');
 let tpfpok=r=>(tf=='all')||(tf=='tp'&&r.tp>0)||(tf=='fp'&&r.fp>0)||(tf=='both'&&r.tp>0&&r.fp>0);
 let f=det.filter(r=>(!ch||r.chrom==ch)&&(!tt.size||tt.has(r.topology_type))&&(!dd.size||dd.has(r.determinacy))&&(!gc.size||gc.has(r.genome_ctx))&&r.n_clusters>=mc&&tpfpok(r)&&(!loh||regLOH(r))&&(!undef||r.undefined)&&(!gene||geneHit(r.region))&&(!q||r.region.includes(q)));
 f.sort(SORT[so]||SORT.coord);if(el('f_sortdir').value=='rev')f.reverse();
 el('cnt').textContent=f.length+' 區';
 el('list').innerHTML=f.slice(0,700).map(r=>`<div class="row" data-i="${det.indexOf(r)}"><b>${r.region}</b> <span class="tag ${TT[r.topology_type]||'t_single'}">${TTLABEL[r.topology_type]||r.topology_type.split('(')[0]}</span><span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span>${regLOH(r)?'<span class="tag" style="background:#9775fa;color:#fff">LOH</span>':''}${geneHit(r.region)?'<span class="tag" style="background:#e64980;color:#fff">🧬癌</span>':''}<br><span class="note">${r.n_sSNV}sSNV·c=${r.n_clusters}·${r.haplotypes}·${r.cn}·TP${r.tp}/FP${r.fp}${r.ambig_nodes>0?'·⚠序未定':''}</span></div>`).join('')+(f.length>700?`<div class="note" style="padding:8px">...前 700（共 ${f.length}）</div>`:'');
 el('list').querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x));
}
// #2 區域內 read 組合矩陣(read 群×sSNV 位點)+ 2^k 可能組合 vs 觀測
function readMatrix(r){
 var pops=r.populations||{},np=r.node_paths||{},glen=(Object.keys(pops)[0]||'').length;if(!glen)return '';
 var vecs=Object.entries(pops).sort((a,b)=>b[1]-a[1]);var tot=vecs.reduce((s,v)=>s+v[1],0)||1;var possible=Math.pow(2,glen);
 var cols='';for(var p=0;p<glen;p++)cols+='<th style="font-size:9px">S'+(p+1)+'</th>';
 var rows=vecs.map(function(kv){var vec=kv[0],n=kv[1],cells='';for(var p=0;p<glen;p++){var ch=vec[p],A=ch=='A';cells+='<td style="background:'+(A?'#f59f00':(ch=='R'?'#e9ecef':'#fff'))+';color:'+(A?'#fff':'#868e96')+';text-align:center;font-weight:'+(A?'700':'400')+'">'+ch+'</td>';}var lab=np[vec]||(vec.indexOf('A')>=0?'—(未定)':'germline');return '<tr><td class="mono" style="color:#1971c2">'+lab+'</td>'+cells+'<td>'+n+'</td><td>'+(100*n/tot).toFixed(0)+'%</td></tr>';}).join('');
 return '<details style="margin:6px 0" open><summary style="cursor:pointer;font-weight:600;font-size:12px">🧬 區域內 read 組合矩陣（read 群 × sSNV 位點；A=橙 R=灰）</summary><div style="font-size:11px;margin:3px 0">可能組合 <b>2^'+glen+'='+possible+'</b> ｜ 觀測到 <b>'+vecs.length+'</b> 種 ｜ 覆蓋 <b>'+(100*vecs.length/possible).toFixed(1)+'%</b>'+(r.truncated?'（截斷:僅前 '+glen+'/'+r.n_sSNV+' 位點）':'')+'</div><table style="font-size:10.5px"><tr><th>lineage</th>'+cols+'<th>reads</th><th>%</th></tr>'+rows+'</table><div class="note" style="margin-top:3px">每列=一種觀測到的 read 基因型組合(=克隆樹一節點,lineage 對應樹);未觀測的 2^k−N 種=無 read 支持。</div></details>';}
// #3 每兩 sSNV 位點 pairwise 關係矩陣(從 populations 推 2×2)
function pw2x2(pops,i,j){var g={RR:0,RA:0,AR:0,AA:0};Object.keys(pops).forEach(function(v){if(v.length<=Math.max(i,j))return;var a=v[i],b=v[j];if((a!='A'&&a!='R')||(b!='A'&&b!='R'))return;g[(a=='A'?'A':'R')+(b=='A'?'A':'R')]+=pops[v];});return g;}
function pwClass(g){var tot=g.RR+g.RA+g.AR+g.AA,eps=Math.max(1,tot*0.02);var ra=g.RA>eps,ar=g.AR>eps,aa=g.AA>eps;if(ra&&ar&&aa)return 'conflict';if(ra&&ar&&!aa)return 'sibling';if((!ra||!ar)&&aa)return 'nested';if(!ra&&!ar)return 'co_linked';return 'partial';}
function pairwiseMatrix(r){var pops=r.populations||{},glen=(Object.keys(pops)[0]||'').length;if(glen<2)return '';
 var CLS={conflict:['衝突(四配子)','#f03e3e'],sibling:['互斥/平行','#f59f00'],nested:['巢狀(祖先→後代)','#1c7ed6'],co_linked:['完全連鎖','#37b24d'],partial:['弱/未定','#adb5bd']};
 var head='<tr><th></th>';for(var j=1;j<glen;j++)head+='<th style="font-size:9px">S'+(j+1)+'</th>';head+='</tr>';var body='';
 for(var i=0;i<glen-1;i++){var tr='<tr><td style="font-size:9px;font-weight:600">S'+(i+1)+'</td>';for(var j=1;j<glen;j++){if(j<=i){tr+='<td></td>';continue;}var g=pw2x2(pops,i,j),C=CLS[pwClass(g)];tr+='<td title="S'+(i+1)+'-S'+(j+1)+' '+C[0]+' ｜ RR'+g.RR+' RA'+g.RA+' AR'+g.AR+' AA'+g.AA+'" style="background:'+C[1]+';width:15px;height:15px;cursor:help;border:1px solid #fff"></td>';}tr+='</tr>';body+=tr;}
 var leg=Object.keys(CLS).map(k=>'<span style="color:'+CLS[k][1]+'">■'+CLS[k][0]+'</span>').join(' ');
 return '<details style="margin:6px 0"><summary style="cursor:pointer;font-weight:600;font-size:12px">🔗 每兩 sSNV 位點 pairwise 關係矩陣（hover 看 2×2）</summary><div class="note" style="margin:3px 0">'+leg+'</div><table style="border-collapse:collapse">'+head+body+'</table><div class="note" style="margin-top:3px">每格=一對位點關係(從 2×2 RR/RA/AR/AA 推;ε=2% 噪聲底);衝突=四配子違反(不能單一樹)、巢狀=一單突變缺、互斥=兩單突變共存從不同時、完全連鎖=只 RR+AA。</div></details>';}
// #4 前後相鄰區關聯(🔴 read-span 不跨區連結→描述性,LOH/CN 可跨區、lineage 不可確認)
function neighborLink(r){var m=r.region.match(/(chr\w+):(\d+)-(\d+)/);if(!m)return '';var chrom=m[1],st=+m[2],en=+m[3];
 var sib=det.filter(x=>x.chrom==chrom).map(function(x){var mm=x.region.match(/:(\d+)-(\d+)/);return {x:x,s:mm?+mm[1]:0,e:mm?+mm[2]:0};}).sort((a,b)=>a.s-b.s);
 var idx=sib.findIndex(o=>o.x.region==r.region);if(idx<0)return '';
 var TTc=t=>(t||'').split('(')[0];
 function cmp(o,gap,side){if(!o)return '<div style="flex:1"><b>'+side+'</b> <span class="note">(同染色體無相鄰區)</span></div>';var x=o.x;var same=[];if(!!regLOH(x)==!!regLOH(r)&&regLOH(r))same.push('🟣LOH');if(x.cn==r.cn&&r.cn!='unknown')same.push('CN='+x.cn);if(TTc(x.topology_type)==TTc(r.topology_type))same.push('拓樸'+TTc(x.topology_type));if(x.haplotypes==r.haplotypes)same.push('HP='+x.haplotypes);var gk=gap>=1e6?(gap/1e6).toFixed(2)+'Mb':(gap/1000).toFixed(0)+'kb';return '<div style="flex:1;min-width:200px"><b>'+side+'</b> <span class="mono" style="cursor:pointer;color:#1971c2" onclick="show('+det.indexOf(x)+')">'+x.region+'</span> <span class="note">('+x.n_sSNV+'sSNV·c='+x.n_clusters+'·'+x.cn+'·gap '+gk+')</span><br>'+(same.length?'<span style="color:#2b8a3e">共享: '+same.join(' / ')+'</span>':'<span class="note">無共享屬性</span>')+'</div>';}
 var prev=idx>0?sib[idx-1]:null,next=idx<sib.length-1?sib[idx+1]:null;
 return '<div style="background:#f4f6ff;border:1px solid #c3d0f5;border-radius:6px;padding:8px 10px;margin:6px 0;font-size:11.5px"><b>↔ 前後相鄰區關聯</b> <span class="note">(依基因組座標;🔴 read-span 不跨區連結 → LOH/CN 段可跨區共享、但 <b>lineage 不可確認連續</b>=描述性非重建)</span><div style="display:flex;gap:14px;margin-top:5px;flex-wrap:wrap">'+cmp(prev,prev?st-prev.e:0,'◀ 前區')+cmp(next,next?next.s-en:0,'後區 ▶')+'</div></div>';}
function show(i,row){el('list').querySelectorAll('.row').forEach(x=>x.classList.remove('sel'));if(row)row.classList.add('sel');let r=det[i];
 let popcount=r.populations;
 let glen=(Object.keys(r.populations||{})[0]||'').length;
 let rdr=(D.rd||{})[r.region]||{};
 let np=r.node_paths||{};
 let cf=fourGamete(r.populations);
 let cyc=cycleCandidates(r);
 let ctr=(D.candtrees||{})[r.region];
 let pt=Object.entries(r.populations).sort((a,b)=>b[1]-a[1]).map(([g,c])=>{let tot=Object.values(r.populations).reduce((a,b)=>a+b,0);return `<tr><td class="mono" style="color:#1971c2;font-weight:600">${np[g]||(g.includes('A')?'—(未定)':'germline')}</td><td class="mono">${g}</td><td>${sLabels(g)}</td><td>${c}</td><td>${(100*c/tot).toFixed(0)}%</td></tr>`}).join('');
 el('detail').innerHTML=`<h3>${r.region} <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type}</span> <span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span></h3>
  <div class="kv"><div class="b">${r.truncated?(glen+'/'+r.n_sSNV+' sSNV(截斷)'):(r.n_sSNV+' sSNV')}</div><div class="b">span ${r.span>=1e6?(r.span/1e6).toFixed(2)+'Mb':(r.span/1000).toFixed(1)+'kb'}</div><div class="b">c=${r.n_clusters} 群</div><div class="b">HP: ${r.haplotypes}</div><div class="b" title="${r.cn=='unknown'?'此樣本未併入 CN census(=未標註,非 CN 不明);6 樣本多為 unknown':'copy-number 狀態:neutral/gain/loss/loh'}">CN: ${r.cn}${r.cn=='unknown'?'(未標註)':''}</div>${regLOH(r)?`<div class="b" style="background:#efe6ff;color:#6741d9" title="longphase-TO LOH.bed(SEQC2-validated)overlap 或 cn=loh">🟣 LOH 區</div>`:''}<div class="b">TP ${r.tp} / FP ${r.fp}</div><div class="b">${r.determinacy}</div>${r.drop_noise_frac>0?`<div class="b">噪聲過濾 ${(r.drop_noise_frac*100).toFixed(0)}%</div>`:''}${r.ambig_nodes>0?`<div class="b" style="background:#fff3bf">⚠ 順序未定 ${r.ambig_nodes}(缺中間群)</div>`:''}${r.truncated?`<div class="b" style="background:#ffe3e3;color:#c92a2a" title="genotype 向量截斷在 8 位點(上游 GCAP=8);此區 ambig/四配子/機率偵測不完整,成環可能為截斷假象">⚠ 截斷 n_sSNV>8(偵測不完整)</div>`:''}</div>
  ${neighborLink(r)}
  ${locusTrack(r)}
  ${(r.determinacy=='incompatible'||r.has_cycle)?`<div style="background:#fff5f5;border:1px solid #ffc9c9;border-radius:6px;padding:7px 10px;margin:6px 0;font-size:11.5px">🔴 <b>成環/衝突成因</b>：${r.cycle_cause||'—'}${r.truncated?' ＋ 截斷(>8 cap)加劇':''}。${(r.cycle_cause||'').indexOf('CN-gain')>=0?'<b>CN-gain multiplicity 假環</b>（同突變多拷貝→pairwise 矛盾）＝<b>非真演化衝突</b>；建議補 CN/mappability mask、不強建樹。':'other-pairwise-cycle：可能真衝突，但須先排除多拷貝/截斷/mapping artifact 才算真四配子衝突。'}<br><span class="note">📊 全樣本 12 成環區：9 為 CN-gain multiplicity 假環 + 3 other；8/12 截斷 → <b>多數非真演化衝突</b>。</span></div>`:''}
  ${r.undefined?`<div style="background:#fff8e6;border:1px solid #ffe0a3;border-radius:6px;padding:7px 10px;margin:6px 0;font-size:11.5px"><b>⚠ 此區分支順序未定/不相容</b> → 下方候選樹為可能位置(等機率·非給答案);甲基裁決見下方統一說明。</div>`:''}
  ${cf.length?`<div style="background:#fff0f6;border:1px solid #f783ac;border-radius:6px;padding:9px;margin:6px 0"><b>⚠ 四配子違反（incompatible）→ 無法成單一樹</b>　錨點 <b>RR=germline</b>（normal 確認 REF）→ <b>AA=雙突變（最遠）</b>；RA／AR 兩單突變並存＝累積順序未定。<table style="margin-top:5px"><tr><th>衝突對</th><th>RR<br>germ根</th><th>RA<br>僅後者</th><th>AR<br>僅前者</th><th>AA<br>最遠</th><th>讀數弱提示</th></tr>${cf.map(c=>`<tr><td class="mono"><b>${c.pair}</b></td><td>${c.g.RR}</td><td>${c.g.RA}</td><td>${c.g.AR}</td><td>${c.g.AA}</td><td class="note">${c.g.AR>c.g.RA?'S'+c.i+' 單突變較多':c.g.RA>c.g.AR?'S'+c.j+' 單突變較多':'兩單突變相當'}（弱·非定論）</td></tr>`).join('')}</table><div class="note" style="margin-top:4px">真實關係非單一樹;可能的固定樹見下方 carousel。</div></div>`:''}
  ${r.n_roots>=2?`<div style="background:#fff4e6;border:1px solid #ffd8a8;border-radius:6px;padding:8px"><b>⚠ 此區跨 H1/H2（${r.n_roots} 棵樹）→ 兩棵分開的 HP 克隆樹（正確）：</b><div class="note" style="margin:3px 0">為何兩棵:somatic 突變散在 H1 與 H2 兩條 germline 單倍型=各自染色體突變(<b>allelic</b>,非 subclone)。${Object.keys(r.populations||{}).filter(g=>g.indexOf('A')>=0).length<=1?`此區 genotype 向量樹 trivial(僅單一 ALT 向量${r.truncated?'+截斷':''})→改以<b>位點層</b>兩棵 HP 樹呈現。`:''}</div>${posTree(r)}</div><details style="margin-top:6px"><summary style="cursor:pointer;color:#868e96;font-size:11.5px">▶ 混合 genotype-向量樹（跨 HP 混合,僅參考）</summary>${tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes)}</details>`:`<b>克隆樹（germline→…；節點=lineage標籤·S-mut-set·reads·%；座標=向量）</b>${tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes)}`}
  ${readMatrix(r)}
  ${pairwiseMatrix(r)}
  ${r.truncated?allPosTree(r):''}
  ${methylVerdict(r)}
  ${(function(){var q=(D.scoring&&D.scoring.queue||[]).find(function(x){return x.region===r.region});return (q&&q.resolution_path)?`<div class="note" style="background:#f1f3f5;border-radius:5px;padding:6px 9px;margin:6px 0">🧭 <b>解法路徑(resolution_path)</b>：${q.resolution_path}${(q.resolution_path.indexOf('VAF')>=0)?`　<span style="color:#a37200">↳ 排序後備鏈:①先用 <b>VAF/CCF(CN-clean)</b> 定先後 → ②VAF tie 才用<b>甲基 L3 後備</b>。但 06-28 cis-control 已測甲基 ρ0.18 不顯著 → 此後備<b>實際乾淨可用≈0</b>(bounded-auxiliary)。</span>`:''}</div>`:'';})()}
  ${cyc?`<div style="background:#f8f0fc;border:1px solid #d0bfff;border-radius:6px;padding:9px;margin:8px 0"><b>🔀 可能的固定樹（左右選·成環區）</b> <span class="note">移除最小衝突邊→列出 <b>${cyc.length}</b> 棵可能的固定樹。<b style="color:#c2255c">⚠ 假環(CN-gain multiplicity)·僅示意·非真演化衝突·🔴甲基無法裁決哪棵更可能</b>。</span><div id="cycbox" style="margin-top:6px">${cycleCard(cyc,0,r)}</div></div>`:''}
  ${(ctr&&ctr.candidate_set&&ctr.candidate_set.length)?(function(){var single=ctr.true_count==1&&!ctr.equiprobable;var bg=single?'#eef9f0':'#f8f0fc',bd=single?'#b2dfc0':'#d0bfff';var badge=single?`<b style="color:#2b8a3e">✅ 唯一確認整樹</b> <span class="note">(true_count=1·resolution=${ctr.resolution||'—'}·中間群被觀測/共event→拓樸唯一,非多候選)</span>`:`<b style="color:#7048e8">⚖ ${ctr.true_count} 種等機率候選（左右滑動）</b> <span class="note">(缺中間群→插虛擬中間節點·${ctr.honest_note}·<b>非給答案</b>·🔴甲基無法裁決)</span>`;return `<div style="background:${bg};border:1px solid ${bd};border-radius:6px;padding:9px;margin:8px 0">${badge}<div id="cttreebox" style="margin-top:6px">${candTreeCard(ctr,0,r)}</div></div>`;})():''}
  <div class="note" style="background:#eef9f0;border:1px solid #b2dfc0;border-radius:5px;padding:6px 9px;margin:6px 0">🔬 <b>read-driven 交叉確認</b>（22-way 平行遍歷 read）：多-ALT read <b>${rdr.rd_multi_alt!=null?rdr.rd_multi_alt:'—'}</b>／distinct combos <b>${rdr.rd_combos!=null?rdr.rd_combos:'—'}</b>（pipeline n_clusters=${r.n_clusters}）／max chain <b>${rdr.rd_max_chain!=null?rdr.rd_max_chain:'—'}</b>${r.truncated?'　🔴 截斷區→read-driven 原始串接（非假樹）':''}</div>
  <div class="note">S1..S${r.truncated?glen:r.n_sSNV}=區內排序 sSNV${r.truncated?`（此區共 ${r.n_sSNV} sSNV,僅前 ${glen} 進向量;樹/標籤只到 S${glen}）`:''}；直系=往下、姊妹=同層分叉；germline 根標 reads·%。tree_shape(pairwise)=${r.tree_shape}。genome_ctx 為近似(±3Mb)。</div>
  <b>細胞群(lineage 標籤 → 向量 → S 突變 → reads → 佔比)</b><table><tr><th>lineage</th><th>向量</th><th>突變(S)</th><th>reads</th><th>佔比</th></tr>${pt}</table>${resolveBlock(r.region)}${geneBlock(r.region)}`;
 window.__cyc=cyc?{cs:cyc,idx:0,r:r}:null;
 window.__ctr=(ctr&&ctr.candidate_set&&ctr.candidate_set.length)?{cs:ctr,idx:0,r:r}:null;
}
function resolveBlock(region){
  var R=(window.__RESOLUTION__||{})[region]; if(!R) return '';
  var kcol={BLOCK:'#f08c00',ORDERED:'#2f9e44',CONFLICT:'#e03131',PARTIAL_ORDER:'#1971c2',AMBIG_NOCOREAD:'#868e96'};
  var klab={BLOCK:'BLOCK 共event（多突變一起·不排序）',ORDERED:'ORDERED 可定序（中間群被觀測）',CONFLICT:'CONFLICT 真衝突（4-gamete）',PARTIAL_ORDER:'PARTIAL 部分定序',AMBIG_NOCOREAD:'真等機率（無 coread）'};
  var pv=R.provenance||{},k=R.n_sSNV;
  var s='<div style="background:#f0f7ff;border:1px solid #a5c8f0;border-radius:6px;padding:9px 12px;margin:8px 0;font-size:12px">';
  s+='<b>🧩 結構解析（gained-pair pairwise 定序 + 區域內所有子read）</b>';
  s+='<div class="note" style="margin:3px 0">provenance：<b>絕對群</b>（單分子跨全部 '+k+' 點）<b>'+(pv.n_absolute_pops||0)+'</b> 群／'+(pv.absolute_reads||0)+' read（絕對比 <b>'+(((pv.abs_frac||0)*100).toFixed(0))+'%</b>）；其餘 pairwise 組合推得</div>';
  var G8=(window.__GT8__||{})[region];
  if(G8){var vc=G8.verdict==='pairwise-tree'?'#2f9e44':'#e03131';
    s+='<div style="margin:5px 0;padding:6px 9px;border-left:4px solid '+vc+';background:#fff">';
    s+='<b style="color:'+vc+'">🌲 >8 全 pairwise 樹（去 8-cap）</b>　'+G8.n_positions+' sSNV → <b>'+G8.n_nodes+'</b> 節點／'+G8.n_tree_edges+' 邊／深 '+G8.max_depth+'（8-cap 只見 '+(G8.vector_cap_nodes!=null?G8.vector_cap_nodes:'—')+' 群）';
    s+='<span class="note"> '+G8.verdict+'（衝突對比 '+G8.conflict_frac+(G8.has_cycle?'·有環':'')+'）·'+G8.provenance+'</span>';
    if(G8.cn_caveat) s+='<div class="note" style="color:#c92a2a;margin-top:3px">🔴 '+G8.cn_caveat+'</div>';
    s+='</div>';}
  (R.edge_resolution||[]).forEach(function(e){
    var oc=(e.order||[]).map(function(o){return o[0]+'→'+o[1];}).join(', ');
    s+='<div style="margin:4px 0;padding:5px 8px;border-left:4px solid '+(kcol[e.klass]||'#ccc')+';background:#fff">';
    s+='<b style="color:'+(kcol[e.klass]||'#333')+'">'+(klab[e.klass]||e.klass)+'</b>　<span class="mono">'+e.parent+'→'+e.child+'</span>';
    s+='<span class="note"> gained '+JSON.stringify(e.gained_pos)+(oc?'　順序 '+oc:'')+(e.tentative?'　⚠tentative(<2%,需放寬層)':'')+(e.vaf_ok===true?'　✅VAF一致':(e.vaf_ok===false?'　🔴VAF不一致':''))+'</span></div>';
  });
  s+='<details style="margin-top:5px"><summary style="cursor:pointer;font-size:11.5px">▶ pairwise 2×2（全對·驅動定序）</summary><table style="font-size:11px;margin-top:4px"><tr><th>對</th><th>RR</th><th>RA</th><th>AR</th><th>AA</th></tr>';
  (R.pairwise||[]).forEach(function(p){var c=p.cells;s+='<tr><td class="mono">S'+(p.i+1)+'-S'+(p.j+1)+'</td><td>'+c.RR+'</td><td>'+c.RA+'</td><td>'+c.AR+'</td><td>'+c.AA+'</td></tr>';});
  s+='</table></details>';
  var sv=R.subread_vectors||{};
  s+='<details style="margin-top:4px" open><summary style="cursor:pointer;font-size:11.5px"><b>▶ 區域內所有子read（'+Object.keys(sv).length+' 種向量·X=未覆蓋）</b>　覆蓋分布 '+JSON.stringify(R.cover_hist||{})+'</summary>';
  s+='<table style="font-size:11px;margin-top:4px"><tr><th>子read 向量</th><th>reads</th><th>覆蓋</th><th>類型</th></tr>';
  Object.entries(sv).slice(0,50).forEach(function(kv){var vec=kv[0],n=kv[1];var nc=k-(vec.split('X').length-1);var typ=nc===k?'✅全跨(絕對)':(nc>=2?'部分('+nc+'點)':'單點');s+='<tr><td class="mono">'+vec+'</td><td>'+n+'</td><td>'+nc+'/'+k+'</td><td>'+typ+'</td></tr>';});
  s+='</table></details>';
  if((R.blind_flags||[]).length) s+='<div class="note" style="color:#a37200;margin-top:4px">⚠ 盲點旗標：'+R.blind_flags.join('、')+'</div>';
  s+='</div>';return s;
}
function geneBlock(region){
  let g=(D.gene||{})[region]; if(!g) return '';
  let cg=Object.entries(g.cancer_genes||{}).map(([n,o])=>`<span style="background:#ffe3e3;border-radius:3px;padding:1px 4px">${n}${o.role?'('+o.role+')':''}${o.tier?' T'+o.tier:''}</span>`).join(' ');
  let dr=Object.entries(g.druggable_genes||{}).map(([n,ds])=>`<span style="background:#d3f9d8;border-radius:3px;padding:1px 4px" title="${ds.join(', ')}">${n}💊</span>`).join(' ');
  let names=(g.protein_coding||g.genes||[]).slice(0,12);
  return `<div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;padding:8px;margin-top:8px">
    <b>🧬 基因註釋(GENCODE+DGIdb${cg?'+COSMIC':''})</b>
    <div class="note">基因(${g.n_genes}): ${names.join(', ')||'(無 protein-coding)'}${g.has_promoter?` · <b style="color:#1971c2">含啟動子</b>(${(g.promoter_genes||[]).slice(0,5).join(',')})`:''}</div>
    ${cg?`<div style="margin-top:4px">🔴 癌症基因(COSMIC): ${cg}</div>`:''}
    ${dr?`<div style="margin-top:4px">💊 可用藥(DGIdb): ${dr}</div>`:'<div class="note" style="margin-top:4px">此區無 DGIdb 可用藥基因</div>'}
  </div>`;
}
['f_chr','f_minc','f_tpfp','f_loh','f_undef','f_gene','f_q','f_sort','f_sortdir'].forEach(id=>{let e=el(id);if(e){e.oninput=render;e.onchange=render}});
render();
// ===== 確認佇列(評分 + 左右判讀) =====
const SC=D.scoring;
el('scoresum').innerHTML=`需確認 <b>${SC.summary.n_need_confirm}</b>/${SC.summary.n_total} 區 · 評分桶 ${JSON.stringify(SC.summary.score_buckets)} · situation ${JSON.stringify(SC.summary.situation_dist)} · <span title="06-28 cis-control 已否決:乾淨可用≈0">曾標需甲基 ${SC.summary.needs_methyl_n}(已否決·非真可用)</span> · 公式: ${SC.summary.score_formula}`;
let Q=SC.queue;
el('q_sit').innerHTML='<option value="">全</option>';[...new Set(Q.map(q=>q.situation))].sort().forEach(s=>{let o=document.createElement('option');o.value=s;o.textContent=s;el('q_sit').appendChild(o)});
const QSORT={score:(a,b)=>a.confidence_score-b.confidence_score,scoreD:(a,b)=>b.confidence_score-a.confidence_score,coord:(a,b)=>a.chrom.localeCompare(b.chrom,undefined,{numeric:true})||a.start-b.start};
const jkey=r=>'topo_judge_'+r;
window.setJ=(r,v)=>{let cur=localStorage.getItem(jkey(r));localStorage.setItem(jkey(r),cur==v?'':v);renderQ()};
window.showByRegion=function(reg){let i=det.findIndex(r=>r.region==reg);if(i>=0){show(i,null);el('detail').scrollIntoView({behavior:'smooth',block:'start'});}else{alert('此佇列區不在拓樸明細中(n_sSNV<2 或未載):'+reg);}};
const scolor=s=>s>=80?'#2b8a3e':s>=60?'#1971c2':s>=40?'#e8590c':'#c92a2a';
function renderQ(){let sit=el('q_sit').value,mo=el('q_methyl').checked,so=el('q_sort').value;
 let f=Q.filter(q=>(!sit||q.situation==sit)&&(!mo||q.needs_methyl));f.sort(QSORT[so]||QSORT.score);
 el('qcnt').textContent=f.length+' 區';
 el('queue').innerHTML=f.slice(0,500).map(q=>{let j=localStorage.getItem(jkey(q.region))||'';
  return `<div class="row" style="display:flex;gap:7px;align-items:center;flex-wrap:wrap">
   <span style="width:140px;cursor:pointer" onclick="showByRegion('${q.region}')" title="點看上方 clone 樹"><b style="color:#1971c2;text-decoration:underline">${q.region}</b></span>
   <span style="width:58px;color:${scolor(q.confidence_score)};font-weight:700" title="confidence 0-100">▮${q.confidence_score}</span>
   <span class="tag ctx_${q.genome_ctx}">${q.genome_ctx}</span>
   <span style="width:112px;font-size:11px">${q.situation}</span>
   <span style="width:96px;font-size:10.5px" class="note">候選 ${q.n_candidates}${q.parsimony_first_rank_prob!=null?` · P1=${q.parsimony_first_rank_prob}`:''}</span>
   <span style="flex:1;min-width:170px;font-size:10.5px" class="note">${q.resolution_path}</span>
   <span style="width:100%;font-size:10px;color:#a33" class="note">🔎 為何: ${q.why_conflict||''}${q.truncated?' ⚠截斷':''}</span>
   <span style="width:100%;font-size:10px;color:#268" class="note">🧬 甲基: ${q.methyl_applicability||''}</span>
   <span style="white-space:nowrap">
     <button onclick="setJ('${q.region}','agree')" style="font-size:11px;background:${j=='agree'?'#d3f9d8':'#fff'}">✓同意rank1</button>
     <button onclick="setJ('${q.region}','alt')" style="font-size:11px;background:${j=='alt'?'#fff3bf':'#fff'}">⇄偏好其他</button>
     <button onclick="setJ('${q.region}','more')" style="font-size:11px;background:${j=='more'?'#ffe3e3':'#fff'}">?需更多資訊</button>
   </span></div>`}).join('')+(f.length>500?`<div class="note" style="padding:8px">...前 500（共 ${f.length}，可篩選縮小）</div>`:'');
}
['q_sort','q_sit','q_methyl'].forEach(id=>{el(id).onchange=renderQ;el(id).oninput=renderQ});
el('q_exp').onclick=()=>{let j={};Q.forEach(q=>{let v=localStorage.getItem(jkey(q.region));if(v)j[q.region]={judgment:v,score:q.confidence_score,situation:q.situation}});
 let b=new Blob([JSON.stringify({n:Object.keys(j).length,judgments:j},null,1)],{type:'application/json'});
 let a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='topology_judgments.json';a.click()};
renderQ();
}  // end bootWS — 每分頁切換時用新樣本資料重跑(function scope,無重宣告衝突)
function selectSample(s){
  window.__DATA__ = window.__SAMPLES__[s];
  document.querySelectorAll('.stab').forEach(t=>t.classList.toggle('active', t.dataset.s===s));
  try{ bootWS(); }catch(e){ ['universe','s_topo','s_clust','s_det','s_root','s_nsnv','unlocatable','list','queue','scoresum'].forEach(id=>{let el2=document.getElementById(id);if(el2)el2.innerHTML='';}); document.getElementById('detail').innerHTML='<b style="color:#c00">⚠ 本樣本('+s+')載入失敗: '+e.message+'。其餘面板已清空避免混樣本誤讀。</b>'; }
}
(function(){
  let names=Object.keys(window.__SAMPLES__||{});
  let bar=document.getElementById('sampletabs');
  let TI={'HCC1395':'凍結基準(SEQC2 truth set)','HCC1395_DORADO':'同細胞株·Dorado 重 basecall(非 SEQC2 truth、census 衍生標籤)'};
  if(bar){ bar.innerHTML = names.map(n=>`<button class="stab" data-s="${n}" onclick="selectSample('${n}')" title="${TI[n]||n+'(census 衍生 TP/FP 標籤弱)'}">${n}</button>`).join(''); }
  if(names.length) selectSample(names[0]);
})();
"""

PROVENANCE_FOOTER = ('<p class="note" style="margin-top:8px;color:#888">'
                     'build_branch: research/subclonal-reconstruction-202606 · '
                     '資料來源(隨樣本)：HCC1395=凍結 @feat/summary-nreadsvalid@5308d9e(SEQC2 truth);其餘 6 樣本=multisample_subclone census 衍生(各分頁宇宙帳本有 caveat,TP/FP 標籤弱)· '
                     '姊妹編號 = 子樹總 read 數遞減（?-1=該 lineage 分支佔比較大，含所有子孫；?-2=較小）· '
                     '甲基 = bounded-auxiliary（見 20260628_cis_control_scope_pilot_verdict_01.md）</p>')

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>多樣本克隆樹拓樸工作站 — {len(SAMPLE_NAMES)} ONT 樣本 sSNV 重建</title><style>{CSS}
.tabs{{display:flex;gap:4px;flex-wrap:wrap;margin:8px 0;border-bottom:2px solid #dee2e6;padding-bottom:0}}
.stab{{padding:7px 14px;border:1px solid #dee2e6;border-bottom:none;border-radius:7px 7px 0 0;background:#f1f3f5;cursor:pointer;font-size:13px;font-weight:600;color:#495057}}
.stab.active{{background:#1971c2;color:#fff}}
</style></head><body><div class="wrap">
<h1>多樣本克隆樹拓樸互動工作站（cluster-first + S/r/m 標籤 + 基因註釋）</h1>
<p class="sub">{len(SAMPLE_NAMES)} ONT 樣本（{", ".join(SAMPLE_NAMES)}）· 每區 genotype 向量→拓樸(perfect-phylogeny+噪聲過濾) · 分頁切換樣本 · S=sSNV/r=read群/m=甲基位點 · 數字由 JSON 注入</p>
<div id="sampletabs" class="tabs"></div>
<div class="zone">📊 整體觀察區（此區塊四張看板全部隨上方分頁樣本變）</div>
<div class="note" style="margin:2px 0 7px;line-height:1.75;background:#f8f9fa;border:1px solid #e9ecef;border-radius:6px;padding:7px 11px">由上而下四張看板：<b>① 記分卡</b>＝此樣本各判讀結果／截斷／LOH／癌基因命中的<b>一眼數量</b>(與評分佇列同口徑)。<b>② 全 sSNV 宇宙帳本</b>＝所有 somatic sSNV 分三桶：<b>linked</b>(有共讀 partner→可建樹)／<b>underpowered</b>(有 partner 無共讀→加深覆蓋可救)／<b>isolated</b>(read-span 內無 partner)。<b>③ HG38 全基因組分布</b>＝每區落在染色體<b>哪個位置</b>(依結果上色)＋<b>LOH 底帶</b>(紫)；toggle 可切 每-sSNV／樹形。<b>④ 四張統計卡</b>＝拓樸型態／群數 c／determinacy／HP 根數(<b>點任一卡</b>看放大圓餅＋類別逐一解釋)。</div>
<details style="margin:2px 0 7px;border:1px solid #ffe0a3;border-radius:6px;background:#fffdf5"><summary style="cursor:pointer;padding:6px 11px;font-weight:600;color:#a37200">⚠ 方法透明度：這張工作站沒顯示什麼 / 已知損失（點開）</summary><div class="note" style="padding:4px 12px 9px;line-height:1.75">
① <b>complete-case（只吃全覆蓋 read）</b>：一條 read 要跨該區全部 sSNV 才算一個 genotype 向量，部分覆蓋的 read <b>靜默丟棄（不 impute）</b>——這是對的（impute 會捏造 linkage，chr17 bug 已證）。代價：HCC1395 有 <b>40.4%(2,887/7,143)</b> 個 n_sSNV≥2 區向量全空被排除單分子拓樸；42 密集區截斷丟 1,132 sSNV。<b>逐區丟棄比例目前未逐區記錄（上游待補）</b>。<br>
② <b>read-span 物理限制</b>：co-read 隨距離衰減（中位 88@&lt;1kb → 11@20-50kb）；A_determined 率隨 span 崩（2-10kb 70% → &gt;50kb 僅 11%）；38.9% 的 sSNV underpowered/isolated（距離&gt;read-span 連不起來）。→ 故本方法是 <b>regional partition（≤read-span）非 genome-wide tree</b>（⭐3 定位）。<br>
③ <b>建樹表達力</b>：row-laminar <b>畫不出 subclone-of-subclone</b>（深分支）；「branched」實為「多根 lineage」（見「拓樸型態」卡）。<br>
④ <b>四配子檢定</b>：C2 修正（07-01）已把非 CN-gain 乾淨違反改判 incompatible（incompatible 12→118 區）；CN-gain multiplicity 假違反仍需 CN-aware genotyping（未做）。<br>
🔴 所有 perfect-phylogeny/IDPP 建構的<b>有效性依賴 infinite-sites</b>，癌症 LOH/CNV（本專案 82-91%）系統違反。完整稽核（Q1-Q4 + 定理來源）→ <b>InterSubMod/docs/methodology/20260701_topology_algorithm_audit_findings_01.md</b>。</div></details>
<div id="scorecard"></div>
<div id="universe"></div>
<details open style="margin:6px 0"><summary style="cursor:pointer;font-weight:700;color:#1971c2;font-size:13.5px;padding:3px 0">🗺️ HG38 全基因組分布（GRCh38 座標 · LOH 底帶 · 結果/樹形/每sSNV；點此摺疊/展開）</summary><div id="ideogram"></div></details>
{GLOSSARY_HTML}
<div class="stats">
<div class="scard"><h4>拓樸型態<span class="more">▸ 點看細節</span></h4><div class="note" style="font-size:10px;margin:-4px 0 3px">read 群在系統樹的<b>形狀</b>:single/linear/branched/star</div><div id="s_topo"></div></div><div class="scard"><h4>群數 c<span class="more">▸ 點看細節</span></h4><div class="note" style="font-size:10px;margin:-4px 0 3px"><b>觀測到的不同 ALT 細胞群數</b>(germline 不計·≤k)</div><div id="s_clust"></div></div>
<div class="scard"><h4>determinacy<span class="more">▸ 點看細節</span></h4><div class="note" style="font-size:10px;margin:-4px 0 3px">樹<b>存在≠能辨識</b>是哪棵:A單分子/B連鎖/C單群/成環</div><div id="s_det"></div></div><div class="scard"><h4>HP 根數<span class="more">▸ 點看細節</span></h4><div class="note" style="font-size:10px;margin:-4px 0 3px">somatic 散在<b>幾條 germline HP</b>;≥2=跨HP(allelic)</div><div id="s_root"></div></div>
<div class="scard" style="cursor:default"><h4>建樹位點分布<span class="more">區×sSNV</span></h4><div class="note" style="font-size:10px;margin:-4px 0 3px">每區有<b>幾個 sSNV</b>(2/3/…/&gt;8)的區數分布</div><div id="s_nsnv"></div></div>
</div>
<div class="note" style="margin:-2px 0 8px">⚠ 分母註:拓樸型態/群數/HP根數＝<b>全區宇宙</b>(含無 genotype 向量區);<b>determinacy 只計有向量的區(n_sSNV≥2)</b>→合計較小、與另三卡分母不同(非錯誤)。兩宇宙對照見「建樹位點分布」卡。</div>
<div class="smbg" id="statmodal" onclick="if(event.target===this)closeStatModal()"><div class="smbox"><button class="x" onclick="closeStatModal()">×</button><div id="statmodal_body"></div></div></div>
<details class="c17" id="chr17wrap"><summary>📚 克隆樹判讀圖鑑（固定教學・範例取自 HCC1395，與上方分頁樣本無關；點開）</summary><div id="teachbody"></div></details>
<details open style="margin:6px 0"><summary style="cursor:pointer;font-weight:700;color:#1971c2;font-size:13.5px;padding:3px 0">📐 拓撲 × 群數 × 樣本 三軸觀察（描述性·L3·partly-artifact；點此摺疊/展開）</summary><div id="threeaxis"></div></details>
<div id="unlocatable"></div>
<div class="zone">🔬 樣本各區域檢視（篩選 → 點選一個區 → 看克隆樹 / 四配子衝突 / 逐區甲基判定）</div>
<div class="ctrl">
chr<select id="f_chr"><option value="">全</option></select>
排序<select id="f_sort"><option value="coord">座標</option><option value="nsnv">複雜度(全sSNV)</option><option value="glen">建樹位點數(genotype)</option><option value="nclust">群數</option><option value="region">region名</option></select><select id="f_sortdir" title="正序/反序"><option value="">↓預設</option><option value="rev">↑反序</option></select>
最少群數<input id="f_minc" type="number" value="0" min="0" max="6" style="width:52px">
TP/FP<select id="f_tpfp"><option value="all">全部</option><option value="tp">只含TP</option><option value="fp">只含FP</option><option value="both">同時TP&amp;FP</option></select>
<label title="LOH 區(longphase-TO LOH.bed overlap 或 cn=loh;7 樣本都用 longphase-TO SEQC2-validated)"><input id="f_loh" type="checkbox">僅 LOH 區</label>
<label title="分支順序未定/不相容;曾標『需甲基輔助』,06-28 cis-control 已否決(乾淨可用≈0)"><input id="f_undef" type="checkbox">僅無法定義(曾標需甲基·已否決)</label>
<label title="該區基因註釋含 COSMIC 癌症基因"><input id="f_gene" type="checkbox">僅癌基因命中</label>
搜尋<input id="f_q" placeholder="chr17:" style="width:120px">
<span id="cnt" class="note"></span>
</div>
<div class="zone" style="margin-top:14px;font-size:12.5px;background:#fff9db;color:#b08900;border-left-color:#ffd43b">☑ 勾選要觀察的狀況（三組可複選；不勾＝該組全納入）</div>
<div class="facets">
 <div class="fcard"><div class="fh">🌳 拓樸型態<span class="fhd">分子累積形狀</span></div><div class="chips" id="cb_topology_type"></div></div>
 <div class="fcard"><div class="fh">🎯 determinacy<span class="fhd">能否唯一辨識</span></div><div class="chips" id="cb_determinacy"></div></div>
 <div class="fcard"><div class="fh">📍 基因體位置<span class="fhd">偽影風險</span></div><div class="chips" id="cb_genome_ctx"></div></div>
</div>
<div class="main"><div class="list" id="list"></div><div class="detail" id="detail"><div class="note">← 左側點選一個區查看克隆樹（或點上方 chr17 worked example）</div></div></div>

<h3 style="margin-top:20px">✓ 候選評分確認佇列（左右選項判讀 + 觀察評分；存瀏覽器 localStorage、可匯出）</h3>
<div id="scoresum" class="note"></div>
<div class="ctrl">
排序<select id="q_sort"><option value="score">評分(低→高,最需關注)</option><option value="scoreD">評分(高→低)</option><option value="coord">座標</option></select>
situation<select id="q_sit"><option value="">全</option></select>
<label title="06-28 cis-control 裁決:這些區甲基乾淨可用≈0,非真能用甲基解"><input id="q_methyl" type="checkbox">曾標需甲基(已否決)</label>
<button id="q_exp">匯出判讀 JSON</button><span id="qcnt" class="note"></span>
</div>
<div class="list" id="queue" style="max-height:62vh"></div>
<p class="note" style="margin-top:12px">⚠ 證據層級：A_determined=單分子向量唯一可辨識(≠對 single-cell 驗證為真)；A_ambiguous=缺中間群順序未定；B_pairwise=拼接非單分子整樹；C_underdetermined=多樹相容。TP/FP=SEQC2 僅觀察不進前處理。genome_ctx 為近似(±3Mb)。甲基不參與拓樸裁決(cis-confounded;06-28 cis-control 已測→bounded-auxiliary,非 resolver)。⭐3 單樣本·regional(≤read-span)非 genome-wide tree·分子共現≠single-cell。</p>
{PROVENANCE_FOOTER}
</div>
<script>window.__SAMPLES__={SAMPLES_JSON};window.__CHR17TREE__={CHR17TREE_JSON};window.__RESOLUTION__={RESOLUTION_JSON};window.__GT8__={GT8_JSON};</script><script>{JS}</script></body></html>"""
with open(MULTI_OUT, "w", encoding="utf-8") as f: f.write(HTML)
print(f"OK wrote {MULTI_OUT} ({len(HTML):,} bytes; samples {SAMPLE_NAMES})")
