#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One-off v6 migration builder: merge graph_v2_flat (30 flat) + graph_ism (ISM deep tree)
into a single complete v3 WBS tree (state/tasks/graph.json) + add 盤點缺漏 nodes.
Pulls real node content from the existing graphs (no re-typing / no fabrication); new nodes
are memory-sourced. Writes via task_graph.dump_graph (one-node-per-line). Throwaway after run."""
import json, os, datetime, importlib.util

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
spec = importlib.util.spec_from_file_location("tg", os.path.join(ROOT, "scripts/tasks/task_graph.py"))
tg = importlib.util.module_from_spec(spec); spec.loader.exec_module(tg)

v2 = json.load(open(os.path.join(ROOT, "state/tasks/graph_v2_flat.json.bak"), encoding="utf-8"))
ism = json.load(open(os.path.join(ROOT, "state/tasks/graph_ism.json"), encoding="utf-8"))
v2n = {n["id"]: n for n in v2["nodes"]}
ismn = {n["id"]: n for n in ism["nodes"]}
TODAY = "2026-06-20"
nodes = []

# confirmed goals (chat confirmation batch 1: T-ISM / T-SL — user 視為確認)
ismn["T-ISM"]["goal"] = "結構×標籤驗證(T-SL)落地 + a-priori 分類可信 + V-2 系統掃描定 R5 → ISM 方法可寫進 ch3 且結果經得起 reviewer"
ismn["T-SL"]["goal"] = "PERMANOVA/ARI/reclassify(已done) + within-HP 選k去循環(C1) + HP-fine組合(C3) + tumor-only重跑(RO) + Stage③(S3) 做完 → 結構↔標籤關聯有 FDR 校準的可信判定；allele≫HP 用 cis-control 分清 cis-ASM vs subclone"
# batch-2 confirmed goals (T-SL 直接子)
ismn["T-SL-PERM"]["goal"] = "（已達）顯著判定不過嚴/鬆，valid 92.9%/Noise 1.5%，dispersion 混淆排除"
ismn["T-SL-ARI"]["goal"] = "（已達）ARI 讓「幾何群是否對齊標籤」可見（median HP0.145/allele0.005）"
ismn["T-SL-RECLS"]["goal"] = "（已達）Noise/Weak→valid 救回 53.9%，orig Strong 降級 FP=0"
ismn["T-SL-C1"]["goal"] = "enable_bootstrap=true + consensus 重抽只留穩定群，去 raw silhouette 過切循環"
ismn["T-SL-C3"]["goal"] = "C3 C++ emit「組合對齊測試欄 + 一群含多標籤 flag」"
ismn["T-SL-RO"]["goal"] = "tumor-only 獨立重跑 + 兩層事實表（只記「tumor顯著/normal不顯著」不判別）"
ismn["T-SL-S3"]["goal"] = "coverage-peak CN 自估 + SEQC2 只驗證後觀察不當輸入"
ismn["T-CSC"]["goal"] = "ISM 完成後把分類結果判讀到 clone/subclone 層（a-priori 條件軸 + 6 樣本外部真值對照）"
# batch-3 confirmed goals (T-S 系列, v2-sourced)
v2n["T-S1"]["goal"] = "a-priori 4-pop 分類定案(ADOPT_WITH_CORRECTIONS)；甲基用於 A 刻畫(B排序illustrative/C drop)；somatic 歸屬待 G-B；within_clean≠subclone(Jaccard 0.123)"
v2n["T-S5"]["goal"] = "收斂 Path B model-based；確認≥2群須對齊 a-priori 軸(CramérV≥0.7+Fisher)非循環；無監督抽象幾群無乾淨解"
v2n["T-S6"]["goal"] = "（已達）對齊=paired 非 tumor-only 完整盤點，verify 53/53"
v2n["T-S2"]["goal"] = "（結案/勿再開）tumor-only 非監督=double-dip NEGATIVE"


def reparent(nid, parent):
    n = dict(v2n[nid]); n.pop("component", None); n["parent"] = parent; return n


def C(cid, title, headline, owner="claude", status="in_progress"):  # container (milestone)
    return {"id": cid, "title": title, "parent": "T-PAPER", "kind": "milestone", "status": status,
            "depends_on": [], "owner": owner, "headline": headline, "added_at": TODAY}


def N(nid, title, parent, kind, status, owner="claude", depends_on=None, headline=None,
      io=None, links=None, missing=None, notes=None, verify=None, goal=None):
    n = {"id": nid, "title": title, "parent": parent, "kind": kind, "status": status,
         "depends_on": depends_on or [], "owner": owner, "added_at": TODAY}
    if headline: n["headline"] = headline
    if goal: n["goal"] = goal
    if io: n["io"] = io
    if links: n["links"] = links
    if missing: n["missing_info"] = missing
    if notes: n["notes"] = notes
    if verify: n["verify"] = verify
    return n


# ---- root (keep ism's T-PAPER) ----
nodes.append(ismn["T-PAPER"])

# ---- 10 main-task containers ----
nodes.append(C("T-ASSET", "資料資產", "6 樣本 BAM/VCF/ISM + normal 甲基（COLO829 缺=主阻塞）", "user"))
nodes.append(C("T-PHASE", "甲基-assist phasing", "甲基救 unphase / haplotag 輔助（V1-V12，單樣本 ⭐3）"))
nodes.append(C("T-ASM", "ASM characterization", "等位基因特異甲基化刻畫（⭐3-4；真實但非 filter）"))
nodes.append(ismn["T-ISM"])  # 完成 ISM 程式 (container, keeps headline/summary)
nodes.append(C("T-METHOD", "方法健全性", "ISM 方法合理性 + 分群驗證 + 已結支撐研究"))
nodes.append(C("T-LOCUS", "位點層驗證", "chr2:18M 等具體位點的 subclone 結構驗證"))
nodes.append(C("T-EXT", "外部驗證", "文獻/工具/真值比對庫 + 引用驗證"))
nodes.append(C("T-NEG", "filter-NEGATIVE 四道脊柱", "甲基非 filter 的 4 道獨立負證（論文 Methods-Negative）", "claude"))
nodes.append(C("T-GATE", "決策 gates（OPEN）", "定論文 Grade / tier / 甲基-subclone 故事的 3 個未決", "claude+user"))
nodes.append(C("T-WRITE", "論文撰寫", "中文碩論章節 / 圖 / 表 / 摘要", "claude+user"))

# ---- ISM deep subtree (all T-SL*) verbatim from graph_ism (keeps claims/problem/verify) ----
for nid, n in ismn.items():
    if nid == "T-SL" or nid.startswith("T-SL-"):
        nodes.append(n)
# clone/subclone judgment (ism T-CSC) under T-ISM
csc = dict(ismn["T-CSC"]); csc["parent"] = "T-ISM"; nodes.append(csc)

# ---- re-parent v2 flat nodes into containers ----
for nid in ("T-A1", "T-A2", "T-M1"):
    nodes.append(reparent(nid, "T-ASSET"))
for nid in ("T-C1", "T-C2", "T-C3", "T-C4", "T-C5", "T-C6"):
    nodes.append(reparent(nid, "T-ASM"))
for nid in ("T-P1", "T-P2"):
    nodes.append(reparent(nid, "T-PHASE"))
for nid in ("T-S1", "T-S2", "T-S5", "T-S6"):   # S3/S4 dropped — covered by ISM subtree RECLS-2 / S3
    nodes.append(reparent(nid, "T-ISM"))
for nid in ("T-V1", "T-V2", "T-V3", "T-V4", "T-V5"):
    nodes.append(reparent(nid, "T-METHOD"))
for nid in ("T-L1", "T-L2"):
    nodes.append(reparent(nid, "T-LOCUS"))
for nid in ("T-E1", "T-E2"):
    nodes.append(reparent(nid, "T-EXT"))

# ---- new nodes: gaps from inventory ----
# assets
nodes.append(N("T-A3", "6 normal 甲基帳號 SPOF 備份", "T-ASSET", "data", "todo", "user",
               missing=["6 normal 全 zhenyu112 帳號 = 單點失效"], headline="6 normal 全 zhenyu112 帳號 → 需備份避免單點失效",
               links={"memory": ["project_subclonal_reconstruction_paper_focus"]}))
# asm cross-sample production
nodes.append(N("T-C-CROSS", "跨樣本 ASM 正式圖表（R6）", "T-ASM", "analysis", "todo", "claude",
               depends_on=["T-C2"], headline="把 6/6 excess-over-null 做成論文正式圖表（現只有計算）",
               io={"inputs": [{"name": "6 樣本 excess 表", "required": True, "ref": "T-C2"}], "outputs": ["R6 跨樣本圖 + Table"]},
               links={"memory": ["project_cross_sample_asm_reproducibility"]}, missing=["待 G-A 跨樣本統計"]))
# external citation
nodes.append(N("T-E3", "citation-verification（投稿前）", "T-EXT", "writing", "todo", "claude",
               headline="每 citation WebSearch+Scholar 驗證才寫入 .bib", links={"memory": ["project_external_validation_library"]}))
# method-soundness 支撐 done（連 memory）
nodes.append(N("T-SUP-CLUST", "clusterability vs coverage/CN（支撐 ⭐3）", "T-METHOD", "analysis", "done", "claude",
               headline="分群力驅動=n_CpG 非 depth/CN；LOH 抑制 6/6", links={"memory": ["project_subclonal_reconstruction_paper_focus"]},
               notes="20260610 報告；併論文 R6"))
nodes.append(N("T-SUP-COPY", "甲基×copy confound（SEQC2，支撐 ⭐3）", "T-METHOD", "analysis", "done", "claude",
               headline="CN vs AUC 無相關(ρ=0.035)；copy 非 driver；副產品可偵 LOH", links={"memory": ["project_asm_cn_confound_pilot"]}))
nodes.append(N("T-SUP-LOCUS", "單位點甲基組合詮釋窮舉（支撐 ⭐3）", "T-METHOD", "analysis", "done", "claude",
               headline="D1-D7×I1-I8 窮舉；最像 subclone 組合 undecidable", links={"memory": ["project_O12_loh_methylation_scenarios"]}))
# NEGATIVE 四道
nodes.append(N("T-NEG-D1", "D1 乾淨 somatic-cis 稀有", "T-NEG", "analysis", "done", "claude",
               headline="HCC1395 816 可測位點僅 1（chr17, perm p=0.001）", links={"memory": ["project_paper_claim_audit_consensus_base_2026_06_12"]}))
nodes.append(N("T-NEG-D2", "D2 甲基=germline-allelic 多重共線", "T-NEG", "analysis", "done", "claude",
               headline="甲基主驅動=haplotype 差非 somatic（normal HP1/HP2 固定）", links={"memory": ["project_methyl_phasing_assist_line"]}))
nodes.append(N("T-NEG-D3", "D3 甲基分群=germline 非 somatic 驅動", "T-NEG", "analysis", "done", "claude",
               headline="無監督分群對齊 germline 載體 85% 非 somatic；tumor-only double-dip", links={"memory": ["project_tumor_only_axis_negative_subclone_classification"]}))
nodes.append(N("T-NEG-D4", "D4 NGroups=phasing 非甲基訊號", "T-NEG", "analysis", "done", "claude",
               headline="NGroups 由 HP-tag 數決定（LabelTest.cpp:265），甲基未參與；HD-4 RESOLVED", links={"memory": ["project_hpfinengroups_subclone_marker"]}))
# 決策 gates (OPEN)
nodes.append(N("T-GATE-HD1", "HD-1：R-SELFREF → 定論文 Grade", "T-GATE", "analysis", "todo", "claude+user",
               headline="phasing 脊柱 by-construction 循環 → 跑 R-SELFREF(~25-50hr) 對照 positive-spine 或降 characterization",
               missing=["未決；影響論文主張強度"], links={"memory": ["project_subclonal_reconstruction_paper_focus"]},
               verify={"how": "R-SELFREF C++ 對照", "ref": "/run-evaluator"}))
nodes.append(N("T-GATE-GA", "G-A：跨 6 樣本 → 定單樣本 ⭐3 vs ⭐4", "T-GATE", "analysis", "todo", "claude+user",
               depends_on=["T-M1"], headline="V1-V12 跨 6 樣本重跑（normal 甲基 5/6 ready）→ 決定 phasing tier",
               missing=["待 COLO829 normal 甲基；G-A 統計未跑"], links={"memory": ["project_methyl_phasing_assist_line"]}))
nodes.append(N("T-GATE-GB", "G-B：within-hap somatic null → 定甲基-subclone 故事", "T-GATE", "analysis", "todo", "claude+user",
               headline="within-hap somatic-vs-baseline 對照（非 germline-het null）→ 甲基-subclone 是 somatic 還 germline",
               missing=["未跑前甲基-subclone 只能寫存在性窄+負"], links={"memory": ["project_apriori_subclone_classification_model"]}))
# 論文撰寫 細化
nodes.append(N("T-W-CH1", "Ch1 緒論", "T-WRITE", "writing", "todo", "claude+user", headline="廣→gap（甲基加值未量化/解析度天花板）→三解析度目標", links={"memory": ["project_thesis_writing_architecture"]}))
nodes.append(N("T-W-CH2", "Ch2 文獻探討", "T-WRITE", "writing", "todo", "claude+user", headline="重建骨幹→ONT 甲基→ASM→演化邊界→分群 4 軸地景→ISM 定位", links={"memory": ["project_ism_vs_external_methylation_tools_comparison"]}))
nodes.append(N("T-W-CH3", "Ch3 材料與方法（骨架已有）", "T-WRITE", "writing", "in_progress", "claude+user", headline="6 cell line×3 癌種→haplotag 6-state→5mCG→ISM 6 核心→null 設計", links={"memory": ["project_thesis_writing_architecture"], "reports": ["InterSubMod/docs/thesis/"]}))
nodes.append(N("T-W-CH4", "Ch4 結果（骨架已有/填數中）", "T-WRITE", "writing", "in_progress", "claude+user", headline="R1 骨幹→R2 germline 強→R3 救援88.5%→R4 subclone 僅存在→R5 LOH 無confound→R6 跨樣本→R7 owned NEGATIVE", links={"memory": ["project_thesis_writing_architecture"]}, missing=["R6 待 G-A；數字待真值"]))
nodes.append(N("T-W-CH5", "Ch5 討論", "T-WRITE", "writing", "todo", "claude+user", headline="最強 takeaway→subclone 天花板 bound→對齊文獻→bounded-claims 專段", links={"memory": ["project_g6_paper_framing_external_corroboration"]}))
nodes.append(N("T-W-CH6", "Ch6 結論", "T-WRITE", "writing", "todo", "claude+user", headline="甲基加了什麼（germline 救援/分支）+ 沒加什麼（subclone 判別/filter）"))
nodes.append(N("T-W-ABS", "摘要（首段點明分工）", "T-WRITE", "writing", "todo", "claude+user", headline="ABT 漏斗；首段 haplotag 骨幹 / 甲基 characterize 有界", missing=["🔴 須點明分工；待 G-A 數字"], links={"memory": ["project_thesis_writing_architecture"]}))
nodes.append(N("T-W-FIG1", "Fig1 方法總覽 schematic ⭐最高優先", "T-WRITE", "writing", "todo", "claude+user", headline="方法總覽示意（spec 已定，實圖未做）", links={"memory": ["project_thesis_writing_architecture"]}))
nodes.append(N("T-W-FIGS", "Fig2-6 其餘圖（含跨樣本 R6）", "T-WRITE", "writing", "todo", "claude+user", headline="Fig2-6 需驗真值；Fig6 跨樣本待 G-A", depends_on=["T-C-CROSS"], missing=["5/6 Fig 物理不存在"]))
nodes.append(N("T-W-TABS", "Table1-3（樣本統計/工具對照/NEGATIVE 摘要）", "T-WRITE", "writing", "todo", "claude+user", headline="Table1 樣本統計 / Table2 工具矩陣 / Table3 NEGATIVE 摘要", missing=["Table1/3 缺"]))
# 整合篇章 + 共識底座 (from v2 W3/W4)
w3 = reparent("T-W3", "T-WRITE"); nodes.append(w3)
w4 = reparent("T-W4", "T-WRITE"); nodes.append(w4)

# ---- v6.1 audit gaps (workflow wabzbk7hf；8 驗證過缺漏) ----
nodes.append(N("T-A4", "驗證 6 樣本 DeepVariant/DeepSomatic 可用性（兩-caller 前置）+ 修 landscape05 CASTLE-COLO829 過度宣稱", "T-ASSET", "data", "todo", "claude+user",
               depends_on=["T-A1"], headline="6/6 canonical 皆無 DeepVariant/DeepSomatic VCF（只 ClairS）→ 兩-caller 判準/Fig2 阻塞",
               goal="確認/補齊第二 caller（ClairS∩DeepVariant）或明標範圍 + 修 landscape05 provenance（COLO829 只在 HCC1395 repo）",
               missing=["6/6 無 DeepVariant；landscape05 CASTLE-COLO829 provenance 錯"], links={"memory": ["project_external_validation_library"]}))
nodes.append(N("T-C-CISEXP", "V-1 乾淨 somatic-cis 真稀有性擴測（normal-anchored × 全合格位點 × 6 樣本）", "T-ASM", "compute", "todo", "claude",
               headline="cis-test 擴到 TAG-C 12,868 + TAG-E 28,254 ×6 樣本算真分母 → 定『稀有』是生物現象 vs under-test artifact",
               goal="全合格位點 normal-anchored cis-test 完成、得真分母、定 headline 措辭 + de-confound R3 + chr17 例外意義",
               missing=["catalog 全位點從未全 cis-tested（稽核校正#2）"], links={"memory": ["project_paper_claim_audit_consensus_base_2026_06_12"]}))
nodes.append(N("T-C7", "V-4 BRCA2 exemplar copy-partition 重驗（illustrative vs partial-cis 口徑）", "T-ASM", "analysis", "todo", "claude",
               depends_on=["T-C1"], headline="BRCA2 疑 copy-confounded(d_within −0.023≪d_copy −0.11) → 37_copy_partition 定有無乾淨 somatic-cis",
               goal="copy-partition null 跑完、裁決 BRCA2 是 illustrative(標caveat,cis錨改chr17) 或部分乾淨 cis；決 copy caveat 能否解除",
               links={"memory": ["project_zar1l_brca2_asm_verification"]}))
nodes.append(N("T-ISM-V2-RECON", "V-2 重建支持位點 pre-reg 判準 + 全基因組×6 掃描 + 頻率/位點 catalog", "T-ISM", "compute", "todo", "claude",
               depends_on=["T-A4"], headline="把 reconstruction 從個案(chr2)升 systematic：pre-reg 判準→genome×6 掃描→位點+頻率+每點假想樹/ARI",
               goal="pre-register『重建支持位點』判準(≥N somatic 兩-caller + 甲基↔突變 ARI 一致 + LOH 脈絡)→全掃→catalog；定 R5 EXEMPLAR climax",
               missing=["依賴 T-A4 第二 caller"], links={"memory": ["project_chr2_18m_subclone_locus_verification"]}))
nodes.append(N("T-METHOD-FDR", "V-7 FDR / 多重檢定校準：跨位點 null + BH-FDR + n_reads 校正", "T-METHOD", "compute", "todo", "claude",
               headline="跨位點 label-shuffle null + BH-FDR + chr17 Bonferroni(×816)/genome-perm + n_reads 回歸校正（投稿前必補）",
               goal="FDR/BH/Bonferroni across loci 算完 + n_reads 校正 epipolymorphism(O11 0.845→0.530)；ch3_methods 補上",
               missing=["投稿前必補 🔴"], links={"memory": ["project_code_methodology_audit_2026_06_10"]}))
nodes.append(N("T-GATE-GC", "G-C：cis vs 突變足跡 ±1-2kb 空間分離（normal-anchored 寬空間控制）", "T-GATE", "analysis", "todo", "claude",
               headline="±1-2kb 系統空間尺度控制 → 決個案能否寫『重建』而非單純 cis（R5 精確度）",
               goal="±1-2kb 空間分離方法跑完，定 chr2/chr17 個案的 cis vs footprint 區辨",
               links={"memory": ["project_apriori_subclone_classification_model"]}))
nodes.append(N("T-GATE-REDLINE", "誠實護欄紅線注入 writing 節點 verify gate + 成稿前跨章節術語掃描", "T-GATE", "writing", "todo", "claude+user",
               headline="6 護欄(含4校正)落地：每 T-W-* 加 verify 引用此 gate + on_fail=拒收/{{待填}}；成稿前禁語掃描",
               goal="14 writing 節點都有 verify gate + 成稿前跑完跨章節禁語掃描（甲基非判別/非重建驅動/守線①-④）",
               missing=["紅線目前只在 memory"], links={"memory": ["project_thesis_writing_architecture"]}))
nodes.append(N("T-W-ISMDEDUP", "去重 graph_ism.json 子樹 + 標 deprecated（唯一真值=graph.json）", "T-WRITE", "analysis", "todo", "claude",
               headline="graph_ism 26 節點全在 graph.json；雙 SoT 矛盾 → 標 deprecated 指回主 graph 或刪",
               goal="graph_ism/TASKS_ism/board_ism 標 deprecated 或刪，消除雙 SoT 矛盾"))

# ---- v6.2 audit-resume gaps (workflow wtb1nzrfa；grep-自驗 9 真新缺漏) ----
nodes.append(N("T-L3", "chr8 全-LOH 第二 exemplar locus 驗證", "T-LOCUS", "analysis", "todo", "claude",
               depends_on=["T-A1"], headline="除 chr2:18M 外的第二旗艦位點（chr8 全-LOH 熱點）→ 加厚 R5 EXEMPLAR",
               goal="chr8 熱點完成 sSNV+甲基分群實測，成為 chr2:18M 之外第二個可寫的 exemplar", links={"memory": ["project_hcc1395_chr8_hotspot"]}))
nodes.append(N("T-C-UNMASK", "LOH-unmask ASM confound（Martin-Trujillo）寫入限制", "T-ASM", "analysis", "todo", "claude",
               headline="LOH 揭露原本沉默的 ASM = 強 confound（Martin-Trujillo）→ 須在限制段明確處理",
               goal="LOH-unmask confound 在 Ch5 限制段明確敘述 + 對 cis 主張的影響界定", links={"memory": ["project_O12_loh_methylation_scenarios"]}))
nodes.append(N("T-METHOD-PSLIMIT", "phase-set 跨 block 無 read/cell 連結（限制）", "T-METHOD", "analysis", "todo", "claude",
               headline="phase-set 跨 block 無 read/cell 連結 → 重建解析度上限的方法限制",
               goal="phase-set 跨 block 限制寫進方法/討論，界定重建解析度天花板", links={"memory": ["project_O12_loh_methylation_scenarios"]}))
nodes.append(N("T-METHOD-TESTS", "Methylation MM/ML parsing 單元測試（零測試補）", "T-METHOD", "compute", "todo", "claude",
               headline="MM/ML 甲基解析目前零單元測試（全碼稽核 3 必修之一）→ 補測試",
               goal="MM/ML parsing 加單元測試覆蓋（含 edge case），ctest 綠", links={"memory": ["project_code_methodology_audit_2026_06_10"]}))
nodes.append(N("T-METHOD-BINVER", "binary_version 一致性檢核", "T-METHOD", "compute", "todo", "claude",
               headline="跨 worktree/結果的 binary 版本一致性（避免混不同 build 口徑）",
               goal="binary_versions 紀錄齊 + 結果標 build SHA，無 merge-blind 漂移", links={"memory": ["project_new_data_integrity_audit_2026_06_17"]}))
nodes.append(N("T-W-DATACODE", "Data & Code availability 聲明", "T-WRITE", "writing", "todo", "claude+user",
               headline="投稿/繳交必備：資料與程式碼可得性聲明", goal="Data & Code availability 段寫完（含 repo/資料存取說明）"))
nodes.append(N("T-W-FRONT", "論文前件 / 系所模板合規（封面/致謝/格式）", "T-WRITE", "writing", "todo", "claude+user",
               headline="碩論前件：封面/書名頁/致謝/目錄/系所格式合規", goal="前件齊全 + 系所模板格式檢查通過"))
nodes.append(N("T-W-SUPP", "Supplementary 材料封裝", "T-WRITE", "writing", "todo", "claude+user",
               headline="補充材料（額外圖表/方法細節/資料表）封裝", goal="Supplementary 封裝完成、正文引用對應"))
nodes.append(N("T-RETIRE", "退役/封存 active.json g1/g6 stale cycle（harness hygiene）", "T-GATE", "analysis", "todo", "claude+user",
               headline="g1/g6 在 active.json stale 18 天（T-C1 DRIFT）→ 走 /conclude-research 正規退役",
               goal="g1/g6 cycle 經 /conclude-research 退役/封存，active.json 無 stale、T-C1 DRIFT 消除", links={"memory": ["project_loop_engineering_harness_review"]}))

# dependency edges: gates / writing flow
for n in nodes:
    if n["id"] == "T-W-CH4":
        n["depends_on"] = ["T-C-CROSS"]  # R6 needs cross-sample

# sanitize: drop dangling depends_on / external-ify dangling io refs (dropped T-S3/T-S4 etc.)
idset = {n["id"] for n in nodes}
for n in nodes:
    n["depends_on"] = [d for d in n.get("depends_on", []) if d in idset]
    io = n.get("io")
    if io:
        for inp in io.get("inputs", []):
            if inp.get("ref") and inp["ref"] not in idset:
                inp["ref"] = None
    pr = n.get("problem")
    if pr and pr.get("fix_child") and pr["fix_child"] not in idset:
        pr["fix_child"] = None

data = {"schema_version": "3.0", "updated_at": TODAY, "focus": "T-SL",
        "note": "完整全樹（v6 migrate）：30節點 flat + ISM 深樹 合併 + 盤點缺漏（3 gate/4 NEGATIVE/章節圖表/支撐研究）。唯一機械真值；TASKS.md/tasks_board.html 自動產生勿手改。status 對 memory 核對非捏造。舊 v2 flat 備份 graph_v2_flat.json.bak。",
        "nodes": nodes}

tg.GRAPH = os.path.join(ROOT, "state/tasks/graph.json")
tg.dump_graph(data, tg.GRAPH)
print(f"[ok] wrote full tree: {len(nodes)} nodes")
