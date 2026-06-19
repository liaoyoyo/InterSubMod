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


def reparent(nid, parent):
    n = dict(v2n[nid]); n.pop("component", None); n["parent"] = parent; return n


def C(cid, title, headline, owner="claude", status="in_progress"):  # container (milestone)
    return {"id": cid, "title": title, "parent": "T-PAPER", "kind": "milestone", "status": status,
            "depends_on": [], "owner": owner, "headline": headline, "added_at": TODAY}


def N(nid, title, parent, kind, status, owner="claude", depends_on=None, headline=None,
      io=None, links=None, missing=None, notes=None, verify=None):
    n = {"id": nid, "title": title, "parent": parent, "kind": kind, "status": status,
         "depends_on": depends_on or [], "owner": owner, "added_at": TODAY}
    if headline: n["headline"] = headline
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
