---
title: COSMIC Cell Lines Project CN 補洞 + 突變交叉驗證 — 結果與對抗驗證
date: 2026-06-29
status: in_progress
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/clp_integration_summary.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/clp_crossval_confound_check.json
evidence_tier: 裁決A=L2(holds,需 SAVANA)；裁決B=L4 假設未證(confound 主導,非獨立驗證)
verify: 4-skeptic adversarial workflow(wf_68b850c5-139) + 主流程獨立再驗 skeptic 數字
---

# CLP CN 補洞 + 突變交叉驗證

## TL;DR
① 接 COSMIC CLP CN 補 6 樣本 cn=unknown → **不足（需 SAVANA）**；② CLP 突變交叉驗證 → **不能宣稱獨立驗證我們的方法**（COSMIC 循環性 + ascertainment bias 主導）。兩裁決經 4-skeptic 對抗 + 主流程獨立再驗。

---

## §1 裁決 A — CLP CN 補 cn=unknown：不足，需 SAVANA（L2 holds）

**真 CN 覆蓋率**（topology `cn` 欄；⚠ 與 §2 突變 region-hit-pct 不同指標，勿混引）：

| 樣本 | CN 填充 | 填充率 | 來源 |
|---|--:|--:|---|
| COLO829 | 0 | 0% | CLP |
| H1437 | 0 | 0% | CLP |
| H2009 | 68 | 1.6% | CLP |
| HCC1937 | 4 | 0.2% | CLP |
| HCC1954 | 36 | 1.8% | CLP |
| HCC1395_DORADO | 176 | 7.4% | CLP |
| *HCC1395（凍結）* | *3885* | *100%* | ***SEQC2（非 CLP）*** |

**裁決**：成立（操作結論：需 SAVANA），但根因改寫：
- ❌ 原措辭「CLP 太稀疏」→ ✅ **範疇錯誤**：CLP 是 gene-level aberration 驗證型 DB，本就不是全基因組 CN 分段器，用它填 WGS CN 是誤用其設計意圖。
- 🔴 **HCC1395 的 100% 來自 SEQC2-like WGS 參考、非 CLP**，不可與其他 6 樣本並列。真正缺口 = 5/6 樣本無等價 WGS CN 參考。
- H2009 cycle：265 區中 68 被填、僅 **17 新重分類**為 CN-gain-multiplicity（51 區先前已註解）。影響有限。
- 全 sSNV 僅 ~62 gain calls 命中、僅 2 個 driver（KRAS/RB1）；CLP 覆蓋在 subclone 複雜度各層一致（不偏好複雜區）。

**修法**：SAVANA 等 segment-based caller 跑我們的 BAM 產全基因組 CN segment。`sm_region_integration.py` 已加 `SM_CN_SPARSE` + CN BED 介面，將來可直接接。

---

## §2 裁決 B — CLP 突變交叉驗證：不能宣稱獨立驗證（L4 confound 主導）

**指標**（`clp_crossval_confound_check.json`）：

| 樣本 | 區含突變% | 癌區命中% | binary 富集 | span 比(mean/med) | 密度比(mut/Mb) |
|---|--:|--:|--:|--:|--:|
| HCC1395 | 4.3 | 6.3 | 1.47 | 1.64 / 1.21 | 0.88 |
| COLO829 | 3.8 | 5.3 | 1.39 | 1.17 / 1.24 | 2.04 |
| H1437 | 4.7 | 8.1 | 1.72 | 2.58 / 1.18 | 0.95 |
| H2009 | 11.9 | 22.0 | 1.85 | 2.53 / 1.26 | 1.17 |
| HCC1937 | 3.8 | 10.0 | 2.63 | 1.05 / 0.98 | 4.69 |
| HCC1954 | 3.7 | 6.3 | 1.70 | 1.03 / 1.02 | 2.55 |
| HCC1395_DORADO | 4.2 | 8.7 | 2.07 | 1.11 / 0.98 | 2.10 |

**裁決：不能宣稱「獨立驗證我們的方法」。** 經對抗 + 獨立再驗的精修結論：

1. **binary hit-rate 是弱指標**（大區天生易命中 ≥1）。改用突變密度（mut/Mb）：**5/7 樣本癌區密度反而更高**（HCC1937 4.69×、HCC1954 2.55×），僅 HCC1395(0.88)/H1437(0.95) 略低。
   - 🔴 **§13.7 校正**：skeptic 宣稱「癌區密度 0.54× 反轉（更稀）」，主流程**獨立再驗未複現**（逐樣本密度比 0.88-4.69，非 0.54）——skeptic 應為 pooled 或單樣本算法，其 headline 反轉數字**不採信**。
2. **🔴 COSMIC 循環性（結構性、無法用 null test 修）**：癌基因註解來源 = COSMIC CGC，突變目錄 = COSMIC CLP，**同源**。「癌基因區富集 COSMIC 突變」部分為**同義反覆**（測 COSMIC 突變是否命中 COSMIC 編目的癌基因）。
3. **ascertainment bias**：COSMIC 偏重定序/編目已知致癌基因 → 兩資料集獨立收斂於同一批 well-sequenced 基因，非「我們區域生物特殊」的證據。
4. **低重疊**：87-96% 區零 COSMIC 命中、僅 3.6-10.7% COSMIC 突變落我們區（平台差異：短讀 coding-focused vs ONT 全基因組 sSNV-共現）。
5. **span mean vs median 衝突已釐清**：mean 比 1.03-2.58（少數大癌區拉高），median 比僅 0.98-1.26（典型癌區只略大）。

**正確口徑**：CLP 突變交叉驗證**最多作 sanity-check**（我們區域確實落在已知變異位點，比率符合平台預期），**不可當作我們重建方法有效性的獨立外部驗證**——因為 confound（循環性 + ascertainment）結構性主導。獨立驗證仍需 single-cell / multi-region 正交（與主軸 ⭐3→⭐4 升級條件一致）。

---

## §3 對抗驗證紀錄（4 skeptic + 主流程再驗）

- workflow `wf_68b850c5-139`：2 裁決 × 2 獨立 skeptic（Explore/haiku）+ 1 綜合（opus）。
- 裁決 A：skeptic 一致同意「需 SAVANA」（措辭分歧不影響操作結論）。
- 裁決 B：2 skeptic 一致 false/high（ascertainment + 循環 + 無 null + 低重疊）。
- **主流程獨立再驗（§13.7）**：對 skeptic 1 的「0.54× 密度反轉」逐樣本重算 → 未複現（0.88-4.69，5/7 >1）。採信其**結構性批評**（循環/bias/null 缺），不採信其**特定反轉數字**。

## §4 未完成 / future work
- **SAVANA 全基因組 CN**（解 §1 cn=unknown，pipeline 級）。
- §2 若要進一步：size+gene-density-matched 隨機區 null + permutation + 已知/新發現區分層——但**循環性結構性無法用 null 消除**，故 §2 核心結論（非獨立驗證）不會因 null 翻轉。
- CLP 資料（CN BED + crossval JSON）已落檔各樣本 workdir，腳本 `clp_cn_bed.py` / `clp_crossval.py` 可重跑。
