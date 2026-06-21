---
title: "口徑修正草稿 — T2 OVERSTATED（甲基只證 1-1/2-1 可分；歸 H3 未證）"
date: 2026-06-22
status: draft
task: T-P2
audience: 論文作者（phasing-assist 章節 / Methods caveat 撰寫材料）
build_branch: docs/limitations
note: >
  本檔是「另寫的口徑修正草稿」，**不修改任何既有 committed 結論檔**
  （VERIFIED_RESULTS.md V11b 為 SoT，本檔只把其 T2 結論精煉成可貼進論文的措辭並標明引用）。
data_sources:
  - docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md
  - docs/methodology/20260622_ism_method_cognition_and_open_questions_01.md
related_memory:
  - project_methyl_phasing_assist_line
provenance_note: >
  所有數字（AUC 0.90 / null 0.732 / 15-18% / 90.4% / 5.5%）逐一引自
  VERIFIED_RESULTS.md（V11/V11b 區段，行內標來源行號），本 session 不重算。
---

# 口徑修正 — T2「甲基歸 H3」OVERSTATED（draft）

> **這段要解決的問題**：早期 V11 初版把 T2 寫成「甲基可把 H3 read 歸到 H1-1/H2-1，AUC 0.90 → 可行」。經 4-agent 對抗審查（workflow wf_b5b391cf-ea4，所有反駁獨立從 per-chr JSON 重現）修正為 **OVERSTATED / NOT-SUPPORTED**。本草稿把修正後的精確口徑整理成可貼進論文/答辯的版本，並標清「已證 vs 未證」的界線，避免論文沿用舊措辭。

---

## P2.1 一句話口徑（headline，論文用此版本）

> 甲基**只證明**了：在**有 germline 真值的 1-1 vs 2-1 read** 上，甲基可區分其所屬的 germline 分支（LOO AUC 0.90，null 0.732）。**未證明**的是：把這個能力外推到 **H3**（依定義無 germline 真值的 read），亦即「用甲基把 H3 歸到 H1-1/H2-1」**尚未被驗證**。

---

## P2.2 已證 vs 未證（精確分界 — 全部引自 VERIFIED_RESULTS.md V11b）

| | 內容 | 數字（來源）| 狀態 |
|--|------|------------|------|
| ✅ **已證** | 甲基在「有 germline 真值的 1-1 vs 2-1 read」上可分其 germline 分支 | LOO AUC **0.90** / null **0.732** / 73% 位點顯著（VERIFIED_RESULTS.md V11 表 + V11b L148）| SUPPORTED |
| ❌ **未證（V11 誤推）** | 把 0.90 外推到 H3（用甲基歸 H3）| H3 **無 AUC、無 null**，rigor_t2 對 H3 **只算 margin**（grep 確認，V11b L149）| NOT-SUPPORTED |
| ❌ **可指派性低** | H3 能與質心共享 ≥5 CpG 被指派者 | 僅 **15–18%**（chr22 異常 54%），>80% 連指派都不行（V11b L149）| 適用性受限 |
| ⚠ **「56% 高信心」是假象** | 該比例是 margin>0.1，**無 null 對照** | 幾何上 ambient read 也產生 margin（V11b L149）| 不可引用為證據 |
| ℹ️ **H3 組成** | no_germline(case a) **90.4%** / inconsistent(case b) **5.5%**（建議 unTag）| VERIFIED_RESULTS.md L123 | 背景事實 |

---

## P2.3 為什麼會 OVERSTATED（機制 + 結構根因，一段話）

機制上，甲基訊號是 **germline-haplotype 層級**（V10 已證 H1 vs H2 甲基真實且非 copy）。因此甲基「分不同 haplotype」強（T1 救 unphase、T2 分 1-1 vs 2-1 germline 分支都成立），但**「歸 H3」需要的是把無 germline 真值的 read 指派到正確分支**——這一步**沒有 ground truth 可驗**（H3 依定義 germline 未知），所以只能量 margin、不能量 AUC/null。把「有真值子集上的 0.90」當成「全 H3 母體可行」就是把**在有真值 read 上量到的強數字外推到無真值母體**（V11b cross-cutting caveat #2）。

---

## P2.4 論文撰寫指引（DO / DON'T）

**DO（可寫）**
- 「甲基可區分**有 germline 真值的** somatic read 屬於 H1 還是 H2 分支（1-1 vs 2-1，AUC 0.90）。」
- 把此結果定位為「甲基攜帶 germline-branch 資訊」的證據（characterize/corroborate）。
- 明寫 tier 上限 **⭐2-3**（單樣本 HCC1395 + 單一 longphase-S pipeline，被矯正的 HP tag 本身是 longphase-S 產物 → T2/T3 self-reference）。

**DON'T（禁寫）**
- 🔴 禁寫「甲基**可把 H3 歸**到 H1-1/H2-1」「H3 reassignment 可行」「56% H3 高信心可指派」。
- 🔴 禁把 0.90 描述成「H3 reassignment 的準確率」。
- 🔴 禁省略「需 H3 外部真值（long-range linkage / trio）才能成立」這條前提。

> **成立條件（若未來要 reopen T2-as-H3）**：需 H3 的**外部真值**（long-range linkage 或 trio phasing），或 margin 的 null 對照。在此之前，「歸 H3」維持 **未證**。

---

## P2.5 與全論文紅線一致

- 甲基 = **characterize / corroborate**（這裡是「佐證 read 帶 germline-branch 資訊」），**不是** driver / reconstruct / filter。
- 此線全程**單樣本 HCC1395 + 單 longphase-S pipeline（self-reference）** → tier 硬上限 ⭐2-3；跨樣本 COLO829 未做。
- 同檔 T3（within-haplotype 亞群 vs 母本）的「可用性」維持 NEGATIVE，**勿與 T2 混談**（T3 存在性有窄翻案 V11c，可用性仍 NEGATIVE；見 VERIFIED_RESULTS.md V11c）。

---

## 待填 / 缺資訊（投稿前補）

- `{{待填}}` 若論文要量化「T2 可寫範圍」，需引用 VERIFIED_RESULTS.md 的 1-1 vs 2-1 位點數（436 位點，L172）與顯著比例——使用前再 Read 回確認。
- COLO829 跨樣本驗證未做 → tier 不可升過 ⭐3。
