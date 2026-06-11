<!--
建立時間: 2026-06-11
狀態: reconcile note (session 間 reconcile — 不是新 hub；指向 SoT + 交接本 session 新資產)
報告類型: paper_focus_session_reconcile_note
受眾: 廖子游 · 並行/後續 session
provenance_note: 不產新數字；指向既有 SoT + 本 session ledger 95–99。
-->
<!-- provenance-verified: 指向既有 foundation/整合篇章 SoT + 本 session ledger 95-99；無新數字。 -->

# 主軸升級 reconcile note — Subclonal Reconstruction（多 session 收斂）

> **⚠ 本檔不是 hub，是 reconcile note。** 2026-06-11 偵測到**並行 session 已建立**兩份 SoT，本檔指向它們並交接本 session 唯一未收錄的新資產。**勿在此檔展開內容**（避免第三個競爭 hub）。

## SoT 指向（新主軸的權威文件 — 已存在，勿重建）

| 文件 | 面向 | 職責 |
|------|------|------|
| `InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md` | **A 面：甲基-phasing-assist**（V1-V12）| 新主軸 **foundation SoT**；6 樣本資產 + G-A~G-E gap |
| `InterSubMod/docs/concepts/2026/06/20260611_Subclonal_Reconstruction_Paper_Focus_整合篇章_01.md` | **B 面：ASM-characterization + NEGATIVE**（本 session 這條線）| 兩線 reconcile + 標題拆解 + 章節映射 + critical path |

> 兩份已**充分涵蓋**本 session 多數工作（catalog / chr17 / BRCA2 退役 / 四道 NEGATIVE / 6/6 excess）。**新主軸的 hub = 上面兩份，不是本檔。**

## 🔑 本 session 唯一未被收錄的新資產（交接給 SoT 整合）

**clusterability vs coverage/CN（2026-06-10，ledger 97/98）** —— 兩份 SoT **皆未含**，但這是支撐主軸最直接的新 mechanistic 證據：

> **甲基 read-clusterability 非 coverage/CN-dosage 假象；LOH 抑制分群 6/6 樣本（OR 0.01–0.107, 3 癌種）= 甲基 read-clustering 需 ≥2 haplotype 的等位甲基差** → 直接證明「**甲基 profile 隨 germline-haplotype/subclone 結構化**」（= 主軸 methylation pillar 的機制地基）。真 driver=n_CpG（非 coverage）。
> 報告 `InterSubMod/docs/experiments/in_progress/2026/06/20260610_methylation_clusterability_vs_coverage_CN_characterization_01.md`；已併 catalog R6 §2.5 + G6 backbone §N2b。

**建議**：把此資產加進**整合篇章 B 面**的 Methods 機制段（「甲基 = haplotype-structured 非技術假象」），它與 foundation §4「甲基=germline-haplotype 層級」**機制上互證**（LOH 抑制分群 ⇔ 需兩個 haplotype ⇔ 甲基是 haplotype-level）。

## 收斂確認（兩 session 各自獨立到同結論 = 強驗證）

- **同一關鍵 gap**：foundation **G-D「真正重建 demo」** = 本 session 盤點的「characterization vs reconstruction」gap = **唯一 paper-level blocker**。
- **同一開放問題**：整合篇章 **G-B「subclone 甲基 = somatic-specific 還是 germline-allelic baseline」**（within-haplotype null 未跑）= 甲基-subclone 貢獻強度的關鍵閘。
- **同一誠實天花板**：single-pipeline ⭐3 / HD-1 phasing circularity / 甲基非 variant-filter。

## ⚠ 並行 session 警示（governance）
3 個 session 同時動同主軸文件 + 共用 git HEAD → 有 concurrent-edit 互撞風險（CURRENT_FOCUS / foundation / 整合篇章本 session **未碰**，留給 owner session）。**建議指定一個 session/文件為 owner**，其餘交接資產後停手。
