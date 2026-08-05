<!--
建立時間: 2026-07-22 Asia/Taipei
目標: 在進行 HCC1395 跨 PS 的 HP 方向敏感度驗證前，明列假設、反例、可證偽條件與放行標準
處理範圍: HCC1395 chr1-chr22 全量 region census；BAM 層只做預先定義的 bounded region probes，未修改 production pipeline
關聯檔案: InterSubMod/research/20260722_hcc1395_ps_block_hp_orientation_audit/scripts/audit_hcc1395_ps_regions.py; InterSubMod/research/20260722_hcc1395_ps_block_hp_orientation_audit/scripts/probe_hcc1395_ps_orientation.py
-->

# HCC1395 跨 PS 與 HP 方向問題 pre-decision audit

用 Assumption Map + Falsifier：**初始判定為 PROBE；若不同 PS 的相對 HP flip 會改變 T、Topo 或 HP-family unit，則 current mixed-PS topology 不得稱為已矯正。**

## 1. 任務分類與研究目標

- HCC1395 `chr1-chr22` mixed-PS census：Task Type B（全 HCC1395 retained regions）。
- BAM/sidecar orientation sensitivity：Task Type A（bounded probe；不能估計 865 區的受影響率）。
- 服務 G1、G2、G4：改善 LongPhase-S／ISM 整合的 phase 語意，並提高跨樣本可重現性。

研究啟動五問：

1. Thread D：否；本題是 genetic read-state／PS block。
2. Thread B 撤回範圍：否。
3. KDE-corrected：不適用。
4. VCF caller AF：不需要；主要證據為 HP／PS／read pattern。
5. 長計算／C++／檔案搬移：否；只做 Python 唯讀 census 與狹窄 BAM region probe。

## 2. 決策問題

> HCC1395 的距離型 region 是否包含多個 PS；若包含，現行只按 HP 第一碼分組的做法，是否依賴不同 PS block 任意的 HP1／HP2 方向？

## 3. 關鍵假設

| ID | 假設 | 風險 | 驗證方式 |
|---|---|---|---|
| A1 | HP1／HP2 只在同一 PS block 內有一致方向 | 高 | KB／phase-set contract；不同 PS 不可直接比較 |
| A2 | HCC1395 sidecar 的 HP／PS 有精確接回原 BAM alignment | 高 | exact identity join；missing/conflict 必須為 0 |
| A3 | 全部 PS 一起翻轉只是全域命名交換 | 中 | 固定一個 anchor，只枚舉其餘 `2^(b-1)` 相對方向 |
| A4 | 無 PS read 無法安全做 block flip | 高 | 另報 missing-PS；probe 中保持原標籤，不宣稱已定向 |
| A5 | mixed-PS 不等於拓撲必錯 | 高 | 每區重跑相對方向；只以結果不變性判斷 |

## 4. 反方解釋與 falsifier

- 反方 1：LongPhase-S 可能已跨 block 提供全域 HP orientation。反駁門檻不是名稱，而是需要外部 anchor 或跨 block signed bridge evidence；僅有不同 PS 的同數字 HP 不足。
- 反方 2：mixed-PS 可能只是 1–2 條零星 read。驗證次要 PS read 數與比例，不能只數 PS label。
- 反方 3：即使方向未知，樹可能對所有 flip 都不變。這是主要 falsifier。

**可證偽條件**：若所有受測 mixed-PS region 在 `2^(b-1)` 個相對方向下，HP-family unit、T、Topo 與 morphology 全不變，則此問題只需保留 QC 警示，不構成實質拓撲風險。

## 5. Step → Verify

1. 程式 contract audit → 驗證：指出 region 分組、HP pooling、PS census、solver input 的 frozen source 行號。
2. HCC1395 全量 census → 驗證：8,222 個 MLHP/current-v5 region key 一對一，PS 三類守恆。
3. Raw read reproduction → 驗證：native read patterns、T、Topo 與 canonical artifact 完全相符。
4. 相對 PS flip → 驗證：固定第一 PS，交換第二 PS 的 HP1／HP2，重跑同一 frozen solver。
5. 判讀 → 驗證：區分「風險存在」「特定 region 敏感」「865 區盛行率尚未估計」。

## 6. 初始評分與 verdict

| 維度 | 分數（20） | 理由 |
|---|---:|---|
| 理論基礎 | 20 | PS block 與相對 orientation 定義明確 |
| 既有觀察 | 15 | 已有 mixed-PS aggregate，但尚無 HCC flip 實測 |
| 機制可測 | 20 | 可對每個 PS block 做明確 HP swap 並重跑 solver |
| 反例控制 | 15 | 可用 CN-neutral、k=2、balanced PS 作乾淨反例 |
| 資源可行 | 20 | `b<=4`，相對方向空間小 |
| **總分** | **90/100** | **PROBE，禁止先宣稱全部 865 區錯誤** |

## 7. Probe 後決策更新

- HCC1395 full census：865/8,222 mixed-PS。
- CN-neutral 反例：相對 flip 使 `T=1→3`、`Topo=1→2`。
- 高對比 panel：12/12 native 重現 canonical；12/12 HP-family unit signature 改變。

因此：

- **GO**：實作 PS-aware unit 與 orientation-invariance gate。
- **NO-GO**：把 current mixed-PS topology 直接稱為「已矯正／confirmed」。
- **尚未判定**：865 區中實際有多少在所有合理方向下改變；需全量 1,802 configurations 才能回答。
