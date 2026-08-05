<!--
建立時間: 2026-07-10 23:40 +08:00
目標: 導航分層重建 v2 端到端程式碼、方法、設定、判斷與驗證稽核產物
處理範圍: active layered reconstruction v2；含 ClairS/LongPhase-S 上游、Python 重建主鏈、InterSubMod C++ 平行支線與報告層
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md
task_type: B Comprehensive validation
framework: Diátaxis + A3/ADR
build_branch: research/subclonal-reconstruction-202606
build_commit: 4fb9e742482b63a660de19a1f1bd07d49d713111
worktree: /big7_disk/liaoyoyo2001/InterSubMod (dirty; audited files include modified/untracked worktree content)
驗證方式: source line audit + 7/7 input preflight + synthetic/golden tests + C++ build/ctest + HTML browser QA
證據等級: 程式語意 L1；輸入可讀性 L1；全量 v2 動態輸出未執行，因此端到端結果僅 L4 pending
-->

# 分層重建 v2：端到端概念流程、當前執行邊界與判斷稽核

> **⚠ 動態驗證範圍聲明**：本目錄完成全主鏈靜態稽核、7/7 輸入 preflight 與輕量測試。動態快照鎖定於 `2026-07-11 00:19:04 +0800`、host `bip7`：immutable run root `20260710_232501_layered_reconstruction_v2` 只記錄 HCC1395／HCC1395_DORADO 的 preflight PASS 與 mlhp START，沒有 sample complete 或 `verification_summary.json`；該host亦未見相符 pipeline process。因此不可把此未完成 run、舊產物或單元測試寫成「v2 全流程已驗證」。

## 主要產物

- `InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/20260710_分層重建端到端流程程式碼判斷稽核_01.md`：文字版 source of truth。
- `InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/20260710_分層重建端到端流程程式碼判斷稽核_01.standalone.html`：獨立 HTML 說明書。
- `InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/input_preflight.json`：7 datasets 輸入路徑、header、index、VCF/CN 基本契約與 checksum。
- `InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/source_inventory.tsv`：主鏈程式碼快照清單與 SHA-256。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/`：三條獨立稽核支線的原始證據筆記。

## 最終狀態語意

| 對象 | 狀態 | 可宣稱範圍 |
|---|---|---|
| 本次程式碼／方法稽核報告 | `reviewed` | 可依行號重查判斷、設定、風險與驗收條件 |
| 7-dataset 輸入 preflight | `7/7 PASS` | 路徑、index、基本 VCF/CN/header 契約可讀；不含 BAM 全內容 checksum 或 HP scope 一致性 |
| v2 合成／golden 測試 | `PASS with warnings` | 小型人工案例符合目前實作；不代表真實資料全量行為 |
| 文件 fresh-reader closure | `99/100; no material ambiguity` | 能獨立重述 current/target、8 P0、Gate-0–4 與 claim ceiling；不是 pipeline validation |
| v2 7-dataset full run | `INCOMPLETE / FINAL NOT FOUND` | run root 已建立但停在 2 datasets 的 mlhp START；不可宣稱 7/7 output、funnel、V1–V7 或跨樣本結果已由本版程式驗證 |
| 生物學 clone confirmation | `OUT OF CLAIM SCOPE` | 目前只能稱 regional mutation-state tree；缺 single-cell／multi-region orthogonal truth |
