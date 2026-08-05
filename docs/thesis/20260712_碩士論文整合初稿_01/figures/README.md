<!--
建立時間: 2026-07-12
目標: 記錄碩士論文整合初稿圖組的用途、資料來源與更新 gate
處理範圍: 本資料夾 6 張原創 SVG；不含舊 Meeting 或 historical layered 圖
關聯檔案: ../20260712_碩士論文_01.md；../20260712_碩士論文_01.html
-->

# 圖組索引與來源

本資料夾的 SVG 均為本次論文整合初稿重新繪製，沒有直接複製舊簡報圖片。Markdown 與 HTML 皆以相對路徑引用。

| 圖 | 用途 | 數據／概念來源 | 更新條件 |
|---|---|---|---|
| `01_research_narrative_map.svg` | 研究主軸與 claim ceiling | `InterSubMod/external_validation/_schema/paper_claims.tsv`；`InterSubMod/docs/CURRENT_FOCUS.md` | claim registry 變更時 |
| `02_layered_evidence_pipeline.svg` | 方法與證據分層 | `InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md`；LongPhase-S KB | input contract 或 HP vocabulary 變更時 |
| `03_claim_boundary.svg` | 工程、技術與生物證據邊界 | `paper_claims.tsv`；solver red-team audits；2026-07-12 raw-VAF/PyClone-VI 報告 | 新正交 truth 或 solver reconciliation 完成時 |
| `04_current_evidence_snapshot.svg` | 2026-07-12 數據快照 | raw-all closeout；legacy solver artifact + red-team audit；raw-VAF/PyClone-VI；clean-v3 terminal receipts | clean-v3／solver stress 修復重跑後必更新 |
| `05_thesis_completion_gates.svg` | 論文收斂與替換計畫 | pre-decision audit + current runtime state | C2/C3 completion gates 完成時 |
| `06_worked_candidate_tree_example.svg` | read matrix → partial subcube → minimum candidates 的可手算示例 | solver formal model；synthetic didactic case | 演算法 contract 改動時 |

> **狀態**：`04_current_evidence_snapshot.svg` 是有日期的 partial snapshot；它記錄 workers 7/7 完成但 verifier terminal FAILED，不得當成 clean layered-v3 終版結果。
