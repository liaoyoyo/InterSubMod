<!--
建立時間: 2026-05-12 01:30
目標: 標明 sessions/ 已遷到 research/igv_sessions/；本目錄只保留 PNG snapshots
-->

# IGV V5 Audit Figures（snapshot 產出目錄）

## 結構

本目錄僅保留 IGV snapshot **PNG 產出檔**：

- `by_HP_4ver/` — V5/V3F/V2b/baseline 四版並列 audit snapshots（SP1/SP2/SP3）
- `by_HP_v6/` — V6 audit snapshots（fast-batch + audit-quality 兩種）
- `v4ver_audit_summary.tsv` / `v5_audit_summary.tsv` — audit 數據摘要

## Sessions 已遷移

⚠ **IGV session XML + annotations + README** 已於 2026-05-12 遷到：
**`InterSubMod/research/igv_sessions/`**

理由：IGV session XML 為跨報告 reusable 研究級配置，與 `research/*` 並列管理更合理；本目錄專注於 audit snapshot PNG 輸出。

## 對應關係

| 此目錄（PNG 輸出） | session XML 來源 |
|---|---|
| `by_HP_4ver/D_SP{1,2,3}_*.png` | `InterSubMod/research/igv_sessions/v5_all_4versions.xml` |
| `by_HP_v6/D_SP{1,2,3}_*_v6.png`（fast-batch） | `by_HP_v6/igv_batch_v6.txt`（自包含 4-BAM） |
| `by_HP_v6/D_SP{1,2,3}_*_v6_audit.png`（audit-quality） | `InterSubMod/research/igv_sessions/v6_germline_absent_audit.xml` |

詳見：`InterSubMod/research/igv_sessions/README.md`
