<!--
建立時間: 2026-05-01 00:00
目標: 說明 v5_purity_compare_with_paired IGV session 的 phase、tag 與 TP/FP 多層級審查用途
處理範圍: IGV session、phase-context VCF、TP/FP VCF、audit marker BED 與 site-level manifest
關聯檔案:
  - docs/reports/pi_reports/2026/04/figures/igv_v5_audit/sessions/v5_purity_compare_with_paired.xml
  - docs/reports/pi_reports/2026/04/figures/igv_v5_audit/sessions/annotations/site_layer_manifest.tsv
-->

# IGV V5 Purity Compare With Paired Session

## 目的

此 session 用於人工確認 HCC1395 LongPhase-TO baseline / V5 在不同 purity 條件下的 haplotag 結果，並與 paired LongPhase-S ground truth 對照。新增的 track 目標是同時檢查三件事：

1. phase graph / phase VCF 在特殊位點附近是否有 anchor 與 `GT|PS` 資訊。
2. BAM read 的 `HP` tag 是否與 paired truth、baseline、V5 一致。
3. 每個特殊位點在 LongPhase-S / LongPhase-TO 的 TP、FP 層級判斷是否一致。

## Track 分層

| Panel | Track | 用途 |
|---|---|---|
| `VariantPanel` | ClairS-TO@0.93 / ClairS-TO@0.6 | 原始 tumor-only candidate VCF |
| `PhaseContextPanel` | LP-S normal phase context | paired normal phase anchor 參考 |
| `PhaseContextPanel` | LP-TO BL/V5 phase context | tumor-only baseline / V5 phase anchor 對照 |
| `TPFPVariantPanel` | LP-S paired TP/FP VCF | paired truth 層級的 TP/FP 參考 |
| `TPFPVariantPanel` | LP-TO PASS TP/FP VCF | tumor-only 層級的 TP/FP 結果 |
| `AuditSiteMarkerPanel` | all colored / TP / FP / V5max / self-phasing markers | 特殊位點定位與分類顏色 |
| `BedPanel` | LOH / GE BED | phase block 與 genome event 背景 |
| BAM panels | paired、baseline、V5、purity 0.6 BAM | 用 `HP` tag 檢查 read-level tag 結果 |

## 顏色規則

| 分類 | 顏色 | RGB |
|---|---|---|
| TP | 藍色 | `37,99,235` |
| FP | 紅色 | `220,38,38` |
| V5max | 紫色 | `124,58,237` |
| Self-phasing | 橘色 | `234,88,12` |

## 判讀順序

1. 先看 `Audit markers: all colored sites` 定位特殊位點與分類。
2. 看 `TPFPVariantPanel`：確認該位點是否出現在 LP-S paired TP/FP 或 LP-TO TP/FP VCF。
3. 看 `PhaseContextPanel`：確認特殊位點附近是否有 phased anchor、`GT` 是否為 phased genotype、`PS` 是否落在合理 phase set。
4. 看 BAM panels：以 `colorByTag=HP` 與 `groupByOption=PHASE` 檢查 paired、baseline、V5 的 read-level haplotag 是否一致。
5. 交叉看 `BedPanel` 的 LOH/GE track，避免把 LOH 或 genome event 背景誤判為 tagging 改善。

## 重要限制

- `LP-TO V5@0.93 phase context` 使用 `pononly_v2b/tumor_phased.vcf` 作為 phase backbone；目前 V5 somatic fallback run 主要產出 tagged BAM，未看到獨立 V5 phase VCF。
- `site_layer_manifest.tsv` 的 `*_phase_n1kb` 是特殊位點正負 1 kb 內 phase-context VCF record 數，代表可用於人工審查的 phase anchor 密度，不等於該位點一定被當作 anchor。
- V5max 與 self-phasing marker 主要是 tag 行為觀察點，部分 marker 沒有 ref/alt 欄位，因此不會在 TP/FP VCF exact-match 欄位標成 `Y`。

## 產物

| 檔案 | 說明 |
|---|---|
| `v5_purity_compare_with_paired.xml` | IGV session 主檔 |
| `annotations/audit_sites_all_colored.bed` | 全部特殊位點合併 colored BED |
| `annotations/audit_sites_TP.bed` | TP 特殊位點 |
| `annotations/audit_sites_FP.bed` | FP 特殊位點 |
| `annotations/audit_sites_V5max.bed` | V5 改變最大位點 |
| `annotations/audit_sites_SelfPhasing.bed` | self-phasing 觀察位點 |
| `annotations/site_layer_manifest.tsv` | 每個特殊位點的 TP/FP VCF 命中與 phase-context anchor 摘要 |
