<!--
建立時間: 2026-07-11 Asia/Taipei
目標: 驗證 LongPhase-S _sc.vcf FILTER=PASS 作為 sSNV tree input 的四層位點守恆與 exact HP/PS join
處理範圍: PARTIAL / HCC1954 chr22 bounded real-data probe；不可作 7-dataset comprehensive evidence
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/20260711_ClairS_LongPhaseS_sSNV位點與ReadTag契約稽核_01.md
-->

# HCC1954 chr22 LongPhase-S PASS Tree Contract Probe

## 結論

此 bounded probe 通過四層 record-key 與 sSNV/read-tag 守恆：raw ClairS 426 records、ClairS PASS/LongPhase-S input 181、LongPhase-S `_sc.vcf` all 181、`_sc.vcf` PASS/tree input 167。14 個 LongPhase-S LowQual records 保留在 site ledger，但不進建樹。

tree input 的 167 個 sSNV 完整分解為 positional singleton 99 + retained 68；31 個 retained groups 共暴露 4,017 alignment-group observations，sidecar exact matches 4,017、missing 0、conflicts 0。

## 輸入

- raw ClairS：`/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/output.vcf.gz`
- ClairS PASS / LongPhase-S input：`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz`
- LongPhase-S no-truth chr22 output（原始 bounded execution）：`/tmp/lps_fifo_chr22_HCC1954/HCC1954_clean_sc.vcf`
- original tumor BAM：`/big8_disk/data/HCC1954/ONT/HCC1954.bam`
- exact sidecar（原始 bounded execution）：`/tmp/lps_fifo_chr22_HCC1954/HCC1954.read_tags.digest.tsv.gz`

原始 `/tmp` 工作區不是持久化證據；本目錄已保存 raw/input/`_sc`-all/`_sc`-PASS 四層 indexed VCF、`mlhp_part_chr22.json` 與 sidecar validation。下列 persistent replay 完全由本目錄 artifact 重建，不依賴 `/tmp` 路徑。

## 執行

1. `bcftools view -r chr22` 建立 raw 與 LongPhase input bounded VCF。
2. `bcftools view -f PASS` 從 `_sc.vcf` 建立 canonical tree VCF。
3. `sm_multilocus_combinations.py chr22` 使用 original BAM + exact sidecar 執行 sSNV census。
4. `build_ssnv_site_ledger.py` 驗證 raw PASS=LPS input、LPS input=`_sc` all、`_sc` PASS=tree input 與完整 disposition。

## 輸出與驗證

- `ssnv_site_ledger_HCC1954_chr22.summary.json`：六項 checks 全為 `true`。
- `ssnv_site_ledger_HCC1954_chr22.tsv.gz`：426/426 raw records，各有唯一 disposition；raw ClairS 與 LongPhase-S output 的 ID/QUAL/FILTER/INFO/FORMAT/sample values 並列保存。
- `ssnv_site_ledger_HCC1954_chr22.tsv.gz.tbi`：可依 genomic interval 直接查詢完整位點結果。
- `ssnv_site_ledger_HCC1954_chr22.persistent.tsv.gz` 與 `.tbi`：由本目錄 artifact 重新產生的現行 schema ledger，另含 per-site `phase_set_counts_json` 與 `phase_set_HP_counts_json`。
- `ssnv_site_ledger_HCC1954_chr22.persistent.summary.json`：輸入與輸出路徑均指向本目錄，七項 checks 全為 `true`，可作持久化 readback。
- `mlhp_part_chr22.json`：scope conservation、exact HP/PS join 全通過。
- `stream_sidecar_validation.json`：producer/capture count、HP vocabulary、record keys 全通過。
- `artifacts.sha256`：本目錄小型證據檔完整性清單。

## 限制

這是 HCC1954 chr22 partial probe，只證明資料流與 producer-to-consumer tag fidelity；不證明 LongPhase-S phase assignment 的生物真實性，也不取代進行中的 7-dataset genome-wide production validation。
