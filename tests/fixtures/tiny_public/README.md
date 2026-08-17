<!--
建立時間：2026-08-13
目標：提供可公開、可重生的 InterSubMod tiny synthetic E2E fixture
處理範圍：1 synthetic contig、1 SNV、12 synthetic reads；DEMO only
關聯檔案：InterSubMod/scripts/handoff/build_tiny_public_fixture.sh；InterSubMod/scripts/handoff/run_tiny_public_e2e.sh
-->

# Tiny public synthetic fixture — DEMO ONLY

> **PARTIAL / DEMO：此 fixture 只驗證軟體執行、I/O contract 與 199-column schema；不可作科學 validation evidence。**

- `reference.fa`：200 bp synthetic contig `chrTiny`。
- `variants.vcf`：單一 `chrTiny:101 A>T` synthetic SNV。
- `tumor.sam`：12 條 synthetic reads；6 條 HP1/REF、6 條 HP2/ALT，全部具有 Dorado-style `MM:Z:C+m?` 與 `ML:B:C` tags。
- BAM、BAI 與 FAI 由 `InterSubMod/scripts/handoff/build_tiny_public_fixture.sh` 在隔離工作目錄重生，Regular Git 不保存 binary。

`clustering/tree.nwk` 是由 read–read methylation distance 產生的 **read dendrogram**。它不是 cellular lineage、clonal ancestry、phylogeny truth、subclone prevalence 或生物正確率的證據。

本 fixture 服務 G4（reproducibility）與 G5（外部可驗證交接），不驗證 G1/G3 的生物學 claim。
