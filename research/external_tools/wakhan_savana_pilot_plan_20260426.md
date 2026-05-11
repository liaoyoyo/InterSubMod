# External CN/SV Pilot 設計 — Wakhan + SAVANA on HCC1395

**日期**：2026-04-26
**狀態**：plan only — pending user approval, not auto-execute
**Scope**：HCC1395 paired_full smoke only（不擴展其他樣本）

## 1. 動機

- Thread D LOH-constrained phasing 已成主軸；ISM 內部用 CovM proxy 推 CN tier
- 需外部 CN caller 驗證 CovM tier 是否與真實 CN call 對齊
- SV caller 用於驗證 phase block boundary 與 SV breakpoint 是否協調
- HCC1395 為 SEQC2 reference 樣本，truth set 較完整（COSMIC + SEQC2 high-confidence consensus）
- 若外部工具與 ISM 內部 proxy 一致 → 強化 Thread D 證據鏈；若不一致 → 暴露 R-EXT-CN 風險，主軸需加 caveat

## 2. 工具選擇與版本

### 2.1 Wakhan（CN caller for long-read tumor）

- **Repo**：`https://github.com/KolmogorovLab/Wakhan`
- **版本（查詢日期 2026-04-26）**：v0.4.2（bioconda 最新；最後更新 2026-02-19）
- **Publication**：medRxiv preprint posted 2025-12-15（"Wakhan: reconstruction of chromosome-scale copy number profiles of tumor genomes with long-read sequencing"）
- **Install 估時**：~30 min（建議 bioconda channel 安裝；含 Severus 依賴）
- **Input**：tumor BAM + 參考基因組 + phased germline VCF + somatic SV VCF（可由 SAVANA 或 Severus 提供）
- **Output**：CN segments BED + per-segment haplotype-specific CN level + tumor purity/ploidy estimate
- **核心特性**：haplotype-specific somatic CN profiling；phase switch 校正；breakpoint-driven segmentation

### 2.2 SAVANA（SV caller for long-read tumor）

- **Repo**：`https://github.com/cortes-ciriano-lab/savana`
- **版本（查詢日期 2026-04-26）**：repo 最後更新 2026-03-03（請於啟動時確認 `git tag` 取最新 stable release）
- **Publication**：Nature Methods 2026（"SAVANA: reliable analysis of somatic structural variants and copy number aberrations using long-read sequencing"）
- **Install 估時**：~30 min（conda 或 source；含 minimap2/winnowmap 依賴）
- **Input**：tumor BAM + matched normal BAM + 參考基因組
- **Output**：SV BEDPE + VCF + SCNA segments + tumor purity/ploidy + allele-specific absolute CN
- **核心特性**：single-haplotype resolution SV/SCNA；經 99 tumor-normal pair 驗證；sensitivity 與 specificity 顯著優於 Severus / Sniffles2 等

## 3. Success Criteria（pre-registered）

寫死在啟動前，不事後調整：

### A. Wakhan vs ISM CovM 對齊

- **A1**：Wakhan CN segment level vs ISM CovM tier Spearman correlation ≥ 0.7（per-region 配對）
- **A2**：80% 以上 ISM region 落在一致 CN class（diploid / loss / gain 三類分類一致率 ≥ 0.80）
- **A3**：tumor purity 估值差異 vs ISM 內部估算 ≤ 10 個百分點

### B. SAVANA vs ISM phase block 協調

- **B1**：SAVANA SV breakpoint 與 ISM phase block boundary 距離 ≤ 500 bp 的比例 ≥ 40%
- **B2**：ISM phase block 跨越 SAVANA SV breakpoint 的比例 ≤ 10%（block 不應跨真實 SV）
- **B3**：LOH 區域內 SV breakpoint 密度 vs 非 LOH 區域 enrichment ratio ≥ 2.0（驗證 LOH-SV 關聯）

**通過條件**：A1+A2+B1+B2 全達標 → Thread D §3 證據鏈強化；任一失敗 → 標記為 caveat 並開 R-EXT-CN risk 卡

## 4. Pilot 範圍

- **樣本**：HCC1395 only（SEQC2 reference）
- **資料**：paired_full mode（已 KDE-corrected；對應 master dataset）
- **Tumor BAM 路徑**：待啟動時由用戶確認（big8 archive 或 big7 working copy）
- **Normal BAM 路徑**：待啟動時由用戶確認
- **參考基因組**：GRCh38（與 ISM pipeline 一致）
- **計算量**：~6-8 hr（SAVANA ~4 hr 含 normal pair；Wakhan ~3 hr 依賴 SAVANA SV VCF）
- **儲存**：~10 GB output（含中間 phased VCF、SV BEDPE、CN BED、purity/ploidy report）
- **Cores**：建議 16-32 threads（SAVANA OpenMP）

## 5. 不執行原因

- 用戶已明示 P2-2「不立即啟動」（2026-04-26）
- 等 Thread D 主軸報告（P0-1）完成後再決定優先序
- HCC1954 case panel（P0-3）優先級高於外部工具
- 啟動前需確認 SEQC2 truth set 路徑與比對協議（避免重複既有 truth comparison）

## 6. 啟動 checklist（未來執行時）

- [ ] 用戶批准（Hard Gate — 不可逆計算 6-8 hr）
- [ ] Wakhan + SAVANA install + smoke test（小 region 確認可跑）
- [ ] HCC1395 tumor + normal BAM 路徑確認（big7 / big8）
- [ ] 計算資源配置（disk 預留 15 GB headroom；cores 16-32）
- [ ] 預設 success criteria（§3）寫死，不事後調整
- [ ] Pre-register pilot output 路徑：`output/synthesis/external_tools/wakhan_savana_HCC1395_<DATE>/`
- [ ] 設定 timeout（每階段 wall-clock ≤ 5 hr，超時暫停回報）

## 7. 與 Thread D 的關係

- **若 §3.A 通過**：Wakhan 確認 CN tier 對齊 ISM CovM → 強化 Thread D §3 CovM-as-CN-proxy 推論
- **若 §3.B 通過**：SAVANA 確認 phase block 與 SV breakpoint 協調 → 補強 LOH-constrained phasing 機制（block boundary 真的有結構基礎）
- **若 §3.A 失敗**：CovM tier ≠ true CN → R-EXT-CN 風險，Thread D §3 需加 caveat「CovM 為 read-density proxy，未經外部 CN caller 直接對齊」
- **若 §3.B 失敗**：phase block 與 SV 不協調 → 機制需重審（可能 phase block 是 noise / local artifact）

## 8. 連結

- Thread D 主軸報告（佔位）：`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`
- AutoResearch direction（佔位）：`InterSubMod/research/autoresearch/research_direction.md`
- Wakhan upstream：`https://github.com/KolmogorovLab/Wakhan`
- SAVANA upstream：`https://github.com/cortes-ciriano-lab/savana`
- SEQC2 truth set 協議：`InterSubMod/docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md`
- LOH-constrained phasing 記憶：`project_loh_constrained_phasing_discovery.md`
