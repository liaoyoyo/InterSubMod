<!--
建立時間: 2026-04-27
目標: 紀錄 0.6 purity 模擬實驗計畫（先紀錄，之後實作）
受眾: 研究團隊（之後執行此實驗時參考）
處理範圍: BAM mix simulation + LongPhase-TO baseline/V5 重跑 + HP tag concordance 分析
狀態: planned_pending_execution
priority: P2（補強 audit suite 但非必要）
relates_to:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md
trigger_question: "V5 tag 是否在 0.6 purity 下更偏向 paired pure tag？baseline vs V5 在 0.6 上的影響？"
estimated_runtime: ~3.5-4 hr 機器時 + ~1 hr 分析撰寫
blockers:
  - ClairS-TO 環境未設定 → caller-level F1 不可重跑（用 workaround：重用既有 VCF）
disk_requirement: ~550 GB（/big7_disk 6 TB 可用，邊緣 OK）
mini_phase_options:
  - Mini Phase 1（~1 hr）：BAM mix + 量化檢查
  - Mini Phase 2（+1.5 hr）：LongPhase 兩版本完整跑
  - Mini Phase 3（+45 min）：分析 + 報告
status_note: "計畫已 planned，等用戶觸發實作。當前 audit suite 已含理論推論章節。"
-->

# 0.6 Purity 模擬實驗計畫（pending execution）

> **狀態**：📝 已紀錄，**未實作**。等用戶觸發後執行。
> **觸發後實作位置**：實際輸出將存於 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` 與對應 `data/`、`figures/` 子目錄。
> **計畫副本**：同步存於 `/bip7_disk/liaoyoyo2001/.claude/plans/streamed-spinning-wilkinson.md`（Claude 內部 plan file）。

---


## Context

**觸發**：PI 問「V5 tag 是否在 0.6 purity 下更偏向 paired pure tag？baseline vs V5 在 0.6 上的影響？」

之前理論推論：
- 0.6 purity 下 baseline self-phasing 強度降低（17.3:1 → 預期 3-5:1）
- V5 PON-only flag 不依賴 purity，行為仍一致
- F1 差異仍噪音範圍

**本任務**：實際做 simulation 驗證上述推論，量化 V5 在低 purity 下對 paired ground truth 的 concordance 改善。

**Explore Phase 1 已釐清**（feasibility）：

| 項目 | 狀態 |
|------|:----:|
| HCC1395 tumor BAM | ✅ 264 GB，無 HP tag，RG 確認 |
| HCC1395BL normal BAM | ✅ 136 GB，含 HP tag + PS（同個體 cell line pair）|
| 磁碟空間 /big7_disk | ⚠ 6 TB 可用（需 ~550 GB，9% 使用率，OK 但邊緣）|
| LongPhase-TO baseline binary | ✅ `longphase-to-baseline`（已分離，commit 8b8c1fd 之前）|
| LongPhase-TO V5 binary | ✅ `longphase-to`（含 --pon-only-phasing flag）|
| PON files × 4 | ✅ 全部存在 |
| **ClairS-TO 環境** | ❌ **無**（容器/工具/腳本都沒設定）|
| Existing ClairS-TO VCF | ✅ `/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz` |
| Single-run 預估時間 | ~45-50 min（baseline log 2895s for full HCC1395）|

---

## 範圍切分（Scope Boundary）

### ✅ 在範圍內（可驗證）

- **Phasing-layer 差異**：baseline vs V5 在 0.6 purity 下的 HP tag concordance to paired
- **Self-phasing 強度量化**：HP1:HP2 bias 在 0.93 vs 0.6 purity 下的對比
- **Per-site 改善**：15 既有位點在 0.6 場景下的 V5 vs baseline 表現
- **Imbalance ratio variation**：0.6 vs 0.93 場景下的 ratio distance

### ❌ **不在範圍內**（受 ClairS-TO blocker 限制）

- **Caller-level F1**：不能在 mixed BAM 上重新 call variants
- **Recall 影響**：低 purity 對 caller TP/FP rate 的衝擊（已知理論：Raw F1 ~0.7166 → ~0.55-0.65，但無法實測）

### Workaround：重用既有 VCF

對 mixed BAM 跑 LongPhase-TO 時，**直接用既有 HCC1395 tumor VCF**（snv.vcf.gz）作 input。這意味：
- ✅ 可驗證 phasing 層 V5 vs baseline 差異（核心問題）
- ❌ 不可驗證 calling 層 F1 變化（次要問題）
- ⚠ Caveat：mixed BAM 中 normal reads 對應的 somatic 位點將顯示低 AF（可能影響 phasing graph anchor 信心）

---

## 執行流程（5 階段）

### 階段 A：製造 0.6 Purity Mixed BAM（30-45 min）

**輸入**：
- Tumor: `/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` (264 GB, no HP)
- Normal: `/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam` (136 GB, has HP—需 strip)

**步驟**：
```bash
WORK=/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation
mkdir -p $WORK

# 1. Strip HP/PS tags from normal BAM (避免污染 V5 重 phasing)
samtools view -h NORMAL.bam | sed 's/HP:i:[0-9]*//g; s/PS:i:[0-9]*//g' | samtools view -bS - > $WORK/normal_stripped.bam

# 2. Random-sample 60% from tumor + 40% from normal
samtools view -s 0.6 -b TUMOR.bam > $WORK/tumor_60.bam     # ~158 GB
samtools view -s 0.4 -b $WORK/normal_stripped.bam > $WORK/normal_40.bam  # ~54 GB

# 3. Merge + sort + index
samtools merge -@ 16 $WORK/mixed_06.bam $WORK/tumor_60.bam $WORK/normal_40.bam
samtools sort -@ 16 -o $WORK/mixed_06_sorted.bam $WORK/mixed_06.bam
samtools index $WORK/mixed_06_sorted.bam
rm $WORK/tumor_60.bam $WORK/normal_40.bam $WORK/normal_stripped.bam $WORK/mixed_06.bam  # 清掉中間檔
```

**驗證**：`samtools view -c $WORK/mixed_06_sorted.bam` 約 ~250M reads（粗估）；`.bai` 索引存在

**輸出**：`mixed_06_sorted.bam` (~213 GB)

### 階段 B：跑 LongPhase-TO Baseline on Mixed BAM（45-50 min）

```bash
cd /big7_disk/liaoyoyo2001/longphase-to-mod
./longphase-to-baseline phase \
  -s /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz \
  -b $WORK/mixed_06_sorted.bam \
  -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
  --ont --caller clairs_to_ssrs \
  --pon-file /big7_disk/liaoyoyo2001/data/PON/1000g-pon.sites.vcf.gz,/big7_disk/liaoyoyo2001/data/PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz \
  --strict-pon-file /big7_disk/liaoyoyo2001/data/PON/dbsnp.b138.non-somatic.sites.vcf.gz,/big7_disk/liaoyoyo2001/data/PON/gnomad.r2.1.af-ge-0.001.sites.vcf.gz \
  --loh --out-ge -t 48 \
  -o $WORK/baseline_06/tumor_phased

# Haplotag step
./longphase-to-baseline haplotag \
  -s $WORK/baseline_06/tumor_phased.vcf \
  -b $WORK/mixed_06_sorted.bam \
  -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 48 -o $WORK/baseline_06/tumor_tagged
samtools index $WORK/baseline_06/tumor_tagged.bam
```

**驗證**：`baseline_06/tumor_tagged.bam` 存在 + 含 HP tags + run.log 記 purity 數值

### 階段 C：跑 LongPhase-TO V5 on Mixed BAM（45-50 min）

```bash
./longphase-to phase \
  ... (同 baseline) \
  --pon-only-phasing \                # ← 唯一差異
  -o $WORK/v5_06/tumor_phased

./longphase-to haplotag \
  ... (同 baseline) \
  -o $WORK/v5_06/tumor_tagged
samtools index $WORK/v5_06/tumor_tagged.bam
```

**驗證**：`v5_06/tumor_tagged.bam` 存在 + 含 HP tags + run.log 記 purity（預期 ≈ 0 由於 PON-only flag）

### 階段 D：HP Tag Concordance 分析（30-45 min）

寫 Python 腳本 `InterSubMod/scripts/analysis/v5_purity06_concordance.py`：

```python
# 對 15 既有位點（同 v5_audit_suite/data/per_site_concordance.tsv）
# 取 reads 從 4 BAM:
#   - baseline_93 (existing): /big7_disk/.../baseline/tumor_tagged.bam
#   - baseline_06 (NEW):     $WORK/baseline_06/tumor_tagged.bam  
#   - v5_93 (existing):      /big7_disk/.../pononly_v5_somatic_fallback/tumor_tagged.bam
#   - v5_06 (NEW):           $WORK/v5_06/tumor_tagged.bam
#   - paired_pure (gt):      /big8_disk/.../HCC1395BL_ONT_5khz_..._tagged.bam
# 
# 計算 4 metric for each (BL_06, V5_06) pair:
#   L1 HP family concordance to paired
#   L2 HP exact concordance to paired
#   L3 HP1:HP2 ratio distance to paired
#   L4 orientation-corrected match
# 
# 比較 0.93 場景與 0.6 場景的 V5 vs baseline 改善幅度
```

**輸出**：
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/purity06_concordance.tsv`
- 4 圖表（per-site 0.93 vs 0.6 對比、ratio scatter、self-phasing bias 量化、win/loss summary）

### 階段 E：撰寫驗證報告（30-45 min）

新增到 audit suite：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md`

結構：
- Section 1: 模擬設計（mix BAM + 重 phasing）
- Section 2: 0.93 vs 0.6 兩場景的 self-phasing 強度對比
- Section 3: V5 vs Baseline 在 0.6 purity 下的 HP tag concordance to paired
- Section 4: Per-site 改善幅度（與 0.93 場景對比）
- Section 5: 推論 vs 實證對比（驗證理論預測是否成立）
- Section 6: F1 caveat（caller-level 不可驗證，已知 phasing-level）
- Section 7: 結論：低 purity 場景下 V5 仍可信？

並更新 `00_INDEX.md` 加入新文件入口。

---

## 關鍵檔案路徑

**輸入**：
- Tumor BAM: `/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`
- Normal BAM: `/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam`
- Existing VCF (workaround): `/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz`
- Reference: `/big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta`
- PON × 4: `/big7_disk/liaoyoyo2001/data/PON/`

**輸出**：
- Mixed BAM: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/mixed_06_sorted.bam`
- Baseline 0.6 tagged: `.../baseline_06/tumor_tagged.bam`
- V5 0.6 tagged: `.../v5_06/tumor_tagged.bam`
- 分析 TSV: `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/purity06_concordance.tsv`
- 圖表: `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/09_purity06/`
- 報告: `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md`

---

## 預期結果（推論驗證點）

| 推論 | 預期實證結果 | 驗證 metric |
|------|-------------|------------|
| baseline 在 0.6 下 self-phasing 弱化 | HP1:HP2 bias 從 17:1 → 3-5:1 | 階段 D Section 2 |
| V5 在 0.6 下行為一致 | bias 接近 1:1 | 階段 D Section 2 |
| V5 vs baseline 視覺差異變小 | 6/15 翻轉 → 預期 2-4/15 | 階段 D Section 4 |
| 全基因組 paired concordance | V5 仍勝出但幅度減小（+13.3pp → 預期 +6-9pp）| 階段 D Section 3 |
| Imbalance ratio 距離 paired | V5 較 baseline 接近 | 階段 D Section 4 |

如果實證結果**符合**推論 → 強化 V5 在不同 purity 場景的可信度
如果**不符**推論 → 揭示新的限制或機制，需修正 audit suite 結論

---

## 風險與緩解

| R-ID | 風險 | 緩解 |
|------|------|------|
| **R1** | ClairS-TO 不可重跑（已知 blocker）| Workaround：重用既有 VCF；報告明確標註「caller-level F1 不可驗證」|
| **R2** | 磁碟邊緣（6 TB → 預計 -550 GB → 5.5 TB）| 中間檔即跑即刪；最後 offload tagged BAM 到 /big8_disk |
| **R3** | LongPhase-TO 跑 mixed BAM 可能 fail（normal reads 含 HP tag）| 階段 A 必先 strip HP/PS tags（已加入）|
| **R4** | Mixed BAM 與 existing VCF 不匹配（VCF 來自 0.93 樣本）| 預期：低 AF somatic 在 mixed BAM 中支持 reads 變少，phasing graph 較弱 → 仍可分析但結果是「VCF-fixed phasing」而非「AF-aware phasing」|
| **R5** | 兩 BAM 的 read group / index 衝突 | samtools merge 預設 unique RGs；先看 header 確認 |
| **R6** | 預估時間實際翻倍（mixed BAM 較大可能慢）| 串行 ~3 hr → ~6 hr，仍可接受 |
| **R7** | 結果與推論不符 | 報告誠實揭露；可能揭示新發現 |

---

## 時間估計

| 階段 | 估計 | 累計 |
|------|:----:|:----:|
| A. Mix BAM | 45 min | 0:45 |
| B. baseline LongPhase-TO | 50 min | 1:35 |
| C. V5 LongPhase-TO | 50 min | 2:25 |
| D. 分析 + 圖 | 45 min | 3:10 |
| E. 報告 | 45 min | 3:55 |
| **總計** | **~3.5-4 hr** | |

可平行：A 完成後 B 與 C 可並行（兩個 IGV 進程競爭 CPU 但不會死）→ 縮短到 ~3 hr。

---

## 階段切分（漸進式執行選項）

如果擔心 4 hr 太長，可分階段確認：

### Mini Phase 1（短，~1 hr）：BAM mix + 量化檢查
- 只做階段 A
- 抽 chr19 局部驗證 mixed BAM 在 4 個 V5max/SP 位點上的 HP tag 分布
- 看 read counts、AF 是否符合 0.6 purity 預期
- → 若 AF 結構合理，繼續完整流程

### Mini Phase 2（中，+1.5 hr）：完整 phasing 兩版本
- 跑 baseline + V5（並行更省時）
- 不寫完整報告，只看 phasing.log + run.log 確認 purity 數值

### Mini Phase 3（短，+45 min）：分析 + 報告
- 寫 09_purity06_simulation.md
- 整合到 audit suite

---

## 與既有 audit suite 的關係

本任務產出**第 9 個文件**加入 v5 audit suite：

| 既有文件 | 補強 |
|---------|------|
| 06_v5_sanity_bug_check.md | 新增「purity 變化下守恆律是否仍成立」驗證 |
| 07_paired_ground_truth_concordance.md | 新增「0.6 purity 場景對比」section |
| 08_synthesis_conclusions.md | 加入新章節「跨 purity 場景結論穩定性」|
| 00_INDEX.md | 加 09_purity06_simulation.md 入口 |

---

## 合理性總評

| 評估維度 | 判定 | 說明 |
|---------|:----:|------|
| 方向合理性 | ✅ | PI 直接問此情境，需實證 |
| 風險可控 | ⚠ | 磁碟邊緣 + 長計算 + 一個已知 blocker |
| 依據充足 | ✅ | BAM/binary/PON 都齊備 |
| 完成度 | ⚠ | caller-level F1 不可驗證（已標註）|
| 可終止性 | ✅ | 三 mini phase 可中斷（Mini Phase 1 後可決定是否繼續）|
| 對既有結論的影響 | 中 | 結果可能修正或強化 audit suite 結論 |

**建議策略**：執行 Mini Phase 1 看 BAM mix 結果，依需要決定是否進 Phase 2-3。

---

## 開始前需用戶決定

| 決策點 | 選項 |
|-------|------|
| 是否啟動 | A. 完整執行（~3.5-4 hr）／ B. 只 Mini Phase 1（~1 hr）／ C. 維持理論推論不執行 |
| Caller-level F1 | A. 接受不可驗證（workaround 用既有 VCF）／ B. 等 ClairS-TO 環境設定後再做 |
| 並行 baseline + V5 | A. 並行（省時但風險高）／ B. 串行（穩但久）|
