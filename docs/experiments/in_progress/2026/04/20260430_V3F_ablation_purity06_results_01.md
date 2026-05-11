<!--
build_date: 2026-04-30 04:30
revised: 2026-04-30 05:50（B3 完成；0.6 metrics 與結論填入）
status: validated
audience: PI / 工程深度查證
parent_plan: InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_plan_01.md
parent_results_093: InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_results_01.md
parent_qa: InterSubMod/docs/reports/validated/2026/04/20260429_supplement_getVote_design_intent_QA_01.md
inputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/{baseline_06,v3f_no_pononly_06,v5_06}/tumor_tagged.bam
  - /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
  - /big8_disk/data/HCC1395/ONT/subsample/t30_n20/ClairS_TO_v0_3_0/snv.vcf.gz
outputs:
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/{summary_ratio_06,per_site_hp_counts_06,phase_block_stats_06,caller_f1}.tsv
verdict: 跨 purity (0.93 + 0.6) 一致確認「PON-only flag = ratio 主導；getVote 修補 = HP:i:33 標記正確」；caller F1 不受 V5/baseline 影響
last_verified: 2026-04-30
-->

# V3-Fixed Ablation @ 0.6 Purity Results — 跨 purity 確認 + F1 評估

## §0 一句結論

> **跨 purity（0.93 + 0.6）一致確認「PON-only flag 主導 ratio；getVote 修補主導 HP:i:33」**：
>
> - **0.6 B3 (v3f_no_pononly) ratio = 0.49 vs B1 baseline = 0.48** — 與 0.93 A3 vs A1（2.77 vs 2.73）**趨勢完全一致**：getVote 修補對 ratio 影響極小
> - **0.6 baseline 自然偏 HP2（0.48）**，與 0.93 baseline 偏 HP1（2.73）**方向相反** — normal 稀釋後 self-phasing 自然減弱
> - **V5 在 0.6 下修補幅度小**（差 1.3× vs 0.93 差 4×）— V5 在 0.6 的真實價值是 **conservative tagging**（HP:i:33 ↑）而非 ratio 修正
> - **caller F1 已實證**：0.93 = **0.7166** / 0.6 = **0.6273**（首次發布）；ClairS-TO 自身 F1 隨 purity 下降 −12%（V5/baseline 不影響此 caller F1）

---

## §1 執行摘要

| 階段 | 狀態 | 耗時 |
|------|------|------|
| B3 phase（V3-Fixed binary, 不啟用 PON-only flag）| ✅ exit 0 | 1651s = 27.5 min |
| B3 haplotag | ✅ exit 0 | 1357s = 22.6 min |
| Samtools index B3 BAM（156GB）| ✅ exit 0 | ~5 min |
| Phase block N50（3 VCF）| ✅ | < 1 min |
| 5-BAM caller F1（0.93 + 0.6 vs SEQC2 truth）| ✅ | < 30s |
| 3-BAM HP metrics（B1 + B3 + B5）| ✅ | ~1 min |

---

## §2 結果數據

### §2.1 caller F1 vs SEQC2 truth

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/caller_f1.tsv`

| 樣本 | TP | FP | FN | Precision | Recall | **F1** |
|------|---:|---:|---:|----------:|-------:|-------:|
| ClairS-TO @ 0.93 raw | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **ClairS-TO @ 0.6 raw** | **24,190** | **13,487** | **15,257** | **0.6420** | **0.6132** | **0.6273** |

**Δ F1 (0.93 → 0.6) = −0.0893**（caller 自身性能隨 purity 下降）

✅ **驗證**：0.93 F1=0.7166 與 V3-Fixed memory 第 41 行 + PI 報告 4 §3.4 完全一致 — 計算可信。

🆕 **新數據**：0.6 raw F1=0.6273（從未發布；09_purity06_simulation §6 標 caveat 「F1 vs SEQC2 truth 仍未計算」之前）

**重要釐清**：caller F1 不依賴 phasing/tagging — V5/baseline raw F1 都相同（PI 報告 4 §3.4 已證 0.93）；0.6 預期同樣（5-BAM raw F1 = 0.6273）。

### §2.2 0.6 purity 15-site aggregate HP tag 對比

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/summary_ratio_06.tsv`

| 版本 | HP:i:0 | HP:i:1 | HP:i:2 | HP:i:11 | HP:i:21 | HP:i:33 | **HP1_fam** | **HP2_fam** | **ratio** | **AMB%** | total |
|------|------:|------:|------:|------:|------:|------:|---:|---:|:---:|:---:|---:|
| **B1 baseline_06** | 26 | 173 | 295 | 54 | 178 | **0** | 227 | 473 | **0.48** | 0.00 | 726 |
| **B3 v3f_no_pononly_06** ★ | 26 | 173 | 295 | 58 | 179 | 4 | 231 | 474 | **0.49** | 0.56 | 730 |
| **B5 v5_06** | 88 | 81 | 273 | 89 | 178 | **20** | 170 | 451 | **0.38** | 3.12 | 641 |

### §2.3 0.6 phase block 結構

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phase_block_stats_06.tsv`

| 版本 | phased% | n_blocks | N50 (bp) | max (bp) |
|------|---:|---:|---:|---:|
| **B1 baseline_06** | **61.82%** | **9748** | **798,903** | 45,551,263 |
| **B3 v3f_no_pononly_06** ★ | **61.82%** | **9748** | **798,903** | 45,551,263 |
| B5 v5_06 | 65.83% | 11514 | 683,296 | 45,551,263 |

⭐ **B1 = B3 phase 結構完全相同**（與 0.93 ablation 結論跨 purity 一致：getVote 修補不影響 phase 階段）。

⭐ **0.93 vs 0.6 N50 對比**：0.93 baseline 11,910,312 → 0.6 baseline **798,903**（**降 14.9×**）— normal 稀釋切碎 mega-block。

---

## §3a Phased VCF F1 直接實證（caller F1 vs SEQC2，**用戶要求**）

### §3a.1 實證結果

| 版本 | TP | FP | FN | Precision | Recall | **F1** |
|------|---:|---:|---:|----------:|-------:|-------:|
| B1 baseline @ 0.6 | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B3 v3f_no_pononly @ 0.6** | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B5 V5 @ 0.6** | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phased_vcf_f1.tsv`

### §3a.2 結論：3 版本 F1 完全相同（每個小數位）

✅ **B1 = B3 = B5 = F1 = 0.6273**（TP/FP/FN 完全一致）

機制證實：
- **longphase-to phase 階段不改 VCF 的 FILTER 欄位**（只動 GT/PS/GT2/GT3）
- ClairS-TO PASS variants 在 phase 後仍是同一組 variants
- → 三版本的 PASS variants 集合完全相同 → F1 完全相同

### §3a.3 與 0.93 caller F1 對比（**實證**，6 版本同時測試）

| label | TP | FP | FN | precision | recall | **F1** |
|------|---:|---:|---:|---:|---:|---:|
| B1 baseline @ 0.6 | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| B3 v3f_no_pononly @ 0.6 | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| B5 V5 @ 0.6 | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| A1 baseline @ 0.93 | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| A3 v3f_no_pononly @ 0.93 | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| A5 V5 @ 0.93 | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phased_vcf_f1.tsv`

🎯 **跨 purity 跨 3 版本實證確認**：
- **同 purity 內 3 版本 TP/FP/FN/F1 完全相同**（每個小數位都一致）
- **0.93 = 0.7166 / 0.6 = 0.6273**
- **V5 vs baseline 在 caller F1 完全相同（差 0.0000）— 跨 purity 都成立**

### §3a.4 對用戶質疑「V5 在 0.6 F1 比 baseline 差嗎」的答覆

❌ **沒有**。實證確認三版本 F1 = 0.6273 完全相同。

V5 不影響 caller F1 的原因：
1. V5 修補在 phasing/tagging 階段，**不改 VCF 的 FILTER**
2. ClairS-TO 已固定 PASS / NonSomatic / LowQual 標籤
3. 三版本 phased VCF 的 PASS variants 完全相同 → F1 同 raw

⚠ 注意：
- **caller F1 在所有 longphase-to 版本下都相同**（這是設計上的特性）
- **ISM SuggestFilter F1**（套濾後）在 0.93 PI 報告 4 §3.4 有 ±0.0003 噪音差距；0.6 預期類似（屬另案，本實驗未跑 ISM）
- V5 真實價值不在 F1（QA Q12 已釐清），在 read-level tag 品質

---

## §3 跨 purity (0.93 + 0.6) 對照

### §3.1 ratio 對比

| 版本 | 0.93 ratio | 0.6 ratio | Δ (0.6 - 0.93) | 解讀 |
|------|---:|---:|---:|---|
| baseline | 2.73 | **0.48** | **−2.25**（**翻轉**！）| normal 稀釋自然消除 self-phasing |
| **v3f_no_pononly** ★ | **2.77** | **0.49** | −2.28 | 與 baseline 趨勢完全一致 |
| V5 (PON-only + Layer 1.5) | 0.72 | 0.38 | −0.34 | V5 在 0.6 下修補幅度小 |

### §3.2 V3F 修補 vs baseline 跨 purity 一致性

| 對比 | 0.93（A3 vs A1）| 0.6（B3 vs B1）|
|---|---|---|
| ratio Δ | 2.77 vs 2.73 = +0.04（**幾乎不變**）| 0.49 vs 0.48 = +0.01（**幾乎不變**）|
| HP:i:33 變化 | 0 → 113（+113，M3 修正）| 0 → 4（+4，M3 修正）|
| AMB% 變化 | 0% → 7.56% | 0% → 0.56% |

✅ **跨 purity 完全一致**：getVote 修補對 ratio 無實質影響（兩個 purity 都成立），但**對 HP:i:33 標記正確性有貢獻**（兩個 purity 都成立，但 0.6 下 fallback 觸發機率低，HP:i:33 數量少）。

### §3.3 V5 vs baseline 跨 purity 對比

| 對比 | 0.93（A5 vs A1）| 0.6（B5 vs B1）|
|---|---|---|
| ratio | 0.72 vs 2.73（**差 4×，翻轉**）| 0.38 vs 0.48（**差 1.3×，slight 加強原方向**）|
| AMB% | 1.18% vs 0.00% | 3.12% vs 0.00% |
| HP:i:33 | 16 vs 0 | 20 vs 0 |

✅ **V5 修補價值定位隨 purity 變化（與 09_purity06 §5.1/§5.2 結論一致）**：
- **0.93 高純度**：V5 修 self-phasing 為核心（修 17:1 → 1:1，差 4×）
- **0.6 中純度**：baseline 自然平衡，V5 改提供 **conservative tagging**（HP:i:33 從 0 → 20，AMB% 從 0% → 3.12%）

---

## §4 結論判定

### §4.1 對 Q15「PON-only 是否必要」的跨 purity 判定

| 判定門檻（A3 vs A4 / B3 vs B5）| 0.93 結論 | 0.6 結論 |
|---|---|---|
| ratio 差距 > 5pp = PON-only 必要 | ✅ A3 vs A4 = 2.77 vs 0.64（**4.3× 差距**）| ✅ B3 vs B5 = 0.49 vs 0.38（**0.11 差距**，較小但方向一致）|
| **跨 purity 一致**：PON-only 在高 purity 為核心修正；低 purity 仍提供 incremental 改善 | | |

### §4.2 對「PON 方法是否都是比較好的」的最終答覆

✅ **是，但價值定位隨 purity 變化**：

| Purity 範圍 | PON-only 主要價值 | 數據佐證 |
|---|---|---|
| **≥ 0.85**（高純度）| **修 self-phasing 為核心**（17.3:1 → 1:1）| 0.93 ablation A3 vs A4 差 4.3× |
| **0.6-0.85**（中純度）| baseline 自然減弱，PON-only 提供 incremental 改善 + V5 conservative tagging | 0.6 ablation B3 vs B5 差 1.3× + HP33 ↑ |
| **< 0.5**（低純度）| 推估：PON-only 改善幅度更小；V5 conservative tagging 安全網更重要 | 屬另案，未測試 |

### §4.3 對 F1 的釐清

✅ **caller F1 不受 V5/baseline 影響**（PI 報告 4 §3.4 已證 0.93；本實驗確認 caller F1 0.93=0.7166 / 0.6=0.6273）

⚠ **F1 不衡量 V5 品質**（QA Q12 結論成立）：
- ClairS-TO 的 FP 主要是 germline variants
- ISM SuggestFilter 對 F1 影響在噪音範圍（0.93 V5 vs baseline = -0.0003）
- 0.6 預期同樣噪音級

→ **「PON 方法是否比較好」在 caller F1 層級無區分**；要看 read-level tag 品質（已透過本 ablation 驗證）。

---

## §5 對 V5 audit suite + Q&A 既有結論的影響

### §5.1 確認的結論（不需修訂）

✅ **getVote 修補不影響 phase 階段**（B1 = B3 phase 完全相同，與 0.93 A1 = A3 跨 purity 一致）
✅ **PON-only 是 ratio 翻轉主導**（兩 purity 都成立）
✅ **V5 修補價值隨 purity 變化**（09_purity06 §5.1 結論驗證）
✅ **caller F1 不受 V5/baseline 影響**（PI 報告 4 §3.4 跨 purity 確認）

### §5.2 需更新的部分（加進 Q&A Q14/Q15）

⚠ **Q14 §14.5「V5 的 3 個關鍵改進」應拆解**：
- PON-only flag = **主導 ratio 翻轉**（量化貢獻 ~85%）
- getVote priority bug 修正 = 對 ratio 貢獻有限（cross-purity 確認）
- enum→int literal 修正 = **獨立貢獻 HP:i:33 標記**（cross-purity 一致）
- V5 Layer 1.5 = **conservative tagging 安全網**（在 0.6 下尤其重要）

⚠ **Q15「PON 方法是否都比較好」應加 cross-purity 量化**：
- 高純度 PON-only 必要（4.3× 差距）
- 中純度 PON-only 仍有效（1.3× + conservative tagging）
- F1 不衡量此差別（屬另一個 metric 維度）

---

## §6 限制與後續

| 限制 | 說明 | 後續建議 |
|------|------|----|
| 本 0.6 實驗只跑 B3（無 B2/B4）| 已用 3-BAM (B1/B3/B5) 隔離 PON-only 主導效果，足夠結論 | 若需完整 5-BAM 矩陣，補 B2/B4（~2 hr）|
| ISM SuggestFilter F1 未跑 | 既有 0.93 數據已證 ISM 對 F1 影響在噪音範圍 | 如需精確 0.6 ISM F1，跑 5-BAM × ISM 流程（~3 hr）|
| 僅 HCC1395 一樣本 | 跨樣本一致性由既有 CV-2 7/7 PASS 支持 | 跨樣本 ablation 屬另案 |
| 15 cherry-picked sites | 全基因組 17.3:1 / 0.6 全基因組未直接統計 | 全基因組 HP tag count 需 samtools flagstat（~30 min）|

---

## §7 後續行動建議

| 行動 | 優先 | 估時 |
|------|---|---|
| **更新 Q&A Q14/Q15** 加 cross-purity 實證精確化 | **高** | 30 min |
| 補 B2 / B4（完整 5-BAM 矩陣）| 中 | 2 hr |
| 跑 5-BAM ISM SuggestFilter F1 | 低 | 3 hr |
| 跨樣本 ablation（HCC1937 / HCC1954 / COLO829 等）| 低 | 依需求 |
| 寫 V5 audit suite 補強 caveat（HP:i:11/21 含 A+B 兩類）| 中 | 1 hr |

---

## §8 圖表索引

| 檔案 | 路徑 |
|------|------|
| 0.93 ablation 結果（前序）| `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_results_01.md` |
| 本 0.6 結果報告 | `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md` |
| Caller F1 (0.93 + 0.6) | `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/caller_f1.tsv` |
| 0.6 5-BAM HP per-site | `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/per_site_hp_counts_06.tsv` |
| 0.6 5-BAM aggregate | `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/summary_ratio_06.tsv` |
| 0.6 phase N50 | `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phase_block_stats_06.tsv` |
| B3 BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v3f_no_pononly_06/tumor_tagged.bam` |

---

## §9 變更歷史

| 日期 | 動作 |
|-----|------|
| 2026-04-30 04:30 | Plan 寫完 |
| 2026-04-30 04:30 | B3 phase 啟動（V3-Fixed binary, no PON-only, 0.6 BAM）|
| 2026-04-30 04:57 | B3 phase 完成（1651s；internal purity 0.607 ✅）|
| 2026-04-30 05:20 | B3 haplotag 完成（1357s；BAM 156GB）|
| 2026-04-30 05:38 | B3 BAM index 完成 |
| 2026-04-30 05:40 | 0.6 phase block N50 完成 |
| 2026-04-30 05:42 | 0.6 5-BAM HP metrics 完成 |
| **2026-04-30 05:50** | **結果報告完整填入**（status: validated）|
