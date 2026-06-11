---
title: "cis-ASM candidate 的 copy-partition confound 三項量化驗證"
date: 2026-06-03
status: in_progress
tier: "⭐3 (L1 decisive decomposition, single-sample)"
task_type: A_pilot
cycle_id: 20260603_copy_partition_3test_quantified
corrects: 20260603_cis_candidate_copy_partition_confound
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/copy_partition_confirm.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/survivor_permutation.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/cis_scan_full.json
scripts: research/tsg_promoter_asm_reviewer/scripts/37_copy_partition_confirm.py,research/tsg_promoter_asm_reviewer/scripts/38_survivor_permutation.py
sample: HCC1395 paired_full (single sample)
partial_flag: "subset — 7/60 T3 可做 copy/allele 拆解（coverage-limited）；4 loci permutation"
---

# cis-ASM candidate 的 copy-partition confound — 三項量化驗證

> ## ⚠ 修正 banner（2026-06-07，源碼審計後）
> **「copy-partition」應讀作「subclone/copy partition」。** longphase-S 源碼審計（`HaplotagStrategy.cpp:505-516`）確認 **HP1-1 = germline 單倍型1 + somatic ALT = 一個 somatic SUBCLONE tag，不是「一個 copy」**。因此本報告的 `d_copy` 是「**somatic-subclone vs germline 甲基差**」（已證非 CN-dosage → copy-identity 或 subclone 表觀），非純 copy。
> **質性結論仍成立**（d_within 經三重證據驗為有效 subclone-control：HP1-1/ref reads 19/19 帶其他 somatic ALT、focal confident REF、甲基歸 subclone 群）——BRCA2 = **subclone/copy 主導**（焦點 cis 小）。⚠ **2026-06-07 capstone（ledger 94）supersede**：原寫「~80% copy → ~81% subclone/copy（0.81 不變）」的**精確 % 不 robust**（加性分解 `d_HP≈d_copy+d_within` 跨位點不閉合：BRCA2 殘差 9%、chr17 達 40%）→ **改質性框架「subclone/copy 主導」，勿引用精確 %（81%/20%）**。chr17 d_within=0.142 allele 主導 = 乾淨 cis（**強化**）。
> **「7/60 測不了」修正**：subclone CONTROL（d_within）限 leaky-tag；但 raw 焦點對比 `d_focal_CLEAN`（germline 相位+實際 allele）**10/11 可算**，最強位點 untestable 是因焦點突變在其 subclone 內 clonal（無 within-subclone REF），非覆蓋偏差 artifact。
> 完整：`InterSubMod/docs/experiments/in_progress/2026/06/20260607_HP11_tag_definition_correction_record_01.md`（ledger 93）。

## 問題

P3 scan 用 **HP-axis**（HP1 germline vs HP1-1 somatic）在 816 個 nonLOH ASM 位點找到 **60 個 T3 cis-candidate**。
對抗 workflow 的 Angle 4 提出：在 **CN-gain（擴增）區**，HP1 與 HP1-1 兩個 tag 可能對應**不同 copy**，
各自帶 copy-specific 甲基化 → HP-axis 把「**copy 之間的差異**」誤當「**somatic allele 的 cis 效應**」。
本驗證量化回答：**60 個 T3 是真 cis，還是 copy-partition artifact？**

## 方法 — 4 種量測拆解 copy vs allele

對每個位點（tumor reads），用 longphase HP-tag × MSA per-read allele（ALT/REF 由 SNV 定義）算 4 個 paired per-CpG dbeta：

| 量測 | 定義 | 捕捉什麼 |
|------|------|---------|
| `d_HP` | HP1 vs HP1-1（all alleles）| 原始訊號 = **copy + allele 混合** |
| `d_copy` | HP1/ref vs HP1-1/ref（**同 REF allele**）| 純 **copy/tag** 差異（allele held constant）|
| `d_within` | HP1-1 內 alt vs ref（**同 tag**）| 純 **somatic-allele cis**（copy held constant）|
| `d_allele` | 全 tag 的 ALT vs REF（pooled）| ALLELE-axis（檢查是否也 copy-confound）|

近似分解：**`d_HP ≈ d_copy + d_within`**（copy 成分 + allele 成分）。

**殘餘 copy 控制（permutation）**：在 somatic tag 內，把 read 的 allele label 隨機重排 2000 次 → null。
觀測 `d_within` 落在 null 之外 = 真 allele-specific；落在 null 之內 = 僅 within-tag/copy 異質性。

## 結果

### Test 1 — within-tag/copy-aware 重量測（7/60 可算）

> 可算 = somatic tag 內同時有 ALT 與 REF reads（各 ≥3/CpG）；其餘 53 個低覆蓋無法拆解。

| locus | d_HP | d_copy | d_within | n_cpg | dominance |
|-------|-----:|-------:|---------:|------:|-----------|
| locus | d_HP | d_copy | d_within | n_within(CpG) | dominance |
|-------|-----:|-------:|---------:|------:|-----------|
| chr17:79991120 | 0.122 | 0.029 | **0.142** | 126 | **allele 主導（真 cis）**† |
| chr18:11741161 | 0.163 | 0.060 | **0.134** | 27 | allele 主導，但 **CREATES_CpG** ⚠ |
| chr20:7482356 | 0.178 | 0.138 | 0.083 | 18 | copy 主導 + 小真 cis |
| chr20:7482362 | 0.178 | 0.137 | 0.086 | 18 | copy 主導 + 小真 cis |
| chr5:6201328 | 0.129 | 0.055 | 0.059 | 18 | copy≈allele |
| chr14:58041434 | 0.125 | 0.059 | 0.058 | 15 | copy≈allele |
| chr13:32315128 (BRCA2) | −0.122 | −0.110 | −0.023 | 197 | copy 主導（d_within/d_HP=19%）|

† chr17：`d_within(0.142) > d_HP(0.122)`（成分大於整體）——HP-axis pooling 稀釋了真 within-tag 效應，佐證加性分解僅近似。

- **COLLAPSED ≈ 5/7**（heuristic 閾值 |d_within|<0.5·|d_HP|，**post-hoc 非預先註冊**；3 位點貼閾值，取 0.45 則 5/7→2/7，閾值敏感）：HP-axis 量級多由 copy 主導。
- **連續比值**：median |d_HP|=0.129 vs |d_within|=0.083 → **ratio-of-medians 64%**；per-locus 典型保留（median-of-ratios）**47%**。真 cis 量級普遍只剩 HP-axis 約一半。
- **⚠ `d_HP ≈ d_copy + d_within` 僅近似**（三量測在不同 read 子集/CpG 集合估計，殘差 6–40%，chr17 達 40%）；非乾淨加法分解，勿據此精確分帳。

### Test 2 — ALLELE-axis 是否也 copy-confound？（locus-dependent，非全稱）

逐一檢查 7 位點 `d_allele` 較接近 `d_copy`（confound）還是 `d_within`（乾淨）：

| | d_allele 追 d_within（乾淨）| d_allele 追 d_copy（confound）|
|---|---|---|
| 位點 | chr17/chr18/chr20×2/chr14（5/7）| **BRCA2 + chr5**（2/7，copy-dominated）|

- **不是全稱**：5/7 位點 ALLELE-axis 其實追真 cis；只有 copy-dominated 位點（BRCA2、chr5）被 copy 污染。
- BRCA2 是**例外**：`d_allele=−0.099` ≈ `d_copy=−0.110`，遠離 `d_within=−0.023`。
- **修正口徑**：copy-dominated 位點 ALLELE-axis 同樣 confounded（BRCA2 屬此）；allele-dominated（chr17）反而追真 cis。故 6/01「BRCA2 ALLELE-axis cis」撤回**限定 BRCA2**，不可泛化成「所有 ALLELE-axis 都不乾淨」。唯一**普遍**乾淨指標仍是 `d_within`。

### Test 3 — per-read 殘餘 copy permutation（4 loci，2000×）

| locus | within-tag obs | null 95% CI | p | 判定 |
|-------|---------------:|-------------|--:|------|
| chr17:79991120 | 0.142 | [−0.062, 0.066] | **0.001** | 真 allele-cis |
| chr18:11741161 | 0.134 | [−0.102, 0.099] | **0.006** | 真 allele-cis |
| chr20:7482356 | 0.083 | [−0.045, 0.047] | **0.001** | 真 allele-cis（小）|
| chr13:32315128 (BRCA2) | −0.023 | [−0.021, 0.022] | **≈0.02** | **marginal**（緊貼 0.05）|

**分級：3 強顯著（p≤0.006）+ BRCA2 marginal**。即使 copy 主導的 chr20，within-tag allele 效應仍強顯著非零。⚠ BRCA2 是唯一 borderline：obs 僅略超 null CI、n=19/19、兩次 perm 給 p=0.022/0.024（MC 噪音，皆<0.05）→ 讀為「弱證據傾向非零、接近偵測下限」，需 COLO829 確認。

## 結論（修正 entry 80「overturned」過強 + 審查再細緻化）

1. **somatic-cis-ASM 為真現象** — 4/4 permutation 顯著（3 強 + BRCA2 marginal）；不是純 artifact。單樣本，需複現。
2. **但 HP-axis 嚴重高估量級** — copy-partition 主導多數位點 d_HP；BRCA2 真 cis 對 HP-axis 僅 **19%**（非「copy 90%」——加性分解不閉合，用比值陳述）。
3. **真 cis 量級** — BRCA2 僅 −0.023（小、marginal）；**chr17:79991120 (0.142) 是 7 個可做 copy-test 的位點中唯一三重乾淨**（copy-control 存活 + permutation 強顯著 + mech=NEUTRAL；自身 survive Bonferroni）。⚠ 但見限制：最強 6 個位點純 ALT tag 無法 copy-test，chr17 非全域最強乾淨。
4. **chr18:11741161 降級** — d_within 大且顯著，但 **mech=CREATES_CpG**（T>C 造 CpG），permutation 排不掉此機械假象（新 CpG read 恆在 ALT 側）→ 「需先排除 CpG-creation」。
5. **ALLELE-axis copy-confound 是 locus-dependent**（非全稱）— 5/7 位點 ALLELE-axis 追真 cis，只有 copy-dominated（BRCA2/chr5）confounded。**普遍正確指標是 `d_within`**。
6. **BRCA2 不是最強 cis 例子** — copy 主導 + 小 marginal 真 cis；chr17 才是乾淨強 cis。

## 對 pipeline / 論文的意涵

- **pipeline 修正**：P3 `cis_test` 加 `d_within`（within-tag-allele）當**主 cis 指標**，`d_HP` 降 screen-only；tier 改以 copy-control 後 d_within + per-locus permutation 判定。
- **量級修正聲明**：引用 BRCA2 cis 必標「真 allele-cis 對 HP-axis 僅 19%（d_within≈−0.023，**marginal** p≈0.02、n=19/19，需 COLO829）」。撤回 6/01「BRCA2 ALLELE-axis cis」**限定 BRCA2**，勿泛化。
- **後續（值得做）**：chr17:79991120 優先 deep-dive（基因/調控、COLO829 複現＝升⭐4 必要）；chr18 先做 CpG-creation 控制；低覆蓋 53/60 判定未定。
- **全域圖像**：與最早全域 NEGATIVE 自洽 — 沒有「方向一致、量級大」的 somatic-cis-ASM；有的是小而真的 cis + 1 個乾淨強位點，多數 candidate 表面量級是 copy 假象。

## 限制

- 單樣本 HCC1395，需 COLO829 複現；BRCA2 殘餘 cis 為 marginal，尤須複現。
- **多重檢定**：60 T3 用未校正 p<0.05 選；Bonferroni(0.05/816) 只剩 10（chr17 p_cis=0.0 自身 survive）。
- **結構覆蓋偏差（最關鍵）**：copy-test 只能在 som tag 有 ALT+REF（leaky tagging）的位點算 → 僅 7/60。最強 6 個位點（chr4 d_cis=0.706、chr20:59439285 ratio 10.9…）som tag 純 ALT(ref=0) → 無法 copy-test，cis-vs-copy 狀態未知。chr17「乾淨」來自 7 便利子集非強訊號集；此法能 demote leaky/弱候選，無法裁決最強乾淨-tag 位點。
- **association 非 causation**：permutation 證「tag 內 ALT-read 甲基≠REF-read」，非「somatic 突變 cis-驅動」（無機制/pre-post/跨樣本）。characterization-only。
- 僅 7/60 可做拆解；53 個低覆蓋無法判定，不能宣稱「全部如何」。
- 存活者生物意義（基因、調控）未查；chr18 帶 CpG-creation 機械 confound 未排除。
- `d_HP ≈ d_copy + d_within` 僅近似（不同 read 子集/CpG 集合，殘差 6–40%），非乾淨加法分解。
- `d_within` permutation 控制 within-tag 異質性，但殘餘 copy 會作為額外變異**併入** within-tag 對比、使 d_within **偏向高估** → 真 cis ≤ d_within（上界）；逼近真值需 per-read copy-phasing（本樣本無直接 copy marker）。

## Review trail（誠信紀錄）

本報告經 4-agent 對抗審查（workflow wf_cf7ec0ed-ca7）+ 主 agent 獨立重算驗證，修正 3 項初版錯誤：
1. **§5 ALLELE-axis 過度推論**（從 BRCA2 n=1 推全稱「雙軸都不乾淨」）→ 重算證實 5/7 位點 ALLELE-axis 追真 cis，改 locus-dependent。
2. **chr18 漏標 CpG-creation**（scan 已記 mech=CREATES_CpG）→ 從「乾淨 cis」降為「待釐清」。
3. **加性分解過度宣稱**（殘差達 40%、chr17 成分>整體）→ 降為近似成分對照；BRCA2「copy 90%」改比值「真 cis/HP=19%」；BRCA2 p 緩和為 marginal（雙跑 0.022/0.024）。
