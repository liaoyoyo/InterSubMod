<!--
建立時間: 2026-06-11
報告類型: 論文主軸基礎文件 (research landscape) — 新主軸取代 G6
任務類型: D handoff — 供其他 session 接手撰寫論文細節
狀態: foundation (整理現有成果 + 資產盤點 + 誠實 gap)；論文 outline/正文由後續 session 完成
data_sources: docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md
-->

# 論文主軸基礎文件：Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing

> **本文件職責**：把本研究線「所有現有成果整理到乾淨、可開啟論文目標」的狀態。
> **不是論文 outline / 正文**（那由後續 session 依此基礎完成）。
> **決策（2026-06-11 用戶確認）**：此題目**取代 G6 phasing 成為新主軸**；G6 LOH-phasing / G1 ASM 既有成果**降為支撐材料**。資料範圍**先盤點現有跨樣本資產再決定 tier**（已盤點，見 §2）。

---

## §1 一句話定位（thesis）

用 **longphase-S somatic haplotagging**（把 tumor read 分到 germline 單倍型 HP1/HP2 與體細胞子單倍型 HP1-1/HP2-1/HP3）＋ **per-read 甲基化 profile**（ONT 原生 5mCG），在 read 層級重建腫瘤的**亞克隆結構**。核心問題：**甲基化能為「以 somatic 變異為骨幹的亞克隆重建」加上什麼？**

> ⚠ **誠實前提（本研究線反覆驗證的鐵律）**：甲基訊號是 **germline-haplotype 層級**（強：分 HP1 vs HP2）、**subclone 層級弱**（分同一 haplotype 內亞群 vs 母本只有「存在性窄訊號」、不可用於 read 救援）。論文的甲基貢獻必須誠實框在這個天花板內。**絕非** variant TP/FP filter（該方向已 concluded DEAD）。

---

## §2 跨樣本資產盤點（2026-06-11 實掃，決定 tier 的關鍵）

**6 個 cell line × 3 癌種，三要素齊全**（tagged BAM + somatic_pass.vcf.gz + 甲基 modified_bases）：

| 樣本 | 癌種 | somatic-haplotag BAM | somatic VCF | 甲基(5mCG) |
|------|------|:---:|:---:|:---:|
| HCC1395 | 乳腺 | ✅ | ✅ | ✅ (100%) |
| HCC1937 | 乳腺 | ✅ | ✅ | ✅ (~70%) |
| HCC1954 | 乳腺 | ✅ | ✅ | ✅ (100%) |
| H1437 | 肺 | ✅ | ✅ | ✅ (100%) |
| H2009 | 肺 | ✅ | ✅ | ✅ (100%) |
| COLO829 | 黑色素瘤 | ✅ | ✅ | ✅ (100%) |
| (HCC1395_DORADO) | 乳腺(另 basecaller) | ✅ | ✅ | ✅ | 可作 basecaller 穩健性對照 |

- 路徑型式：`big7.../canonical/{sample}/paired_full/{date}_..._complete_matrix/longphase_s/{sample}_tagged.bam` (+ `somatic_pass.vcf.gz`)。
- normal 甲基 BAM（**2026-06-12 samtools MM-tag 實測**，見 `InterSubMod/docs/data_specs/20260612_external_data_dependencies_01.md`）：**5/6 樣本 matched-normal 有甲基**（HCC1395 5mC+5hmC · HCC1937/HCC1954/H1437/H2009 5mC）→ V10 跨樣本可跑；**只 COLO829 R10 normal 無 MM tag**（查 ONT_PAO 或重 basecall）。⚠ 6 normal 全 `zhenyu112` 帳號 + 4/6 symlink = SPOF。
- **意涵**：**單樣本 ⭐3 的硬上限可被打破** → 跨 6 樣本 × 3 癌種驗證可達 ⭐4。但全線**單一 pipeline（longphase-S，HP tag 自我參照）** 仍是 tier 風險（需正交對照）。

---

## §3 現有成果（V1-V12 map 進論文骨架；全數字 SoT = VERIFIED_RESULTS）

> 真值單一來源：`InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md`

| 論文元件 | 現有證據 | 狀態 |
|---------|---------|------|
| **背景：longphase-S tag cascade** | 已讀通原始碼（HaplotagStrategy.cpp:452 + SomaticHaplotagProcess.cpp:461，門檻0.6）；6-state 語意明確 | ✅ 完整 |
| **甲基真實性（非 copy）** | V10：matched normal copy-clean AUC 0.979 ≥ tumor 0.866，6/6 染色體；depth-matched 否證 P-06；imprinting GNAS 正控 Δβ=0.49 | ✅ 決定性（單樣本對） |
| **tag 分布 / 救援 pool** | V1：unphase 45.84%(13.8M) / HP1-1+HP2-1 ~4% / HP3 0.25%(76,738) | ✅ |
| **甲基帶 haplotype 資訊** | V6：held-out 外推 88.5%（可分位點內）；V12：但全 unphase 僅 ~6% 可嘗試（94% 無本地錨點）| ✅ 但適用性窄 |
| **T1 unphase→H1/H2** | SUPPORTED+caveats：誠實 headline V6 0.885；LOH 區未測；abstain gate chicken-egg | ⭐3 |
| **T2 H3→H1-1/H2-1** | OVERSTATED：只證有真值的 1-1/2-1 可分(0.90)；歸實際 H3 未驗(無真值/僅15-18%可指派) | ⭐2-3 |
| **T3 拆 plain H1 是否亞群** | 存在性**窄翻案**(local-allele farCpG AUC 0.85 vs 噪音0.50)；可用性 NEGATIVE(ambiguous lean<0.5) | ⭐3 |
| **subclone 甲基存在性（V11c 關鍵）** | 亞群 ALT vs 母本 REF 甲基可分（非 copy、非 somatic-CpG、無覆蓋 confound）→ **甲基確帶 subclone-discriminating 訊號** | ⭐3，論文核心候選 |

---

## §4 誠實 verdict 與必守紅線（對外引用必遵）

1. **甲基 = germline-haplotype 層級**（V10）；分不同 haplotype 強、within-haplotype 弱。
2. **T3 = 存在性窄翻、可用性 NEGATIVE**：勿無條件説「甲基能拆亞群」。亞群甲基「可分」≠「可救 ambiguous read」（lean 8/8 <0.5）。
3. **T2 勿宣稱「可歸 H3」**（只證 1-1/2-1 可分；H3 無真值）。
4. **絕非 variant filter**（DEAD）；是 subclonal reconstruction / characterization。
5. **cohesion ≠ cis**：farCpG ±100bp 太窄，subclone-program vs 突變 cis 足跡未分。
6. **跨樣本是 tier 解鎖關鍵**，但 single-pipeline 自我參照仍需正交對照。

---

## §5 開放 gap（決定能否升 ⭐4 / 可發表性）

| Gap | 需要什麼 | 優先 |
|-----|---------|------|
| **G-A 多樣本驗證** | 把 V10/V11c 跑遍樣本（**normal 甲基實測 5/6 ready**，§2；COLO829 缺甲基 normal 需補）→ 現象跨癌種復現？| P0（tier 解鎖；5 樣本即可衝 ⭐4）|
| **G-B within-haplotype somatic-vs-baseline 對照** | T3/subclone 甲基的正確 null（非 germline-het null，後者是跨-hap 錯對照）| P0（reopen T3 前置）|
| **G-C subclone-program vs 突變 cis 足跡** | 更寬空間控制（排除突變 ±1-2kb）；normal-anchored cis-test | P1 |
| **G-D 真正「重建」demo** | 從 somatic haplotag + 甲基實際輸出亞克隆結構（不只 read 分類）| P1（論文 figure 1 候選）|
| **G-E 正交 pipeline 對照** | 非 longphase-S 的第二定相/calling，破自我參照 | P2 |

---

## §6 交接給後續 session（其他 session 完成論文細節）

**本 session 已完成（整理放緩）**：
- ✅ 所有現有成果整理（V1-V12 + 誠實 verdict，SoT 維護完整）
- ✅ 跨樣本資產盤點（6 樣本 × 3 癌種三要素齊全）
- ✅ 新主軸設定（本文件 + CURRENT_FOCUS 更新）；G6/G1 降支撐
- ✅ 開放 gap 與優先級

**後續 session 應做（論文目標細節）**：
1. **跑 G-A**：把 V10（matched normal not-copy）+ V11c（local-allele subclone 存在性）跨 6 樣本重現 → 現象跨癌種？（資產已備，可直接用 §2 路徑）
2. **論文 outline**：用 `/narrative-frame` 或 `/structured-tech-report` 產 title/abstract 骨架 + figure 清單
3. **補 G-B/G-C 對照**（reopen T3 存在性前置）
4. **正式 park G6/G1 cycle**（state/cycles）+ `/pivot-direction` 記錄轉向理由
5. **G-D 重建 demo**（論文 figure 1）

**關鍵入口**：
- 真值 SoT：`InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md`（V1-V12）
- 腳本：同目錄（allele_asm / methyl_assist / rigor_t1-3 / t3_local_allele / unphase_inventory 等，皆 read-only 可重跑、可改 --bam 跑其他樣本）
- memory：`project_methyl_phasing_assist_line`（研究線總結 + 紅線）
- 機制比喻白話版：`InterSubMod/docs/experiments/in_progress/2026/05/20260611_methyl_assist_plain_explainer_01.standalone.html`
