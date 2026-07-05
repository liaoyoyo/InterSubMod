<!--
建立時間: 2026-07-05
類型: 工具生產執行結果 (C-production execution result)
狀態: in_progress（CN 生產完成；下游 cn-fill + topology 交並行 session）
主軸: Subclonal reconstruction — CN 輔助層（補 topology 的 cn=unknown）
data_sources: /big7_disk/liaoyoyo2001/cnv_sv_work/{COLO829,H1437,H2009,HCC1954,HCC1937}/savana_wgs/cna_normalhet/*_fitted_purity_ploidy.tsv, /big7_disk/liaoyoyo2001/cnv_sv_work/*/savana_wgs/*_segmented_absolute_copy_number.tsv, /big7_disk/liaoyoyo2001/cnv_sv_work/*/savana_wgs/*_smcnbed.bed
關聯 memory: project_savana_cna_only_pipeline_and_colo829_purity
-->

# SAVANA CN 生產 — 6 canonical 樣本補 cn=unknown（cna-only）

> **版本** v1.0 (2026-07-05) · **狀態：CN 生產 5/5 完成**（HCC1395 已有 SEQC2）· 下游 cn-fill 交並行 session
> 證據 tier：**L1** = 第一手磁碟 SAVANA 輸出重現讀回

---

## 0. 一句話結論

對 5 個 canonical ONT 樣本跑 **`savana cna`（CN-only，跳過記憶體殺手 SV `run`）** 產全基因組 CN，補 topology 工作站的 `cn=unknown`。**5/5 完成**；purity 用 raw BAF 自洽性裁決 → **3 可用**（H1437/H2009/HCC1954）、**2 SAVANA mis-fit**（COLO829/HCC1937 → 維持 cn=unknown）。連同 HCC1395（SEQC2）= **4/7 樣本有可靠 CN**。

---

## 1. 方法（L1）

| 項目 | 值 |
|---|---|
| 工具 | SAVANA 1.3.7（conda env `cnvtools` @bip7）|
| 模式 | **`savana cna` only**（🔴 不跑 `savana run`：其 SV breakpoint clustering 在癌症 WGS ONT 上 held 數百萬 breakpoints → 400-488G/樣本逼 OOM；`cna` 獨立於 run，`-bp` 選配）|
| tumor | **raw big8 BAM**（非 big7 tagged BAM — tagged 的自訂 HP tag `1-1`/`2-1`/`3` 撞 `breakpoints.py:490` KeyError）|
| normal | big8 matched normal（concordance 驗 normal_het 0.996-0.999 全 MATCHED）|
| het SNP | `-v germline_phased_merged.vcf.gz`（matched-normal，**非** `-vg 1000g`）|
| 資源 | cna 輕 ~15G/樣本（windowed allele counting）；2-3 並發 mem 8-18G 安全 |
| runner | `cnv_sv_work/run_savana_5s/run_savana_cna.sh`（repo 外 compute infra）|

---

## 2. CN 生產結果（L1，第一手讀回）

| 樣本 | purity/ploidy | CN 段數 | BAF_p90 | SM_CNBED intervals | 裁決 |
|---|---|---|---|---|---|
| H1437 | **0.95** / 2.9 | 970 | 0.994 | 931 | ✅ 可用 |
| H2009 | **0.95** / 2.22 | 721 | 0.993 | 551 | ✅ 可用 |
| HCC1954 | 0.66 / 3.21 | 843 | **0.840** | 776 | ✅ 可用（自洽）|
| COLO829 | 0.53 / 2.26 | 806 | 0.994 | —（未轉）| ❌ mis-fit |
| HCC1937 | 0.62 / 3.13 | 657 | 0.995 | 628（標 unreliable）| ❌ mis-fit |

---

## 3. 🔴 purity 裁決 — raw BAF 自洽性（purity-無關判準）

原始 meanBAF 是 purity-**無關**觀測量（fitted copyNumber 是 purity-依賴、不可用來反推 purity，循環）。判準：

- **BAF 達 ~0.99（近 1.0）** = LOH 區單 allele 近完全 → 只在 **purity ~1.0** 或**極高擴增（CN~40-60）** 才可能。
- H1437/H2009：BAF 0.99 + purity 0.95 → **一致**（純腫瘤，BAF 達 1.0）✅
- HCC1954：BAF **0.84** + purity 0.66 → **一致**（p=0.66 時 copy-neutral LOH 的 BAF 上限 ~0.83）✅ → 這批 ONT 樣本 tumor fraction 本較低，非 mis-fit。
- COLO829：BAF 0.99 但 purity 0.53，且 **near-diploid（ploidy 2.26）無高擴增** → BAF 0.99 需 CN~40（荒謬）→ **矛盾 = mis-fit**，真 purity ~1.0。
- HCC1937：BAF 0.99 但 purity 0.62，**238 段 CN<3.5 帶 BAF>0.90**（p=0.62 需 CN~60 荒謬）→ **矛盾 = mis-fit**，真 ~1.0（near-triploid+WGD 讓 SAVANA 挑退化解 0.62/3.13）。

**mis-fit 未能修**：`--min_cellularity 0.85` 被 SAVANA 無視（仍回 0.53）；`--overrule_cellularity 1.0` 改了 CN 計算但 grid 仍 0.53；**無外部 COLO829/HCC1937 CN 真值裁決** → 3 次 re-fit（各 ~5-10h read counting）無定論 → 停止再投入（⭐3 characterization 輔助）。COLO829+HCC1937 = **cn=unknown**。

---

## 4. CN 覆蓋（全 7 canonical）

- **4/7 有可靠 CN**：HCC1395（SEQC2 benchmark，既有）+ H1437 + H2009 + HCC1954（SAVANA）
- **3/7 cn=unknown**：COLO829、HCC1937（SAVANA purity mis-fit）+ HCC1395_DORADO（共用 HCC1395）

---

## 5. 交付 + 下游（handoff）

- **已 commit（8aedc86）**：`scripts/analysis/savana_to_smcnbed.py`（SAVANA CN→gain/loss/loh BED，gain-priority 決策樹複現 HCC1395 SEQC2 baseline；LOH vs SEQC2 10kb Jaccard 0.931）+ `scripts/analysis/rerun_cn_integration_from_savana.sh`（convert + 重跑 sm_region_integration 填 cn，SM_CN_SPARSE=0）。
- **🔵 下游 cn-fill + topology 重建 = 並行 session 領域**：topology code（07-04 incompatible 重分類 / 07-05 complexity）正被其活躍改動 → **不單方跑其 pipeline**。3 可用樣本的 SM_CNBED + driver 已就緒，待協調執行。

---

## 6. Provenance / caveat

- **L1**：所有 purity/ploidy/段數/BAF 皆本機 SAVANA 輸出第一手讀回（§2 data_sources 路徑）。
- **caveat**：① COLO829/HCC1937 purity mis-fit（BAF 矛盾）→ 其 CN 不可用、cn=unknown。② HCC1954 0.66 為自洽較低 purity（非 mis-fit）但仍 flag 提醒。③ SM_CNBED 閾值 gain_cn≥2.5/loss_cn≤1.5/loh_minor<0.5（feasibility doc 驗證，diploid-anchored，aneuploid 下 gain/loss 為方向參考、LOH 最穩健）。④ cnv_sv_work/ 在 repo 外（大檔 compute 區）。
