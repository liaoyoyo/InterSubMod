<!--
建立時間: 2026-07-15 16:33 +0800
目標: 在 workstation 加入七資料集拓撲比例比較與 HCC1395 跨技術相似性驗證前，鎖定 grain、effect-size 與 claim ceiling。
處理範圍: current canonical v5；GRCh38 chr1-22；7 datasets／6 biological samples；HCC1395 vs HCC1395_DORADO。
cycle_id: cycle_20260715-1633-sample-topology-comparison
topic: sample_topology_comparison
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology/current_v5_read_af_topology.index.json
  - InterSubMod/research/20260714_hcc1395_unknown_k_clone_state_consistency/results/unknown_k_consistency.json
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/topology_pair_analysis.md
-->

# Pre-Decision Audit：七資料集拓撲比較與 HCC1395 技術配對

> **Verdict**：current-v5 描述性比較與跨技術一致性 dashboard **GO**；把組成接近或同區 agreement 解讀成 clone truth／biological replication **NO-GO**。
> **任務類型**：B Comprehensive validation；全 chr1–22 × 7 datasets，不以單一樣本或 historical layered-v2 代替 current-v5 母體。

## §0 🟤 Cynefin Domain Gate

- **Domain**：Complicated。
- **Test**：同一 sidecar 的互斥類別、region key、TVD、Cohen's κ 與 chromosome-block bootstrap 都能重複得到可預測結果；生物 clone 解讀仍不可由 bulk regional topology 單獨辨識。
- **處置**：固定 current-v5 schema、事先鎖定門檻、同時報 aggregate composition 與 exact-region concordance。

## §1 🟢 Observation Completeness Checklist

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 7/7 sidecar `all_checks_pass=true`，每維度回加 `W_primary` | ✓ | L1 ⭐⭐⭐⭐⭐ | current-v5 sidecar index |
| 7 datasets 實際代表 6 biological samples | ✓ | L1 ⭐⭐⭐⭐⭐ | manifest／HCC1395 5kHz vs DORADO 設計文件 |
| HCC pair retained-site Jaccard=84.87%，exact allele-set region Jaccard=67.31% | ✓ | L1 ⭐⭐⭐⭐⭐ | 20260714 unknown-K current-v5 artifact |
| HCC pair exact-region determinacy agreement=75.90%，pair-complete shape agreement=63.78% | ✓ | L1 ⭐⭐⭐⭐⭐ | 同上 |
| 新 sidecar 三維度 per-sample composition 與 21 pair distance | ✓ | L1 ⭐⭐⭐⭐⭐ | sample_topology_comparison artifact／receipt |
| 新 sidecar exact-coordinate HCC pair 三維 agreement／κ／block CI | ✓ | L1 ⭐⭐⭐⭐⭐ | 同 artifact＋chromosome TSV |
| Desktop/mobile comparison dashboard | ✓ | L1 ⭐⭐⭐⭐⭐ | focused 30/30＋full regression 219/219 |

## §2 🔵 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | compositional TVD、matched agreement、κ 與 block bootstrap 定義清楚 |
| 觀察支撐 | 20 | 7/7 current sidecars + 既有 current-v5 HCC pair evidence |
| 機制清晰度 | 20 | dataset region → mutually exclusive label → composition／exact join → comparison |
| 反例風險 | 10 | coverage、basecaller、region construction、CN/LOH 與 dominant-class prevalence 仍會影響結果 |
| 所需資源 | 10 | 完整分析、HTML 與 QA 約 1–6h；不重跑 BAM／solver |
| **TOTAL** | **80 / 100** | **GO** |

**Falsifier observable**：若 HCC pair 只在 marginal 比例接近、但 exact-region overlap 或 matched agreement 低，便不能稱區域結構可重現；若三維度 HCC profile distance 不在 21 pairs 前列，也不能稱它是最相似 pair。

**Reality-test 三個反例觀察**：

1. dominant unresolved／direct 類別可能抬高 raw agreement，因此必報 κ、balanced agreement 或 confusion，不只報百分比。
2. 7-dataset micro aggregate 會重複計入 HCC1395，因此另報 6-biological-sample macro，先平均 HCC technical pair 再跨 biological ID 平均。
3. exact-coordinate region 相同仍不保證 allele set、HP orientation 或 tree signature相同；dashboard 必連到 current-v5 exact allele-set與 shape 證據。

## §3 🟡 Assumption Map

| Assumption | Importance | Known | Quadrant／處置 |
|---|---|---|---|
| 三個 sidecar label 維度皆是同一 `W_primary` region grain | HIGH | known | (1) 每樣本、每維度守恆 fail closed |
| HCC exact `region` key 能一對一 join | HIGH | unknown | (2) 先驗證 key uniqueness、intersection、coverage |
| profile composition 可代表樣本的全部 biology | HIGH | false | (1) 明示只代表 current pipeline 的 regional topology profile |
| HCC兩列是 biological replicates | HIGH | false | (1) 固定標成同一 biological sample 的 technical/cross-platform pair |
| 單一綜合距離可取代每維度結果 | HIGH | false | (1) 平均距離只作 navigation，三維度仍分列 |
| 逐 region 可視為 IID | HIGH | false | (1) 以 chromosome block bootstrap 表 uncertainty |

## §4 🟠 Quick Pilot 與事先鎖定判讀

1. 讀取 7 sidecars，重算每類 count／`W_primary`。
2. 以三維 normalized composition 計算 TVD；平均三個 TVD 僅命名為 `profile_mean_tvd`，不可稱 biological distance。
3. HCC exact-coordinate primary-both 母體計算 raw agreement、Cohen's κ、confusion；以 22 autosome block bootstrap 產生 95% CI。

**事先鎖定的描述性門檻（operational，不是領域 universal cutoff）**：

- TVD：`≤0.05` 接近、`>0.05–0.10` 中度差異、`>0.10` 明顯差異。
- exact-region Jaccard：`≥0.80` 高、`0.60–<0.80` 中度、`<0.60` 有限。
- matched raw agreement：`≥0.80` 高、`0.65–<0.80` 中度、`<0.65` 有限；必與 κ／confusion 一起讀。
- κ：`≥0.60` substantial、`0.40–<0.60` moderate、`<0.40` limited；只作描述性分帶。
- HCC profile rank：在 21 dataset pairs 的 top 3 才稱「相對接近」；否則只能分維度敘述。

**Checkpoint**：所有守恆／key uniqueness／join denominators 通過才進 HTML；任一核心檢查失敗則 dashboard fail closed。

## §5 🔴 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| 7-sample current-v5 comparison artifact | HIGH | 1–2h | P0 |
| HCC sidecar exact-coordinate three-layer agreement + block CI | HIGH | 1–2h | P0 |
| Heatmap／sample profile／HCC pair audit UI | HIGH | 1–2h | P0 |
| Allele-specific CN、purity、single-cell clone truth | HIGH for biology；不影響描述性比較 | >6h | P2／out of scope |

## §6 🟣 Evidence Conflict Scan

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| HCC pair current-v5 site backbone high、regional shape only moderate | L1 | constraint；預期多層答案而非單一相似／不相似 | 20260714 unknown-K report |
| historical v2 coarse agreement 69.39%、κ=0.497，strict tree-set更低 | L1 engineering snapshot | supporting pattern；版本不可混算 | 20260712 topology pair analysis |
| true clone K／unique clone tree 不可由 current bulk outputs識別 | L2 | claim ceiling | 20260714 unknown-K report |
| gene/drug strata未高於 matched background | L1 | unrelated negative；禁止拿 annotation 當 topology truth | 20260712 gene/drug validation |

**Conflict count**：2 個直接限制；反例風險維持 10，不升 20。

## Red-team Gate

- **Failure mode 1**：用 HCC pair 的同生物來源關係倒推它必須最接近，形成 circular validation。防線：門檻與 rank 事先鎖定，結果可 fail。
- **Failure mode 2**：把 50,215 個 dataset-regions 當獨立 biological observations。防線：dataset／biological-sample 雙口徑，CI 以 chromosome block。
- **Failure mode 3**：把 morphology 或 read-AF first 解讀為 confirmed clone。防線：UI 與 artifact 都保留 regional mutation-state／descriptive ranking claim ceiling。
- **Falsifier 可否證**：可；TVD、rank、Jaccard、agreement、κ 與 CI 都能支持或反駁指定層次。

## §7 ⚫ Decision Path

- **TOTAL**：80 / 100
- **Verdict**：**GO for current-v5 descriptive comparison + HCC technical reproducibility dashboard**。
- **Claim ceiling**：比較的是 current pipeline 下的 regional mutation-state topology profile；不是 biological clone distance、accuracy ground truth 或獨立 replication。
- **Next action**：implementation notes → reproducible artifact → index integration → Playwright QA。
- **Decision lock**：Y；門檻不因觀察結果改動。

## Provenance Footer

- **Commit hash**：`b95f19e7b0ffbbe9322ad91df7858d8c696a4036`
- **Build time**：2026-07-15 16:33 +0800
- **Skill**：`/pre-decision-audit` v0.1
- **Cycle**：`InterSubMod/state/cycles/cycle_20260715-1633-sample-topology-comparison/`
