<!--
建立時間: 2026-07-15T01:46:40+08:00
目標: 驗證 truth-labeled FP sSNV 的 focal-ALT reads 無監督多群能否解讀為高可能 subclone 或 linear topology
處理範圍: 7 datasets / 6 biological samples；最新 LongPhase-S recalibrated PASS 與既有 canonical InterSubMod complete matrices
cycle_id: cycle_20260715-0146-single-fp-alt-multicluster
topic: single_fp_alt_multicluster_subclone_validation
status: verdict_PROBE_claim_GO_validation
audit_version: 0.1
關聯檔案:
  - InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/README.md
  - InterSubMod/docs/experiments/in_progress/2026/06/20260622_cluster_redesign_genomewide_tpfp_01.md
  - InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/20260714_跨HP_clone狀態反推可行性與觀察驗證紀錄_01.md
-->

# Pre-Decision Audit: 單一 FP focal-ALT 多群是否可視為 subclone

> **雙判定**：把「多群」直接升級成高可能 subclone 為 **PROBE（30/100）**；執行全量可否證驗證為 **GO**。
> **研究目標**：G3 / G4 / G5；Task Type B（全基因組、全 7 datasets，不以 subset 代替 final）。

## Frontmatter

- **Topic**: `single_fp_alt_multicluster_subclone_validation`
- **Triggered by**: new research + cross-sample validation + high-impact biological interpretation
- **Last updated**: 2026-07-15T01:46:40+08:00
- **Cycle ref**: `InterSubMod/state/cycles/cycle_20260715-0146-single-fp-alt-multicluster/`（本輪不註冊 state；研究工作區保存 audit）

## §0 Cynefin Domain Gate

- **Domain**: **Complex**。
- **Test**: 相同切群在不同 coverage、CpG density、HP、CN/LOH 與 caller truth class 下未重複產生可預測的 biological clone 結論。
- **Decision**: probe-first；先以 permutation-null、穩定度與反例對照判定「群是否真實」，再談 subclone 歸因。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 7 datasets 的 all-tumor-read silhouette `cansplit` 已量化 | ✓ | L1 | `InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/README.md` |
| 舊 silhouette 會強制 k>=2，純噪音 FP 74.7-90.7% | ✓ | L1 | `InterSubMod/docs/methodology/20260622_phylo_v41_cpp_internalization_audit_01.md` |
| HCC1395 REAL_NOVEL 在 FP 不少於 TP，不能稱 somatic-specific | ✓ | L1/L3 | `InterSubMod/docs/experiments/in_progress/2026/06/20260622_cluster_redesign_genomewide_tpfp_01.md` |
| focal-ALT-only、7-dataset、null-validated 比例 | □ | 尚未產生 | 本輪 P0 |
| 最新 LongPhase-S PASS 與正確 external HP/PS sidecar exact join | □ | 尚未產生 | 本輪 P0 |
| 代表個案由原始 read/CpG 重算及混雜審計 | □ | 尚未產生 | 本輪 P1 |
| single-cell / clone isolate / multi-region orthogonal clone truth | ✗ | 不可用 | 本資料限制 |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 10 | subclone 可伴隨 epigenetic heterogeneity，但並非唯一生成機制 |
| 觀察支撐 | 10 | 有單樣本與 all-read 跨樣本背景，尚無 focal-ALT-only 全量結果 |
| 機制清晰度 | 0 | 單一 truth-FP allele 與 methyl cluster 無法建立 cellular lineage 身分 |
| 反例風險 | 0 | 強制切群、cis-ASM、HP、CN/LOH、CpG/coverage、strand/quality 均可生成多群 |
| 所需資源 | 10 | 最新 sidecar bounded extraction + 全量 null analysis 約 1-6 小時 |
| **TOTAL** | **30 / 100** | **PROBE claim；GO falsification study** |

**Falsifier observable**：若「focal-ALT 多群 = 高可能 subclone」錯誤，則 null-validated 多群會在 FP 對照中同樣或更常見、在 permutation/stability gate 後大量消失、可由 HP/strand/CpG/coverage/CN 解釋，且沒有第二個 genetic event 或正交 clone truth 支持。

**Reality-test 三反例**：

1. permutation 後仍會被 silhouette 強制切成 k>=2，但 phylo-v4.1 應判單群。
2. focal-ALT methyl groups 若主要跟 HP、strand 或 read quality 對齊，屬歸因混雜而非 clone confirmation。
3. 一個 sSNV 只有 ROOT->ALT 的 trivial mutation-state edge；多個 methyl groups 不增加可識別的 mutation nodes。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| `alt_support=ALT` 是焦點 sSNV 本身，不是 `HP=1-1/2-1` carrier proxy | HIGH | KNOWN | verify in schema |
| 最新 LongPhase-S tags 可由 coordinate sidecar exact join，不使用舊 BAM tag 代替 | HIGH | KNOWN contract / 尚待本輪實跑 | MUST validate |
| truth-labeled FP 仍可代表真 somatic clone | HIGH | FALSE as default | claim blocker |
| null-validated methyl groups 等價 cellular clones | HIGH | UNKNOWN | MUST falsify |
| coverage/CpG/HP/strand/CN 不解釋分群 | HIGH | UNKNOWN | MUST validate |
| 單一 sSNV 可識別 linear vs branching evolution | HIGH | FALSE | mathematical blocker |

## §4 Quick Pilot

1. 讀 HCC1395 canonical FP 的 `reads.tsv + methylation.csv + BERNOULLI/matrix.csv`。
2. 僅保留 tumor `alt_support=ALT`，比較 silhouette 與 phylo-v4.1 null verdict。
3. 將 locus 對到最新 LongPhase-S recalibrated PASS，檢查 sidecar exact join 與 read/CpG 重算。

**Checkpoint**：只有在至少 3 biological samples 出現穩定 null-validated focal-ALT 多群，且有獨立 genetic co-occurrence / normal anchor 排除替代解釋時，才允許升級為 subclone candidate；欠缺 orthogonal clone truth 時不得稱 confirmed/high-probability subclone。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| 最新 LongPhase-S PASS truth-FP split | HIGH | 0.5 h | P0 |
| bounded raw BAM + latest HP/PS sidecar materialization | HIGH | 1-3 h | P0 |
| focal-ALT phylo-v4.1 + null/stability + confound census | HIGH | 1-3 h | P0 |
| matched TP control with identical ALT-read strata | HIGH | 1-3 h | P1 |
| single-cell / clone isolate truth | HIGH | external | P2 / current unavailable |

## §6 Evidence Conflict Scan

Repo root 無 `MEMORY.md`（2026-07-15 實查）；已改查 validated/in-progress SoT 與研究紀錄。

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| `cansplit` 與 clean label alignment Jaccard 0.091-0.161 | L1 | conflict | `InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/README.md` |
| HCC1395 FP within-HP multi-group 76.7%，不是乾淨 subclone signature | L1/L3 | conflict | `InterSubMod/docs/experiments/in_progress/2026/06/20260619_tp_fp_structure_association_HCC1395/20260619_tp_fp_structure_label_association_HCC1395_01.md` |
| REAL_NOVEL TP 10.7% vs FP 16.3% | L1 | conflict | `InterSubMod/docs/experiments/in_progress/2026/06/20260622_cluster_redesign_genomewide_tpfp_01.md` |
| 最新 layered L3 為 not_evaluated，methyl 只能 residual corroboration | L1/L3 | dependent | `InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/20260714_跨HP_clone狀態反推可行性與觀察驗證紀錄_01.md` |

## §7 Decision Path

- **Claim score**: 30/100。
- **Claim verdict**: **PROBE**；禁止先驗採用「高可能 subclone」或「linear evolution」。
- **Validation verdict**: **GO**；完整量化能直接關閉錯誤推論或界定可接受的 candidate tier。
- **Decision lock**: high-probability/confirmed subclone 必須具備至少一項獨立 genetic lineage evidence，且排除 HP/technical/CN/normal-cis 主要替代解釋。

## Provenance

- **Commit hash at audit**: `6067568637088838a9f518955e41d222f057e4f1`
- **Branch**: `research/subclonal-reconstruction-202606`
- **Skill**: `/pre-decision-audit` v0.1
- **Topic dir**: `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/`
