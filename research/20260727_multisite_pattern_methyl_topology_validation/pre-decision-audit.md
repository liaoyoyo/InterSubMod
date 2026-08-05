<!--
建立時間: 2026-07-27 09:37
目標: 多 sSNV pattern × methyl association 全量實作前的快速證據稽核
處理範圍: 7 technical datasets / chr1-22 / exact PS × exact raw HP
status: verdict_GO_conditional
audit_version: 0.1
關聯檔案:
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/00_INDEX.md
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md
  - InterSubMod/docs/reports/research_landscape/04_暫停判定與重評估.md
-->

# Pre-Decision Audit：多 sSNV pattern × 甲基關聯

> **Verdict**：`GO_ANALYTIC_CENSUS / NO-GO_LINEAGE_UPGRADE`。只有硬 gate
> 全部落地時才進正式統計。

## §0 Cynefin domain

**Complex → probe-informed comprehensive validation**。read identity與 Bernoulli
計算可重複，但 raw HP 自我參照、read geometry、shared CpG/read及 CN/cis
epiallele交互作用不能靠單一 best practice排除。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 7 datasets × 22 chromosomes strict sparse read calls已完整 | ✓ | L1 | `docs/CURRENT_FOCUS.md:52-60` |
| sparse row保留 `hp_raw`, `phase_set`, R/A calls與 qname SHA-256 | ✓ | L1 | 154 files schema audit |
| LongPhase-S raw HP為9態字串，`1`與`1-1`語意不同 | ✓ | L1 | KB `05_tools/longphase-s.md` §輸出 HP tag格式 |
| exact raw-HP full4 candidate可枚舉 | ✓ | L2 | preflight：family collapse 46/13；exact raw HP 27/10（n≥5/8） |
| user chr22舊圖的 methyl matrix可按 qname重綁 | ✓ | L1 | archive/canonical methyl SHA相同；reads HP標籤不同 |
| 全 7 datasets exact raw-HP methyl effect正式結果 | □ | — | 本 cycle輸出 |
| 跨 region methyl correlation能證 lineage | ✗ | L1 | O13 shared-read confound NEGATIVE |
| 單位點/cluster可直接證 biological subclone | ✗ | L1 | KB `05_tools/intersubmod.md` §生物語意邊界 |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | 同 molecule多 marker pattern與同 read methyl profile grain一致 |
| 觀察支撐 | 10 | candidate存在；全量 methyl effect尚未重算 |
| 機制清晰度 | 20 | exact PS/raw HP → pattern → common-CpG Bernoulli → association DAG清楚 |
| 反例風險 | 10 | 有既有 NEGATIVE；靠硬 gate降級但不能歸零 |
| 所需資源 | 10 | 全量 join/stat/report >6h；輸入 probe已完成 |
| **TOTAL** | **70 / 100** | **Conditional GO** |

**Falsifier observable**：若 frozen eligible family中沒有任何 unit同時通過
BY、effect、restricted PERMANOVA、PERMDISP、共同 CpG、等 N/rarefaction與
geometry gates，則「可評估範圍內存在明顯穩健 pattern-conditioned methyl loci」
被否證。

**三個反例觀察**：

1. effect在 family-collapse顯著、exact raw HP消失；
2. effect在等 N/common-CpG rarefaction或marker masking後消失；
3. PERMANOVA顯著但 PERMDISP、strand/RG或read geometry也顯著。

## §3 Assumption map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| sparse qname digest可唯一對上 methyl read_name | HIGH | 部分已知 | MUST validate |
| 每個 marker tile有可追溯 methyl artifact | HIGH | KNOWN | 20260715 all-sSNV receipt + 本輪逐 marker hash |
| exact raw HP可消除 family mixing | HIGH | KNOWN | verify |
| exact raw HP不造成 collider/self-reference | HIGH | UNKNOWN | claim ceiling + sensitivity |
| same unit-level CpG basis可建立 | HIGH | UNKNOWN | MUST validate |
| HCC1395/DORADO是兩個 biological replicates | HIGH | KNOWN false | 禁止 |

## §4 Pilot checkpoint

1. HCC1395 chr22 user example：methyl矩陣 SHA與 qname join可重現；舊
   `reads.tsv` HP不得採用。後續找到與 July topology 同一 VCF universe 的
   20260715 all-sSNV v2 artifacts，469,849/469,849 位點均有
   reads/methyl/BERNOULLI 三件組；primary 不再依賴 March canonical 交集。
2. full4 preflight：exact raw HP將 family 46降為27個 n≥5候選，證明 primary
   grain會實際改變分析母群。
3. H2009 chr5舊 full4 null (`p=0.325`) 保留為 locked negative canary。

三項支持進 full census，但不支持 lineage claim。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| per-marker methyl artifact-source catalog與hash | HIGH | 1–3h | P0 |
| 154-file qname collision/immutable metadata audit | HIGH | 1–3h | P0 |
| restricted PERMANOVA/PERMDISP與R² API | HIGH | 1–6h | P0 |
| marker-created CpG、CN/LOH與geometry sensitivity | HIGH | 1–6h | P1 |
| second phaser/marker-held-out haplotag | HIGH for causality | separate cycle | P2 |

## §6 Evidence conflict scan

- `MEMORY.md` 不存在；改查 `docs/CURRENT_FOCUS.md`、`docs/experiments/INDEX.md`
  與 append-only `research/autoresearch/evidence_ledger.jsonl`。
- O13 cross-region correlation為 shared-read confound NEGATIVE：本題只能做
  unit內檢定，並強制等 N/common-CpG sensitivity。
- TO HP-based methyl特徵多項 NEGATIVE：本題只用 paired LongPhase-S，仍保留
  tag-dependence/collider警告。
- `20260726_within_hp_near_distal`把 `1`/`1-1` collapse，且 TP不優於 FP；
  本題是 exact raw-HP正式重驗，不得把新 finding歸因 somatic lineage。
- single-locus methyl cluster不是 subclone、strict W/edge不是 ancestry：
  與任何 lineage升級直接衝突，故該 claim固定 NO-GO。

## §7 Red-team gate與決策

獨立紅隊列出六項最強 failure modes：raw HP self-reference、read geometry、
O13型 N/shared-CpG、local cis/CN/density、overlap pseudo-replication、UI
把 halo誤畫成方向。下列 hard gates全數鎖定：

- exact PS × exact raw HP primary；family collapse只作 sensitivity；
- qname join一對一、零 collision/ambiguous；
- unit-level共同 CpG basis、invalid=NA；
- restricted strand × read_group permutation，不可交換即 NOT_EVALUABLE；
- 等 N、CpG rarefaction、5/8/10 reads與marker masking sensitivity；
- PERMANOVA R² + PERMDISP + BY；
- X只作 subcube；methyl永不改 topology；
- HCC1395 technical pair只算 biological n=1。

**Decision**：Conditional GO，decision lock=Y。可以發布
`robust association / local-cis-compatible / tag-dependent / confounded /
not-evaluable`；禁止發布 `subclone / ancestry / cause / true lineage`。

## Provenance

- Commit: `387a101e6a3292e0d7f230ba8d20271c7434972a`
- KB last verified: LongPhase-S 2026-07-11；BAM/InterSubMod 2026-07-12/11
- Audit time: 2026-07-27 09:37 +0800
