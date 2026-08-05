<!--
建立時間: 2026-07-10 22:44
目標: layered reconstruction v2 大改前 pre-decision audit
處理範圍: 核心 schema、7 dataset 重跑、驗證與 claim 邊界
cycle_id: 20260710-2244-layered-reconstruction-v2
topic: 20260710_layered_reconstruction_v2
status: verdict_GO
audit_version: 0.1
build_branch: master
build_commit: 4fb9e742482b63a660de19a1f1bd07d49d713111
worktree: /big7_disk/liaoyoyo2001/InterSubMod
-->

# Pre-Decision Audit: Layered Reconstruction v2

> **Verdict: GO (70/100)**。P0 根因已由 7 dataset census 與程式路徑交叉重現；本輪是修正 canonical semantics，不是重啟已 NO-GO 的 methylation filter。

## §0 Cynefin Domain Gate

- **Domain**: Complicated。
- **Test**: 相同的 denominator／schema 稽核可重複得到一致結果。
- **Rationale**: root-only、CN fallback、verification sampling 與 tree storage 都有明確 deterministic code path；外部 single-cell confirmation 屬另一個 data-availability gate。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---|---|
| 88,358 混合 out-of-scope、cap-excluded、singleton 與 unsupported | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/CURRENT_FOCUS.md` 2026-07-10 funnel entry |
| root-only 單位進入 determined 與 multi-HP 分母 | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py:136` |
| 無 CN bed 時 `cn_state()` 回 neutral | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py:45` |
| runner 設 `SM_VERIFY_EVERY=20` | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh:46` |
| exact shape 與 CCF 使用 stored first 32 | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py:149`；`ccf_and_cn_multisample.py:124` |
| 18,103 非-capped units full V4/V5、4,000 stress 均 0 mismatch | ✓ | L2 ⭐⭐⭐⭐ | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/full_v4v5_verification.py` |
| 7 dataset 修正版 canonical outputs | □ | — | 本 cycle 待產生 |
| single-cell／multi-region orthogonal truth | □ | — | repo 尚無對應輸入資料 |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | lineage、reference control、missingness 與 verification completeness 定義明確 |
| 觀察支撐 | 20 | 7 datasets 同方向，且 HCC1395 有精確重算 |
| 機制清晰度 | 20 | 每個錯誤均可定位至單一 deterministic branch |
| 反例風險 | 10 | legacy consumer、候選樹輸出大小與 CN bed 語意仍需 regression |
| 所需資源 | 0 | 全量 BAM 重讀與 7 dataset full verification 預估超過 6 小時 |
| **TOTAL** | **70/100** | **GO** |

**Falsifier observable**：修正版輸出若仍存在 reference-only／HP3 進 primary denominator、missing CN 產 candidate、non-capped eligible unit 的 V4/V5 skipped 卻標 all-pass、funnel 不守恆，或完整枚舉單位的 analysis candidate count 小於 `n_trees`，則本次修正失敗。Capped units 必須另列 not-applicable，禁止進 exact claim。

**Reality-test 三反例**：

1. 保留 legacy `is_lineage` 的 consumer 可能默默沿用錯分母，因此需要新欄位與 regression 同時更新。
2. 全存候選樹可能讓 JSON／UI 成本不可接受，因此 exact analysis 與 display storage 必須分離且明示完整性。
3. CN bed 若不是「未列區域即 neutral」的完整 segmentation，neutral 判定仍可能是假 missingness；manifest 必須記 source semantics。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| 只有含 ALT 的 HP1/HP2 family 是 primary mutation lineage | HIGH | KNOWN | (1) verify quickly |
| HP3 不是第三 parental lineage | HIGH | KNOWN | (1) verify quickly |
| COLO829/HCC1937 無可靠 CN 時不可推 neutral | HIGH | KNOWN | (1) verify quickly |
| SAVANA smcnbed 未命中區域代表 neutral | HIGH | UNKNOWN | (2) validate first |
| 完整 candidate trees 可在合理儲存成本內保存 | MED | UNKNOWN | (2) validate before publish artifact |
| single-cell truth 可由現有 bulk reads替代 | HIGH | KNOWN FALSE | (1) claim ceiling |

右上象限 gate：先檢查 CN bed coverage／label schema；量測完整候選樹輸出大小與 runtime，再決定 display storage 上限。分析本身不得截斷。

## §4 Quick Pilot

GO 不要求 probe，但實作採 HCC1395 synthetic／existing JSON smoke gate：

1. 以 synthetic root-only、HP3、missing-CN cases 跑 unit tests。
2. 以 HCC1395 既有 mlhp 重建 layered JSON，不重讀 BAM。
3. 檢查 exact-zero thresholds、shape completeness、runtime 與檔案大小；通過後才跑 7 datasets。

中止條件：任何 schema threshold 非零，或單樣本輸出／runtime 超過既有資源兩倍，先修正再 generalize。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| canonical driver 未輸出 reference/control 與 primary lineage 分層 | HIGH | 2-4h | P0 |
| CN missingness 未進 upstream schema | HIGH | 1-2h | P0 |
| runner 缺 manifest、stage status、完整 part 檢查與 checksums | HIGH | 2-4h | P0 |
| shape/read-AF analysis 截斷於 32 trees | HIGH | 2-4h | P1 |
| MAPQ/BaseQ/MINREAD/MAX_SNV sensitivity | MED | 4-8h | P1 |
| COLO829/HCC1937 外部 CN truth | MED | external | P2 |
| single-cell／multi-region orthogonal data | HIGH for clone confirmation | external | P2 |

## §6 Evidence Conflict Scan

`MEMORY.md` 不存在於目前 worktree；以 live `CURRENT_FOCUS`、validated reports 與 code 為替代 conflict sources。

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| methylation TP/FP filter space exhausted | concluded | support：L3 必須 bounded auxiliary | `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md` |
| 7-sample funnel／full V4V5 已旁路驗證 | L2 | support，但不代表 core schema 已修 | `InterSubMod/docs/CURRENT_FOCUS.md` 2026-07-10 |
| CCF winner-clean 不是獨立驗證 | L2 | support：降級為 read-AF ordering heuristic | `InterSubMod/docs/CURRENT_FOCUS.md` 2026-07-10 |
| CURRENT_FOCUS 稱 P0 已修，但 core scripts 保留舊 branch | L1 | conflict：文件狀態超前於執行語意 | 本 audit §1 code references |

**Conflict count: 1**，故反例風險維度只給 10 分。

## §7 Decision Path

- **TOTAL**: 70/100。
- **Verdict**: **GO**。
- **Next action**: core semantics → HCC1395 smoke → 7 dataset rerun → independent verify → docs/ledger sync。
- **Decision lock**: Y。只有新資料、新方法或 upstream 定義改變才可讓 root-only／HP3／missing-CN 重新進 primary claim。

### Independent red team

- Failure mode 1：只新增欄位但舊下游仍讀 legacy denominator。對策：更新 region view、verifier、reports，並用 forbidden-count gate 阻擋。
- Failure mode 2：為求完整而全存 trees，造成 output 過大。對策：analysis 一律全候選；display storage 可設 cap，但必有 exact summary、digest 與 `analysis_complete=true`。
- Failure mode 3：把 read allele fraction 改名後仍暗示真 CCF。對策：schema 名稱、文件 claim 與 UI 都使用 `read_af_ordering`；只有 purity/CN-corrected 後才可稱 CCF。
- 與既有 NO-GO 重疊：否。本輪不重啟 methylation filtering，只把 L3 限制寫進 contract。

## Provenance

- Commit: `4fb9e742482b63a660de19a1f1bd07d49d713111`
- Cycle: `InterSubMod/state/cycles/20260710-2244-layered-reconstruction-v2/`
- Audit JSON: `InterSubMod/state/cycles/20260710-2244-layered-reconstruction-v2/audit.json`
- Skill: `/pre-decision-audit` v0.1

## §8 Gate-1 addendum：clean production tagging 與 full-run launch（2026-07-11 04:02 +08:00）

> **Verdict: PROBE（60/100，紅隊後由 GO 降級）**。無 truth VCF/BED 的 LongPhase-S production tagging 本身有明確機制；但 7 份實體 tagged BAM 預估約 1.94 TB，正式 big7 只剩 932 GB，因此外部 AI 改採原 BAM＋完整 alignment-level HP/PS sidecar。此 consumer join 尚未驗證，必須 probe-first。

### §8.1 Cynefin 與觀察

- **Domain**: Complex（production LongPhase-S=Complicated；sidecar consumer integration=Complex）。
- `[O-L1]` 6/7 historical tagged BAM 的 `@PG` 含 `--truth-bed`；ClairS PASS universe 則為 genome-wide。
- `[O-L1]` 外部 session 於 `2026-07-11 03:56:31 +08:00` 啟動 run `20260711_longphase_s_production_sidecars_v1`；HCC1395 與 HCC1395_DORADO 命令的 truth 旗標為空，`-q 20 --tagSupplementary` 明示。
- `[O-L1]` 7 個 tumor BAM 合計 `1,939,360,284,288` bytes；`/big7_disk` 剩 932 GB，無法安全存放完整副本。
- `[F-L1]` current layered consumer 仍直接從 tagged BAM pileup 讀 base+HP/PS；未支援原 BAM＋sidecar join。

### §8.2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | HP/PS 是 BAM tag，可以 exact alignment identity 外置 |
| 觀察支撐 | 20 | truth-scope mismatch、實際 commands、disk 容量與進程均已直接量測 |
| 機制清晰度 | 10 | capture schema 存在，但 consumer exact join/衝突 policy 未實作驗證 |
| 反例風險 | 10 | supplementary、duplicate identity、tabix boundary、sidecar與 BAM drift 均可改變結果 |
| 所需資源 | 0 | 7-dataset full tagging 明顯 >6h |
| **TOTAL** | **60/100** | 原始門檻可 GO，但紅隊因 exact join falsifier 未通過降為 **PROBE** |

**Falsifier observable**：同一 bounded region 中，歷史／臨時 clean tagged BAM 直讀的 alignment identity→HP/PS，若與「原 BAM＋sidecar exact join」在 eligible pileup observations、R/A/X vectors、family counts 或 output digest 任一不同，sidecar 不得成為 production input。

### §8.3 Reality tests 與假設地圖

1. FIFO capture 若中斷，LongPhase-S 或 capture 可能單邊成功；兩邊 exit、row count、tag count、VCF key conservation 必須全過。
2. supplementary 或重複 QNAME 若只以 qname join，可能把錯誤 alignment 的 HP/PS 套入。
3. sidecar 若未鎖原 BAM header/content identity，即使格式正確仍可對到錯版 BAM。

HIGH×UNKNOWN 必先驗：exact alignment key、supplementary conflict policy、sidecar tabix query boundary、原 BAM＋sidecar 與 clean tagged BAM 的 output digest 對等性。

### §8.4 Probe checkpoint

1. 等 HCC1395/HCC1395_DORADO streaming capture 至少一份完成，驗 `truth_flags_absent`、all input VCF keys preserved、alignment/tag counts exact。
2. 以 chr1 或一個預註冊 region 建立臨時 clean tagged BAM 對照，比較 direct BAM vs sidecar join 的 observation/family digest。
3. 故意篡改 sidecar hash、BAM identity、HP tag 與 manifest scope，validator/runner 必須 non-zero exit，且不建立正式 run root。

**Checkpoint**：上述 exact-equivalence 與四類 adversarial tests 全過 → GO Gate-1；任一 digest 不同 → 停在 PROBE，不啟動 layered full run。

### §8.5 紅隊結論

- 最強反方 1：sidecar 是「省空間的新輸入模式」，不是 BAM 的自明等價物。
- 最強反方 2：current pipeline 的 qname overwrite 本來就有 supplementary ambiguity；sidecar join 不可把此 legacy bug 偽裝成對等。
- 最強反方 3：外部 AI 已先啟動 2-sample 長跑，但「正在跑」不能取代 preflight 或 consumer acceptance。

**Decision lock**: N；目前只核准 sidecar capture 作 safe evidence probe，未核准將它宣稱為 production tagged-BAM replacement。

## §9 Gate-3 addendum：raw ClairS all-site LongPhase-S contract（2026-07-11 06:30 +08:00）

> **Verdict: PROBE（70/100，紅隊後由 GO 降級）**。canonical tree 使用 LongPhase-S `_sc.vcf FILTER=PASS` 仍合理；但若 LongPhase-S input 先裁成 ClairS PASS，原生 `non-PASS → PASS` rescue branch 永遠不可達，不能回答「所有 ClairS 位點是否完整讀取」。現行 PASS-only producer 已在 0/7 完成、無 `_SUCCESS` 時中止並保留為 failed artifact。

### §9.1 Cynefin 與 observation completeness

- **Domain**: Complex，強制 probe-first。raw-all 是官方功能，但本 repo 有具體 crash history，不能直接視為 routine rerun。
- `[O-L1]` `HaplotagVcfParser.cpp:589-598` 明確實作 input non-PASS→PASS rescue 與 input PASS→LowQual remove。
- `[O-L1]` LongPhase-S `docs/somatic_haplotag.md` 說 output 基於 input tumor VCF 重標 PASS/LowQual，其他欄位保持。
- `[O-L1]` 7/7 raw ClairS callsets 全為 biallelic SNV；raw PASS 只是其子集，非 PASS 共 70,179 records。
- `[O-L1]` repo wrapper `scripts/pipeline/steps/01_longphase_s.sh:124-125` 記錄曾因 LowQual `chrX:72880028 QUAL=0` crash；目前 HCC1395 raw VCF 仍含該 record。
- `[O-L1]` 被中止的 PASS-only run 只啟動 HCC1395/HCC1395_DORADO，0/7 terminal PASS，`_FAILED` 存在且 `_SUCCESS` 不存在。

### §9.2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | all-site input 才能完整啟用 LongPhase-S 雙向 recalibration |
| 觀察支撐 | 10 | source/官方文件與 raw records 已核對，但現行 binary raw-all 尚未成功執行 |
| 機制清晰度 | 20 | raw all→normalized all→`_sc` all→`_sc PASS`→tree，transition matrix 可逐 key 驗證 |
| 反例風險 | 0 | 已有具體 crash history；QUAL=0、GT=1/1、零 ref support 或 header typing 都可能觸發 |
| 所需資源 | 20 | chr22 與單一 chrX locus probe 可在一小時內完成 |
| **TOTAL** | **70/100** | 原門檻為 GO；紅隊因已知 crash falsifier 未解而降為 **PROBE** |

**Falsifier observable**：normalized raw-all bounded input 若 current binary non-zero/crash、input record-key multiset不等於 `_sc` all、FILTER transition無法逐 key唯一對應、sidecar capture不守恆，或 known chrX crash locus仍使執行失敗，則不得推廣到 7 datasets。

### §9.3 Assumption map 與 reality tests

| Assumption | Importance | Known | Gate |
|---|---|---|---|
| raw non-PASS 必須進 LongPhase-S 才可能 rescue | HIGH | KNOWN | source verified |
| header-only normalization足以避免歷史 crash | HIGH | UNKNOWN | MUST probe |
| raw-all 不會讓 tag assignment／purity機制失真 | HIGH | UNKNOWN | compare PASS-only bounded output |
| `_sc` 能對 raw-all 每個 key只輸出一次 | HIGH | UNKNOWN | exact multiset + uniqueness gate |

Reality tests：
1. HCC1954 chr22 normalized raw-all 必須 426 input→426 `_sc` all，另報四格 FILTER transition matrix與 HP/PS exact join。
2. HCC1395 `chrX:72880028` targeted probe 必須 current binary正常退出；若需排除此位點，便不再是 all-site contract，應明示 NO-GO 或修 upstream/tool。
3. raw-all 與 PASS-only bounded結果需比較 shared PASS site FILTER、tag census及新增 rescued sites；差異不得被寫成獨立 truth improvement。

### §9.4 Decision path

- **Verdict**: **PROBE**。
- **Next action**: 先做 HCC1954 chr22 raw-all；通過後做 known-crash chrX targeted probe；兩者通過才重建 raw-all manifest/producer/receipt/closeout contracts並跑 7 datasets。
- **Decision lock**: N。7-sample production 在兩個 bounded probes 前禁止啟動。
- **最強反方**: raw-all 理論上完整但 current binary 對 caller LowQual 極端狀態不健壯；強行「全部餵入」可能把方法完整性換成 runtime failure，必須由實際 probe 而非文件裁決。

### §9.5 Probe checkpoint result（2026-07-11 07:00 +08:00）

- HCC1954 chr22、original binary：426→426，rescue 53、remove 14，217,023 sidecar rows exact。
- HCC1395 chrX:72880028、original binary：如預註冊反例，`tumorPosReadCorrBaseHP` missing，exit 1。
- copied-source最小 patch：只把 low-confidence/no-eligible-read的 fatal改為 warning+no-op；target locus 1→1、LowQual→LowQual、64 rows全 HP=.、exit 0。
- HCC1954 chr22 patched regression：sidecar bytes、426 VCF record lines、purity與四格 transition matrix皆與 original binary完全相同。
- HCC1395 whole chrX patched：35,823→35,823；LowQual→PASS 959、PASS→LowQual 1,616；2 no-read warnings；1,090,570 sidecar rows、unknown/conflict=0；artifact hashes全PASS。

**Re-score for full launch**：theory 20 + observation 10 + mechanism 20 + counter-risk 10 + full-run resource 0 = **60/100**。紅隊通過，Verdict升為 **GO_WITH_FAIL_CLOSED**：可實作7-dataset raw-all producer與啟動全量；任一新 fatal、key loss、payload drift、unknown/conflict或hash gate仍須停止且不得建立 `_SUCCESS`。
