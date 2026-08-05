<!--
建立時間: 2026-07-19 05:38 +08:00
目標: LongLineage C++/HTSlib 全流程重建進入實作前的 pre-decision audit
處理範圍: P0-P8；engineering GO 與 scientific phase gates
cycle_id: cycle_20260719-longlineage-cpp-rebuild
topic: longlineage_cpp_rebuild
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/00_INDEX.md
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/implementation-notes.md
-->

# Pre-Decision Audit: LongLineage C++／HTSlib 全流程重建

> **Verdict**: `GO` 建立獨立 repo 與逐 phase 驗證；P3-P7 不繼承 GO，必須各自通過 parity／scale gate。

## §0 Cynefin Domain Gate

- **Domain**: Complex。
- **Test**: 相同方法尚未曾以獨立 C++ 全量重建；M2 v2 real pilot曾 timeout。
- **決定**: 採 phase-bounded safe-to-fail probes；禁止由 P0/P1 scaffold 推論 P3-P7 已可用。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 7/7 latest HP/PS sidecar closeout PASS，164,253,537 mapped alignments | ✓ | L1 | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/raw_all_receipt_closeout.json` |
| topology Python authority bytes符合鎖定 SHA | ✓ | L1 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`，2026-07-19 `sha256sum` |
| M1/M2 census與 source-attested receipts存在 | ✓ | L1/L2 | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/` |
| exact-preserving remediation 341/341 tests PASS | ✓ | L1 | `InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/20260718_Hypercube_exact_preserving修補與驗證結果_01.md` |
| v2 M2 real pilot 8h timeout、partial gzip、正式 NO-GO | ✗ | L1 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md` |
| 正式 full co-occurrence C++/Python authority | □ | — | 尚不存在；本 cycle 要建立第一份雙 C++ authority |
| 7-dataset C++ 24/40 worker parity | □ | — | P7 待執行 |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | typed contracts、deterministic kernels、independent validation有明確工程基礎 |
| 觀察支撐 | 10 | 7-dataset upstream與Python synthetic evidence完整，但沒有 C++ port evidence |
| 機制清晰度 | 20 | production/evaluation/presentation與 producer/validator 邊界已鎖 |
| 反例風險 | 10 | 已知 RNG、MM/ML、identity、solver scale反例，仍待 C++ 實測 |
| 所需資源 | 0 | 明顯超過 6 小時且 P7 為長計算 |
| **TOTAL** | **60 / 100** | **GO，附 phase locks** |

**Falsifier observable**：若「可忠實 C++ 重建」為錯，frozen vectors會出現 read-set、partition、
permutation exceedance、status、objective、family digest或parent count任一非零 mismatch；該 phase
必須保持 `BLOCKED`，不得寫 validated receipt。

**Reality-test 三反例**：

1. PCG64/RNG consumption、tie-break或浮點次序讓 decision/status漂移。
2. raw BAM projection在 RG-only duplicate以外出現 typed aux差異，exact join必須 fail closed。
3. topology family在 cap/deadline前未列完；只能標 incomplete，不能產生 winner。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| frozen source與goldens足以定義decision parity | HIGH | PARTIAL | (2) MUST validate |
| local HTSlib可供開發；production由1.18 container鎖定 | HIGH | UNKNOWN | (2) MUST validate |
| latest sidecar可對所有production reads exact join | HIGH | upstream已驗、C++未驗 | (2) MUST validate |
| direct HiGHS pin可重建SciPy route | HIGH | UNKNOWN | (2) MUST validate |
| Python presentation可完全只讀chart-ready資料 | HIGH | KNOWN by design | (1) |
| 效能足以完成7-dataset | HIGH | UNKNOWN | (2) MUST validate |

## §4 Quick Pilot / Phase Gate

1. 建立 synthetic BAM/VCF/sidecar＋malformed fixtures。
2. 以 C++ preflight、typed parser、validator執行 positive/negative tests。
3. 以 1/2/4 workers及兩種chunk順序產生相同 semantic SHA。

Checkpoint：

- schema、negative gates、semantic SHA全 PASS → P0-P2 GO。
- 任一 scientific parity mismatch → 對應 P3-P5保持 PROBE/BLOCKED。
- 任何 truth leakage或 forged PASS可通過 → 全 cycle NO-GO，先修治理。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| C++ M1 frozen vectors與PCG64 consumption port | HIGH | >6h | P0 |
| formal exact K×2 test與global FDR parity | HIGH | >6h | P0 |
| direct HiGHS pinned integration | HIGH | >6h | P1 |
| full co-occurrence first authority | HIGH | 長計算 | P1 |
| 7-dataset 24/40 workers雙跑 | HIGH | 長計算 | P1 |
| orthogonal cellular lineage truth | HIGH | external | P2 / claim ceiling |
| repo根 `MEMORY.md` 不存在，與既有 cold-start說明不一致 | MED | <1h | P0（新 repo修正） |

## §6 Evidence Conflict Scan

- `InterSubMod/MEMORY.md`：實際不存在；此 governance gap列入新 repo cold-start test。
- `find InterSubMod/docs/reports/validated -name "*NEGATIVE*"`：檔名掃描未提供本題可取代 current
  M2 NO-GO 的正式正面證據。

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| M2 frozen v2 pilot timeout／partial gzip／NO-GO | L1 | conflict：禁止宣稱scale已解決 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md` |
| exact-preserving remediation只 PASS_FOR_PROBE | L1 | dependent：C++ port需重驗 | `InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/20260718_Hypercube_exact_preserving修補與驗證結果_01.md` |
| cellular clone/ancestry claim仍 NO-GO | L1/L2 | support：採mutation-state family ceiling | `InterSubMod/docs/experiments/INDEX.md` |

## Red-team gate

1. **最強反方 A**：12k+ lines Python scientific authority的 C++ port容易在隱性 RNG/tie/precedence上漂移。
   防線是逐decision frozen vectors與「任何 status漂移即 FAIL」。
2. **最強反方 B**：建立漂亮 repo與CI可能製造「已完成」錯覺。
   防線是 phase-state enum、release gate與 production kernel未驗證時 hard stop。
3. **最強反方 C**：full scope資源超標重演 v2 timeout。
   防線是 resource pilot、queue byte cap、deadline不產winner、24/40雙跑前置 gate。

Falsifier可觀察且能阻止 release；未與已 concluded-dead FP-filter方向重疊。紅隊通過，但要求 phase locks。

## §7 Decision

- **Score**: 60/100。
- **Verdict**: **GO** 啟動工程與逐 phase full verification。
- **Decision lock**:
  - production禁止 truth。
  - Python禁止科學計算。
  - 未通過 P3/P4/P5 parity不得進 P7。
  - cap/deadline/incomplete不得產winner。
  - 正式主張不得超過 lineage-compatible mutation-state families。

## Provenance

- InterSubMod commit at audit: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`（worktree dirty，僅作基準）。
- Audit timestamp: `2026-07-19T05:38:30+08:00`。
- Skill: `/pre-decision-audit` v0.1。
- Machine record: `InterSubMod/state/cycles/cycle_20260719-longlineage-cpp-rebuild/audit.json`。

