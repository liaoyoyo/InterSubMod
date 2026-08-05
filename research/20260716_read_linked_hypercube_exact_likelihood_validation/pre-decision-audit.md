<!--
建立時間: 2026-07-16 06:05 +08:00
目標: 在全量實作前稽核 read-linked region、partial-read symbolic groups、k>8 exact/certified solver 與 read-pattern likelihood 的可行性
處理範圍: chr1-22、7 datasets / 6 biological samples、LongPhase-S recalibrated FILTER=PASS canonical v5
cycle_id: cycle_20260716-0618-read-linked-hypercube-exact-likelihood
topic: 20260716_read_linked_hypercube_exact_likelihood_validation
status: verdict_PROBE_M2_PS_contract
audit_version: 0.3
關聯檔案:
  - InterSubMod/research/20260712_read_weighted_hypercube_tree/pre-decision-audit.md
  - InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
-->

# Pre-Decision Audit：read-linked Hypercube Steiner 全量驗證

> **2026-07-16 07:45 amendment — M2 full 暫回 PROBE**：獨立 red-team 證實 schema 1.1 雖保存 `phase_set`，但 HP component 與 ranking 未依 PS 分層；不同 PS block 的 HP1/HP2 orientation 可翻轉，不能合併為同一 haplotype。HCC1954 chr22 舊 pilot 實際有 50 個 PS，且 k>1 HP units 已出現跨 PS 混合。M0 與 symbolic solver 驗證仍有效；M2 full 必須先升級為 `HP family × known PS × read-linked component`、完成 pilot 與 receipt QA，才可重新升 GO。正在跑的 schema 1.1 chr22 v2只作診斷，不可作 final。

> **原始 verdict（已由上方 amendment 取代）：GO（90/100）**。計算方法 pilot 已通過，但這個GO未涵蓋後來發現的PS-block orientation問題；M2 full目前不得依此舊verdict啟動。  
> **服務目標**：G3（read-level evidence）、G4（跨樣本 reproducibility）、G5（可外部驗證的 solver/data contract）。

## §0 任務與 Cynefin Gate

- **Task Type**：B — Comprehensive validation。
- **固定 scope**：chr1-22、7 datasets / 6 biological samples，不可用 subset 取代最終交付。
- **Domain**：Complex。
- **判斷**：相同 read-pattern evidence 是否能穩定選出相同 vertex set 尚未在多樣本反覆驗證；CN、error、partial masks 與相同 vertex-set/different-edge 的不可辨識性會交互作用。因此先做 safe-to-fail pilot，再升 full cycle。

### 研究啟動五問

1. **Thread D？** 是 read-level / HP evidence 的延伸，但不重新宣稱甲基化可救 phasing。
2. **Thread B 撤回範圍？** 不在 production FP filter；不得把甲基化或 VAF 排序包裝成 truth。
3. **KDE-corrected？** 本輪樹輸入不依賴 KDE；若使用 CN/CovM，只接受已標 provenance 的 current data。
4. **VCF caller AF？** joint read likelihood 由 lossless molecule table計算；caller AF只作獨立/敏感度欄，不與同 reads重複計分。
5. **長計算／C++／搬移／NO-GO？** 是長計算與新 spec；不改 C++，不刪搬檔，先 Python research implementation與 exact pilot。

## §1 已知觀察

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| canonical v5 為 chr1-22、7 datasets、LongPhase-S recalibrated PASS | ✓ | L1 | `InterSubMod/docs/CURRENT_FOCUS.md` |
| 469,849 sSNV；51,815 multi-locus positional groups | ✓ | L1 | canonical v5 aggregate / current focus |
| 5,854 groups k>8；cap 排除 225,268 sites | ✓ | L1 | canonical v5 group census重算 |
| 現行 k>8 只保留 densest contiguous 8，非分塊或 overlapping windows | ✓ | L1 | frozen `sm_multilocus_combinations.py` code audit |
| 現行 partial pattern是 subcube coverage constraint，count不進 objective | ✓ | L1 | frozen solver code audit |
| 現行 solver輸出第一個 minimum-hidden layer的所有 feasible node sets與 parent choices | ✓ | L1 | frozen solver code audit + synthetic execution |
| read count可安全直接視為 evolution-edge support | ✗ | L2 | `InterSubMod/research/20260712_read_weighted_hypercube_tree/pre-decision-audit.md` |
| lossless molecule pattern table可由現行 aggregate完整還原 | ✗ | L1 | 現行 `MINREAD` censorship、full/partial重疊、counts丟失 code audit |
| symbolic partial constraint與 explicit compatible-state enumeration在 k<=8 等價 | ✓ | L1 | 本輪 9,840 patterns／2,015,538 state checks，0 mismatch |
| k=9-12 exact MILP可在本機求解且回報 bound/gap | ✓ | L1 | 本輪 12/12 synthetic cases PASS；最大 4,096 variables、gap=0 |
| 現行 solver與新MILP在可比較k<=8案例的minimum hidden／vertex set一致 | ✓ | L1 | 66/66 checks一致；9 legacy-capped案例依契約SKIP |
| snapshot likelihood不會虛假區分same-vertex/different-edge | ✓ | L1 | score差 `1.42e-14`；數值容差內為0 |
| read-linked region可在7 datasets全量重建且穩定 | □ | L5 | 需 BAM/sidecar全量執行與 sensitivity |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | group coverage、directed monotone hypercube、mixture likelihood均可形式化 |
| 觀察支撐 | 20 | canonical、exhaustive symbolic、legacy cross-check與k9-12 pilot均有SHA-bound收據 |
| 機制清晰度 | 20 | candidate generation、node-state likelihood與edge non-identifiability可分層 |
| 反例風險 | 10 | 三項預註冊反例均PASS；CN與bulk snapshot claim ceiling仍保留 |
| 所需資源 | 20 | pilot約9秒且輸出小；full run採排程與壓縮輸出，不與既有Task-B爭用 |
| **TOTAL** | **90 / 100** | **GO（工程驗證）；科學主張仍受claim ceiling限制** |

### Falsifier observable

若方法不成立，會看到至少一項：

1. symbolic partial solver與現行 explicit 2^u subcube在 k<=8 的 optimum/node-set不一致；
2. 相同 vertex set、只改 parent edges，likelihood分數不同；
3. k=9-12 solver無法提供一致的 incumbent/bound/gap或重跑結果不穩；
4. read-linked grouping對 bridge threshold 1/2/3/5極度敏感，且held-out/replicate無穩定區間；
5. 由 joint pattern table重算的 per-site REF/ALT denominator不能與來源 evidence精確對帳。

### Reality-test 三反例

1. `AAA` 100 reads仍不能區分3個mutation的6種order；若模型唯一化，FAIL。
2. `RRA + AAA + RAX`的 `RAX`只支持 `{RAR,RAA}` 群，不直接支持 `RAR->RAA`；若edge得到read count，FAIL。
3. 兩棵樹 vertex set相同但parent不同時，snapshot likelihood必須完全同分；若不同，FAIL。

## §3 Assumption Map

| Assumption | Importance | Known? | Quadrant / action |
|---|---|---|---|
| 一個 alignment segment等於一個獨立 molecule | High | No | 必須先定義QNAME/primary/supplementary去重 |
| X是未觀察而非REF/ALT/OTHER | High | No | 新 schema 必須拆 `uncovered/deletion/refskip/other` |
| 50 kb能代表read linkage | High | No | 改用cut bridge support；50 kb僅作sensitivity |
| min hidden是結構第一層objective | High | Yes | 保留 lexicographic first objective |
| joint read likelihood後可再加同read VAF | High | Yes=False | 禁止雙計數；VAF作posterior predictive diagnostic |
| snapshot reads可辨識parent edge | High | Yes=False | 相同vertex set全部標 EDGE_NONIDENTIFIABLE |
| SciPy/HiGHS status足以稱數學證明 | High | Partial | 只稱 solver-certified within tolerance；保存bound/gap/status |
| CN/LOH不影響raw VAF ancestry | High | No | CN-clean sensitivity或ABSTAIN_CN |

## §4 Quick Pilot（先過才升 GO）

1. **Symbolic equivalence**  
   輸入：現行 frozen solver golden cases + seeded k<=8 cases。  
   輸出：explicit vs symbolic compatible-set、minimum hidden、optimal vertex-set digest。  
   → 驗證：所有cases 0 mismatch。
2. **MILP cross-check**  
   輸入：現行 complete non-capped k<=8單位抽樣 + synthetic k=9-12。  
   輸出：objective、UB、LB、gap、status、runtime。  
   → 驗證：k<=8 objective/optimal-set一致；k=9-12至少每k 3個case `gap<=1e-9`。
3. **Likelihood identifiability**  
   輸入：synthetic full/partial pattern counts。  
   → 驗證：same-V/different-E score差=0；true vertex-set在adequate-depth simulation不劣於錯誤set；低depth允許tie/abstain。
4. **Resource gate**  
   → 驗證：不與既有469,849-site Task-B競爭；big7可用空間需保留安全餘裕，輸出採壓縮TSV/JSON且有預估上限。

### Pilot checkpoint

- 1-3全PASS且輸出估算可控：升 **GO for full engineering validation**。
- 任一 exact/edge-invariance FAIL：維持 PROBE並修spec，不啟動全量。
- lossless欄位缺失：允許結構層GO；likelihood科學結論維持BLOCKED。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---:|---:|---:|
| lossless `(molecule,HP,PS,mask,pattern,BQ,MAPQ)` table | High | >6h full | P0 |
| read-linked cut bridge全量counts | High | >6h full | P0 |
| exact solver golden cross-check | High | <1h pilot | P0 |
| k>12 exact-safe separator decomposition | High | 2-6h implementation + full runtime | P0 |
| CN/purity/multiplicity independent gate | High | data-dependent | P1 |
| external lineage/single-cell truth | High | unavailable | P2 / claim ceiling |
| full proof logging solver | Medium | environment unavailable | P2；先用HiGHS bound/gap |

## §6 Conflict Scan

| Prior conclusion | Relation | 處理 |
|---|---|---|
| raw read count / uniform inflow作edge weight為NO-GO | direct conflict | 不實作為主方法，只保留negative control |
| joint likelihood後再加同reads VAF會double count | dependent constraint | VAF只做derived diagnostic或read-disjoint evidence |
| bulk HP marginals不能唯一配對cellular clone | claim conflict | 結果只稱regional mutation-state candidate |
| methylation-only TP/FP與subclone claim多輪NEGATIVE | claim conflict | methylation只作orthogonal follow-up，不作tree edge evidence |
| canonical v5 capped units是不完整候選 | support | 新結果分開 `OPTIMAL_FULL_SCOPE` / `INCOMPLETE/CAPPED` |

## §7 Red Team 與決策

### 最強反方

1. **模型精確但問題錯了**：exact optimum只證明固定objective的最優，不能證明 biological ancestry。
2. **資料被censor後再精算**：若仍用MINREAD後aggregate，likelihood只是thresholded alignment-exposure分數，不是lossless molecule likelihood。
3. **相同state、不同edge不可辨識**：任何用snapshot VAF/read counts選唯一parent的做法都可能製造假確定性。

### Verdict

- **目前**：GO 90/100（2026-07-16 06:18 +08:00）。
- **Pilot實證**：symbolic exhaustive 2,015,538 checks零差異；legacy/MILP 66個可比較vertex-set checks零差異；k=9-12 12/12 exact cases PASS；likelihood三個negative/positive controls PASS。
- **允許**：新增獨立research scripts/tests、全7-dataset lossless extraction與read-linked sensitivity；必須保存terminal receipt與solver status/gap。
- **暫停**：canonical replacement、把bulk snapshot winner稱作真實clone／parent edge、在同一read likelihood後重複加入VAF。
- **Decision lock**：raw-edge-count主方法永久NO-GO；若有new longitudinal/single-cell/lineage data才可重開edge-identifiability。

### Pilot execution receipt

- **Input**：canonical v5 frozen solver、7 datasets的MLHP JSON、synthetic k=1-12 fixtures。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_exact_symbolic_pilot.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot --canonical-per-dataset 5`
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot/pilot_receipt.json`，SHA256 `9939a9aeb89758d17657876089faac20f0e4fb248dfdb2091b61fea54149d8b0`。
- **Actual**：exit=0、`all_pass=true`；7/7 datasets present；完整數值如上。

## Step -> Verify 執行契約

1. 凍結source/input schema → 驗證：manifest path/hash、sample數=7、chr=1-22。
2. 建symbolic group與MILP pilot → 驗證：exact cross-check、edge-invariance、gap/status。
3. 建lossless extractor → 驗證：QNAME去重與sitewise counts exact reconciliation。
4. 建read-linked regions → 驗證：threshold 1/2/3/5、held-out/replicate stability、50kb差異表。
5. 全量run → 驗證：每sample terminal receipt、無missing/duplicate、UB/LB/gap與abstain原因完整。
6. 報告 → 驗證：所有比例有total/relative denominator，HTML與source tables精確一致。

## Provenance

- Commit：`fd4c7c2cfbe8ffc0081cb7addba0313507379676`
- Canonical input：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`
- Runtime：Python + SciPy 1.13.1 (`scipy.optimize.milp` available)
- Current resource caveat：48 cores正在執行另一個全sSNV Task-B；big7 99% used（717G available at preflight）。
