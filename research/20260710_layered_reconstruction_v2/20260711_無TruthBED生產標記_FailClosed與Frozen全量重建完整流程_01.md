<!--
建立時間: 2026-07-11 05:04 +08:00
目標: 完整說明無 truth VCF/BED 的 LongPhase-S production tagging、HP/PS sidecar、fail-closed validator、frozen provenance 與 7-dataset layered full run 的判斷、設定、命令、驗證與限制
處理範圍: ⚠ PARTIAL EXECUTION — 方法與程式稽核涵蓋全流程；active producer 尚為 0/7 complete，新的 layered 7-dataset run 尚未 launch
完整 scope 未驗證: 7 production sidecars × chr1-22 runtime；35/35 layered splits；final v3 receipt/_SUCCESS
補強計劃: producer 7/7 → v3 frozen lock → resource gate → 7-dataset full run → verifier/evidence ledger → 更新本檔與 HTML
關聯檔案:
  - InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/04_clean_tagging_launch_audit.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/05_fail_closed_wiring_audit.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/06_launch_resource_redteam.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/07_layered_v3_lifecycle_contract.md
task_type: B Comprehensive validation
goals: G4,G5（次要 G1,G3）
framework: Hybrid-PI-Report L0-L3 + A3/ADR + Diátaxis；launch decision 使用 WRAP + Pre-Mortem
cycle_id: 20260710-2244-layered-reconstruction-v2
build_branch: research/subclonal-reconstruction-202606
build_commit: 4fb9e742482b63a660de19a1f1bd07d49d713111
worktree: /big7_disk/liaoyoyo2001/InterSubMod
data_version: active production run 20260711_longphase_s_production_sidecars_v1
status: in_progress / NO-LAUNCH layered full run
證據等級: L2 ⭐⭐⭐⭐（多路 code/audit/synthetic + bounded real probe；production 7/7 pending）
-->

# 無 Truth BED Production Tagging、Fail-Closed 與 Frozen 全量重建完整流程

> ⚠ **PARTIAL EXECUTION / NO-LAUNCH**：完整方法、角色、程式與 gate 已審查；截至本版截點，production sidecar 仍為 0/7 complete，新的 layered full run 未啟動。本文件不把 `START`、檔案存在或單一 `pass:true` 寫成完成。

## 0. 目前正確動作是等待 producer，不能重複啟動

現有 LongPhase-S production run 已在執行 HCC1395 與 HCC1395_DORADO。命令沒有 `--truth-vcf` 或 `--truth-bed`，但保留 LongPhase-S production default internal filtering，並輸出 recalibrated VCF。第二個相同 producer 或重疊 layered full run 都是 NO-GO；原因是 I/O 已飽和且 downstream artifacts 尚未完成。

| Gate | 截點狀態 | Verdict | 證據 |
|---|---:|---|---|
| production truth flags | active commands 0 個 truth flags | 方向正確；完成度 pending | ⭐⭐⭐⭐ L2 process observation |
| production completion | 0/7 PASS | ✗ 不可 launch layered | ⭐⭐⭐⭐ L2 run-status readback |
| sidecar join method | synthetic 3/3；COLO829 bounded 35/35 | ✓ bounded PASS | ⭐⭐⭐⭐ L2 |
| strict v3 input contract | 17/17 tests；receipt normalizer 8/8 | ✓ component PASS | ⭐⭐⭐⭐⭐ L1 synthetic mechanism |
| atomic lifecycle | 13/13 tests | ✓ component PASS | ⭐⭐⭐⭐⭐ L1 synthetic mechanism |
| frozen runner / verifier | runner 12/12；verifier 16/16 | ✓ integrated component PASS；real inputs pending | ⭐⭐⭐⭐⭐ L1 synthetic mechanism |
| physical 7 BAM capacity | 1.674 TiB payload vs 932 GiB free | ✗ HOLD | ⭐⭐⭐⭐ L2 measured size/capacity |

Effect size / CI：本表是工程 invariant 與容量 gate，不是效果估計；Cohen ribbon 與 confidence interval 不適用。所有數字以 exact count、bytes 或 checksum驗證。

## 1. 任務範圍固定為 7 datasets、6 biological samples、chr1–22

- Datasets：HCC1395、HCC1395_DORADO、COLO829、H1437、H2009、HCC1937、HCC1954。
- Biological samples：6；HCC1395 與 HCC1395_DORADO 是同一細胞株的兩個資料處理版本。
- Primary genome scope：chr1–22；chrX/chrY只可列 out-of-scope census。
- Task type：B Comprehensive validation；subset不能冒充 final evidence。
- Research goals：G4 reproducibility、G5 external auditability；HP/read-level evidence次要服務 G1/G3。

## 2. 四種 VCF 與兩種 BAM 角色必須分開命名

| Role | 定義 | 是否進 tree | 不可混淆點 |
|---|---|---:|---|
| `caller_raw_vcf` | raw ClairS 全 records | 否；進 lossless site ledger | 可能含 non-PASS |
| `longphase_input_vcf` | ClairS FILTER=PASS，供 LongPhase-S | 否；是 producer input | 不等於 recalibrated PASS |
| `longphase_recalibrated_all_vcf` | LongPhase-S `_sc.vcf` 全 keys，含 PASS/LowQual | 否；作完整 recalibration ledger | key-preserving不等於tag-all |
| `longphase_recalibrated_pass_vcf` / tree input | `_sc.vcf` FILTER=PASS | 是；canonical primary sSNV universe | 2026-07-11 04:50 最新決策 |
| raw tumor BAM | base/quality/coordinates authority | 是 | embedded historical HP/PS必 ignore |
| production tag sidecar | HP/PS authority | 是 | 不是 persistent BAM，不含 base/quality |

### 2.1 Production policy 的精確語意

目前命令沒有 `--disableFilter`，所以 policy 名稱固定為 `production_default_filter`：

- ClairS PASS 全部交給 LongPhase-S parser。
- LongPhase-S 仍執行內部 somatic filtering；高信心位點參與 tagging。
- `--output-somatic-vcf` 保存所有 input record keys並標 recalibrated FILTER。
- 只有 recalibrated FILTER=PASS 進 canonical tree backbone。
- 若未來要「所有 ClairS PASS 都作 tagging seed」，那是另一個 `caller_pass_tag_all` role，必須加 `--disableFilter` 並另跑；不可事後改名。

## 3. 大型 workspace 只沿索引與 manifest 導航

本任務禁止預設對 `/big7_disk`、`/big8_disk` 做 broad recursive search。查詢順序固定：

1. `InterSubMod/docs/CURRENT_FOCUS.md` 與研究 `00_INDEX.md`。
2. active manifest、`run_context.json`、`upstream_dependency.tsv`、run-root status/hash ledger。
3. 依 manifest 精確開啟 producer/consumer/validator檔。
4. 只在單檔用 bounded `rg pattern file` 或精確 `stat/sha256sum/jq`。

這個策略同時避免 NFS metadata storm、archive誤判為current、以及同一資料在不同run被混用。

## 4. Historical tagged BAM 不能作 clean production authority

歷史 7 個 tagged BAM 中，6/7 LongPhase-S `@PG` 含 `--truth-bed`，7/7含 historical truth VCF。Source audit顯示 benchmark BED會在 `tagRead` 前移除BED外 somatic variants；因此 historical embedded HP/PS只適合作 engineering baseline，不能證明 genome-wide clean tagging。

因果鏈：

```text
truth BED/VCF flags
  → 可參與 somatic tagging 的 variant universe 改變
  → read HP/PS assignment 改變
  → per-read R/A/X 與 germline-family evidence 改變
  → regional tree candidate/support 改變
```

主要 confounders：LongPhase internal filter policy、raw BAM identity、sidecar join collision、mutable manifest/source、I/O overlap。Mediator是HP/PS evidence；本報告不把 run completion status 當 biological outcome，也不對 collider作因果調整。

## 5. Production LongPhase-S 命令的每個參數都有固定角色

Canonical command shape：

```bash
longphase-s somatic_haplotag \
  -s <germline_phased.vcf.gz> \
  -b <normal.bam> \
  --tumor-snv-file <ClairS.PASS.vcf.gz> \
  --tumor-bam-file <raw_tumor.bam> \
  -r <GRCh38.fasta> \
  -t 12 --tagSupplementary -q 20 \
  --output-somatic-vcf \
  -o <production_prefix>
```

| Option | Value/origin | 基本判斷 |
|---|---|---|
| `-s` | explicit | phased germline scaffold；需VCF index與reference dictionary一致 |
| `-b` | explicit | matched normal alignment |
| `--tumor-snv-file` | explicit | ClairS paired FILTER=PASS；不是truth |
| `--tumor-bam-file` | explicit | raw tumor alignment payload |
| `-r` | explicit | GRCh38 no-alt reference |
| `-t` | 12 explicit | 每sample threads；外層最多2 samples並行 |
| `--tagSupplementary` | true explicit | supplementary alignment也輸出HP/PS |
| `-q` | 20 explicit | MAPQ基本門檻 |
| `-p` | 0.6 default | read assignment threshold；必在effective-options鎖定default origin |
| `--output-somatic-vcf` | true explicit | 產recalibrated `_sc.vcf` |
| `--truth-vcf/--truth-bed` | forbidden/absent | 任一 token存在即fail |
| `--disableFilter` | absent | internal filtering enabled；policy=`production_default_filter` |

## 6. FIFO sidecar 是容量受限下的可稽核替代，不是實體 BAM

7 個 historical tagged BAM合計約 1.674 TiB（1.841 TB decimal）；以20%安全餘裕，physical plan需要約2.3 TB。big7只餘932 GiB，因此不能安全保留7份新BAM/BAI。

既有 producer將 LongPhase輸出的 BAM寫入named pipe；capture consumer逐mapped alignment保存：

```text
CHROM START0 END0 QNAME FLAG MAPQ CIGAR_B2 HP PS
```

保存raw TSV、bgzip TSV與tabix index；BAM payload不落地。預估全7 datasets約20–40 GiB，60 GiB為hard planning cap。

### 6.1 coordinate_join_v1 的安全條件

Alignment key：`QNAME + RNAME + START0 + END0 + FLAG + BLAKE2b8(CIGAR)`。

此v1 key不含SEQ/QUAL digest，故只能在以下全部成立時使用：

1. Producer raw tumor BAM以 `storage_identity_v1` 綁定，不能換版。
2. Sidecar與index做full SHA-256。
3. 全域 duplicate rows=0、conflict=0、malformed=0、unknown HP=0。
4. Runtime每split exact matches=alignment exposures，missing/conflict=0。
5. Embedded BAM HP/PS fallback明確 forbidden。
6. 報告標示 lower assurance，不稱為 `alignment_identity_v2` 或 full payload identity。

### 6.2 已完成的等價證據

- Synthetic 3/3：同QNAME primary/supplementary分離、缺列檢出、衝突列檢出。
- COLO829 bounded real probe：chr1兩個sSNV、35 alignment exposures；direct與sidecar calls/HP/PS/family digest完全相同。
- Combined payload SHA-256：`b466e584ab820448769bd1aede31434021f3c6ce1762788afc7363c0fe2d71d7`（兩路相同）。

限制：real probe使用historical tags，只證join方法，不證production clean tags；真實supplementary由synthetic fixture補，不是full-scope evidence。

## 7. Fail-closed 的定義是「任何未知都停止」

至少要拒絕：

- 0 datasets造成 `all([])=true`。
- 不是exact 7 datasets / 6 biological IDs、duplicate sample、extra sample dir。
- schema 2.0/2.1 mixed mode或任一sidecar欄缺失。
- truth option token、producer argv/input/hash不符。
- sidecar/index/validation hash漂移或validation subject不綁BAM。
- specified tabix index未實際使用。
- source/env/manifest在preflight後改變。
- worker任一non-zero、heartbeat stale、PID reuse、非法state transition。
- output manifest、五parts、funnel、site ledger、region/layered invariant任一不守恆。

Exit 0只能在所有明確條件成立時出現；「看不到錯誤」不是PASS。

## 8. Frozen provenance 將意圖、觀察、執行與成功分層

Trust chain：

```text
reviewed source manifest v3
  → strict validator重算artifact/producer/sidecar identity
  → atomic frozen input lock
  → source bytes bundle + environment lock
  → launch receipt（互綁所有digest）
  → workers只讀frozen lock並執行bundle bytes
  → verifier獨立重算outputs/post-input identity
  → SUCCEEDED state
  → _SUCCESS最後原子建立
```

### 8.1 `storage_identity_v1` 不是 full BAM hash

大型BAM若未做full SHA，必明示lower assurance：realpath、device/inode、size、mtime/ctime、first/middle/last fixed chunks、header/ref dictionary、BAI full SHA與quickcheck。它無法排除刻意維持metadata與sampled chunks的中段mutation，因此不可寫成full content hash。

### 8.2 Lifecycle component 已驗證的機制

`InterSubMod/scripts/layered_v3_lifecycle.py` 的13/13 tests涵蓋parent flock、preflight no-publish、source byte tamper、unknown `SM_*`、heartbeat stale/PID reuse、child fail-fast、真實SIGTERM與 `_SUCCESS` last。

Strict manifest/lock contract另有17/17 tests，涵蓋zero/duplicate/non-comprehensive、truth/disableFilter option、effective-option origin、BAM/storage mutation、subject swap、global duplicate override、incomplete native validation、measured-CN semantics freeze與錯誤tree backbone。Post-run producer receipt normalizer另有8/8 tests。Active producer真實readback於未完成時exit 3 `E_INDEX_MISMATCH`，只寫 `valid_lock_written:false` failure report，未建立manifest。

Receipt normalizer後續補入真實 symlink regression，現為8/8 tests：producer inventory依GNU `stat -c`比logical-path `lstat`，同時由`storage_identity_v1`另鎖target realpath/metadata/chunks與BAI；避免HCC1395_DORADO的80/79-byte symlink被拿去和91.36/250.15 GB target size誤比。

### 8.3 Verifier component 已驗證的機制

`InterSubMod/scripts/verify_layered_v3.py` 的16/16 tests涵蓋empty/duplicate/extra、source/env/input/output mutation、five-part與read-unsupported守恆、layered/region/site-ledger recomputation、tree/input role防交換、false payload-v2 claim與premature `_SUCCESS`。`run_layered_v3.py` 的12/12 tests另涵蓋跨run-ID global flock、production resource-policy freeze、LongPhase衝突、NFS mountstats、baseline期間新process TOCTOU與pre-bundle failure。Cross-component adapter已整合；real 7-dataset artifacts仍未具備，所以component-ready不等於launch-ready。

### 8.4 實際程式碼 inventory 與 authority

| 程式 / schema | 角色 | Authority / status |
|---|---|---|
| `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_longphase_production_sidecars.sh` | 7-dataset LongPhase producer/FIFO orchestration | active run-start hash authority；執行中不可修改 |
| `.../capture_longphase_tagged_bam_sidecar.py` | BAM stream → 9-column HP/PS sidecar | active producer source；coordinate_join_v1 |
| `.../validate_streamed_longphase_sidecar.py` | native capture/log/VCF validation | 原始 evidence，不可被post-run normalizer覆寫 |
| `.../prepare_clean_layered_manifest_v3.py` | completed producer → reviewed source manifest v3 | 16-case contract suite的一部分 |
| `.../validate_layered_v3_inputs.py` | strict artifact/producer/method validation → frozen lock | fail-closed authority；valid/failure outputs分離 |
| `.../schemas/layered_input_manifest_v3.schema.json` | source intent schema | exact 7/6、roles、producer/read-tag/CN contract |
| `.../schemas/layered_input_lock_v1.schema.json` | observed frozen input schema | worker/verifier唯一input authority |
| `.../schemas/longphase_production_capture_receipt_v1.schema.json` | post-run producer evidence schema | schema frozen；normalizer 8/8 PASS |
| `InterSubMod/scripts/layered_v3_lifecycle.py` | flock/staging/source/env/state/heartbeat/trap/_SUCCESS | 13/13 component PASS |
| `InterSubMod/scripts/run_layered_v3.py` | lifecycle + validator + scientific worker + verifier adapter | 12/12 integrated component PASS；real launch gated |
| `InterSubMod/scripts/verify_layered_v3.py` | post-input/output independent verifier | 16/16 integrated component PASS |
| `.../sm_linkage_genomewide.py` | pileup、exact sidecar join、HP/PS diagnostics | scientific consumer core |
| `.../sm_multilocus_combinations.py` | multi-locus R/A/X/family groups與five-part outputs | frozen bundle worker |
| `.../layered_tree_reconstruction.py` | V1–V7 candidate tree/solver | frozen bundle worker |
| `.../build_region_view.py` | region-level determinacy/funnel | frozen bundle worker |
| `.../build_ssnv_site_ledger.py` | raw/input/recal/tree four-role lossless ledger | frozen bundle worker |

表中 `...` 皆展開自 `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/`；frozen source bundle必保存實際bytes與full path mapping，不能只保存這張人讀表。

## 9. Layered full run 的設定不可依隱含預設漂移

| Parameter | Locked value | 方法意義 |
|---|---:|---|
| genome | chr1–22 | comprehensive autosomal scope |
| datasets / biological | 7 / 6 | exact cardinality；non-empty/unique |
| parallel samples | ≤2 | 避免NFS疊加 |
| chromosome parts | 5 | 22 autosomes固定partition，無漏無重 |
| `MINREAD` | 3 | multi-locus基本read支持門檻 |
| `MAX_SNV` | 8 | 每group分析上限 |
| `TIER_R` | 50,000 bp | 相鄰sSNV group距離 |
| `MAPQ_MIN` | 20 | consumer alignment gate，與sidecar綁定 |
| `BASEQ_MIN` | 0 | pileup base gate；明示而非library default |
| `VERIFY_EVERY` | 1 | 每個eligible unit跑V1–V7 |
| `ANALYSIS_TREE_CAP` | 0 | 不截斷scientific candidate set |
| `DISPLAY_TREE_CAP` | 32 | 只限制展示，不得反饋分析 |

固定5-part partition：

1. chr1,6,11,16,21
2. chr2,7,12,17,22
3. chr3,8,13,18
4. chr4,9,14,19
5. chr5,10,15,20

## 10. Launch gate 同時看方法、資料、資源與排他性

Conditional GO 必須同時成立：

1. Producer `verification_summary.json` exact 7/7 pass；所有truth flags absent。
2. 每sample sidecar unique、hash、BAM/producer input binding與recalibrated VCF flow全PASS。
3. v3 manifest/lock/source/env/receipt cross-component tests PASS。
4. COLO829 bounded join方法與synthetic edge cases PASS。
5. 無另一LongPhase/layered同類process；parent flock取得。
6. launch前300-second baseline：iowait <20%、big8 NFS READ <80 decimal MB/s、big7 free ≥500 GiB、RAM available ≥128 GiB、load/logical CPU ≤1.25。
7. `SM_PARALLEL_SAMPLES≤2`，且production與layered不得重疊。

任一不成立即NO-LAUNCH；不可用「快完成」或人工目測override。

### 10.1 Reviewed handoff supervisor：現在只等待，不代表full run已啟動

為避免14–24小時producer期間由人工作業漏掉gate，新增one-shot supervisor：`InterSubMod/scripts/continue_layered_v3_after_producer.py`。其28/28 synthetic tests涵蓋exact 30秒polling、status append/order/regression、雙supervisor flock、source/tool/PATH drift、partial closeout/receipts、no-overwrite與foreground runner failure。

Reviewed plan為`InterSubMod/research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_reviewed_v1.json`，byte SHA-256=`7e0ab871ee2fa15e772d17da75ebe0e836cc7dcc6c2726e7b136a3b305e3da6a`，permissions=`0444`；5 inputs與25 executable/source/schema/tool pins均由主代理重算PASS。PID `929656`於`2026-07-11 06:16:43 +0800`以JSON authorization、CLI execute flag與out-of-band SHA三重gate啟動並持有唯一handoff flock。

目前handoff workspace與`20260711_layered_reconstruction_v3_full_no_truth_v1` run root皆不存在，故狀態是**active wait / NO-LAUNCH**。只有7/7+ALL PASS、producer/LongPhase/capture quiescent、closeout與7 receipts逐一readback後，才會建立v3 manifest並以前景方式呼叫frozen runner；第二個supervisor會以`E_HANDOFF_LOCKED`拒絕。限制：等待期stdout/stderr連到parent Codex pipe，尚無persistent monitor log；正式stage開始後才在handoff workspace保存command/log/PID/exit receipts。

## 11. Step → Verify 執行鏈

1. 等既有production producer完成  
   → 驗證：7/7 sample各有PASS status、sidecar/VCF/index/validation/hash，無START-only。
2. Freeze producer provenance  
   → 驗證：run-start code 4/4 readback、LongPhase binary SHA、command argv與完整inputs/outputs綁定。
3. 建source manifest v3與strict frozen lock  
   → 驗證：exact 7/6、roles正確、sidecar/BAM/producer subjects重算相等；failure不產valid lock。
4. 建source bundle/environment lock與launch receipt  
   → 驗證：一byte mutation、unknown `SM_*`、missing import均non-zero；正式run root尚不存在。
5. 通過resource/exclusivity gate後publish run root  
   → 驗證：parent flock、state `READY→RUNNING`、heartbeat與child registry可讀。
6. 跑7 datasets × 5 parts  
   → 驗證：35/35 split exact=exposure、missing/conflict=0；任何child fail取消並wait siblings。
7. 跑layered/region/site ledger與v3 verifier  
   → 驗證：五part/funnel/V1–V7/site rows/hash/post-input identity全部重算PASS。
8. Finalize  
   → 驗證：state先SUCCEEDED，`_SUCCESS`最後建立；evidence ledger與本報告/HTML同步。

## 12. 失敗、重跑與復原規則

- 不刪除partial output；保留staging/FAILED作forensic evidence。
- 沒有 `_SUCCESS` 就不是成功，禁止resume-by-file-existence。
- 同run ID已存在即fail；retry使用新 `_retryNN` run ID。
- Producer或worker沒有可靠checkpoint時，不把partial BAM/TSV當可續跑點。
- 若sidecar equivalence或global uniqueness失敗，停止使用sidecar；回退需要≥2.3 TB正式容量的physical tagged-BAM方案，或另設計可審核分段策略。
- 若production_default_filter不是研究問題所需，建立新role/manifest重跑，不覆寫現有artifact語意。

## 13. 目前已完成與仍缺的證據

| Evidence item | Status | Claim tier |
|---|---|---|
| KB/source control-flow audit | PASS | ⭐⭐⭐⭐⭐ L1 code fact |
| 7 input quickcheck/index/dictionary | PASS | ⭐⭐⭐⭐ L2 observed inputs |
| production command truth flags | active two commands absent | ⭐⭐⭐⭐ L2 process observation |
| coordinate_join synthetic | 3/3 PASS | ⭐⭐⭐⭐⭐ L1 mechanism fixture |
| COLO829 bounded join | 35/35 PASS | ⭐⭐⭐⭐ L2 bounded real data |
| v2 external fail-closed runner | 4/4 PASS | ⭐⭐⭐⭐⭐ L1 synthetic mechanism |
| v3 strict input contract | 17/17 PASS；normalizer 8/8；active incomplete root exit 3、無valid lock | ⭐⭐⭐⭐⭐ L1 mechanism + L2 state observation |
| v3 lifecycle | 13/13 PASS | ⭐⭐⭐⭐⭐ L1 synthetic mechanism |
| v3 runner / verifier | runner 12/12；verifier 16/16；cross-component integrated | ⭐⭐⭐⭐⭐ L1 synthetic mechanism |
| handoff supervisor | 28/28 PASS；reviewed plan 30/30 pins；PID 929656 active wait | ⭐⭐⭐⭐⭐ L1 mechanism + L2 process observation |
| production outputs | 0/7 complete at initial cutoff | ⚠ not yet evidence |
| new layered full run | not launched | ⚠ not yet evidence |

## 14. Evidence index

- `InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md` — scope與remediation gates。
- `InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md` — PROBE verdict與falsifier。
- `InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md` — live decisions/deviations。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/04_clean_tagging_launch_audit.md` — producer/source/容量。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/05_fail_closed_wiring_audit.md` — strict schema與adversarial matrix。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/06_launch_resource_redteam.md` — host/resource/exclusivity。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/07_layered_v3_lifecycle_contract.md` — lifecycle component。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/08_layered_v3_runner_integration.md` — frozen runner、resource/NFS/process gate與cross-component整合。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/09_layered_v3_handoff_supervisor.md` — reviewed plan、30 pins、wait/handoff lifecycle與28-case tests。
- `InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/equivalence_probe.json` — bounded real receipt。
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1/` — active producer output root。

## Scope Limitation

- 已驗證：方法 control flow、角色、capacity、bounded join、synthetic fail-closed/lifecycle/verifier components。
- 未驗證：active producer 7/7 final artifacts、production global sidecar-to-BAM binding、300秒實機resource baseline、35/35 real splits與final outputs。
- 推論限制：不得把本版寫成production validation complete；不得把coordinate_join_v1稱為payload-v2；不得把historical truth-tagged BAM當clean authority。
- 補強動作：依§11順序完成，不因時間或檔案存在跳 gate。

## Provenance

- Cycle：`20260710-2244-layered-reconstruction-v2`
- Commit：`4fb9e742482b63a660de19a1f1bd07d49d713111`
- Branch：`research/subclonal-reconstruction-202606`（11 worktrees並行；current只相對本worktree）
- Generated：2026-07-11 05:04 +08:00
- Narrative pre-registration：`InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md` 2026-07-11 04:57條目。
