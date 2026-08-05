<!--
建立時間：2026-07-11 04:23 Asia/Taipei
目標：建立最終 HTML 流程說明書的 claim-to-evidence map、流程節點、判斷/設定 registry 與發布 gate
處理範圍：InterSubMod/research/20260710_layered_reconstruction_v2 的 pre-decision、implementation notes、audit_notes；指定 production/sidecar/layered scripts；active production run root 的 manifest/params/code hash/status
關聯檔案：InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md；InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md
限制：本文件為唯讀 evidence mapping；未啟動計算、未修改 live runner/consumer；active run 與 working-tree scripts 均為 moving state
任務類型：B — Comprehensive validation；本 evidence map 不代表 full validation 已完成
-->

# 最終 HTML 流程說明書 Evidence Map

> **用 SCQA + Claim–Evidence–Status：production sidecar 已在執行，但最終 HTML 目前只能標「PROBE / full launch NO-GO」；主要阻擋是 7/7 尚未完成、sidecar 與 clean tagged BAM 尚未完成直接等價驗證、fail-closed/frozen provenance 尚未閉合，以及「filtering-enabled sidecar」和「caller-PASS-preserving tag-only」契約衝突。**（影響：高；信心：高）

## 0. 讀者先看：狀態與證據詞彙

### 0.1 觀察截點

- Evidence cutoff：`2026-07-11T04:21:21+08:00`；audit directory 再確認於 `04:22:27+08:00`。
- Active production root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1/`。
- 截點狀態：`HCC1395`、`HCC1395_DORADO` 各只有 `production_tagging START`；**0/7 可由 root-level status 證明 PASS**。
- `audit_notes/05_fail_closed_wiring_audit.md` 在本截點尚不存在；本圖只引用實際存在的 `01`、`02`、`03`、`04`、`06`。
- Working-tree scripts 可能持續被其他 agent 修改；本文件列出的 SHA-256 是取證快照，不自動等於 active run 啟動後所有 sample 實際執行的 source identity。

### 0.2 Status vocabulary

| Status | 在本報告中的精確意思 |
|---|---|
| `CONFIRMED` | 可直接由指定 code/artifact 在截點讀得；不代表 biological truth。 |
| `IN_PROGRESS` | 已有 active process/root 或 START evidence，但完成條件尚未成立。 |
| `PROBE` | 只核准 bounded equivalence/adversarial validation；不可擴成 7-dataset claim。 |
| `PENDING` | 有設計或程式入口，但必要 artifact 尚未存在或未驗收。 |
| `NO-GO` | 在列出的 falsifier/gate 通過前，禁止重複 launch、full launch 或 promotion。 |
| `SUPERSEDED` | 舊 artifact/claim 可作歷史證據，不可作 current canonical evidence。 |

### 0.3 Confidence vocabulary

| Confidence | 使用條件 |
|---|---|
| 高 | 直接 code/artifact、機器輸出或兩份獨立稽核一致。 |
| 中 | 直接 code 加工程推論，或 point-in-time resource observation。 |
| 低 | 設計假設、未執行的 acceptance、缺少直接對照 artifact。 |

## 1. 最終 HTML 的 canonical 流程圖節點

```text
[N0 Scope/contract]
  7 datasets / 6 biological samples
  LongPhase: full input traversal；Layered: chr1-22
        |
        v
[N1 Historical-confound audit] ----------------------------- CONFIRMED
  6/7 historical tagged BAM truth-BED-conditioned
  caller PASS != LongPhase recalibrated PASS != truth
        |
        v
[N2 Production input manifest] ------------------------------ CONFIRMED snapshot
  raw ClairS + caller PASS + germline phased VCF
  original tumor/normal BAM + reference
        |
        v
[N3 Immutable production root/preflight] -------------------- PARTIAL
  copies manifest + hashes 3 scripts + params/status
  but workers still read mutable source manifest
        |
        v
[N4 LongPhase-S somatic_haplotag] --------------------------- IN_PROGRESS
  no --truth-vcf / --truth-bed
  -q20 --tagSupplementary -t12
  FILTERING ENABLED + --output-somatic-vcf
        |
        +--------------------------> recalibrated all/pass VCF
        |
        v
[N5 FIFO BAM stream -> alignment HP/PS sidecar] ------------- IN_PROGRESS
  BAM payload not persisted；no BAI deliverable
        |
        v
[N6 Per-sample streamed validator] -------------------------- PENDING 7/7
  truth log / counts / HP vocabulary / exact conflicts / VCF keys
        |
        v
[N7 bgzip+tabix + output hashes] ---------------------------- PENDING 7/7
        |
        v
[N8 production verification_summary 7/7] ------------------- PENDING / HARD GATE
        |
        v
[N9 schema-2.1 clean layered manifest] ---------------------- PENDING
  original BAM + exact HP/PS sidecar + raw/caller/recal VCF roles
        |
        +--> [N10 fail-closed input preflight] -------------- CODE PARTIAL
        |
        +--> [N11 direct-BAM vs sidecar bounded equivalence] - PROBE
        |          exact observations/vectors/families/digest
        |
        v
[N12 7-dataset layered full run] ---------------------------- NO-GO NOW
  2 samples parallel × 5 chr shards；chr1-22
        |
        v
[N13 MLHP -> solver -> region -> site ledger -> verifier] --- PENDING NEW RUN
        |
        v
[N14 canonical HTML/report builder] ------------------------- PENDING / NO canonical builder
```

### 1.1 節點必附的 evidence anchors

| Node | Primary evidence | Secondary audit | Current status |
|---|---|---|---|
| N0 | active `input_manifest.json`; active `params.json` | `03_verification_risk_audit.md` §1.1、§8 | CONFIRMED scope；「genome-wide」用詞仍需拆層 |
| N1 | `01_upstream_core_audit.md` §4–§5 | `04_clean_tagging_launch_audit.md` §1–§2 | CONFIRMED historical confound |
| N2 | active `input_manifest.json` | `04_clean_tagging_launch_audit.md` §3 | CONFIRMED 7/6 input snapshot |
| N3 | `run_longphase_production_sidecars.sh`; active `code.sha256` | `06_launch_resource_redteam.md` §4 | PARTIAL frozen provenance |
| N4 | `run_longphase_production_sidecars.sh` command | `04_clean_tagging_launch_audit.md` §2.3、§11 | IN_PROGRESS; role mismatch unresolved |
| N5 | `capture_longphase_tagged_bam_sidecar.py` | `implementation-notes.md` 04:08、`06` §5.1 | IN_PROGRESS; not persistent BAM |
| N6 | `validate_streamed_longphase_sidecar.py` | `06` §7 acceptance | Validator exists; 7/7 results absent |
| N7–N8 | producer runner tail + active `run_status.tsv` | `06` §7 | PENDING |
| N9 | `prepare_clean_layered_manifest.py` | `implementation-notes.md` 04:02 | Code exists; output cannot exist before 7/7 |
| N10 | `validate_layered_v2_inputs.py`; layered runner | `02_layered_chain_audit.md` P0-04/P0-05 | PARTIAL; not fully fail-closed |
| N11 | `pre-decision-audit.md` §8.4 | `implementation-notes.md` 04:08；`06` §7 | PROBE; full equivalence pending |
| N12–N13 | `run_layered_7samples_newbb.sh` | `02` §5 S1–S4；`03` §9 | NO-GO until all gates |
| N14 | `02` §5 S7 | `03` §3、§9 | No current canonical end-to-end HTML builder |

## 2. Claim → evidence → confidence → status

### 2.1 Scope、input role 與 production semantics

| ID | Claim 可寫法 | Evidence path/section | Confidence | Status / caveat |
|---|---|---|---|---|
| C01 | 本矩陣有 7 datasets、6 biological samples；HCC1395/DORADO 是同一 biological ID。 | active `input_manifest.json`; `03_verification_risk_audit.md` V13 | 高 | CONFIRMED；禁止寫 7 independent samples。 |
| C02 | Layered primary reconstruction scope 是 chr1–22；不是含 X/Y 的全基因組。 | `run_layered_7samples_newbb.sh` `SPLITS`; `03` §1.1 | 高 | CONFIRMED。LongPhase input traversal 與 layered scope需分開寫。 |
| C03 | 6/7 historical tagged BAM 曾受 truth BED 改變 tagging universe；COLO829 無 truth BED。 | `01_upstream_core_audit.md` §5.2 | 高 | CONFIRMED historical blocker；不能沿用舊 BAM 當 clean evidence。 |
| C04 | ClairS raw、ClairS caller PASS、LongPhase recalibrated all/pass、truth set 是不同 role。 | `implementation-notes.md` 04:12；`01` §4；active manifest contract | 高 | CONFIRMED；禁止用 `somatic_pass` 一詞跨 role。 |
| C05 | Primary operational backbone 是 paired ClairS FILTER=PASS biallelic SNV，不是 biological truth。 | active manifest; `01` §10；`03` §8 | 高 | CONFIRMED wording gate。 |
| C06 | Active production command不帶 `--truth-vcf/--truth-bed`，使用 `-q 20 --tagSupplementary -t 12`。 | producer runner command; active `params.json`; active status detail | 高 | CONFIRMED command design / START；每 sample完成後仍須 log readback。 |
| C07 | Active production command仍啟用 LongPhase variant filtering，並帶 `--output-somatic-vcf`。 | producer runner；`04_clean_tagging_launch_audit.md` §2.3、§11 | 高 | CONFIRMED；這不是 audit §7 定義的 caller-PASS-preserving tag-only。 |
| C08 | `--disableFilter` 才能讓所有已供應 caller-PASS variants保持可參與 tagging；是否採用是高影響方法決策。 | `04_clean_tagging_launch_audit.md` §0.1、§2.3、§7.1 | 高 | **UNRESOLVED CONTRACT**；最終 HTML 不得模糊成「all caller PASS tags」。 |
| C09 | Active output 是 original BAM＋HP/PS sidecar＋recalibrated VCF，不是 persistent tagged BAM/BAI。 | producer/capture scripts；active manifest `output_mode`; `04` §11 | 高 | CONFIRMED deviation；使用者原始「tagged BAM」deliverable尚未達成。 |
| C10 | 7 個 persistent tagged BAM 的容量需求超出當時 big7空間；sidecar是容量驅動的替代 probe。 | `pre-decision-audit.md` §8；`04` §10；`06` §5 | 高 | CONFIRMED at audited point；resource數字需標 observation time。 |

### 2.2 Producer、capture 與 streamed validator

| ID | Claim 可寫法 | Evidence path/section | Confidence | Status / caveat |
|---|---|---|---|---|
| P01 | Production manifest由 layered manifest與 historical tagged BAM `@PG`/run context回推出原 BAM、normal、reference、germline與caller PASS。 | `prepare_longphase_production_manifest.py` | 高 | CONFIRMED code behavior；是派生 provenance，不是外部 registry oracle。 |
| P02 | Manifest固定7/6、record counts、tag vocabulary、MAPQ20、supplementary on、無truth flags。 | active `input_manifest.json` | 高 | CONFIRMED snapshot；未固定 LongPhase binary hash/version。 |
| P03 | Producer拒絕既有 run root、建立 samples、複製 manifest、寫 params/status/code hashes。 | `run_longphase_production_sidecars.sh` top/prelude | 高 | PARTIAL immutability；見 P04/F02。 |
| P04 | Producer worker的 `value()`與 sample list仍讀原始 `$MANIFEST`，不是 frozen run-root copy。 | producer runner `value()`/`mapfile` | 高 | FAIL-CLOSED GAP；run中 manifest漂移可跨 sample生效。 |
| P05 | FIFO只作 streamed BAM transport；capture略過 unmapped，將每個 mapped alignment寫9欄。 | capture script | 高 | CONFIRMED code；完成結果待7/7。 |
| P06 | Sidecar identity欄是 `chrom,start,end,QNAME,FLAG,MAPQ,CIGAR_B2`，另存HP/PS；CIGAR digest為BLAKE2b 8 bytes。 | capture script | 高 | CONFIRMED schema；consumer exact key是否含相同所有欄需另外明列。 |
| P07 | Capture要求 coordinate order、不接受未知HP，且 mapped row數需等於 total-unmapped。 | capture script summary/pass | 高 | CONFIRMED check definition；未證 full content與原BAM相同版本。 |
| P08 | Stream validator檢 truth log空、benchmark removal absent、parser/input VCF count、alignment/tag/row counts、malformed/unknown/conflict、recalibrated VCF key conservation。 | `validate_streamed_longphase_sidecar.py` `checks` | 高 | CONFIRMED implementation；不是 BAM↔sidecar直接等價 oracle。 |
| P09 | Recalibrated VCF「preserves all input keys」只驗 key multiset，不代表 FILTER/classification 保持不變。 | validator `vcf_keys`; implementation notes 04:12 | 高 | CONFIRMED distinction。 |
| P10 | Validator沒有 sidecar SHA、自身 source BAM identity/header/@PG binding、producer binary hash或manifest contract hash。 | validator全文；active `code.sha256` | 高 | OPEN P0 provenance gap。 |
| P11 | Per-sample validation後才bgzip/tabix、建立all/pass recalibrated VCF與output hashes。 | producer runner post-validation block | 高 | CONFIRMED order；但沒有 atomic `ACCEPTED` publish。 |
| P12 | Producer沒有 trap/heartbeat/sample lock/atomic status；partial root可能只停在START。 | producer runner全文；`01` U08；`06` §7 | 高 | OPEN durability gap。 |
| P13 | Root aggregate只有7份 `sidecar_validation.json`且全pass才產7/7 PASS gate。 | producer runner tail | 高 | CONFIRMED designed gate；active root尚未達成。 |
| P14 | 截點 active root只有2 START、0 aggregate verify PASS。 | active `run_status.tsv` | 高 | IN_PROGRESS；任何「7/7 clean tags完成」claim為假。 |

### 2.3 Clean handoff、layered preflight 與 launch

| ID | Claim 可寫法 | Evidence path/section | Confidence | Status / caveat |
|---|---|---|---|---|
| H01 | Clean handoff manifest將 original tumor BAM與sidecar分角色，並保留historical truth-restricted BAM欄。 | `prepare_clean_layered_manifest.py` | 高 | CONFIRMED code design。 |
| H02 | Handoff只有所有 required files存在且每sample validation.pass=true才可寫schema2.1。 | prepare-clean script | 高 | PARTIAL gate；未驗validation檔與sidecar/hash/BAM身份關聯。 |
| H03 | Handoff沒有先驗 dataset sample-set完全一致、duplicate sample、production root ID/hash與aggregate 7/7 summary。 | prepare-clean script全文 | 高 | OPEN fail-closed gap。 |
| H04 | Layered runner接受schema2.0/2.1，但實際 `run_sample` 無條件要求 sidecar/recal/validation欄。 | layered runner preflight/run_sample | 高 | CONTRACT AMBIGUITY；schema2.0 acceptance與required fields不一致。 |
| H05 | Layered runner先建立正式run root、copy manifest、寫status/code hash，之後才執行Python input validator。 | layered runner順序 | 高 | FAIL-CLOSED GAP；invalid input仍留下formal-looking root。 |
| H06 | Layered runner的worker與final verifier仍讀原始 `$MANIFEST`；不是 `$RUN_ROOT/input_manifest.json`。 | `manifest_value`; final verifier invocation | 高 | OPEN P0；違反frozen manifest consumption。 |
| H07 | Input validator檢BAM header/reference/BAI、PASS-only ClairS VCF、CN labels、sidecar index fetch chr22、validation pass/truth flags。 | `validate_layered_v2_inputs.py` | 高 | CONFIRMED checks。 |
| H08 | Input validator不hash sidecar本體、不核 producer root/code/input hash、不驗 schema2.1 必有 sidecar、不驗exact identity contract。 | input validator全文 | 高 | OPEN P0。 |
| H09 | MLHP每sample分5個互斥chr shard，預設最多2 samples，參數MINREAD3/MAX8/gap50kb/MAPQ20/BQ0。 | layered runner env/splits/params | 高 | CONFIRMED design；full sidecar run未測。 |
| H10 | 每part gate要求scope守恆、sidecar missing/conflicts=0、exact matches=alignment exposures。 | layered runner post-MLHP jq | 高 | GOOD fail gate；source只用`!= BAM_HP_PS`負面判斷，仍應改成exact mode/path/digest。 |
| H11 | Site ledger覆蓋raw ClairS records，primary solver仍只用caller PASS sSNV。 | layered runner `build_ssnv_site_ledger.py`; implementation notes 04:12 | 高 | CONFIRMED intended separation；新full artifact pending。 |
| H12 | Layered `code.sha256`沒有hash runner自身、manifest-preparer、source bundle/patch或environment lock。 | layered runner code-hash block | 高 | OPEN frozen provenance gap。 |
| H13 | Layered run沒有trap/heartbeat/atomic stage output；status多process append也無lock。 | layered runner全文；`02` P1-15 | 高 | OPEN durability gap。 |
| H14 | 現在不得啟動sidecar-aware 7-dataset full run。 | `pre-decision-audit.md` §8；`06` §2/§7 | 高 | HARD NO-GO until N8/N10/N11 gates pass。 |

### 2.4 Core method、verification 與 claim ceiling

| ID | Claim 可寫法 | Evidence path/section | Confidence | Status / caveat |
|---|---|---|---|---|
| M01 | Region grouping是相鄰sSNV gap≤50kb的connected component；total span可>50kb。 | `02_layered_chain_audit.md` S1/decision D02；implementation lore | 高 | CONFIRMED；禁止寫region≤50kb。 |
| M02 | group<2為singleton；>8取最小span連續8點，tie取最左。 | `02` S1/D03–D04 | 高 | CONFIRMED；被cap點要進funnel。 |
| M03 | pileup預設MAPQ≥20、BASEQ≥0；R/A/X由exact REF/ALT/其他或未覆蓋決定。 | `02` S1/D05–D08 | 高 | CONFIRMED；BQ0 sensitivity未完成。 |
| M04 | HP `1*→family1`,`2*→family2`,`3→family3`,其他含HP4→none。 | `02` S1/D09 | 高 | CONFIRMED current behavior；HP4 role是P1未決。 |
| M05 | full/partial genotype需count≥3；pooled retention可能通過但各family均不足3。 | `02` S1/P0-02 | 高 | CONFIRMED counterexample；retained→unit守恆仍有P0。 |
| M06 | Solver以最少hidden nodes找root-connected unit-flip tree；hidden≤4、每level budget150k。 | `02` S2/D14–D15 | 高 | CONFIRMED model；capped不是exact result。 |
| M07 | Analysis candidate cap=0（全候選），display cap=32；display不得反過來改分析。 | layered runner params；implementation notes | 高 | CONFIRMED intended contract；HTML舊builder違反。 |
| M08 | Primary lineage只含mutation-bearing HP1/HP2；root-only、HP3、none不進primary denominator。 | implementation notes；`02` S2 | 高 | ACCEPTED method contract；新full artifact pending。 |
| M09 | V1–V7是implementation invariant suite；V5驗count一致，不是獨立tree-set completeness oracle。 | `03_verification_risk_audit.md` §4/R02 | 高 | CONFIRMED limitation；禁止稱formal independent oracle。 |
| M10 | Missing CN應為unavailable；neutral只在declared segmentation contract下成立。 | implementation notes；`02` S3 | 高 | ACCEPTED；SAVANA uncovered semantics仍待證。 |
| M11 | recurrence、capped、ambiguous與CN verdict要分軸；單一region label可能隱藏capped。 | `02` P0-03 | 高 | OPEN P0 display/schema gap。 |
| M12 | read-AF只能稱ordering heuristic；不是CCF/posterior truth。 | implementation notes；`02` S5；`03` R07 | 高 | CONFIRMED claim ceiling。 |
| M13 | methyl L3目前只能not-evaluated/bounded auxiliary，禁止rank/confirm。 | implementation notes；`02` S6 | 高 | CONFIRMED canonical restriction。 |
| M14 | 現有HTML/workstation不相容current v2且含hardcoded/舊schema claim。 | `02` S7/P0-06/P0-07；`03` §3 | 高 | SUPERSEDED / NO-GO as validation evidence。 |
| M15 | 沒有single-cell/multi-region orthogonal truth，最多稱model-based regional mutation-state inference。 | implementation notes；`03` §9.3 | 高 | Biological validation未完成。 |

### 2.5 Frozen provenance、可重現性與執行狀態

| ID | Claim 可寫法 | Evidence path/section | Confidence | Status / caveat |
|---|---|---|---|---|
| F01 | Active production root保存input manifest、params、status與producer/capture/validator hashes。 | active root四檔 | 高 | CONFIRMED partial provenance。 |
| F02 | Production `code.sha256`不包含LongPhase binary/source commit、Python/pysam/bgzip/tabix/bcftools版本、host/env、handoff preparer。 | active `code.sha256`; producer runner | 高 | NOT FROZEN ENOUGH。 |
| F03 | Active run hash鎖定producer script SHA=`06d6…`, capture=`c28f…`, validator=`42c4…`。 | active `code.sha256` | 高 | CONFIRMED launch hash；須在每sample前/後readback才防live source drift。 |
| F04 | 本取證時 layered runner SHA=`d8a803…`、input validator=`6e7ee4…`、prepare-clean=`a8ce0f…`。 | 2026-07-11 04:21前本輪`sha256sum` | 高 | Snapshot only；未綁定任何新layered run。 |
| F05 | Hash可偵測內容但不能重建uncommitted source；active run root沒有source/patch bundle。 | `03` §1；active root listing | 高 | OPEN reproducibility gap。 |
| F06 | Current source manifest copy存在，但producer/layered仍可在後續sample讀到mutable original。 | producer/layered runners | 高 | OPEN execution-identity gap。 |
| F07 | 大BAM未做full-content hash；多處只靠path/stat/header/BAI。 | `04` §3.2/§12；input validator | 高 | Declared limitation；不能寫full input content frozen。 |
| F08 | 目前沒有新的layered run ID、7/7 completion、verification/output hashes可引用。 | active root status；`06` verdict | 高 | PENDING；最終HTML數值區必保持空/待完成。 |

## 3. 最終 HTML 必須逐項解釋的最基本判斷與設定

| Topic | 必寫精確內容 | 不可只寫 |
|---|---|---|
| Production role | filtering-enabled/no-truth/recalibrated sidecar，或disableFilter/tag-only persistent BAM；二選一並寫理由 | clean tagged data |
| Truth | flags absent、benchmark removal absent、truth只可另作evaluation | truth-free所以是真值 |
| Dataset n | 7 datasets / 6 biological samples / HCC pipeline replicate | 7 samples |
| Scope | LongPhase traversal scope與layered chr1–22分開 | whole genome |
| Caller universe | paired ClairS PASS operational backbone | true somatic variants |
| Sidecar identity | exact key各欄、coordinate convention、CIGAR digest、duplicate/conflict policy | qname join |
| Alignment policy | unmapped/secondary/supplementary/duplicate/QC-fail/orphan的納入規則 | all reads |
| LongPhase params | `-q20`,`-p`是否明示、supplementary、filtering、tumor DNA fraction、threads | defaults |
| VCF roles | raw、caller PASS、recal all、recal PASS、truth set | somatic VCF |
| Grouping | adjacent gap 50kb、singleton、densest8/tie-left | 50kb region |
| Read calls | MAPQ20、BQ0、R/A/X、deletion/refskip、support3 | supported reads |
| HP families | 1*/2*/3/4/none/extended tag policy | phased/unphased |
| Solver | minimal hidden、unit-flip、hidden cap/budget、capped fallback | enumerate all trees |
| Verification | V1–V7各自檢查、eligible分母、skipped/n/a、V5 limitation | all pass |
| CN | source、coordinate system、coverage、anchor、missing、LOH0.7 heuristic | CN-adjusted |
| read-AF | formula/threshold/temperature、無purity/CN correction | CCF posterior |
| L3 methyl | bounded/not evaluated/prohibited uses | methyl confirms tree |
| Provenance | run ID、frozen manifest/code/env/input/output hashes、dirty/patch state | reproducible |
| Completion | START、PASS、ACCEPTED、aggregate7/7的不同層級 | finished |

## 4. 已知模糊、矛盾與 NO-GO register

| ID | Ambiguity / contradiction | Why it matters | Required closure | Current gate |
|---|---|---|---|---|
| A01 | 使用者要persistent production tagged BAM；active run只保留sidecar。 | deliverable與可審核性不同。 | 正式接受sidecar替代，或取得≥2.3TB合法空間產BAM/BAI。 | NO-GO宣稱BAM完成 |
| A02 | active sidecar filtering enabled；audit tag-only contract要求`--disableFilter`。 | 哪些caller-PASS variants參與HP tagging可能改變。 | 方法owner明確選定角色；若改契約要新run ID，不能改live run。 | HARD decision gate |
| A03 | active manifest寫「all record keys preserved」容易被誤讀成all variants被tagging採用。 | key conservation≠tagging inclusion。 | HTML拆成VCF-key gate與tagging-universe gate。 | Wording NO-GO |
| A04 | sidecar exact join smoke有正面證據，但direct clean BAM vs sidecar full observation/family digest尚未完成。 | sidecar不是自明等價。 | 預註冊bounded direct-vs-sidecar digest全等。 | PROBE only |
| A05 | production validator不bind sidecar到原BAM/header/@PG/binary。 | 可拿對格式但錯版本artifact通過。 | source identity/hash/header contract + adversarial swap test。 | NO-GO full launch |
| A06 | producer與layered worker讀mutable manifest。 | frozen copy不一定是真正consumed input。 | validation前freeze；之後所有讀取只指run-root copy；hash readback。 | NO-GO provenance |
| A07 | layered run root在validator前建立。 | fail input會留下formal-looking root。 | preflight在root外或staging root；PASS後atomic publish。 | NO-GO fail-closed claim |
| A08 | input validator未hash sidecar本體或驗production aggregate。 | individual pass檔可被替換/混run。 | sidecar/index/validation/output SHA與root aggregate綁定。 | NO-GO |
| A09 | schema2.0/2.1入口與runner實際required fields不一致。 | schema acceptance語意模糊。 | schema2.1 sidecar模式獨立、required欄hard gate。 | NO-GO |
| A10 | 目前只有2 START、0/7 completion。 | 不能進clean handoff/full run。 | production `verification_summary`: 7/7 pass且hash readback。 | HARD NO-GO |
| A11 | 相同production job已在執行且I/O曾高壓。 | duplicate launch可能損壞資源與證據。 | attach-only；無同類process且resource baseline過門檻。 | NO-GO additional launch |
| A12 | 新layered consumer/source仍在moving working tree。 | code hash不能重建且可能跨sample漂移。 | 停止改動、commit或保存source/patch bundle、啟動前hash。 | NO-GO full launch |
| A13 | retained pooled support可無family unit。 | funnel PASS仍可能漏tree/region denominator。 | zero-family branch或retention改法＋fixture/invariant。 | P0 method |
| A14 | recurrence label可隱藏capped；HTML用display prefix做science。 | UI可改變/隱藏結論。 | multi-axis flags；HTML只讀all-candidate core evidence。 | NO-GO canonical HTML |
| A15 | V5 count一致被誤稱tree-set completeness。 | 截斷仍可能PASS。 | 更名＋第二實作canonical set oracle/proof。 | Claim ceiling |
| A16 | threshold sensitivity未完成。 | MAPQ/BQ/support/gap/cap可能驅動結論。 | prereg sensitivity matrix。 | 方法robustness pending |
| A17 | 無orthogonal clone truth。 | internal consistency不能確認true clone/history。 | 外部single-cell/multi-region資料，否則維持claim ceiling。 | Biological NO-GO |

## 5. Promotion gates：HTML 只能依序升級

### Gate P — Production artifacts

1. 7/7 sample有validation PASS、truth flags absent、input/output VCF keys missing/extra=0、unknownHP=0、exact conflicts=0。
2. active `verification_summary.json` 顯示dataset_count=7、n_pass=7、all_pass=true。
3. input/code/output hashes readback全過；source manifest與真正consumed frozen copy相同。
4. 明確決定current filtering-enabled sidecar是否就是production contract；若不是，現run只能標experiment/probe。

### Gate E — Exact equivalence and adversarial fail-closed

1. 同一bounded region direct clean tagged BAM vs original BAM＋sidecar：eligible observations、R/A/X vectors、family counts、digest全等。
2. 篡改sidecar hash、交換sample BAM、刪index、改truth flag、製造duplicate conflicting tag：validator/runner全部nonzero。
3. Invalid input不得建立可誤認正式的run root/ACCEPTED/verification PASS。

### Gate L — Layered full launch

1. schema2.1 7/7 sample-set完全相等，所有sidecar及produceraggregate hashes綁定。
2. COLO829 bounded probe 5 parts皆`exact=exposures, missing=0, conflicts=0`，且source mode/path exact match。
3. no duplicate LongPhase/layered process；resource baseline符合`06` §7。
4. frozen source bundle/commit/patch、manifest、params、environment與runner自身hash全部保存。

### Gate R — Reportable engineering result

1. 新run 7/7 `complete PASS`，每sample stage/status/output hash齊全，final verification summary/TSV/hash存在。
2. non-capped eligible units V1–V7 zero fail/zero skipped；capped只列n/a。
3. funnel、retained→family/region、candidate completeness、denominator、CN missingness有獨立重算或明確限制。
4. raw ClairS site ledger、caller PASS solver universe、recalibrated output三個role不混用。

### Gate H — Canonical HTML

1. Builder只讀immutable run root；任何missing sample/schema/hash即nonzero。
2. 所有數字由run-scoped JSON/TSV注入；禁止hardcode舊HCC counts/percentages。
3. 頁面同時顯示aggregate、canonical sample、extreme outlier、well-explained sample四層。
4. 頁首顯示status ribbon；footer顯示run ID、scope、commit/patch、manifest/code/output hashes、claim ceiling。
5. 舊HTML標`SUPERSEDED / NOT VALIDATION EVIDENCE`，不可與新數字拼接。

## 6. 最終 HTML 建議章節與 evidence source mapping

| HTML section | Must answer | Machine/human evidence |
|---|---|---|
| 1 Executive verdict | 現在可否full launch/發布？ | active status + Gate table |
| 2 Scope/data roles | 7/6、chr scope、VCF/BAM/sidecar roles | input manifest + params |
| 3 Historical root cause | truth BED怎麼改tagging universe？ | `01` §4–§5 + `04` §2.3 |
| 4 Production command | 每個flag/default如何改方法？ | producer script + frozen argv/log（待） |
| 5 Sidecar schema | 每alignment怎麼唯一join？ | capture script + consumer diagnostics |
| 6 Fail-closed validator | 每個錯誤在哪一層停止？ | streamed/input validators + adversarial outputs（待） |
| 7 Frozen provenance | 真正consumed code/input/env是哪版？ | run manifests/hashes/source bundle（部分待） |
| 8 Layered Stage1 | grouping/read/family/support/CN | `02` S1 + new MLHP JSON（待） |
| 9 Solver/verification | objective/cap/V1–V7/limitations | `02` S2 + `03` §4 + new outputs（待） |
| 10 Full-run results | 7/7 runtime/resource/funnel/class | new immutable run（尚無） |
| 11 Robustness | parameter/caller/CN/adversarial sensitivity | sensitivity matrix（尚無） |
| 12 Limitations | internal vs biological validation | `03` §9.3 + implementation notes |
| 13 Reproduction | exact command、paths、hash checks | frozen argv/run bundle（尚未閉合） |

## 7. Evidence inventory與取證快照

### 7.1 Human-readable evidence

- `InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md`：原GO與2026-07-11 PROBE addendum。
- `InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md`：accepted/proposed decisions與claim ceiling。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/01_upstream_core_audit.md`：upstream truth/VCF role/core contract。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/02_layered_chain_audit.md`：完整call chain、decision registry、P0/P1。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/03_verification_risk_audit.md`：claim-to-source、V1–V7、risk/ambiguity/release gates。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/04_clean_tagging_launch_audit.md`：persistent tag-only contract、capacity與active sidecar角色差異。
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/06_launch_resource_redteam.md`：duplicate launch/resource/consumer thresholds。

### 7.2 Code snapshot

| Script | SHA-256 at evidence read | Role |
|---|---|---|
| `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_longphase_production_sidecars.sh` | `06d6d7be754977dfe9e66cc37ceee2a1b690a4275fc2fb4d7f97e67737346888` | producer runner |
| `.../capture_longphase_tagged_bam_sidecar.py` | `c28f9acedfd9af002afdf499456c36f717fb0adc7aeaf61840c2e8facaf486db` | FIFO capture |
| `.../validate_streamed_longphase_sidecar.py` | `42c41f3574234b5f65d518e368a27a1c272f3c7b86b0235f3ea5375bfa8da6b7` | per-sample validator |
| `.../prepare_longphase_production_manifest.py` | `2a284474c0c44968daec2cffe92b4c85ea394034fa47a7ff4b5b4889f87c6f40` | production manifest derivation |
| `.../prepare_clean_layered_manifest.py` | `a8ce0fc4674382f67d1f66d745527c7b8cd537cb23db33f730a984760994a8d2` | schema2.1 handoff |
| `.../run_layered_7samples_newbb.sh` | `d8a8037f71f9185333ab9e4c311fdb5ee532099b616815c4364619478e63ef2a` | layered orchestrator |
| `.../validate_layered_v2_inputs.py` | `6e7ee45b0a5861761f678d67050649c9e4f139dfa955e80725caee34ae000500` | layered input preflight |

> `...` 只在本表表示共同前綴 `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/`；最終 HTML source manifest應展開完整路徑，不使用省略號。

### 7.3 Active run evidence

| Artifact | What it proves | What it does not prove |
|---|---|---|
| `input_manifest.json` | 7/6 input paths/roles/counts與declared contract | 真正每sample consumed identity未漂移 |
| `params.json` | threads=12、parallel=2、truth_flags=false、MAPQ20、supplementary、sidecar mode | binary/env/filter behavior完整設定 |
| `code.sha256` | manifest＋producer/capture/validator launch hashes | LongPhase binary、env、source bundle、future source readback |
| `run_status.tsv` | 兩sample在03:56:31記START | process health、completion、PASS、7/7 aggregate |

## 8. 當前一句發布判定

**可以發布的是「稽核與進行中狀態說明」，不能發布的是「production clean tagged BAM已完成、fail-closed/frozen provenance已閉合、或新的7-dataset full layered validation已完成」。** 下一個合法里程碑是既有production run 7/7完整驗收，再做direct-BAM vs sidecar bounded equivalence與fault-injection；兩者皆通過後，才可能把full launch從`NO-GO`升為conditional `GO`。
