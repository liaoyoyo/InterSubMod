<!--
建立時間: 2026-07-10 Asia/Taipei
目標: 完整稽核 layered reconstruction v2 的上游資料產生、核心執行鏈、CLI/default、I/O、錯誤與 fallback
處理範圍: ClairS paired PASS、LongPhase-S canonical BAM/VCF、7-dataset manifest、multi-locus read census、tree solver、region view、final verifier、InterSubMod C++ 舊/平行鏈、Shell/Python wrappers
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md; InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md
任務類型: B — Comprehensive validation
服務目標: G2 LongPhase-S + ISM 協同契約; G3 read-level epigenetic/HP 語意; G4 多樣本 reproducibility; G5 外部可驗證
快照分支: research/subclonal-reconstruction-202606
快照 HEAD: 4fb9e742482b63a660de19a1f1bd07d49d713111
working tree: dirty; v2 核心 scripts 多為 modified/untracked，本文稽核的是 2026-07-10 當下 workspace，不等於只 checkout HEAD
狀態: upstream_core_audit_complete; end-to-end_v2_not_complete
限制: 未執行全基因組重建、未執行完整 C++ ctest；只做 read-only/light validation 與 /tmp preflight
-->

# Layered reconstruction v2：上游與核心執行鏈完整稽核

用 SCQA + call-chain audit + claim-to-evidence matrix：先回答目前能不能視為「完整、無模糊、可發布流程」，再逐層列出真正的程式入口、判斷、設定、輸入輸出、錯誤與修正 gate。

> **TL;DR — 目前不能把 7-dataset v2 稱為完成的 comprehensive validation。最嚴重問題不是 solver，而是上游契約不一致：6/7 canonical tagged BAM 在 LongPhase-S haplotag 前被 `--truth-bed` 限制，v2 卻與 genome-wide ClairS PASS VCF 配對；此外 `somatic_pass.vcf.gz` 是 LongPhase-S 輸入，不是 `_sc.vcf` 重校準輸出。現有 preflight 與 verifier 都沒有檢查這兩件事。**（影響：高；信心：高）`[O-L1]`

## 0. 稽核判斷、證據語彙與發布 gate

### 0.1 判斷

- **[O-L1] P0 blocker 1 — benchmark BED 進入 haplotagging 路徑**：LongPhase-S v1.0.0 原始碼會在 `_sc.vcf` 寫出之後、tag reads 之前，刪除 truth BED 外的 TUMOR/TRUTH variant state；7 個 manifest BAM 中有 6 個帶 `--truth-bed`，只有 COLO829 例外。
- **[O-L1] P0 blocker 2 — VCF 角色被混稱**：HCC1395 `somatic_pass.vcf.gz` 有 113,997 個 ClairS paired PASS SNV，是 LongPhase-S 的輸入；真正輸出 `HCC1395_tagged_sc.vcf` 是 113,997 records，其中 108,530 PASS、5,467 LowQual。文件多處把兩者寫成同一物。
- **[F-L1] runner preflight 不完整**：正式 runner 只查檔案存在；可驗 ClairS header、PASS-only、VCF index、BAM header/ref dictionary 的 `validate_layered_v2_inputs.py` 沒被 runner 呼叫。
- **[O-L1] 7/7 input files 本身通過手動 preflight**：本稽核手動執行 validator，7/7 pass；這證明目前路徑可讀與 VCF 為 ClairS PASS，不代表正式 runner 已強制此契約，也不檢查 truth-bed-induced BAM scope。
- **[F-L1] v2 C++ 關係必須重寫清楚**：目前 v2 runner 不呼叫 `build/bin/inter_sub_mod`，也不消費 `significance_summary.csv`；其 L3 明寫 `not_evaluated/bounded_auxiliary`。InterSubMod C++ 是 canonical/legacy methylation benchmark 的平行下游，不是本輪 L0-L2 樹重建執行鏈。
- **[O-L1] 小型驗證通過**：Shell `bash -n`、5 個 v2 regression tests、8 個 solver golden cases通過；CTest 只列出 234 個 C++ tests，Python v2 tests 未納入 CTest/CI。
- **[U-L5] full-run gate 未滿足**：稽核時找不到完成的 `verification_summary.json`；`20260710_232501_layered_reconstruction_v2` 只停在兩個 dataset 的 `mlhp START`，沒有活程序、stage FAIL 或 completion 紀錄。

### 0.2 Claim tag

| Tag | 定義 |
|---|---|
| `[F-L1]` | 直接讀版本化／workspace 程式碼或機器產物得到的事實。 |
| `[O-L1]` | 本稽核實際執行命令得到的觀察。 |
| `[I-L2]` | 至少兩個 L1 事實支持的工程／方法推論。 |
| `[U-L5]` | 尚未執行或現有資料不能回答。 |
| L2-L4 | repo 內分析、文件或間接證據；不能冒充直接 source/artifact。 |

### 0.3 發布 gate

在以下四項全部完成前，總報告應標示 **BLOCKED FOR COMPREHENSIVE VALIDATION / PARTIAL ENGINEERING AUDIT**：

1. `[I-L2]` 以**不帶 truth-vcf/truth-bed 的 production LongPhase-S tagging run**重產 7 個 tagged BAM，或以實驗證明現有 truth-bed selection 不改 v2 主結果；建議前者。
2. `[F-L1]` 將 ClairS PASS input 與 LongPhase-S recalibrated output 拆成不同、不可覆寫的 canonical 檔名與 manifest role。
3. `[F-L1]` runner 內建完整 preflight、BAM `@PG` contract、group→unit→region conservation、interruption trap。
4. `[U-L5]` fresh 7-dataset full run 產生 7/7 completion、checksums、verification summary、H8 sensitivity matrix。

## 1. 真正存在的兩條流程，不可混寫

### 1.1 Canonical/legacy upstream + C++ benchmark chain

```text
ClairS paired output.vcf.gz
  └─ bcftools view -f PASS
       └─ somatic_pass.vcf.gz             # LongPhase-S INPUT
            ├─ LongPhase-S somatic_haplotag
            │    ├─ *_tagged.bam          # HP-tagged BAM
            │    ├─ *_tagged_sc.vcf       # recalibrated OUTPUT; PASS/LowQual
            │    ├─ *_purity.out
            │    └─ *_somatic_haplotag.metrics
            └─ benchmark split (_sc PASS ∩ truth)
                 ├─ filtered_snv_tp.vcf.gz
                 └─ filtered_snv_fp.vcf.gz
                      └─ InterSubMod C++ TP/FP methylation analysis
                           ├─ per-region reads/methylation/distance/tree files
                           └─ significance_summary.csv
```

- `[F-L1]` orchestrator：`InterSubMod/scripts/pipeline/run_benchmark.sh:159-246` 呼叫 LongPhase-S step；`:274-298` 才呼叫 C++ step。
- `[F-L1]` LongPhase-S step：`InterSubMod/scripts/pipeline/steps/01_longphase_s.sh:124-153` 建立 PASS input；`:155-185` 執行 LongPhase-S；`:217-291` 由 `_sc.vcf` 建 TP/FP。
- `[F-L1]` C++ wrapper：`InterSubMod/scripts/pipeline/steps/02_intersubmod.sh:93-129` 建共同 CLI，`:144-188` 先 TP 後 FP，僅檢查程序 rc 與 summary 存在。

### 1.2 Layered reconstruction v2 actual chain

```text
layered_v2_input_manifest.json
  ├─ tumor_bam = historical LongPhase-S *_tagged.bam
  ├─ somatic_vcf = historical ClairS paired PASS somatic_pass.vcf.gz
  └─ CN BED / integer-CN BED (optional)
       ↓
run_layered_7samples_newbb.sh
  └─ 7 datasets, sample parallelism=2; each sample has 5 chromosome processes
       ↓
sm_multilocus_combinations.py
  ├─ autosomal PASS-SNV positional components
  ├─ densest contiguous MAX_SNV=8 cap
  ├─ tumor read REF/ALT vector census
  └─ split by LongPhase-S HP family
       ↓
layered_tree_reconstruction.py
  ├─ tree_enumeration_solver.py
  ├─ L0 family role
  ├─ L1 all-minimum tree enumeration + V1-V7
  ├─ L2 post-tree CN annotation
  └─ L3 = not_evaluated/bounded_auxiliary
       ↓
build_region_view.py
  └─ region-centric JSON + six-branch funnel
       ↓
verify_layered_v2.py
  └─ cross-file invariant summary
```

- `[F-L1]` actual entry：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh:53-145`。
- `[F-L1]` **沒有**任何 `inter_sub_mod`、`significance_summary.csv` 或 methylation JSON call；`layered_tree_reconstruction.py:329-348` 將 L3 固定為 `not_evaluated`。
- `[I-L2]` 因此「InterSubMod C++ 已驗證 layered v2 tree」是不成立的敘述；正確說法是 v2 使用歷史 LongPhase-S tagged BAM，另有 repo-local Python reconstruction stack。

## 2. 端到端 call chain 與每步判斷

### 2.1 `run_layered_7samples_newbb.sh`

| 步驟 | Source | 判斷／設定 | 輸出／失敗 |
|---|---|---|---|
| bootstrap | `InterSubMod/.../run_layered_7samples_newbb.sh:5-31` | Python 可執行、manifest schema=2.0、run root 必須不存在 | 建 immutable run root、copy manifest、status header；不會 resume。 |
| sample input | `:53-92` | BAM、`.bai`、VCF、宣告 CN files 存在 | `input_files.tsv`；SHA 只含 VCF、BAI、VCF index、CN，**不含 BAM content**。 |
| ML/HP | `:94-121` | 每 sample 5 splits；`MINREAD=3, MAX_SNV=8, TIER_R=50k, MAPQ=20, BASEQ=0` | 每 part 只查 schema + local funnel conservation。 |
| layered | `:123-130` | `VERIFY_EVERY=1, ANALYSIS_TREE_CAP=0, DISPLAY_TREE_CAP=32` | 只查 `all_eligible_V1V7_pass=true`。 |
| region view | `:132-143` | join ML raw evidence、重數 VCF universe | 只查 six-branch conservation；hash outputs。 |
| aggregate verify | `:172-180` | `xargs -P 2` 跑 7 samples，再跑 final verifier | `verification_summary.{json,tsv}`；僅全成功才寫 ALL verify PASS。 |

`[F-L1]` 未使用變數：top-level `SPLITS` (`:19-20`) 與 `PROFILE` (`:17`) 沒有實際 consumer；每 sample 又重建 local `splits`。

### 2.2 VCF universe 與 positional grouping

- `[F-L1]` `sm_multilocus_combinations.py:28-43` 用 `pysam.TabixFile.fetch(chrom)`，只保留 biallelic one-base SNV。
- `[F-L1]` 程式**不讀 FILTER 欄**，而是無條件在 tuple 寫 `PASS` (`:36-40`)；正式 runner 沒先跑 PASS-only validator，因此 mixed-filter VCF 會被當 PASS。
- `[F-L1]` Tabix chromosome/index錯誤被 `except (ValueError, OSError): return []` 靜默吞掉 (`:41-42`)；錯誤延遲到 region-view upstream count mismatch，local part 仍可能產合法的空 JSON。
- `[F-L1]` `sm_linkage_genomewide.py:145-156` 以**相鄰** sSNV gap >50 kb 切 component；total span 無上限。
- `[F-L1]` group >8 時，選 span 最小的連續 8 個 (`sm_multilocus_combinations.py:108-114`)；被排除 sSNV 記入 cap branch，不進 read census/solver。

### 2.3 Per-read allele census

- `[F-L1]` `sm_linkage_genomewide.py:83-114` 對整個 group span做一次 pileup，REF/ALT/OTHER 按 read QNAME聚合。
- `[F-L1]` MAPQ gate 是 20 (`:104-106`)；base quality由 runner傳 0，故 Q0 base 也可進樹底料。
- `[F-L1]` installed pysam `stepper="samtools"` 預設排 unmapped/secondary/QCfail/duplicate，**不排 supplementary**；程式也沒加 `flag_filter=BAM_FSUPPLEMENTARY`。
- `[F-L1]` 同一 QNAME 的 primary/supplementary/重複 alignment 會共用 `calls[rn]`、`hp[rn]` (`:107-113`)；重疊位置最後一次覆寫，跨 segment 則可能拼成一個 vector。
- `[O-L1]` HCC1395 `chr17:7,550,000-7,560,000` 小窗有 203 unique QNAME、35 個多 alignment QNAME、2 個 supplementary QNAME、1 個 primary+supplementary QNAME；這證明該分支不是純理論，但全基因組影響未量化。

### 2.4 HP family、read support 與保留

- `[F-L1]` `sm_multilocus_combinations.py:69-83`：`1*→family1`、`2*→family2`、`3→family3`、其他/缺值→`none`。
- `[F-L1]` 因使用 `startswith`，`1-1/1-2` 與 `2-1/2-2` 均正確 collapse 到 germline family；LongPhase-S 完整 enum 仍有 `HP:Z:4`，會落入 `none`，沒有 explicit H4 auxiliary。
- `[F-L1]` 每 read 產 R/A/X vector (`:129-147`)；full vector 與 partial subread 都計數。
- `[F-L1]` pooled full/partial 支持達 MINREAD=3 即保留 group (`:148-154`)；per-family 再各自做 MINREAD (`:155-161`)。
- `[I-L2]` **denominator hole**：三個不同 family 各 1 read 可能令 pooled pattern=3、group retained，但每 family 都 <3，於是 `populations_by_hp/subread_groups_by_hp` 全空。layered driver把該 group計入 `L0.regions_total`，卻產 0 unit；region view只從 detail建 `byreg`，因此完全不顯示該 region。現有 verifier沒有 `retained_groups == regions_with_or_without_units` 檢查。

### 2.5 CN 判斷

- `[F-L1]` CN file不存在→`None`，每 region→`unavailable` (`sm_multilocus_combinations.py:46-66`)。
- `[F-L1]` 有 segmentation 但該位置沒 interval→`neutral`；這依賴 manifest 的「non-neutral intervals + unlisted=neutral」契約。
- `[F-L1]` lookup 以 VCF 1-based position做 `s <= pos <= e` (`:63-65`)；BED通常是 0-based half-open，interval 兩端有 off-by-one 風險。
- `[F-L1]` solver driver 的 integer-CN lookup也用 inclusive (`layered_tree_reconstruction.py:38-41`)。
- `[F-L1]` `layered_tree_reconstruction.py:23-24` standalone default 指向 HCC1395 SEQC2 integer CN；正式 runner對無值傳空字串是安全的，但手動執行其他 sample 若忘記 env，會錯套 HCC1395 CN。

### 2.6 Tree solver

- `[F-L1]` full genotype `R/A`→ALT bit set；partial `X`→subcube (`tree_enumeration_solver.py:25-42`)。
- `[F-L1]` 覆蓋公理：所有 full node必存在、每 partial subcube至少命中一個 tree node (`:66-73`)。
- `[F-L1]` parent edge只能加一個 ALT bit；逐 hidden-node level搜尋最小 closed node set (`:98-165`)。
- `[F-L1]` hard limits：`extra_cap=4`、`per_level_budget=150,000`；超出改 greedy closure並標 `capped`，不是 exact answer。
- `[F-L1]` tree數以每 node可選 unit predecessor數乘積後跨 feasible node sets求和 (`:166-203`)；`ANALYSIS_TREE_CAP=0` 才生成完整候選集合。
- `[F-L1]` classification (`:238-253`)：有 recurrence-free minimal solution時依 tree count分 determined/ambiguous；全部 minimal solution需 recurrence才叫 recurrence_required；capped不允許 determined。
- `[F-L1]` V1-V3驗每棵 stored/generated tree；V4重算 minimal hidden；V5只重算 analytical count；V6/V7檢 flag/no-overclaim (`:255-429`)。

### 2.7 Driver role 與 denominator

- `[F-L1]` `layered_tree_reconstruction.py:115-210` 每 region×family產 unit。
- `[F-L1]` primary lineage = mutation-bearing HP1/HP2；reference-only = HP1/HP2/HP3但沒有 ALT；HP3 mutation-bearing = auxiliary；none = unphased auxiliary (`:137-152`)。
- `[F-L1]` L2 只在 mutation-bearing recurrence_required、non-capped後執行 (`:155-161`)；CN不改 L1 tree set。
- `[F-L1]` full/noncapped verification 必須 zero skipped/failed才 `full_pass` (`:162-174`)。
- `[F-L1]` display只存前32，但 analysis在 cap=0時用全部 trees算 exact shapes/digest (`:175-208`)。
- `[F-L1]` aggregate `eligible` 包含 primary、reference、H3、none 所有 non-capped units (`:268-300`)；`all_eligible_V1V7_pass` 不是只驗 primary lineage。

### 2.8 Region view 與 final verifier

- `[F-L1]` `build_region_view.py:20-33` 讀 raw ML JSON；任何 JSON exception被靜默 `continue`。
- `[F-L1]` region只由 `detail` 建 (`:47-52`)；因此 0-unit retained group消失。
- `[F-L1]` raw join找不到時，positions/coverage/populations會是 `None`，仍可輸出成功 (`:58-95`)。
- `[F-L1]` `underdetermined` 在 region-level mapping被放入預設 `A`，最後顯示 `has_ambiguous` (`:35-45`)；ambiguous與underdetermined語意被混合。
- `[F-L1]` six-branch conservation會重數 VCF SNV (`:123-160`)；它能抓染色體漏載，抓不到 retained group在region view消失，也不檢 FILTER。
- `[F-L1]` `verify_layered_v2.py:28-119` 檢角色、verification status、analysis completeness、funnel、missing CN、hp multiplicity；不重讀 BAM/VCF evidence、不核對 raw join completeness、不核對 stored checksum、不核對 BAM `@PG`/truth-bed。

## 3. 7-dataset input contract：手動 preflight 實測

### 3.1 執行

輸入：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json`

命令：

```bash
PYTHONDONTWRITEBYTECODE=1 /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_layered_v2_inputs.py \
  --manifest InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json \
  --output /tmp/20260710_layered_v2_preflight_audit.json
```

實際輸出：

```text
INPUT PREFLIGHT: 7/7 pass -> /tmp/20260710_layered_v2_preflight_audit.json
  HCC1395: PASS
  HCC1395_DORADO: PASS
  COLO829: PASS
  H1437: PASS
  H2009: PASS
  HCC1937: PASS
  HCC1954: PASS
```

`[F-L1]` validator 檢查：VCF gzip scan/SNV counts/filter set/source/version/cmdline (`validate_layered_v2_inputs.py:30-69`)；BAM header/reference dictionary/BAI hash (`:72-86`)；CN state (`:89-100`)；PASS-only與ClairS gate (`:103-161`)。

`[F-L1]` 正式 runner沒有任何 `validate_layered_v2_inputs.py` 引用；此 7/7 是 audit 手動執行，不能寫成 runner built-in preflight。

### 3.2 Preflight summary

| Dataset | biological_id | SNV total | chr1-22 | out-of-scope | FILTER | NAF=0 fraction | caller | CN source |
|---|---|---:|---:|---:|---|---:|---|---|
| HCC1395 | HCC1395 | 113,997 | 80,234 | 33,763 | PASS | 0.968078 | ClairS 0.4.0 | SEQC2 |
| HCC1395_DORADO | HCC1395 | 112,387 | 79,120 | 33,267 | PASS | 0.984616 | ClairS 0.4.0 | SEQC2 same cell line |
| COLO829 | COLO829 | 38,196 | 36,585 | 1,611 | PASS | 0.975469 | ClairS 0.4.1 | unavailable |
| H1437 | H1437 | 75,578 | 73,243 | 2,335 | PASS | 0.792797 | ClairS 0.4.1 | SAVANA |
| H2009 | H2009 | 157,405 | 150,370 | 7,035 | PASS | 0.988171 | ClairS 0.4.1 | SAVANA |
| HCC1937 | HCC1937 | 49,548 | 15,915 | 33,633 | PASS | 0.874667 | ClairS 0.4.1 | unavailable |
| HCC1954 | HCC1954 | 20,969 | 19,743 | 1,226 | PASS | 0.971768 | ClairS 0.4.1 | SAVANA |

- `[F-L1]` 7 datasets = 6 biological samples；HCC1395與DORADO不可當獨立生物 replicate。
- `[F-L1]` NAF=0只是一個 FORMAT descriptor；validator沒有把非零NAF判fail。
- `[I-L2]` H1437的NAF=0只有79.3%，直接反證「paired PASS 等同 NAF=0」；paired caller可在存在低量normal ALT evidence時仍判PASS。

## 4. P0：`somatic_pass` 與 `_sc.vcf` 不是同一物

### 4.1 Source-level role

| Artifact | 真正角色 | 產生點 | HCC1395實測 |
|---|---|---|---:|
| `somatic_pass.vcf.gz` | ClairS paired output經 `bcftools view -f PASS` 的 **LongPhase-S input** | `InterSubMod/scripts/pipeline/steps/01_longphase_s.sh:124-150` | 113,997 PASS SNV |
| `HCC1395_tagged_sc.vcf` | LongPhase-S recalibrated **output**，保留PASS/LowQual | `01_longphase_s.sh:155-168`; LongPhase-S docs `somatic_haplotag.md:120-139` | 108,530 PASS + 5,467 LowQual |
| `filtered_snv_tp/fp.vcf.gz` | `_sc.vcf` PASS與truth交集/差集；legacy C++ benchmark label | `01_longphase_s.sh:217-287` | 非v2 input |

HCC1395 canonical log直接證據：

- `[O-L1]` `/big7_disk/.../HCC1395/.../longphase_s/longphase_s.log:3-18`：tumor SNP input=`somatic_pass.vcf.gz`；somatic calling output=`HCC1395_tagged_sc.vcf`。
- `[O-L1]` 同 log `:41-53`：Tumor SNP count=113,997；somatic variant count(Flag)=108,530；之後寫 output VCF。
- `[O-L1]` `_sc.vcf` FILTER tally command輸出：`108530 PASS`, `5467 LowQual`。
- `[O-L1]` `somatic_pass.vcf.gz` header含 `##source=ClairS`, `##clairs_version=0.4.0`, paired `--tumor_bam_fn + --normal_bam_fn`，沒有 `##longphase_s_version`；`_sc.vcf` 多出 `##longphase_s_version=1.0.0` 與 LongPhase commandline。

### 4.2 Current wrapper filename collision

- `[F-L1]` `01_longphase_s.sh:126` 將 prefilter input命名為 `somatic_pass.vcf.gz`。
- `[F-L1]` 同 script `:235-246` 又把 `_sc.vcf` PASS寫回**同一路徑**。
- `[F-L1]` `:289-291` 最後刪掉該路徑及index。
- `[I-L2]` 因此 current wrapper完成後既不能保留原始 ClairS PASS input，也不能讓 v2 manifest穩定引用該檔；現有 canonical `somatic_pass.vcf.gz` 是歷史/另一步驟保留物，current wrapper本身無法以相同結果重產。

### 4.3 Required rename

建議 canonical role固定為：

```text
clairs_paired.pass.input.vcf.gz
longphase_s.recalibrated.all.vcf.gz      # PASS + LowQual
longphase_s.recalibrated.pass.vcf.gz
benchmark.filtered_snv_{tp,fp,fn}.vcf.gz
```

任何 manifest欄位必須寫 `role` 與 `producer`，不能只寫模糊的 `somatic_vcf`。

## 5. P0：6/7 tagged BAM 受 truth BED 限制

### 5.1 LongPhase-S source evidence

- `[F-L1]` `/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp:80-87` 先寫 `_sc.vcf`。
- `[F-L1]` 同檔 `:90-101` 接著在 tagRead前呼叫 `somaticBenchmark.removeVariantsOutBedRegion()`。
- `[F-L1]` `/big8_disk/.../SomaticBenchmark.cpp:508-543`：BED外若無NORMAL state則整個position erase；有NORMAL則移除TUMOR/TRUTH state、保留NORMAL。
- `[F-L1]` LongPhase-S KB也明載 `--truth-bed` 是「僅使用此區域內variants進行 tagging與評估」：`/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md:126-142`。

### 5.2 7 BAM `@PG` 與 log核對

命令：

```bash
samtools view -H <manifest.tumor_bam> | inspect @PG PN:longphase-s CL
rg -n 'removing tumor & truth somatic variants outside bed regions' <longphase_s.log>
```

| Dataset | `@PG` truth VCF | `@PG` truth BED | canonical log removal line | 結論 |
|---|---|---|---|---|
| HCC1395 | SEQC2 HC sSNV | SEQC2 HC BED | `longphase_s.log:53` | truth-BED restricted |
| HCC1395_DORADO | SEQC2 HC sSNV | SEQC2 HC BED | `longphase_s.log:53` | truth-BED restricted |
| COLO829 | NYGC truth | **NONE** | no removal line | 唯一無BED restriction |
| H1437 | orthogonal-tools truth | orthogonal-tools BED | `longphase_s.log:53` | truth-BED restricted |
| H2009 | orthogonal-tools truth | orthogonal-tools BED | `longphase_s.log:53` | truth-BED restricted |
| HCC1937 | orthogonal-tools truth | orthogonal-tools BED | `longphase_s.log:53` | truth-BED restricted |
| HCC1954 | orthogonal-tools truth | orthogonal-tools BED | `longphase_s.log:53` | truth-BED restricted |

完整 BAM header commandline以各 manifest `tumor_bam` 的 `@PG CL` 為 artifact evidence；HCC1395另可見 `/big7_disk/.../longphase_s.log:103` 完整命令。

### 5.3 HCC1395 scope magnitude

輸入：

- `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz`
- `/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`

命令與結果：

```bash
bcftools view -H -R /big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed \
  /big7_disk/.../longphase_s/somatic_pass.vcf.gz | wc -l
# 31295
```

- `[O-L1]` 只有 31,295/113,997=27.5% ClairS PASS records在SEQC2 HC BED；約82,702在BED外。
- `[I-L2]` BED外 read仍可能因normal germline SNP得到HP1/HP2，但somatic-derived HP1-1/2-1/3/4/1-2/2-2 evidence會被移除或降為germline/unphased。
- `[I-L2]` v2 primary tree仍直接用read allele，可在BED外產樹；問題是 family assignment/support、H3/none比例與read retention受benchmark BED選擇影響，而COLO829的規則不同。
- `[I-L2]` 這是 benchmark-induced feature selection／跨樣本契約不一致；不能以「只拿truth做metrics」帶過。

### 5.4 Required fix and acceptance test

1. production tagging完全移除 `--truth-vcf`、`--truth-bed`；benchmark metrics另開只讀評估步驟。
2. 7 BAM `@PG` 必須 machine-check：無 truth flags、相同LongPhase version、MAPQ、percentage threshold、supplementary policy。
3. HCC1395至少比較 old-vs-clean BAM：HP vocabulary/count、assigned rate、retained region、primary units、multi-HP、determinacy、recurrence；任一主比例≥5pp需標backbone-sensitive，≥10pp否證H7/H8。

## 6. CLI/default registry

### 6.1 v2 runner environment

| Variable | Default | 實際效果／限制 | Source |
|---|---:|---|---|
| `SM_PYTHON` | `/bip7.../python3` | 必須可執行 | runner `:6,22` |
| `SM_INPUT_MANIFEST` | bundled v2 JSON | schema只查2.0與dataset_count | `:7,23-24` |
| `SM_RUN_ID` | timestamp | run root immutable | `:8,25-28` |
| `SM_PARALLEL_SAMPLES` | 2 | 最多2 samples×5 parts=10 Python/BAM readers | `:10,98-108,173-174` |
| `SM_VERIFY_EVERY` | 1 | >1會產partial verification，runner layered gate會fail；實質上comprehensive只能1 | `:11,125,128` |
| `SM_ANALYSIS_TREE_CAP` | 0 | >0可使candidate set incomplete；final verifier會fail；comprehensive只能0 | `:12,125` |
| `SM_DISPLAY_TREE_CAP` | 32 | 僅UI storage | `:13,126` |
| `SM_MINREAD` | 3 | per pooled/per-family pattern支持門檻 | `:102,150` |
| `SM_MAX_SNV` | 8 | densest contiguous cap | `:102,150` |
| `SM_TIER_R` | 50,000 | adjacent gap，不是total span | `:103,150` |
| `SM_MAPQ_MIN` | 20 | per-read mapping quality | `:103,151` |
| `SM_BASEQ_MIN` | 0 | 接受Q0 allele，需sensitivity | `:104,151` |

### 6.2 LongPhase-S actual canonical settings

`[O-L1]` HCC1395 log `:21-34`：threads=24、tag region=all、MAPQ=20、percentage threshold=0.6、tag supplementary=enabled、purity=auto、variant filtering=enabled、write calling VCF=enabled。

`[F-L1]` tool defaults若直接執行則不同：MAPQ=1、threads=1、supplementary=false、threshold=0.6；見 `/big8_disk/.../docs/somatic_haplotag.md:57-96`。報告必寫**effective command**，不可只寫tool default。

### 6.3 InterSubMod C++ CLI/default

| Parameter | Config default | CLI effective/default | 問題 |
|---|---:|---|---|
| threads | 16 (`Config.hpp:43`) | help寫1 (`ArgParser.hpp:57`)，未給值實際仍16 | 文件與執行值矛盾 |
| distance metric | BERNOULLI (`Config.hpp:40`) | parser local default NHD並clear Config (`ArgParser.hpp:86-90,160-180`) | direct Config與CLI不同 |
| window | 1000 | 1000 | 一致；legacy wrapper覆寫5000 |
| MAPQ/read length/baseQ | 20/1000/20 | CLI可改 | v2 Python為20/無read-length/0，不同方法 |
| methyl high/low | 0.8/0.2 | CLI可改 | 一致 |
| min site/common CpG | 5/3 | common可改；site coverage無CLI | hidden hardcode |
| PMD gating | true | 無CLI、且pmd path無CLI | 行為不可由help完整控制 |
| clustering | on, UPGMA, min reads10 | 無CLI | hidden hardcode |
| expected coverage | 0=KDE | CLI可改 | legacy C++ only |
| germline HP only | false | flag | legacy C++ only |

- `[O-L1]` `./build/bin/inter_sub_mod --help`印help但exit code=1；`main.cpp:15-17` 把help與parse error都映成1。
- `[F-L1]` `Config::validate()` 對tumor index缺失只warning (`Config.cpp:34-38`)，但`BamReader` constructor會throw (`BamReader.cpp:31-37`)；normal index連warning都沒有。
- `[F-L1]` VCF loader只收PASS biallelic SNP並取first sample FORMAT.AF (`SomaticSnv.cpp:121-176`)；目前ClairS是single SAMPLE，契約成立，但未驗sample identity。
- `[F-L1]` per-region exception被記為failed result (`RegionProcessor.cpp:1375-1382`)；main不檢failed count，仍return0 (`main.cpp:63-81`)。
- `[F-L1]` summary開檔失敗只log後return (`RegionProcessor.cpp:1385-1394`)；main仍可能return0。
- `[F-L1]` summary只輸出success+significance-computed rows (`:1475-1479`)；wrapper只看檔案存在，可能把大量partial failure當成功。
- `[F-L1]` per-region artifacts在significance前先寫 (`RegionProcessor.cpp:836-849`)；下游失敗會留下看似完成的region files。

## 7. I/O schema registry

### 7.1 Inputs

| Role | Required fields/contract | 現有檢查 | 缺口 |
|---|---|---|---|
| tumor BAM | coordinate-sorted/indexed; GRCh38_no_alt; LongPhase HP string; primary/supp policy | exists `.bai`; manual validator hashes header/ref dict | runner不檢@PG、truth-bed、HP vocabulary、BAM content hash |
| somatic VCF | bgzip+index; ClairS paired; PASS-only; biallelic SNV; chr naming | runner只exists；manual validator較完整 | loader不讀FILTER、Tabix error靜默 |
| CN BED | chrom,start,end,state | exists+recognized state | coordinate convention未宣告/驗證 |
| integer CN BED | chrom,start,end,number | exists | coordinate convention、source compatibility未驗 |
| integration JSON | optional HCC1395 only | runner不validate | 缺檔時region view靜默少一欄 |

### 7.2 ML/HP part JSON

`[F-L1]` producer `sm_multilocus_combinations.py:210-222`：

- `schema_version`, `params`, `input_funnel`
- aggregate all/clean/neutral
- each group：coordinates、positions、cap counts、full/partial populations、per-family populations/subreads/coverage、reads_by_hp、CN、dominant HP summary。

`[F-L1]` runner只驗 `schema_version` 與 `input_funnel.check_scope_conservation`；沒有JSON schema、field type、group uniqueness、chrom split completeness。

### 7.3 Layered JSON

`[F-L1]` producer `layered_tree_reconstruction.py:329-348`：analysis contract、params、merged funnel、L0/L1/L2/L3、detail units。

每 unit含 role、class、tree counts、display trees、exact digest/shape、CN verdict、V1-V7 status、trace。完整候選樹可能很大；output以普通`json.dump(open(...))`非atomic write。

### 7.4 Region view JSON

`[F-L1]` producer `build_region_view.py:171-176`：sample、contract、L3、census、regions；每region內嵌lineages及raw observed coverage/populations/subreads。

`[F-L1]` 缺raw join不會fail；0-unit retained region不會出現；`region_determinacy=has_ambiguous`同時包含underdetermined。

### 7.5 Verification outputs

`[F-L1]` `verify_layered_v2.py:132-158` 產`verification_summary.json`與同stem TSV；final pass是所有內部checks pass，並非生物truth validation。

## 8. Error、fallback 與模糊行為矩陣

| ID | 場景 | 現行行為 | 風險 | 應改 |
|---|---|---|---|---|
| E01 | VCF index/chrom fetch error | return empty list | 延遲、難診斷 | fatal + sample/chrom/error |
| E02 | mixed FILTER VCF | FILTER不讀，全部當PASS | universe污染 | loader assert FILTER=PASS |
| E03 | BAM/VCF truth-bed scope不一致 | 不檢 | benchmark-selection confound | @PG contract gate |
| E04 | pooled support但0 per-family unit | retained、之後region消失 | denominator leak | explicit `no_family_supported` branch |
| E05 | raw ML JSON parse fail | silently continue | region evidence變None仍success | fatal + exact file |
| E06 | solver level太大 | greedy fallback + capped | 非exact |可保留，但須coverage/cap報告 |
| E07 | analysis cap>0 |候選不完整；初步V1-V7仍可能pass | config語意誤導 | comprehensive mode拒絕>0 |
| E08 | verify_every>1 | partial pass；runner再fail |可配置但不可成功 | comprehensive mode固定1 |
| E09 | CN absent | unavailable | 正確 | verifier維持 |
| E10 | CN interval miss | neutral |需依manifest contract | validate coverage semantics |
| E11 | runner被kill | status停在START，無trap | run state不可信 | EXIT/TERM trap寫ABORTED |
| E12 | concurrent status append |多process直接`>>`同檔 |小行通常但無明確鎖 | flock或per-sample status |
| E13 | Python JSON write interruption | truncated output |後續可能parse fail；無atomic | temp+fsync+rename |
| E14 | C++ region failures | log後main可能0 | wrapper假成功 | fail count/threshold→nonzero |
| E15 | C++ summary write fail | log+return，main可能0 |無summary或stale file | throw/propagate |

## 9. Findings register

| ID | Severity | Finding | Evidence | Fix / acceptance |
|---|---|---|---|---|
| U01 | **P0 blocker** | 6/7 tagged BAM受truth BED限制，VCF卻是genome-wide；COLO規則不同 | §5 source/log/@PG `[O-L1]` | clean re-haplotag 7/7；@PG no truth flags；rerun sensitivity |
| U02 | **P0 blocker** | `somatic_pass`被誤稱LongPhase output；input 113,997 vs recalibrated PASS108,530 | §4 `[O-L1]` | role-specific filenames+manifest producer |
| U03 | **HIGH** | 強preflight未接正式runner | runner `:22-24,67-74`; validator全檔 `[F-L1]` | runner先執行validator並保存lock |
| U04 | **HIGH scientific** | 文件把ClairS PASS/NAF=0寫成「真somatic/真值」 | `CURRENT_FOCUS.md:33`; data-model spec `:19-23,80,150-151`; `verify_newbb_numbers.py:27-54` | 改「operational caller universe」；NAF只作descriptor |
| U05 | **HIGH** | pooled-supported、per-family unsupported group可從region denominator消失 | ML `:148-161`; layered `:248-266`; view `:47-52` | group/unit/region conservation test |
| U06 | **HIGH method** | tree allele底料預設BASEQ=0且無sensitivity | runner `:103-105`; pileup `:96` | prereg BaseQ10/20 grid；報flip率 |
| U07 | **HIGH provenance** | current pipeline wrapper重用並刪`somatic_pass`，無法重產manifest artifact | `01_longphase_s.sh:124-150,235-246,289-291` | 永不覆寫/刪 canonical input；immutable artifacts |
| U08 | **HIGH observability** | active v2 run停在START、無程序、無FAIL/completion；runner無trap | run root status；runner全文 | trap + heartbeat + final state |
| U09 | **HIGH release** | v2核心與tests多為modified/untracked，HEAD無法重建；full run未完成 | git status；`00_INDEX.md:46-53` | commit code/manifest/test後再run |
| U10 | **MED** | HP4落none；完整9-state vocabulary沒preflight；prevalence未量化 | KB `longphase-s.md:197-211`; family code `:69-83` | explicit HP enum + unknown counter/gate |
| U11 | **MED** | supplementary alignment按QNAME merge/覆寫 | pileup code + pysam docs + small-window observation | define molecule policy；conflict counter/test |
| U12 | **MED** | BED 0/1-based inclusive boundary不明 | CN lookup code | normalize BED half-open，unit tests at boundaries |
| U13 | **MED** | standalone driver默認HCC1395 integer CN | layered `:23-24` | default empty；require explicit source |
| U14 | **MED** | region view raw parse error靜默、underdetermined併ambiguous | view `:20-45` | fail-fast；separate category |
| U15 | **MED** | C++ CLI default/help與Config不一致；help rc=1 | §6.3 | single source defaults；help rc=0 |
| U16 | **MED** | C++ partial failure/summary failure不傳exit code | main/RegionProcessor | completeness threshold+nonzero |
| U17 | **MED docs** | README/QUICKSTART只描述legacy C++ per-SNV流程，不含v2 | `README.md:33`; `QUICKSTART.md:33-105` |新增v2 manual；清楚標兩條chain |
| U18 | **MED stale verifier** | `verify_newbb_numbers.py`讀舊pilot/multisample outputs、glob第一個canonical run；missing sample不令allok false | script `:20-25,59-76` | manifest/run-root only；missing=failure |
| U19 | **MED stale census** | `funnel_census_7samples.py`硬寫僅HCC1395有CN、讀legacy outputs，與v2 5/7 CN manifest矛盾 | script `:17-33,132-140`; v2 manifest | retire/supersede；only current run summary |
| U20 | **LOW** | Python tests有unclosed-file ResourceWarning、JSON非atomic | test output；open/json.dump calls | context managers+atomic writer |

## 10. NAF、PASS 與「truth」的正確說法

- `[O-L1]` HCC1395第一個PASS record `chr1:94893 A>G` 的 `NAF=0.0789`；caller仍判PASS。
- `[F-L1]` HCC1395有110,358/113,997 NAF=0，仍有3,639 NAF非0；H1437 NAF=0比例只有79.3%。
- `[I-L2]` paired caller整合normal evidence後的PASS是**caller判定**，不是每筆normal必須零ALT，更不是正交truth。
- `[F-L1]` `verify_newbb_numbers.py:42` 更有一個細節：缺NAF時default字串是`"0"`，會把missing NAF算成NAF=0；manual validator則default `nan`，兩者口徑不同。
- `[F-L1]` 正確報告用語應為：

> 「Primary operational backbone is the PASS SNV universe emitted by paired ClairS. It uses matched-normal evidence but is not an orthogonal biological truth set. FORMAT.NAF=0 is reported as a QC descriptor, not a truth criterion.」

## 11. Build/test correspondence

### 11.1 Commands and outputs

```bash
bash -n InterSubMod/docs/.../run_layered_7samples_newbb.sh \
  InterSubMod/scripts/pipeline/run_benchmark.sh \
  InterSubMod/scripts/pipeline/steps/{00_prepare_germline,01_longphase_s,02_intersubmod}.sh
# exit 0, no output
```

```bash
PYTHONDONTWRITEBYTECODE=1 /bip7.../python3 \
  InterSubMod/docs/.../test_layered_reconstruction_v2.py
# Ran 5 tests in 0.079s — OK
# ResourceWarning: unclosed ml.json / layered.json
```

```bash
PYTHONDONTWRITEBYTECODE=1 /bip7.../python3 \
  InterSubMod/docs/.../tree_enumeration_solver.py
# 8/8 golden cases PASS; V1-V7 all checkmarks
```

```bash
ctest -N --test-dir InterSubMod/build
# Total Tests: 234
```

- `[F-L1]` `InterSubMod/CMakeLists.txt:143-178` 只註冊C++ GoogleTest；v2 Python tests與runner contract tests不在CTest。
- `[F-L1]` `InterSubMod/tests/test_snv_loading.cpp` 是獨立main且未列入CMake；不構成自動VCF ingestion regression。
- `[U-L5]` 本稽核未跑234個CTest，只確認discovery；不能寫「C++ full test pass」。

### 11.2 Missing tests

1. truth-bed-containing BAM `@PG` must fail production manifest。
2. `somatic_pass` vs `_sc` role/header/filter count。
3. mixed FILTER VCF and broken/missing Tabix chromosome。
4. pooled-supported but zero-family-unit group。
5. HP4/1-2/2-2/unknown vocabulary。
6. supplementary overlap with conflicting allele/HP。
7. BASEQ 0/10/20、MAPQ20/30、MINREAD3/4、MAX_SNV8/6 sensitivity。
8. CN BED start/end exact boundary。
9. raw ML JSON corrupt/missing join。
10. SIGTERM/interruption status transitions and atomic JSON。
11. C++ all-regions-failed/summary-unwritable must return nonzero。

## 12. Step → Verify remediation plan

1. **重產clean LongPhase-S BAM** → 驗證：7/7 `@PG`無truth flags；HP enum/count與BAM checksum保存。
2. **拆分VCF artifact roles** → 驗證：input header無LongPhase；recalibrated header有LongPhase；counts 113,997 vs108,530/5,467可重算且manifest role明示。
3. **wire preflight** → 驗證：故意混入LowQual、缺index、truth-bed BAM、錯reference dictionary都在任何BAM pileup前exit nonzero。
4. **修group conservation** → 驗證：`retained_groups = regions_with_units + no_family_supported_regions`，7/7 exact equality。
5. **定義pileup quality/molecule policy** → 驗證：BASEQ sensitivity與supplementary conflict counters進verification summary。
6. **修runner durability** → 驗證：SIGTERM測試後status=`ABORTED`、無半成品被當PASS、可用新run ID安全重跑。
7. **fresh full run** → 驗證：7 datasets/6 biological samples、5 ML parts each、7/7 final pass、checksums readback、所有H1-H9 thresholds有值。
8. **同步README/method spec/HTML** → 驗證：全文搜尋不再把`somatic_pass`稱`_sc`、不再把NAF=0/PASS稱truth、C++與v2 call chain分開。

## 13. Unverified boundaries

- `[U-L5]` 未全BAM掃描HP4/1-2/2-2與supplementary conflict prevalence。
- `[U-L5]` 未量化truth-bed restricted vs clean tagged BAM對7-dataset主結果的差異；目前只有source-level因果與HCC1395 scope magnitude。
- `[U-L5]` 未完成7-dataset v2 reconstruction；無法宣告runtime、memory、all-pass或H7/H8 robustness。
- `[U-L5]` 未跑完整C++ test suite或production BAM C++ smoke。
- `[U-L5]` 未以single-cell/multi-region truth驗證biological clone；claim ceiling仍是regional mutation-state inference。
- `[U-L5]` external second caller只在部分dataset存在；不可宣告cross-sample caller robustness。

## 14. 稽核命令與實際輸出摘要

| Command | Input | Output | Actual excerpt |
|---|---|---|---|
| manual v2 preflight | manifest + 7 BAM/VCF/CN | `/tmp/20260710_layered_v2_preflight_audit.json` | `7/7 pass` |
| `_sc` FILTER tally | HCC1395 tagged_sc.vcf | stdout | `108530 PASS; 5467 LowQual` |
| `somatic_pass` count/header | HCC1395 somatic_pass | stdout | `113997`; source ClairS0.4.0 paired |
| HC region count | HCC1395 somatic_pass + SEQC2 BED | stdout | `31295` |
| BAM `@PG` audit | 7 manifest BAMs | stdout | 6 truth-bed; COLO none |
| log removal audit | 7 canonical logs | stdout | 6×line53 removal; COLO no line |
| HP small-window sample | HCC1395 tagged BAM chr1/chr17 windows | stdout | common `1,2,1-1,2-1`; supplementary branch present |
| Shell syntax | 5 wrappers | stdout | rc0 |
| v2 tests | test_layered_reconstruction_v2.py | stdout | 5/5 OK + ResourceWarning |
| solver golden | tree_enumeration_solver.py | stdout | 8/8 ALL PASS |
| CTest discovery | build | stdout | 234 tests listed |
| C++ help | build/bin/inter_sub_mod | `/tmp/20260710_inter_sub_mod_help.txt` | help printed, rc1 |

## 15. 一句結論

**目前核心solver的小型案例可運作，但上游BAM受truth BED限制、VCF角色混稱、正式preflight未接線與full run未完成，使整體仍不能被稱為無模糊且可外部重現的完整流程；首要修正是重產clean tagged BAM並鎖死input/output contract，而不是先美化HTML。** `[I-L2]`
