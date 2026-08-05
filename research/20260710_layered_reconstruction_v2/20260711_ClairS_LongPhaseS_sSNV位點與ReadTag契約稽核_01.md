<!--
建立時間: 2026-07-11 04:10 Asia/Taipei
目標: 確認 raw ClairS 全位點、LongPhase-S input/output record keys、tagging scope 與 sSNV HP/PS consumer 是否完整
處理範圍: 7 datasets（6 biological samples）；normalized paired ClairS raw-all → LongPhase-S bidirectional recalibration → `_sc` PASS tree input + exact HP/PS sidecar；歷史 PASS-only run另作 scope audit
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/clairs_longphase_ssnv_contract_audit.json; InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/audit_clairs_longphase_ssnv_contract.py
任務類型: B — Comprehensive validation
服務目標: G2/G3/G4/G5
狀態: raw_all_production_authorized_pending_start
證據等級: L1 artifact/code audit；兩個 bounded probes與 patched regression已過，7-dataset raw-all production/closeout仍未完成
-->

# ClairS → LongPhase-S → sSNV 位點與 Read Tag 契約稽核

用 SCQA + Step→Verify：先直接回答「是否全部讀取與正確使用」，再分開呈現 caller record、tagging scope、sSNV disposition 與 read-tag consumer。

> **TL;DR（2026-07-11 07:00 現行）**：完整 LongPhase-S recalibration universe 已更正為 **normalized paired ClairS raw-all biallelic sSNVs**；同一 run 的 `_sc.vcf FILTER=PASS` 是 canonical tree input，ClairS PASS 只作 sensitivity/audit。PASS-only producer 在 0/7 完成時中止。HCC1954 chr22 與 patched HCC1395 whole-chrX probes 已通過，但 7/7 raw-all production、closeout、U0–U7 與 fresh downstream run 均尚未完成。歷史 tagged BAM 另只有 COLO829 是 genome-wide scope，其餘 6/7 受 truth BED 限制。**（影響：高；信心：高）**

> **2026-07-11 07:00 方法更正（優先於下文歷史 PASS-input 描述）**：LongPhase-S 原生可把 input non-PASS rescue 成 PASS；因此只餵 ClairS PASS 雖有 record continuity，卻不是完整 all-site recalibration。0/7 完成的 PASS-only producer已中止且無 `_SUCCESS`。HCC1954 chr22 normalized raw-all probe為 426→426、rescue 53、remove 14、sidecar 217,023 rows exact；original binary 在 HCC1395 chrX:72880028 no-read LowQual site exit 1，copied-source最小 patch使該位點保留 LowQual且不中止，並在 HCC1954 chr22 得到 byte/record/purity/transition-identical regression。patched HCC1395 whole chrX 已完成 35,823→35,823、LowQual→PASS 959、PASS→LowQual 1,616、2 zero-read warnings、1,090,570 sidecar rows、unknown/conflict=0。pre-decision audit重評為 `GO_WITH_FAIL_CLOSED`，意義是**允許啟動** 7-dataset raw-all producer，不是 production 已完成。

## 1. 直接判斷

1. **完整 production contract 是否應把 raw ClairS 全部 biallelic sSNV records 送 LongPhase-S？是。** 現行 producer會 header-normalize paired ClairS raw-all，不先做 `-f PASS`；但 7/7 full production尚未完成。歷史 wrapper 先過濾 PASS 的行為已淘汰。
2. **歷史 ClairS PASS 位點是否完整送入當時的 LongPhase-S？是。** 7/7 historical datasets 的 raw-PASS record multiset 與 `somatic_pass.vcf.gz` 完全相等，missing=0、extra=0；這只證 subset continuity。
3. **LongPhase-S 能否在 raw-all 上保留所有 records並雙向改 FILTER？bounded probes 是。** HCC1954 chr22與patched HCC1395 whole chrX均 input=`_sc all`；7/7仍待驗。
4. **LongPhase-S `_sc.vcf` 在 historical PASS-only run 是否保留所有 input keys？是。** 7/7 historical input 與 recalibrated-all record multiset相等；但這不能涵蓋 non-PASS→PASS rescue branch。
5. **歷史 tagged BAM 是否對所有 genome-wide PASS 位點做一致的 somatic tagging？否。** 6/7 command 帶 `--truth-bed`，LongPhase-S 在 tag reads 前移除 BED 外 tumor state；只有 COLO829 無此限制。
6. **舊 sSNV consumer 是否使用完整正確 tag 契約？否。** 舊版以 QNAME 合併 alignment、HP4 落入 `none`、PS 不落輸出、raw fine HP 不保存。這四項已在 clean consumer 修正並通過 chr22 smoke，7-dataset clean rerun 尚在執行。
7. **演化樹應使用哪一層 sSNV？** production 主樹使用同一 normalized raw-all run 的 LongPhase-S `_sc.vcf FILTER=PASS`；ClairS PASS 只作 sensitivity/audit subset。此決策由使用者於 2026-07-11 明確確認。

### 1.1 建樹輸入的合理性判斷

| 候選輸入 | 角色 | 判斷 |
|---|---|---|
| normalized paired ClairS raw-all | 完整 LongPhase-S recalibration universe／lossless ledger | 不可直接建主樹；含 caller non-PASS，須由 LongPhase-S 重校準 |
| ClairS `FILTER=PASS` | historical continuity + backbone sensitivity/audit subset | 不再作完整 LongPhase-S input；否則 non-PASS→PASS rescue branch不可達 |
| LongPhase-S `_sc.vcf` all | 完整 recalibrated output／filter audit | 不可直接建主樹；含 LongPhase-S LowQual |
| LongPhase-S `_sc.vcf FILTER=PASS` | integrated LongPhase-S sSNV backbone | **canonical 主樹輸入** |

這個選擇的前提是研究主張為「ClairS 候選經 LongPhase-S recalibration 與 HP/PS tagging 後的 read-level 演化重建」。同一個 no-truth、normalized-raw-all LongPhase-S run 必須同時提供 `_sc.vcf` 與 exact alignment sidecar，不能把另一個受 truth BED 限制、PASS-only 或不同參數 run 的 tags 接上來。

限制是 LongPhase-S FILTER 與 read co-occurrence evidence 都來自同批 reads，存在 selection coupling；因此 canonical 結果是 operational integrated pipeline，不是獨立 caller／獨立 phasing confirmation。正式結果必須另跑 7-dataset ClairS-PASS backbone sensitivity，並報告 retained-site、primary-unit 與 topology digest 的差異，不能只報 canonical tree。

`_sc FILTER=PASS` 定義的是 canonical sSNV input universe，不保證每個位點都形成一條樹邊。每筆位點仍會依 scope、連通 component、MAX_SNV、read support 與 retained-group 狀態得到唯一 disposition。現階段可支持的產物是 bulk long-read 的局部 read-supported mutation candidate trees；沒有 single-cell／multi-region orthogonal truth 時，不可升格成已確認的全域 clone phylogeny。

## 2. 歷史 PASS-only 全 dataset record 與 scope 結果（只作 continuity/scope audit）

| Dataset | raw ClairS | historical raw PASS=LPS input | historical LPS `_sc` all | historical `_sc` PASS | PASS key完整 | truth BED內 | truth BED外 | production scope |
|---|---:|---:|---:|---:|---|---:|---:|---|
| HCC1395 | 134,122 | 113,997 | 113,997 | 108,530 | PASS | 31,295 | 82,702 | FAIL |
| HCC1395_DORADO | 123,240 | 112,387 | 112,387 | 107,986 | PASS | 30,574 | 81,813 | FAIL |
| COLO829 | 43,014 | 38,196 | 38,196 | 37,458 | PASS | n/a | n/a | PASS |
| H1437 | 83,261 | 75,578 | 75,578 | 74,498 | PASS | 68,250 | 7,328 | FAIL |
| H2009 | 168,638 | 157,405 | 157,405 | 155,326 | PASS | 133,948 | 23,457 | FAIL |
| HCC1937 | 59,161 | 49,548 | 49,548 | 46,500 | PASS | 13,171 | 36,377 | FAIL |
| HCC1954 | 26,823 | 20,969 | 20,969 | 20,433 | PASS | 18,101 | 2,868 | FAIL |

機械來源：`InterSubMod/research/20260710_layered_reconstruction_v2/clairs_longphase_ssnv_contract_audit.tsv`。

7/7 datasets 的 raw ClairS、raw PASS、LongPhase-S input 與 `_sc.vcf` 中，`records = biallelic_snvs`；因此上表沒有把 indel 或 multiallelic record 混入 sSNV 分母。這項判斷來自同一份 audit JSON 的 `counts.records` 與 `counts.biallelic_snvs` 逐層相等，不是由工具名稱推定。

### 2.1 「全部位點」與 autosomal tree scope

歷史 LongPhase-S run 接收 caller PASS VCF（含 chrX/chrY）；現行 production contract 則接收 normalized raw-all（同樣含 chrX/chrY），而 canonical tree scope 預註冊為 chr1–22。raw ClairS 所有 contigs 仍逐筆寫入 ledger；進到 LongPhase-S 但不屬 chr1–22 的 sSNV 會明示為 `out_of_scope_non_autosomal`，不會靜默遺失或誤塞入 diploid-autosome topology。下表仍是 historical PASS subset census，不是尚未完成的 raw-all 7/7 count。

| Dataset | raw non-autosomal | ClairS PASS/LPS-input non-autosomal |
|---|---:|---:|
| HCC1395 | 35,829 | 33,763 |
| HCC1395_DORADO | 35,388 | 33,267 |
| COLO829 | 1,894 | 1,611 |
| H1437 | 2,877 | 2,335 |
| H2009 | 7,459 | 7,035 |
| HCC1937 | 35,965 | 33,633 |
| HCC1954 | 1,513 | 1,226 |

此排除是分析 scope，而不是 LongPhase-S ingestion loss。chrX/chrY 若要進樹，需另定性別、ploidy、PAR、LOH/CN 與 HP-family 解讀；在這些契約未建立前直接混入會比明示排除更不可靠。

機械來源：同一份 `clairs_longphase_ssnv_contract_audit.json/.tsv` 的 `non_autosomal_records` 欄位；7/7 的 raw-PASS non-autosomal、LongPhase-S input non-autosomal 與 `_sc`-all non-autosomal counts 相等。

## 3. VCF 角色不能再混稱

```text
raw ClairS output.vcf.gz
  └─ header/sample normalization（record keys/payload守恆）
       └─ normalized paired ClairS raw-all（LongPhase-S完整 input）
            └─ LongPhase-S recalibrated-all _sc.vcf
                 ├─ PASS → primary sSNV tree input
                 └─ LowQual → excluded_by_longphase_filter（仍留 ledger）

ClairS FILTER=PASS → historical continuity + backbone sensitivity/audit only
```

`somatic_pass.vcf.gz` 是歷史 ClairS PASS input；`*_tagged_sc.vcf` 才是 LongPhase-S output。兩者 record keys 相等不代表角色相同，也不能用歷史 key equality 取代 raw-all non-PASS rescue/remove audit。

## 4. Read tag 完整語彙與 clean consumer

LongPhase-S L1 程式碼定義 9-state `HP:Z:`：`./1/2/3/4/1-1/2-1/1-2/2-2`。

| raw HP | clean sSNV role | primary denominator |
|---|---|---|
| `1/1-1/1-2` | parental family 1；raw fine tag另存 | 只有 mutation-bearing unit 才進 |
| `2/2-1/2-2` | parental family 2；raw fine tag另存 | 只有 mutation-bearing unit 才進 |
| `3` | unresolved H3 auxiliary | 否 |
| `4` | shared H4 auxiliary | 否 |
| `.`/missing | unphased | 否 |

clean alignment identity：`QNAME + chrom + start + end + FLAG + BLAKE2b8(CIGAR)`。HP、PS、MAPQ 與 primary/supplementary status 由 FIFO tagged-BAM stream逐 alignment 保存，不以 QNAME 最後一筆覆寫。**HP 進 family stratification；PS 保存為 phase-block QC，不作 topology edge 或 lineage label。**

本輪「tag 正確」的可證明範圍是 producer-to-consumer fidelity：LongPhase-S 寫出的 HP/PS 能逐 alignment 無遺失、無衝突地被 sSNV consumer讀取；九種 HP 語彙不被錯誤遺失，且實際進分析的是 HP family。PS 不承擔樹拓樸，避免把 phase-set 標籤誤當獨立共現證據。這不是 LongPhase-S 生物相位 assignment 的獨立 truth validation；後者仍需 orthogonal phasing／single-cell 證據。

HCC1954 chr22 bounded probe 的 31 個 retained regions 中，30 個只有 1 個 PS、1 個含 2 個 PS（chr22:46,039,334/46,055,696）；因此正式全量輸出新增 per-region `phase_set_counts`、`phase_set_HP_counts` 與跨樣本 mixed-PS census。證據：`InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_HCC1954_chr22_longphase_pass_contract_v1/phase_set_policy_probe.json`。

HCC1954 chr22 smoke 實際結果：

```text
LongPhase parser input       = 20,969 PASS SNP
streamed mapped alignments   = 217,023
LongPhase tagged alignments  = 124,449
sidecar tagged alignments    = 124,449
unknown HP                   = 0
exact identity conflicts     = 0
sSNV alignment exposures     = 4,017  (_sc PASS tree input)
sidecar exact matches        = 4,017
sidecar missing/conflict     = 0/0
```

## 5. 所有位點輸出

`build_ssnv_site_ledger.py` 對 raw ClairS 每筆 record 保留 ID/REF/ALT/QUAL/FILTER/INFO/FORMAT，並給唯一 disposition：

- `excluded_before_longphase_nonPASS`
- `excluded_by_longphase_filter`
- `out_of_scope_non_autosomal`
- `positional_singleton`
- `max_snv_excluded`
- `read_unsupported`
- `retained`

retained 位點另帶 pooled/family REF-ALT read counts、raw HP census、HP-with-PS census、unique PS 數、實際 PS count、HP×PS count 與 LongPhase recalibrated FILTER。

HCC1954 chr22 **historical PASS-only** smoke ledger：raw ClairS 426、ClairS PASS/LPS input 181、`_sc.vcf` all 181、`_sc.vcf` PASS/tree input 167；disposition 為 pre-LPS non-PASS 245、LongPhase LowQual 14、singleton 99、retained 68。這只驗 consumer/ledger branch。現行 raw-all probe另為 426→426，LowQual→PASS 53、PASS→LowQual 14；完整證據在 `20260711_HCC1954_chr22_raw_all_longphase_contract_v1/`。

## 6. Step → Verify

1. raw-all production LongPhase-S 移除 truth flags → 驗證：7/7 logs 的 truth VCF/BED 空白，無 benchmark-removal line。
2. normalized raw-all input → 驗證：raw biallelic keys/payload=input，input=`_sc` all，四格 FILTER transition逐 key守恆。
3. FIFO 擷取 exact HP/PS sidecar → 驗證：execution alignment/tagged counts 與 sidecar逐筆相等；unknown/conflict=0。
4. 原始 BAM + sidecar 做 sSNV → 驗證：每 split exact matches=exposures、missing=0、conflict=0。
5. `_sc` PASS=tree input；ClairS PASS另跑 sensitivity → 驗證：retained-site、primary-unit、topology digest分開報。
6. clean chr1-22 × 7 rerun → 驗證：7/7 producer receipt、U0–U7與fresh verifier全 PASS 後才恢復 comprehensive claim。

## 7. 現況

- 歷史資料稽核：完成，僅 1/7 production-scope PASS。
- PASS-only 7-dataset producer：在 0/7 terminal completion中止，改名 `20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1`；`_FAILED` 存在、`_SUCCESS` 不存在、無殘留 producer process。
- HCC1954 chr22 raw-all、patched regression與 HCC1395 whole-chrX raw-all probes：完成，全部 bounded hard gates PASS。
- pre-decision verdict：`GO_WITH_FAIL_CLOSED`（允許啟動）；**7-dataset raw-all production目前未啟動/未完成，無 production `_SUCCESS`**。
- 08:32 execution update：raw-all producer/receipt/layered-v3 canonical+sensitivity gate **84/84 PASS**，跨版本 gate **45/45 PASS**。另已凍結 7 個既有 paired-full LongPhase-S tagged BAM 的 sampled-content identity baseline（set SHA-256=`ce6c63d42e3f334d6847a1a6d3e46ead165b59a03197acb098319be5c67fcf90`）；正式 run 後須 exact-match 才能證明未覆蓋。producer 仍在 resource queue，最近 load1=62.95、iowait=64%，target run root 尚未建立，不能寫成已啟動或已完成。
- 舊 layered run：保留作 upstream-mismatched engineering baseline，不作 comprehensive evidence。

## 8. 可重跑命令

```bash
python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/audit_clairs_longphase_ssnv_contract.py \
  --manifest InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json \
  --output InterSubMod/research/20260710_layered_reconstruction_v2/clairs_longphase_ssnv_contract_audit.json
```

```bash
# 授權但尚未執行完成；runner會拒絕既有 immutable RUN_ROOT並 fail-closed。
SM_RUN_ID=20260711_longphase_s_raw_all_production_sidecars_v2 \
SM_PARALLEL_SAMPLES=2 SM_LPS_THREADS=12 \
bash InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_longphase_raw_all_production_sidecars.sh
```

## 9. 來源

- KB：`/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md`，HP 9-state、`--truth-bed` tagging scope、`--tagSupplementary`。
- LongPhase-S source：`SomaticHaplotagProcess.cpp:80-101`，先寫 `_sc.vcf`、再移除 BED 外 tumor state、最後 tag reads。
- InterSubMod producer：`InterSubMod/scripts/pipeline/steps/01_longphase_s.sh`。
- InterSubMod consumer：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_linkage_genomewide.py` 與 `sm_multilocus_combinations.py`。

## 10. 2026-07-11 07:00 checkpoint 與下一 gate

- Pre-decision audit：`research/20260710_layered_reconstruction_v2/pre-decision-audit.md` §9.5，re-score **60/100**，verdict=`GO_WITH_FAIL_CLOSED`。
- Patch scope：只把 low-confidence/no-eligible-read fatal改為 warning+no-op；保留 LowQual；尚非 upstream-reviewed release。
- Bounded evidence：HCC1954 chr22 original 426→426；patched regression byte/record/purity/transition identical；HCC1395 whole chrX patched 35,823→35,823、2 zero-read warnings、sidecar 1,090,570 rows、unknown/conflict=0。
- 下一個可改 claim status 的事件不是「runner啟動」，而是 7/7 immutable raw-all producer receipts + `_SUCCESS`、U0–U7與 fresh downstream verifier 全部通過。完成前 `ISM-M03/R03` 維持 `pending_rerun`。
