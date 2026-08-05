<!--
建立時間: 2026-07-19T21:58:26+08:00
目標: 以目前 worktree 的 machine-readable artifacts 獨立重算 latent molecular substructure 比例、分母與 claim ceiling
處理範圍: Task Type B；7 datasets / 6 biological samples / chr1-22 / 469,849 LongPhase-S recalibrated FILTER=PASS biallelic sSNV dataset-sites
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/20260719_latent_molecular_substructure證據分層與比例_01.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_site_audit.tsv.gz
-->

# Latent substructure 比例與 claim ceiling 獨立 agent 稽核

> **Verdict: APPROVE。** `469,849 / 459,928 / 102,842`、all-site M2
> `919 / 1,867`、positional singleton `50,432` 及 singleton
> `30 PASS + 18 FAIL + 5,913 NOT_EVALUABLE + 44,471 NOT_RUN` 均由 frozen
> TSV/JSON 獨立重算一致，比例與分母口徑正確。此 APPROVE 只涵蓋本文件所列
> 數字與受限 claim，不代表正式 Task-B downstream/release 完成，也不批准
> biological prevalence、confirmed cellular subclone 或 linear ancestry claim。

## 1. 稽核框架、scope 與前提

使用「來源綁定 -> raw row 重算 -> 分母矩陣 -> 分層偏差 -> claim ceiling」稽核
矩陣。本任務是 **(B) Comprehensive validation**，服務 G2/G3/G4；統計單位是
`(dataset, chrom, pos, ref, alt)`，不是 cell、clone、patient 或獨立 biological
sample。

| 啟動問題 | 判定 |
|---|---|
| Thread D read-level epigenetic 相關？ | 是；M1/M2 都是 focal-ALT carrier reads 的甲基分群證據。 |
| Thread B 撤回範圍？ | 未重新開啟已撤回的 FP-filter claim；truth 分層只作 post-hoc 稽核。 |
| KDE-corrected？ | ⚠ 本次必讀 artifacts 與 M1/M2 gate 未宣告 KDE-correction contract；本稽核只核對 frozen row、gate 與分母，不把 KDE-corrected 當已證明前提。 |
| 需要 VCF caller AF？ | 不需要。`caller_af` 雖存在 site TSV，但不進 positional、M1 或 M2 分母判定。 |
| 長計算 / C++ / 搬移 / NO-GO gate？ | 執行全量唯讀 Python/TSV replay；未跑 C++、未跑 producer、未刪除或搬移檔案。 |

關鍵假設：

1. `all_ssnv_site_results.tsv.gz` 的 469,849 rows 是本稽核的 formal site
   universe；manifest、summary 與 primary audit 對同一 SHA-256 綁定。
2. `M1 stable` 由 row 欄位獨立套用
   `m1_evaluable AND coarse_ng>=2 AND NOT unstable AND
   modal_assignment_ari_min>=0.8`。
3. `M2-determinate = eligible + evaluable_ineligible`；`NOT_EVALUABLE`
   不是 biological negative，`NOT_RUN` 也不是 FAIL。
4. `truth_label` 是 variant benchmark stratum，不是 clone truth 或 lineage
   truth。

## 2. 實際輸入與版本綁定

題示的
`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/positional_singleton_summary.v4.json`
在目前 worktree **不存在**。目前可讀且有 atomic marker 綁定的 singleton
machine-readable source 是：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_audit_summary.json`

它的 directory/schema 是 `v1_source_attested / 1.0.0`，但 source chain 記錄其
生成時驗證的是 v4 source authority；不能把不存在的檔名當版本證據。

| 輸入 | SHA-256 | 用途 |
|---|---|---|
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/20260719_latent_molecular_substructure證據分層與比例_01.md` | `2d3dbbc92d3902147ae53a535913e0d8ee54846db1680611e373bf52a330e3e2` | 待審 claim |
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json` | `fae5380ebf7c468e1e712b9dda8ea473bc8fcfd30a36f5794769aa50c8ca91b2` | 7-dataset universe / truth / branch contract |
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json` | `4a0f411b9ce39128d07e5d15abfc5a5c999797b8a3f1eb33886185e830031daf` | all-site M2 reference recount |
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md` | `da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af` | claim ceiling |
| `.../all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_summary.json` | `454949afd1081f7e9475a3893d5ef0bde37c7e0d2b553fd231c8d8ae9fdf1f80` | all-site summary |
| `.../all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz` | `a8871af3a8c3955bf31aec5eeef0c93aca0683f52cf6d6f1e06fbbb713324f74` | 469,849-row獨立重算 |
| `.../all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_stable_multigroup_read_assignments.jsonl.gz` | `82af076d8ce70c66c7f75c4ed4bdeacb7d4c5c43d0905859303b121df483a4ba` | 102,842 M1 assignment proofs / M2 replay |
| `.../positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_audit_summary.json` | `6d4128be0535f2e16b1382cf038c558054ac42d0d7bde75ab4854b7a5be7aedc` | singleton summary |
| `.../positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_site_audit.tsv.gz` | `2d5d24790918fca34a32c313b9965f1de8c186c031e31d1a643febba721d90ce` | 50,432-row singleton重算 |

`_SUCCESS.json` 對 singleton summary 與 site TSV 記錄的 digest，與本次
`sha256sum` 完全一致。

### Source drift

singleton artifact 內記錄的 producer auditor SHA-256 是 `5b259c7b...23695f`；
目前 worktree 同一路徑
`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/audit_positional_singleton_methyl_multigroup.py`
為 `d4263f51...ab65`，且該檔目前是 untracked，無法由 artifact 所列 git HEAD
`0ee2fa1b...` 用 `git show` 重建。相對地：

- current independent M2 auditor=`3e63aad0...63d9`，與 artifact 綁定值一致。
- current M2 production gate=`f3d93238...3b79`，與 M2 recount reference 一致。
- current all-site analyzer=`390228ce...0b2`，與 v10 pinned analyzer identity一致。

因此，frozen artifacts 的數值 identity 是完整的，但「用 current singleton
auditor 精確重播 artifact producer」不是成立的 provenance claim。本次以自寫
唯讀重算避開此 drift；它是 residual reproducibility limitation，不是數值
mismatch。

## 3. 唯讀命令與實際輸出

主要輸入路徑如上；所有重算輸出只到 stdout，沒有中間檔。執行命令：

```bash
sha256sum \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json \
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz \
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_site_audit.tsv.gz

gzip -cd \
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz \
  | sed -n '1,3p'

python3 -I -B - <<'PY'
# Standard-library csv/gzip replay:
# 1. stream all 469,849 rows;
# 2. recompute M1 from coarse_ng/unstable/modal_assignment_ari_min;
# 3. sort each (dataset,chrom), split only when adjacent gap > 50,000;
# 4. compare recomputed component_size=1 keys with singleton TSV;
# 5. count singleton M1/M2 by dataset and truth;
# 6. import logic-independent audit_independent_m2_gate_recount.py and call
#    recount(site_results, stable_assignments) in memory.
PY
```

核心 stdout：

```text
all_site sites=469849 evaluable=459928 m1_stable=102842
all_site_m2 eligible=919 evaluable_ineligible=948
all_site_m2 not_evaluable_axis_indeterminate=100974 group_count_gt10=1
positional singleton_count=50432 keyset_symmetric_difference=0
singleton m1_evaluable=48347 m1_flagged=5961
singleton M2 PASS=30 FAIL=18 NOT_EVALUABLE=5913 NOT_RUN=44471
mismatches duplicate_keys=0 m1_formula=0 positional_membership=0
```

## 4. 全 universe 分母矩陣

| 指標 | 分子 / 分母 | 重算比例 | 分母口徑 | 判定 |
|---|---:|---:|---|---|
| focal-ALT evaluable | 459,928 / 469,849 | 97.888471% | all-universe | 一致 |
| M1 stable yield | 102,842 / 469,849 | 21.888309% | all-universe | 一致；不是 prevalence |
| M1 stable conditional | 102,842 / 459,928 | 22.360456% | M1-evaluable | 一致 |
| M2-determinate | 1,867 / 102,842 | 1.815406% | M1-flagged | 一致；`919+948` |
| M2 PASS conditional | 919 / 1,867 | 49.223353% | M2-determinate | 一致；高度選擇分母 |
| M2 PASS / M1 | 919 / 102,842 | 0.893604% | M1-flagged | 一致 |
| M2 PASS / all | 919 / 469,849 | 0.195595% | all-universe | 一致；operational yield |

守恆：

```text
102,842 M1
= 919 M2 PASS
 + 948 M2 determinate FAIL
 + 100,974 axis-indeterminate
 + 1 group-count > 10
```

`919/1,867` 不能替代 `919/102,842` 或 `919/469,849`；四個比例回答不同問題。

## 5. Positional singleton 獨立重建

程式定義與本次重算皆為：

1. 依 `(dataset, chrom)` 分組。
2. 依 `(pos, ref, alt)` 排序。
3. 相鄰位置 gap `<=50,000 bp` 留在同一 component；只有 gap `>50,000`
   才斷開，因此是 transitive component。
4. `component_size=1` 才是 positional singleton。

469,849 rows 重建得到 50,432 keys，與 singleton TSV 的 50,432 keys 完全相同；
雙向 set difference 都是 0，minimum finite singleton nearest gap=`50,003 bp`。

| 指標 | 分子 / 分母 | 重算比例 | 分母口徑 |
|---|---:|---:|---|
| singleton / all | 50,432 / 469,849 | 10.733661% | all-universe |
| M1-evaluable / singleton | 48,347 / 50,432 | 95.865720% | all singleton |
| M1 FLAGGED / singleton | 5,961 / 50,432 | 11.819876% | all singleton |
| M1 FLAGGED / evaluable | 5,961 / 48,347 | 12.329617% | singleton M1-evaluable |
| M2-determinate / FLAGGED | 48 / 5,961 | 0.805234% | singleton M1-flagged |
| M2 PASS / FLAGGED | 30 / 5,961 | 0.503271% | singleton M1-flagged |
| M2 PASS / singleton | 30 / 50,432 | 0.059486% | all singleton |
| M2 PASS / determinate | 30 / 48 | 62.500000% | singleton M2-determinate |

完整守恆：

```text
50,432
= 44,471 M1 NOT_FLAGGED / M2 NOT_RUN
 + 5,913 M1 FLAGGED / M2 NOT_EVALUABLE
 + 18 M2 FAIL
 + 30 M2 PASS
```

### Positional singleton 不等於 read-sharing degree-zero

此判定只使用位置與 50 kb component contract，未建立 read-to-site bipartite
graph。兩者不等價：

- ONT read 可能跨越 `>50 kb` 並同時覆蓋另一個 sSNV，所以 positional
  singleton 不保證 read-sharing degree=0。
- 反之，兩個位置落在同一 `<=50 kb` component，也不保證有任何一條 read
  同時覆蓋兩者。
- actual degree-zero 必須由正式 read identity/cooccurrence table 對每個 focal
  site 枚舉共同 reads 後重算。

因此 `50,432` 只能稱 positional singleton；G1/G2/formal R1 在此 supplemental
artifact 中全為 `NOT_RUN`。

## 6. Dataset 與 truth 分層

### 6.1 All-site dataset

| Dataset | All | M1-evaluable | M1 | M1/evaluable | M2 determinate | M2 PASS | determinate/M1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| COLO829 | 37,788 | 35,484 | 3,579 | 10.0862% | 0 | 0 | 0% |
| H1437 | 77,080 | 74,961 | 10,187 | 13.5897% | 41 | 16 | 0.4025% |
| H2009 | 154,465 | 152,594 | 54,644 | 35.8101% | 1,603 | 816 | 2.9335% |
| HCC1395 | 79,687 | 78,629 | 12,838 | 16.3273% | 37 | 4 | 0.2882% |
| HCC1395_DORADO | 79,739 | 78,637 | 14,789 | 18.8067% | 68 | 13 | 0.4598% |
| HCC1937 | 18,690 | 17,886 | 1,938 | 10.8353% | 86 | 50 | 4.4376% |
| HCC1954 | 22,400 | 21,737 | 4,867 | 22.3904% | 32 | 20 | 0.6575% |
| **Total** | **469,849** | **459,928** | **102,842** | **22.3605%** | **1,867** | **919** | **1.8154%** |

來源稿的 all/evaluable/M1 per-dataset table與 raw TSV 完全一致。新增的 M2
分層揭示 selection bias：H2009 單一 dataset 佔 M2-determinate
`1,603/1,867=85.8597%`，也佔 PASS `816/919=88.7922%`。因此 pooled
`49.2234%` 幾乎是 H2009 加權結果，不是 6 個 biological samples 的平均。

### 6.2 All-site truth

| Truth | All | M1-evaluable | Evaluable/all | M1 | M1/evaluable |
|---|---:|---:|---:|---:|---:|
| TP | 335,296 | 330,896 | 98.6877% | 74,312 | 22.4578% |
| FP | 7,745 | 4,967 | 64.1317% | 582 | 11.7173% |
| UNASSESSED | 126,808 | 124,065 | 97.8369% | 27,948 | 22.5269% |

FP 的 evaluability 明顯低於 TP/UNASSESSED，不能把 truth-stratified M1 比例直接
解讀為生物差異或 specificity。

### 6.3 Singleton dataset/truth

來源稿的 7-dataset singleton、M1-evaluable、M1 FLAGGED、M2 PASS 數全部與
singleton TSV 一致。singleton/all 在 dataset 間由 H2009 `1.8470%` 到
HCC1937 `45.3130%`，顯示該分支強烈受 callset site density 與 dataset 結構
影響，不是可跨樣本平均的 biological frequency。

| Truth | Singleton | M1-evaluable | M1 FLAGGED | M2 P / F / NE / NR | Evaluable/singleton |
|---|---:|---:|---:|---:|---:|
| TP | 45,193 | 44,171 | 5,494 | 30 / 16 / 5,448 / 39,699 | 97.7386% |
| FP | 1,084 | 514 | 52 | 0 / 1 / 51 / 1,032 | 47.4170% |
| UNASSESSED | 4,155 | 3,662 | 415 | 0 / 1 / 414 / 3,740 | 88.1348% |

全部 30 M2 PASS 是 TP 的數值陳述正確，但 FP 的 M2-determinate 只有 `n=1`；
這不能估 specificity，也不能說 M2 已區分 TP/FP。TP/FP 是 variant truth，
更不是 subclone truth。

### 6.4 技術資料集重現

HCC1395 與 HCC1395_DORADO 是同一 biological sample。由 singleton TSV
重新建立 biological locus set：

| Tier | Intersection | Union | Jaccard |
|---|---:|---:|---:|
| positional singleton | 7,484 | 9,116 | 82.0974% |
| M1 FLAGGED | 407 | 1,289 | 31.5749% |
| M2 PASS | 0 | 4 | 0% |

來源稿的 `82.10% / 31.57% / M2 intersection=0` 一致；目前沒有 locus-level
strict M2 technical replication。

## 7. Claim ceiling

| Claim | Verdict | 理由 |
|---|---|---|
| operational M1 stable multigroup yield | 可稱 | 明確限定 dataset-site 與 evaluable/all 分母。 |
| M2 measured-axis residual read-level epigenetic partition | 可稱 | 919 all-site、30 singleton；仍可能有未測 confound。 |
| latent molecular substructure compatible/candidate-generating signal | 可稱但須受限 | 只能表示與該模型相容；不可變成 cellular identity。 |
| biological prevalence | **不可稱** | 分母是選擇後的 dataset-sites，不是 cells、clones、patients 或獨立樣本。 |
| confirmed cellular subclone | **不可稱** | 缺至少第二個獨立 genetic marker 的同 read co-membership、normal/REF、CN/purity/CCF 與正交 cellular evidence。 |
| linear ancestry / lineage order | **不可稱** | 單一 focal ALT + methyl groups 無法區分 within-clone state、branching 或 linear order。 |
| positional singleton = actual read-sharing degree-zero | **不可稱** | 前者是位置 component，後者需要 read-sharing graph。 |

`confirmed cellular subclone=0` 與 `linear ancestry=0` 只能解讀為「目前零筆
達到確認門檻」，不能解讀為真實 biological prevalence 是 0%。來源稿若使用
「0% confirmed」，必須保留這個語義。

## 8. Mismatch 與 selection-bias table

| 項目 | 來源 claim | 獨立重算 / 檢查 | 判定 |
|---|---|---|---|
| formal universe | 7 datasets / 469,849 | 7 / 469,849；duplicate key=0 | MATCH |
| focal-ALT evaluable | 459,928 | 459,928 | MATCH |
| M1 stable | 102,842 | 102,842；formula mismatch=0 | MATCH |
| all-site M2 | 1,867 determinate / 919 PASS | 948 FAIL + 919 PASS | MATCH |
| positional singleton | 50,432 | 50,432；key-set symmetric difference=0 | MATCH |
| singleton M1/M2 | 48,347 / 5,961；30/18/5,913/44,471 | 完全相同且守恆 | MATCH |
| per-dataset | 來源稿兩表 | raw TSV 全部一致 | MATCH |
| singleton truth | TP/FP/UNASSESSED 三列 | raw TSV 全部一致 | MATCH |
| HCC1395 technical Jaccard | 82.10% / 31.57% / 0% | 82.0974% / 31.5749% / 0% | MATCH |
| 題示 singleton source path | `results/positional_singleton_summary.v4.json` | 不存在；current 是 output workspace `v1_source_attested/...summary.json` | RESOLVED PATH MISMATCH |
| singleton producer source | artifact=`5b259c...` | current worktree=`d4263f...`，且 untracked | PROVENANCE DRIFT；數值非 mismatch |
| `0 confirmed` | 零筆確認 | 不是 biological absence/prevalence=0 | SEMANTIC LIMIT |
| pooled M2 ratio | 919/1,867=49.2234% | H2009 佔 determinate 85.86% / PASS 88.79% | SELECTION BIAS |
| singleton truth comparison | 30 PASS 全 TP | FP determinate n=1、FP evaluability 47.42% | SELECTION BIAS |

## 9. 結論與 residual limitations

**APPROVE** `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/20260719_latent_molecular_substructure證據分層與比例_01.md`
中的核心數字、分母區分、positional-singleton 限定與「不支持 confirmed subclone /
linear ancestry」結論。

Residual limitations：

1. M2-determinate 是極小且 dataset-skewed 的選擇子集；`919/1,867` 與
   `30/48` 都不能外推。
2. 7 datasets 只有 6 biological samples，且 HCC1395 技術對沒有共同 M2 PASS
   locus。
3. truth missingness 差異大；所有 singleton PASS 為 TP 不足以估 specificity。
4. M2 只處理八個 measured axes，未排除 CN、purity、cis-ASM、未測 technical
   或 biological state。
5. actual read-sharing degree-zero、正式 G1/G2/R1、matched-normal、CN/CCF
   與 cellular/lineage truth 仍未由此 supplemental 分支證明。
6. singleton artifact 的精確 producer auditor source 未在 current worktree
   保持相同 SHA；後續正式 release 若要宣稱 exact code replay，需另有
   content-addressed frozen source 證據。

本稽核未執行 producer，未建立除本 review 外的任何檔案。
