<!--
建立時間：2026-07-19
目標：記錄 compressed VAF ranking P1 fail-closed 修補與 focused oracle。
處理範圍：Task A exploratory pilot；synthetic m<=5；09b hard 未執行。
關聯檔案：receipt.json、receipt.json.sha256、../../../scripts/compressed_vaf_rank_probe.py、../../../tests/test_compressed_vaf_rank_probe.py
-->

# Compressed VAF ranking P1 fail-closed 驗證

> **PARTIAL。** 本輪只證明小型結構DP與full-enumeration current-likelihood
> differential oracle；普通浮點UB剪枝只可作diagnostic，不是machine certificate。
> 09b hard未執行，沒有hard-case或production speedup結論。

## Assertion–Evidence

**Assertion：float-pruned結果不再發布authoritative winner/ties。**

- Core固定`numerical_bound_certified=false`。
- 只要`pruned_candidate_count>0`，就回
  `INCOMPLETE_UNCERTIFIED_FLOAT_BOUND_PRUNING`、
  `ranking_complete=false`、`exact_publishable=false`、
  `best_log_likelihood=null`與空winner IDs。
- 只有`evaluated_candidate_count==logical_family_count`、`pruned=0`且tie輸出
  完整時，才可發布current machine endpoint結果。

**Evidence：m=5共有120 candidates。**

- Full evaluation：120/120實評、0 prune；與獨立parent-vector enumeration加
  exhaustive current scorer的best score、24-way tie IDs完全相同。
- Float diagnostic：實評24、剪96；因沒有directed-rounding certificate，
  正確fail closed，不發布best/ties。

## Focused validation

命令：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest -v \
research/20260718_solver_methyl_edge_probe/tests/test_compressed_vaf_rank_probe.py
```

實際結果：14/14 PASS，5.953s，exit 0。

測試涵蓋：

- m=1…5 independent parent-vector family oracle。
- branching、mandatory full、mixed R/A/X、gapped active bits。
- branch family互斥、count守恆、每個completion包含於該branch
  `possible_vertices`。
- candidate cap、tie cap、strict control types。
- LL/gap NaN、Inf與負gap拒絕。
- union/rowwise普通float pruning fail closed。
- source snapshot含`hypercube_exact.py`並可偵測hash drift。

## Receipt

輸出：

- `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v2/p1_hardened_r2/receipt.json`
- `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v2/p1_hardened_r2/receipt.json.sha256`

Runner wall=`5.764620s`；verdict=
`PASS_SMALL_ORACLE_HARD_NOT_RUN_P1_GATE`。兩檔mode=`0444`、nlink=`1`，
`sha256sum -c` PASS。Receipt綁定current M2 ranker、`hypercube_exact.py`、
prototype、counter、R3 validator hashes及NumPy `1.23.5`；寫出使用
`allow_nan=false`。

## 複雜度與claim ceiling

- DP count：`O(3^m poly(m))`。
- 目前會先materialize全部top-level branches與其`possible_vertices`，最壞約
  `O(m*3^m)`；不是recursive B&B。
- 沒有directed-rounding certificate時，machine-complete ranking仍需
  `Theta(F*Fit)`全候選實評。
- 若T個候選同分，完整顯式tie輸出仍有`Omega(T)`下界。
- Core deadline不涵蓋DP、全部branch/possible-vertex建構與完整單次optimizer，
  不是hard wall bound；因此hard固定`NOT_RUN_P1_REVIEW_GATE`。

## P2／未決

- selected hard unit尚未固定唯一ID，也尚未逐欄驗unit與case structural payload
  完全相等；在補完前hard必須保持NOT_RUN。
- 現有source/input pre/post hash可偵測一般drift，但尚未升級為same-FD、
  inode/O_NOFOLLOW等更強TOCTOU證明。
- 若要讓剪枝也machine-publish，需真正outward-directed rounding或interval
  likelihood certificate；本輪未擴張實作。
