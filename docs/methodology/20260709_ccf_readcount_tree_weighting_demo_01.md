<!--
建立時間: 2026-07-09
類型: 小規模 demo — read-count(CCF pigeonhole)給等機率樹加權;使用者提案驗證
狀態: demo 驗證(有效但 reach 有限);可一鍵重驗
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_region_integration.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/ccf_tree_weighting_demo.json
provenance: 數字由 scripts/ccf_tree_weighting_demo.py 於 2026-07-09 對 sm_region_integration.json populations(含計數)重算。partial flag: HCC1395 / ambiguous 子集。
-->

# read-count CCF 加權給等機率樹 — 小規模 demo(使用者提案)

> **提案**:把 enumerate_min_trees 的「N 棵等機率樹」用每群 read 數(→per-mutation CCF)加權 → 加權 posterior + 較可能者。約束不變(RRR 根 + 有向 unit-flip),只用 read count(遺傳、非甲基、非循環)。原理:pigeonhole/crossing rule 祖先突變 CCF ≥ 後代突變 CCF。

## §1 方法
- 每突變 CCF = Σ(含該突變的基因型 read 數) / total。
- tree score = Σ_edge Σ_(anc∈parent)(CCF[anc] − CCF[acquired]);softmax(TEMP=0.05)→ posterior。
- violation = 後代 CCF > 祖先 CCF + margin(0.05)= 明確違反 pigeonhole。

## §2 結果(HCC1395,1528 ambiguous 區,n_trees≥2 非 capped)〔L1〕
| 指標 | 值 |
|---|---|
| **tie 被破(top posterior≥0.6)** | **35(2.3%)** |
| 保持等權重(<0.6) | 1493(97.7%) |
| **winner pigeonhole-clean(0 violation)** | **35/35 = 100%** |
| CN-gain ambiguous(read≠CCF 須標記) | 824(54%) |

範例(CN-clean 破成功):chr11:22267628 CCF 0.17/0.52 → posterior 1.0(高頻當祖先);chr1:64969029 CCF 0.47/0.25 → 1.0。

## §3 判斷:有效合理〔L1 + L4〕
- ✅ **有效**:確實把等機率變加權 posterior;**贏家 35/35 全 pigeonhole-consistent**(生物有據);純 read count(非甲基、非循環);約束不變。
- ⚠ **reach 有限 = 2.3%,是資料本質非方法弱**:97.7% 等機率是**真頻率對稱**(單一雙突變、中間群未觀測 → 兩突變 CCF 相等)→ read 破不了、**也不該破**(對稱硬排=造假)。read 只破「中間群被觀測 + 單突變頻率不對稱」的那批。
- 🔴 **CN 前提**:54% ambiguous 是 CN-gain,read≠CCF(multiplicity)→ 必須排除/標記(否則把 Gap B 多拷貝假象當演化)。
- **soft/hard(TEMP)**只調贏家 sharpness,不改「能破哪些」(由 CCF 是否對稱決定,TEMP-independent)。

## §4 裁決
> **有效、合理、可大規模應用,但是誠實的小幅精修**(破 ~2.3% ambiguous、贏家永遠生物一致、CN-gated),非整批收斂。大部分等機率**正確地保持等機率**(真對稱)= 誠實優點。

**落地建議(若要)**:在 `tree_enumeration_solver` 枚舉後加 CCF-pigeonhole 加權層(soft beta-binomial、CN-gated、對稱保持等權重),輸出每棵樹 posterior。⚠ full-only demo 是 reach 下界;含 partial(col_coverage)的完整 reach 需 mlhp 輸入,為 follow-up。

## §5 一鍵重驗
```bash
cd InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 ccf_tree_weighting_demo.py
```
只讀不寫源碼,重算 tie-broken 2.3% / winner-clean 100% / CN-gain 54%。2026-07-09 實跑確認。

**關聯**:layered solver `tree_enumeration_solver.py`;形式化 `20260704_formal_problem_statement...`;VAF pigeonhole 是 §8 允許輔助通道(甲基不是);gap 文件 `20260709_subcube_recovery_gap_quantification_01.md`。
