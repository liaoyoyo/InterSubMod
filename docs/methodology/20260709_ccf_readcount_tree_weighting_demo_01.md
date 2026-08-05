<!--
建立時間: 2026-07-09
類型: 歷史小規模 demo — family-specific read-AF ordering（舊稱 CCF）
狀態: superseded_by_layered_v2_read_af_ordering
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_region_integration.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/ccf_tree_weighting_demo.json
provenance: 數字由 scripts/ccf_tree_weighting_demo.py 於 2026-07-09 對 sm_region_integration.json populations(含計數)重算。partial flag: HCC1395 / ambiguous 子集。
-->

# Family-specific read-AF ordering — 歷史小規模 demo

> 🔴 **2026-07-10 語意更正**：本頁舊稱「CCF」的量是未經 purity／CN/multiplicity 校正的 family-specific ALT read fraction，**不是真 cancer cell fraction**；winner-clean 也是同一 pigeonhole 規則下的自洽度，非獨立驗證。新的權威入口為 `read_af_tree_ordering_multisample.py`：只納入 mutation-bearing HP1/HP2、對所有候選樹重枚舉，並報 temperature／posterior／margin sensitivity。本頁保留作歷史 demo，不再支撐 current claim。

> **歷史提案**:把 enumerate_min_trees 的「N 棵等機率樹」用每群 read 數（per-mutation read AF）加權。約束不變（RRR 根 + 有向 unit-flip），只用 read count。原理：pigeonhole/crossing heuristic 期待祖先突變 read AF ≥ 後代突變 read AF。

## §1 方法
- 每突變 read AF = Σ(含該突變的基因型 read 數) / total。
- tree score = Σ_edge Σ_(anc∈parent)(readAF[anc] − readAF[acquired]);softmax(TEMP=0.05)→ posterior。
- violation = 後代 read AF > 祖先 read AF + margin(0.05)。

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

## §4 裁決(v1 full-only;§6 有含 partial 完整版)
> **有效、合理、非循環、約束不變、winner 生物一致** —— v1 full-only reach 2.3%(下界);**§6 含 partial 的真實 reach = 66.8%(CN-clean 可信部分 22.3%)**,值得 CN-gated 落地。

## §6 全面觀察 v2 — 含 partial reads(col_coverage CCF)〔L1〕
`ccf_tree_weighting_full_observe.py`:用 mlhp `col_coverage_by_hp`(每位點 nALT,含 partial reads=更多 read)算 CCF,對 layered solver 的 full+partial ambiguous 樹加權。HCC1395 **5959 ambiguous 單位**:

| 指標 | 值 |
|---|---|
| **tie 被破(≥0.6)** | **3980 = 66.8%**(🔺 vs full-only 2.3%) |
| winner pigeonhole-clean | **3956/3980 = 99.4%** |
| top posterior ≥0.9(強贏家) | 2219 |

**🔴 按 CN 分類(決策關鍵)**:
| ambiguity × CN | 總 | 破 | 破% | 可信 |
|---|---|---|---|---|
| structure \| **CN-clean** | 1960 | **1331** | **67.9%** | ✅ read≈CCF |
| structure \| **CN-gain** | 3999 | 2649 | 66.2% | 🔴 read≠CCF(multiplicity)不可信 |

- **「更多 read → 更能定」驗證**:col_coverage(含 partial)深度遠高於 full-only → reach 2.3%→66.8%。
- **可信 reach = CN-clean 1331 / 全 5959 = 22.3%**(CN-clean ambiguous 67.9% 乾淨解出)。
- 全部 ambiguous 為 `structure(多完成)` 型(partial 帶進多完成歧義)。

**§7 決策裁決**:✅ **值得落地,CN-gated** —— CN-clean 區乾淨解出 67.9% ambiguity(winner 99.4% 生物一致)= 明顯可觀測收斂;🔴 CN-gain(67%)保持誠實(不排序/標 multiplicity-confounded)。落地=`tree_enumeration_solver` 枚舉後加 soft CCF-pigeonhole 加權層(CN-gated、對稱保持等權重、beta-binomial)。一鍵重驗 `ccf_tree_weighting_full_observe.py`。

## §5 一鍵重驗
```bash
cd InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 ccf_tree_weighting_demo.py
```
只讀不寫源碼,重算 tie-broken 2.3% / winner-clean 100% / CN-gain 54%。2026-07-09 實跑確認。

**關聯**:layered solver `tree_enumeration_solver.py`;形式化 `20260704_formal_problem_statement...`;VAF pigeonhole 是 §8 允許輔助通道(甲基不是);gap 文件 `20260709_subcube_recovery_gap_quantification_01.md`。
