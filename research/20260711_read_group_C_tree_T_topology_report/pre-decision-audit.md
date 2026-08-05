<!--
建立時間: 2026-07-11 16:48 +0800
目標: 校正 read-group C、候選樹 T、拓撲形狀與 VAF ranking 後，建立 HCC1395 範例及全 7 datasets HTML 數據報告
處理範圍: chr1-22 × 7 datasets；historical layered-v2 engineering snapshot；clean layered-v3 尚未完成
cycle_id: 20260710-2244-layered-reconstruction-v2 (derivative report correction)
topic: read_group_C_tree_T_topology_report
status: verdict_GO_descriptive; verdict_PROBE_independent_VAF_confirmation
audit_version: 0.1
build_branch: research/subclonal-reconstruction-202606
build_commit: 6067568637088838a9f518955e41d222f057e4f1
worktree: /big7_disk/liaoyoyo2001/InterSubMod
data_sources: /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2, InterSubMod/research/20260710_layered_reconstruction_v2/input_manifest_lock.json
關聯檔案:
  - InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md
  - InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/20260710_分層重建端到端流程程式碼判斷稽核_01.md
-->

# Pre-Decision Audit：read-group C、tree T 與 topology census

> **Verdict**：描述性工程報告 **GO**；將 read-AF/VAF 寫成「獨立確認真樹」則維持 **PROBE**。

## 2026-07-11 21:04 extension：0.50 門檻與並列第一 census

- **研究問題**：依使用者定案，以嚴格 `top weight >0.50` 比較 0.60，並判斷原 0.60 以下的 14,497 regions 可拆成多少 unique top、tied top/same Topo、tied top/different Topo。
- **Task type / scope**：B comprehensive validation；chr1-22 × 7 dataset rows，historical layered-v2 engineering snapshot。
- **假設 H1**：嚴格 `>0.50` 可保證單棵候選取得過半 relative weight；`.50/.50` 不通過並保留為並列，仍須與 `n_top` 聯合呈現。
- **假設 H2**：最高 score 等價集合的 `n_top` 與其 canonical shape 數可由 frozen raw groups + complete solver candidate set 決定性重算。
- **替代方案**：A=`>=0.60` 原 gate；B=`>=0.50 + n_top/n_top_shapes`；C=只報 rank/gap，不設 priority gate。傾向先比較 A/B，不在結果前替換 canonical gate。
- **Tie 定義**：主分析由整數 REF/ALT counts 建立 `Fraction` exact read-AF 與 exact raw ordering score，以有理數完全相等定義 co-top；另以 absolute score gap `1e-12`、`1e-9`、`1e-6` 作 near-tie 敏感度。softmax weight 只作排序強度，不是 calibrated probability。
- **Region 規則**：各 primary HP unit 分別取 top set；region top-tree combinations 為各 unit top-set 的 Cartesian product，數量可相乘，但不相乘 VAF weights。joint Topo 以帶 family 的 canonical-shape tuple 計數。
- **成功條件**：重現既有 0.60 region counts；candidate/shape count 全部吻合；14,497、29,053、39,885 三層守恆；三個 tolerance 的 headline class 穩定或差異被完整揭露。
- **失敗條件**：candidate mismatch、shape mismatch、分母不守恆、或 tie 結果對數值 tolerance 高度不穩定。
- **Verdict**：**GO for descriptive recomputation**；任何「0.50 是生物真值機率門檻」或「VAF 確認真樹」維持 **NO-GO**。

## §0 Cynefin front-gate

- **Domain**：Complicated（可重現 census）＋ Complex（VAF 對真實演化樹的生物解讀）。
- **理由**：raw group、候選樹與 AHU shape 可決定性重算；VAF 排序受 CN、purity、multiplicity 與 read-depth 影響，不具相同輸入必得真樹的保證。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 7 datasets historical mlhp/layered raw JSON 存在 | ✓ | L2 ⭐⭐⭐⭐ | `InterSubMod/research/20260710_layered_reconstruction_v2/input_manifest_lock.json` |
| `C` 可由 `populations_by_hp[HP]` 中含 A 的 full-span read groups重算 | ✓ | L1 ⭐⭐⭐⭐⭐（工程語義） | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py` |
| `T` 可由 non-capped 完整 candidate set 與 `n_trees` 驗證 | ✓ | L1 ⭐⭐⭐⭐⭐（工程語義） | historical `layered_reconstruction_*.json` |
| 拓撲 shape 可用 rooted unordered AHU form canonicalize | ✓ | L1 ⭐⭐⭐⭐⭐（同構計數） | `layered_tree_reconstruction.py:_canon_shape` |
| 7 datasets 每個 ambiguous primary unit 都有 family-specific read-AF | ✗ | L2 ⭐⭐⭐⭐ | `/tmp/20260711_read_af_tree_ordering.json`；各樣本有 3–603 units missing |
| clean raw-all → layered-v3 7/7 完成 | □ | U | active producer 尚未 aggregate success |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | Boolean hypercube、rooted group-Steiner、AHU 同構均定義清楚 |
| 觀察支撐 | 10 | historical 7/7 可算；clean-v3 未完成 |
| 機制清晰度 | 20 | C→T→shape→VAF ordering 的資料依賴可逐欄追蹤 |
| 反例風險 | 10 | pooled-vs-HP、CN/VAF confound、partial-only C=0 均需顯式處理 |
| 所需資源 | 10 | 全候選 VAF/shape census 約 1–6h 含 QA |
| **TOTAL** | **70 / 100** | **GO for descriptive report** |

**Falsifier observable**：若定義或實作錯誤，會出現 C bins 不回到 primary-unit/region 母數、`T` 重枚舉與 `n_trees` 不一致、canonical shape 數與 `n_distinct_shapes_exact` 不一致，或 VAF ranking candidate mismatch >0。

**Reality-test**：

1. 雙 HP 的 pooled C 若小於 HP-specific C 總和，證明 pooled 會合併家族，不可當主 C。
2. `C=0` 仍可有 partial-supported tree，證明 0 不是「沒有樹」。
3. CN-altered winner 若占多數，VAF winner 不得升格為獨立確認。

## §3 Assumption map

| Assumption | Importance | Known | 處置 |
|---|---|---|---|
| C 不含 virtual ROOT、hidden Steiner node 或 all-reference full group | HIGH | known by user definition + solver contract | 主報告固定；另列 all-full group 作 audit |
| 雙 HP 的 C 應保留 `(C_HP1,C_HP2)` 並可報和，不相乘 | HIGH | known | protected decision |
| T_joint 為各 primary HP candidate-tree 數乘積 | HIGH | known | 取代舊誤名 C_region |
| family read-AF 可選真實樹 | HIGH | unknown | 只允許 relative ordering；不寫 confirmation |
| historical snapshot 可代表 clean-v3 | HIGH | false | partial/scientific-NO-GO ribbon |

## §4 HCC1395 probe gate

1. 讀 HCC1395 mlhp/layered JSON → 驗證 C/T/shape 欄位可連接。
2. 重算 HCC1395 C、T、Topo、hidden、VAF → 驗證所有母數守恆且 candidate mismatch=0。
3. 通過後擴全 7 datasets；任一樣本守恆失敗則停止 HTML 發佈。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| clean layered-v3 7/7 | HIGH：限制正式 scientific claim | upstream long run | P0 external |
| 1,005 ambiguous units 缺 family read-AF（7 datasets 加總） | MED：VAF ranking coverage不全 | 需查零 coverage/coordinate | P1 |
| purity/CN/multiplicity corrected VAF likelihood | HIGH：不能叫獨立 confirmation/CCF | >6h | P1 |
| single-cell/multi-region truth | HIGH：不能確認完整癌症樹 | external data | P2 |

## §6 Conflict scan / red team

| Prior conclusion | Relation | 來源 |
|---|---|---|
| first-32 display prefix ranking 不可信 | support：必須全候選重枚舉 | `InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/20260710_分層重建端到端流程程式碼判斷稽核_01.md` |
| read-AF ordering 是 heuristic，不是 CCF/independent validation | conflict with「VAF確認真樹」措辭 | 同上 |
| current L3 canonical status = not_evaluated | dependent：本報告的新重算要另列 derivative provenance | `InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md` |
| repo root `MEMORY.md` 不存在 | gap | `rg --files -uu` 無結果 |

**最強反方**：C 在 pooled 與 HP grain 混用會改變分布；VAF 祖先序在 CN gain/LOH 或 subclonal multiplicity 下可能反轉。兩者均以雙口徑表與 claim ceiling 處理，故描述性報告可 GO。

## §7 Decision path

- **GO**：HCC1395 probe → 7-dataset deterministic census → partial HTML report。
- **PROBE**：任何「VAF confirmed true tree」「confirmed clone/subclone」主張。
- **Decision lock**：本報告鎖定 descriptive engineering scope；clean-v3 或 corrected VAF likelihood 完成時才重評。

## Provenance

- Commit：`6067568637088838a9f518955e41d222f057e4f1`
- Build time：2026-07-11 16:48 +0800
- Topic：`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/`
