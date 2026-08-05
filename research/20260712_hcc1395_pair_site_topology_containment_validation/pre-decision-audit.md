<!--
建立時間: 2026-07-12T04:54:04+08:00
目標: 在逐位點與 HP forest 層驗證 HCC1395 兩技術來源的 topology equality / induced-substructure / conflict
處理範圍: chr1-22；5,720 exact-coordinate complete-both regions；historical layered-v2 engineering snapshot
cycle_id: 20260712-0454-hcc1395-site-topology-containment
topic: hcc1395_pair_site_topology_containment_validation
status: verdict_GO_engineering_probe_scientific_no_go
audit_version: 0.1
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_units.tsv
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_candidates.tsv.gz
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/
-->

# Pre-Decision Audit：HCC1395 逐位點 topology containment

> **Verdict：GO（70/100）僅限 fail-closed engineering probe；scientific proof 仍 NO-GO。** 可全量比較同一 cell line 的兩個 cross-source technical datasets，但不能把相容候選、VAF 排序或癌症註記升格成 biological tree truth。

## §0 Task / Goal / Cynefin

- Task type：**B Comprehensive validation**；全 chr1-22、全 5,720 exact complete-both pairs，不使用 subset 代替 headline。
- 服務目標：G4（跨資料來源一致性）為主，G5（外部可驗證證據鏈）與 G3（read-level constraints）為輔。
- Domain：**Complex**。候選樹多重性、HP assignment、site inventory 與 VAF/CN/coverage 共同影響可觀察結果，需 safe-to-fail probe 與 matched null。
- Thread D：間接相關（read-level pattern），但不是 methylation endpoint。
- Thread B：不在撤回方向內。
- KDE-corrected：不適用。
- VCF caller AF：本輪不用 merged AF；VAF 層使用 same-HP per-site read-AF，明示未校正 CN/purity/multiplicity。
- 長計算／C++／搬移：Python 全量重枚舉；不改 C++，不搬移既有檔案。

## §1 問題、假說與 falsifier

**問題**：相同 exact genomic region 中，兩側 read-pattern-compatible candidate trees 在共同 sSNV 的祖先／平行關係是否一致？若位點集合不同，較小樹是否等於較大樹刪除 private-site events 後的 induced topology？

**假說 H1**：控制 shared-site k、HP components 與候選集合規模後，同一 HCC1395 cell line 的兩技術來源會有高於條件隨機基準的逐位點 structure agreement，且 conflict rate 有限。

**Falsifier**：若 observed compatibility 不高於 matched null、只在 `k=2` 或大型 candidate sets 成立、VAF 與 read 層方向翻轉，或差異主要是 HP split/merge，則不能支持 cross-source structural reproducibility。

## §2 已知觀察與資料品質

| Observation | 結果 | Tier / 限制 |
|---|---:|---|
| Exact complete-both 母群 | 5,720 regions | L2 historical layered-v2 |
| Region site set 完全相同 | 5,534（96.75%） | caller loci；不是 HP-specific universe |
| Strict subset site sets | 185（3.23%） | 子結構規則的直接適用母群 |
| Partial-overlap site sets | 1（0.02%） | 只能稱 shared-core |
| Shared site allele identity | 15,713/15,713 PASS | frozen VCF CHROM/POS/REF/ALT |
| HP unit count 不同 | 703（12.29%） | 不可由 HP swap 消除 |
| Shared `k=2` | 3,121（54.56%） | 只能辨方向 vs parallel |
| Full candidate materialization | 97,663 unit trees | 全 unit non-capped、solver 可重枚舉 |
| VAF comparison-count varies | HCC 71.93%；DORADO 87.27% ranked units | exact-sum 候選間尺度風險 |
| Biological topology truth | 缺 | 無 single-cell / colony / synthetic truth |

## §3 Operational definitions（決策鎖）

1. 位點 ID 使用 frozen VCF 的 `(chrom,pos,REF,ALT)` 驗證；樹解碼依 `mlhp.positions` 對應 bit index。
2. 每個 `region×HP` 的 mutation-bearing loci 取 frozen solver `universe`，不是整個 region positions。
3. Non-recurrent candidate 在 shared loci 上以每對位點 `{i→j, j→i, parallel}` 關係矩陣作 canonical induced-topology signature；private-site edge 收縮。
4. Shared-site recurrence fail-closed，另列 `recurrent_projection`；不壓成單一 event。
5. `H_` 只表示 hidden/partial-supported state。Primary 忽略顯示前綴但保留 event path；strict sensitivity 才比較 hidden/observed support mask。
6. Ordered HP 是 primary；整個 region 只允許一次 HP1↔HP2 swap sensitivity。HP count mismatch 與含 shared mutation 的 unmatched component 另列，不能默認 compatible。
7. Read endpoint 是 **full read-pattern-compatible candidate set**，不是 read truth。候選投影集合必報 exact、proper subset（方向）、partial overlap、disjoint、Jaccard、雙向 coverage；`any-compatible` 只作 optimistic sensitivity。
8. VAF endpoint 使用 exact Fraction score argmax、保留 ties、fixed T1 直接加入、缺 VAF fail-closed；另報 comparison-count-constant subset。稱 heuristic top set，不稱 posterior 或確認樹。

## §4 Assumption map / Red team

| Assumption / risk | Importance | Known | Action |
|---|---|---|---|
| Shared locus 在兩 VCF 是同 allele | HIGH | KNOWN | 15,713/15,713 QA；輸出 provenance |
| 任一 candidate match 可代表重現 | HIGH | FALSE | 禁作 headline；候選數越多越易偶合 |
| HP1/HP2 label 可直接跨來源對齊 | HIGH | UNKNOWN | ordered + whole-region swap；分列 assignment conflict |
| `k=2` 能代表複雜 branching | HIGH | FALSE | k=2/3/4/5+ 分層 |
| VAF exact-sum 在所有 candidate 間同尺度 | HIGH | OFTEN FALSE | constant-count subset + normalized sensitivity |
| 兩來源是 biological replicates | HIGH | FALSE/UNKNOWN | 只稱 same-cell-line cross-source technical data |
| 癌症基因／藥物 annotation 是 topology truth | HIGH | FALSE | 只作 context，不作 validation endpoint |

**Reality tests**：

1. Direct-heavy／small-k 演算法 prior 即使錯配 region 也會有高 compatibility。
2. `|Q_A||Q_B|` 很大時，existential match 接近必然。
3. 一側 1 HP、另一側 2 HP 若忽略 unmatched component，會把 split/merge 誤算成子樹。
4. VAF 與候選樹均來自同批 reads，不能形成獨立證據鏈。

## §5 Step → Verify

1. 重枚舉兩側全部 primary HP candidates → 驗證：unit `n_trees`、non-capped、candidate ID 100% 對上 frozen outputs。
2. 建 HP-specific genomic universe 與 shared-site projection → 驗證：golden cases涵蓋 private event contraction、direct/parallel conflict、hidden relaxed equality、recurrence fail-closed。
3. 分別比較 read full set 與 VAF top set → 驗證：5,720 row conservation；分類 exact+subset+overlap+disjoint+not-evaluable=5,720。
4. 跑 complexity-matched conditional null → 驗證：同 chromosome、shared k、HP count、candidate-size bin；observed 與 null 用同一 HP-swap 規則。
5. 做 chromosome-block bootstrap／leave-one-chromosome-out → 驗證：95% CI 與 22 個 autosome 方向，不把 region 當完全獨立。
6. 產出 HTML／TSV／JSON → 驗證：兩次重跑 deterministic、所有整數可回加、HTML desktop/mobile browser QA。

## §6 Gate 與 claim ceiling

工程 GO checkpoint（非領域標準）：evaluable coverage ≥70%；observed−conditional-null ≥10 pp；block-bootstrap 95% CI lower >0；conflict 與 complexity strata 全揭露。未通過仍交付負面結果，不調門檻。

若通過，最高只可寫：

> 同一 HCC1395 cell line 的兩個 cross-source technical datasets，在共同 sSNV 上重現高於條件隨機基準的部分 regional mutation-state constraints／induced substructures，支持部分技術再現性。

不得寫：真 clone tree、biological accuracy、獨立 biological replicate、VAF 已確認祖先方向、癌症基因／藥物註記證明方法有效。

## §7 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | unit-flip event、shared-site induced topology 可正式定義 |
| 觀察支撐 | 20 | 5,720 全量配對、shared allele QA、完整候選可重枚舉 |
| 機制清晰度 | 15 | read constraints、HP forest、VAF heuristic 可拆層 |
| 反例風險 | 5 | 已知 small-k / candidate multiplicity / HP mismatch；null 尚未跑 |
| 所需資源 | 10 | 約 10 萬 candidate trees，Python 可行 |
| **TOTAL** | **70/100** | **GO engineering probe；scientific proof NO-GO** |

## Provenance Footer

- Commit hash：`6067568637088838a9f518955e41d222f057e4f1`
- Build time：2026-07-12T04:54:04+08:00
- Skill：`/pre-decision-audit` v0.1

