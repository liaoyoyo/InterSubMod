---
title: 拓撲結構三軸觀察（群數 c × 樹拓撲 × 樣本）— 統計分析 + 對抗驗證
date: 2026-07-02
status: analysis
build_branch: research/subclonal-reconstruction-202606
evidence_grade: L3 / ⭐3（descriptive characterization, single-bulk, partly-artifact）
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json, /big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{sample}/topology_per_region.json, docs/methodology/_assets/20260702_topology_three_axis/three_axis.json, docs/methodology/_assets/20260702_topology_three_axis/canonical_shapes.json
method: 7 樣本 C2/C3 修正後同版本;main-loop cross-tab + scipy 卡方;3-agent workflow(獨立重算 + confound skeptic + 收斂)獨立驗證
---

# 拓撲三軸觀察 — c × 拓撲 × 樣本

## 🔴 一句話裁決
數字**逐位可重現**（獨立 scipy），跨樣本差異**真實存在但只是部分生物訊號（`partly_artifact`）**；證據級 **L3/⭐3 = 描述性 characterization，非 subclone 判別**。兩個「branched-heavy」離群樣本（COLO829、HCC1954）**成因不同、不可合併敘述**。

## TL;DR
- **c=1 → single 100%（機械）**;c=2 是唯一有 branched-vs-linear 自由度的格。
- **c=2 branched% 隨樣本顯著變**：COLO829 87.0 / HCC1954 78.6 vs 其他 5 樣本緊聚 53–63%。卡方 chi2=525.9、CramérV=**0.227**(small–medium)。
- 🔴 **但**：branched% 是**read-未驗證的幾何上界**;CN 是決定性未控 confound;p 值因 pseudoreplication 無意義（真 n=7）;COLO829 的高 branched 主要是**低 coread artifact**。

## 📊 數據總表（canonical 樹形組成 + 一致性 + read 驗證；依 2平行% 排序）
| 樣本 | 總區 | 單群% | **2平行%** | 2直系% | c≥3% | 形狀種類 | 一致性(成環/c>k+1) | c=2 read-confirmed% |
|---|--:|--:|--:|--:|--:|--:|:--:|--:|
| HCC1954 | 1979 | 41 | **44** | 12 | 3 | 8 | 0 / 0 ✅ | 52.6 |
| HCC1937 | 1636 | 36 | **35** | 25 | 3 | 8 | 0 / 0 ✅ | 67.7 |
| COLO829 | 3786 | 67 | **28** | 4 | 0 | 4 | 0 / 0 ✅ | 🔴 39.3 |
| HCC1395_DORADO | 2379 | 48 | **28** | 21 | 2 | 7 | 0 / 0 ✅ | 60.1 |
| H1437 | 4768 | 53 | **27** | 16 | 3 | 8 | 0 / 0 ✅ | 58.7 |
| H2009 | 4243 | 43 | **26** | 21 | 6 | 10 | 0 / 0 ✅ | ⚠66.7 |
| HCC1395 | 3885 | 51 | **25** | 19 | 3 | 9 | 0 / 0 ✅ | 57.6 |

**欄位**：單群=germline→1突變(single·`(())`);2平行=兩條平行 lineage 各1突變(`(()())`,幾何上界);2直系=一條累積2突變(`((()))`);c≥3=較複雜;形狀種類=此樣本不同 canonical 樹形數;一致性=成環(無效樹)/perfect-phylogeny c>k+1 違反(全 0=100% 有效);read-confirmed%=c=2 有 read 驗證結構的比例(對照 2平行 幾何上界)。

**三組比較**：
- **偏平行(2平行 ≫ 2直系)**：HCC1954(44 vs 12)、COLO829(28 vs 4)、H1437(27 vs 16) — 🔴 但 COLO829 read-confirmed 僅 39%(低 coread artifact)、HCC1954 undetermined。
- **均衡(2平行 ≈ 2直系)**：HCC1937(35 vs 25)、DORADO(28 vs 21)、H2009(26 vs 21)、HCC1395(25 vs 19)。
- **結構最簡**：COLO829(67% 單群、僅 4 種形狀、c≥3=0);**最複雜**：H2009(10 種形狀、c≥3 6%,⚠ 資料品質旗標)、HCC1937(單群僅 36%)。
- **跨樣本效應**：c=2 branched-vs-linear CramérV=**0.227**(small-med);**不報 p 值**(pseudoreplication·真 n=7)。
> 🔴 **一致性欄全 0 = 可放心當定論**（11 種形狀全涵蓋、100% 有效樹）;**但 2平行% 的跨樣本差異仍受 branched=幾何上界 / CN未控 / coread 三 confound**（§②③④）— 統計正確 ≠ 生物差異證據。

## 前置：7 樣本已同版本（C2/C3 修正後）
本分析前先把 6 樣本重生為 C2/C3 修正版（讀已存 06-29 sm_linkage,不讀 BAM,cn 保 unknown）。incompatible 全上升（signature 正確）：HCC1395 12→118 / COLO829 1→15 / H1437 65→263 / H2009 265→812 / HCC1937 0→44 / HCC1954 1→27 / DORADO 4→55。

## ① 三軸主結果（每樣本 c=2：branched% 與 read-confirmed% 並列）
| 樣本 | n | c=1(single)% | c=2 n | **c=2 branched%**(幾何上界) | **c=2 confirmed%**(tree_shape read-驗證) | 缺口 |
|---|--:|--:|--:|--:|--:|--:|
| COLO829 | 3786 | 67% | 1225 | **87.0** | **39.3** | 🔴 47.7 |
| HCC1954 | 1979 | 41% | 1095 | **78.6** | 52.6 | 26.0 |
| H1437 | 4768 | 53% | 2053 | 62.3 | 58.7 | 3.6 |
| HCC1937 | 1636 | 36% | 985 | 57.7 | 67.7 | −10 |
| HCC1395 | 3885 | 51% | 1744 | 56.7 | 57.6 | −0.9 |
| DORADO | 2379 | 48% | 1184 | 56.5 | 60.1 | −3.6 |
| H2009 | 4243 | 43% | 2093 | 53.0 | 66.7 | −13.7 |

→ **COLO829 branched% 87% 但 read-confirmed 僅 39.3%（缺口 47.7,全樣本最大）** = 幾何斷言但 read 未驗證。

## ①b Canonical 樹形 + 一致性驗證（2026-07-02 補：正確且完整的結構統計）
用**樹同構標準式（AHU canonical form，只看形狀不看標籤）**統計所有觀察到的樹,回答「是否能正確、一致地統計所有結構拓撲」：
- 🔴 **全 7 樣本合計只有 11 種不同 canonical 樹形**（非 n^n;perfect-phylogeny + c≤k+1 把空間壓到小有限集）。全部淺（深度多為 1-2,深度≥3 極稀少）。
- 🔴 **一致性 = 100%**：全 7 樣本**成環(無效樹) 0 / c>k+1 違反 0** → 所有建出的樹都是**有效 acyclic perfect-phylogeny**。
- 11 種形狀（合計數,樹同構式·標籤）：`(())` 單群 single 11490 / `(()())` 2平行 lineage 6542 / `((()))` 2直系 linear 3688 / `((())())` 2平行+1直系 513 / `(()()())` 3平行 star 112 / `(((())))` 3直系 chain 31 / 其餘 c=4/5/6 各 1-19。
- 每樣本 canonical 組成（top）：COLO829 67% 單群/28% 2平行;HCC1954 44% 2平行/41% 單群;HCC1937 36% 單群/35% 2平行/25% 2直系;H2009 43%/26%/21%。→ 結構組成隨樣本變（同 c=2 branched% 訊號,但 caveat 同下）。
- 資料:`_assets/20260702_topology_three_axis/canonical_shapes.json`（per-sample shape counts + descriptors + consistency,皆重算 grep-able）。
> 結論:**可以正確且完整地統計**（11 種形狀全涵蓋、0 未分類）、**且 100% 一致**（0 無效樹）。但跨樣本組成差異仍受下述 confound（branched=幾何上界、CN、coread）— 統計正確 ≠ 生物差異證據。

## ② topology_type vs tree_shape — 兩套不同計算，該用哪個
- **topology_type='branched'（topology_analysis.py:110-111,從 genotype 向量樹 `has_branch`）= 確定性幾何陳述**：2 條被 call 的 ALT 向量彼此**非嵌套**→ 必各自掛 ROOT 分岔（非 a→b 線性）。含義 = 「≥2 條非嵌套 / 平行 / 多根相容譜系」,**不需 read 共現證據**、完整可重現。
- **tree_shape（sm_region_integration,從 pairwise 巢狀/姊妹邊）= 較保守**：需 read 共現確認 nested/sibling 邊,證據弱→`no_confirmed_structure`。
- **兩軸回答不同問題**：branched% = 幾何**上界**;confirmed% = read 驗證**下界**。發散是 by construction 非 bug。
- 🔴 **報告/顯示用法**：以 branched% 為描述性主軸（完整可重現）**但務必永遠與 confirmed% 並列**;措辭「≥2 條非嵌套 ALT 向量（平行/多根相容）」,**禁**「confirmed parallel subclones / branched biology」。

## ③ Confound 裁決（3-agent 對抗驗證）
| Confound | 裁決 | 關鍵證據 |
|---|---|---|
| 數字重現性 | 非 confound | 獨立 scipy 全小數點級吻合 |
| sSNV 密度 | **排除** | c=2 mean n_sSNV：COLO829 2.35 / HCC1954 2.13 ≈ 非離群(2.09–2.41) |
| Basecaller | **排除** | HCC1395 56.7% vs DORADO 56.5%（同樣本兩 basecaller,近乎相同）|
| **Coverage/coread** | 🔴 **COLO829 主導(多為 artifact)/ HCC1954 否決** | COLO829 medCoread=16(最低)、51% c=2 群<6 reads、branched 中僅 4.7% 成對確認;Pearson(medCoread,confirmed)=+0.867。HCC1954 medCoread=40、confirmed 38.8% 健康 |
| **CN/LOH** | 🔴 **material 且未控(決定性)** | HCC1395(唯一有 CN)c=2 branched%：neutral **75.4** / gain 67.6 / loss 42.9 / **LOH 27.5**（擺 47.9 點）。**6/7 樣本 cn=unknown → 無法 CN-adjust** |
| Pseudoreplication | 確認 | chi2 把 10230 區當獨立,實則同樣本共享 purity/coverage/CN → **p=2.2e-110 不可解讀**,真 n=7 |
| c↔topology 機械相依 | 真實但已由 within-c=2 控制 | c=1 強制 single;c=2 內唯一自由度=branched-vs-linear,回流 coread confound |
| TP/FP 標籤 | 弱 | 僅 HCC1395 SEQC2 真值;HCC1954 fp=0、COLO829 fp=255,6 樣本 census-derived |
| H2009 資料品質 | 旗標 | branched 區 mean 26.6 sSNV vs linear 5.6 + incompatible 96 → noise/FP,生物結論前單獨審查 |

## ④ 兩離群的分解（務必分開講）
- **COLO829 87%** = **主要低 coread artifact**（medCoread=16、51% 無法成 link、僅 4.7% 成對確認、confirmed 僅 39.3%）。
- **HCC1954 78.6%** = **非 coread artifact,但 UNDETERMINED**（coread 充足、confirmed 52.6%;可能真分岔、也可能 CN-gain multiplicity,因 cn=unknown 無法排除）。

## ⑤ 需再驗證才能確認
- **取得 6 樣本 CN/purity/coverage**（尤其 COLO829/HCC1954）→ 唯一能 CN-adjust branched% 的路。
- **HCC1954 leg**：需 CN 或更強 read 驗證,區分「真平行」vs「CN-gain 多重性」。
- **TP/FP-conditioned 拓樸 claim**：需 SEQC2 級真值。
- **H2009**：資料品質存疑,單獨審查前不宜納生物結論。

## 結論（不 overclaim 底線）
這是「≥2 條非嵌套 ALT 向量在部分樣本較多」的**可重現描述性觀察**;其跨樣本差異**部分為 artifact（COLO829 coread）、部分受未控 CN 混淆**,**不能**作為平行 subclone 存在或跨樣本生物差異的證據。**適合放進 HTML 當 descriptive 數據觀察,但須帶下述 caveat。**

## HTML 顯示 caveat（Phase C 落地依據）
1. 措辭「c=2 非嵌套 ALT 向量比例（幾何上界）」,禁「confirmed parallel subclones」。
2. branched% **永遠與 confirmed% 並列**。
3. 呈現 per-sample 表 + CramérV=0.227,**不**呈現 chi2 p 值（或明標 pseudoreplication·真 n=7）。
4. CN caveat 醒目 + 附 HCC1395 CN 分層（27.5%↔75.4%）作 confound 存在證明。
5. 兩離群分開（COLO829=低 coread artifact / HCC1954=undetermined）。
6. H2009 資料品質旗標。
7. 頂部 ribbon：characterization / single-bulk / weak census labels / no CN control / L3-⭐3。
8. 已排除項（sSNV 密度、basecaller）也寫,示已對抗檢驗。

## 附錄:定義與計算方法（每指標精確定義 + 公式 + 源碼位置，可重現）

### A. 資料單位與上游欄位
| 名詞 | 定義 | 來源 |
|---|---|---|
| **region（區）** | 最大單分子連鎖區（≤50kb 的 somatic-sSNV 鏈）;分析只取 `n_sSNV≥2` 的區 | sm_region_integration.py:19,141 |
| **n_sSNV (k)** | 該區的 somatic sSNV 位點數 | 上游 census |
| **populations** | 該區「genotype 向量 → read 數」dict;向量 = 每個 sSNV 位點的 R/A 字串（如 `ARR`）;**最多 8 位點**（>8 截斷,truncated=true） | sm_multilocus_combinations.py（MAX_SNV=8）|
| **nested_edge** | 一對 sSNV 的 ALT-集合呈 ⊂（祖先→後代）| sm_region_integration.py:build_tree |
| **sibling_pair** | 一對 sSNV 共現但互不為子集（平行分支）| 同上 |
| **co_linked** | 一對 sSNV **總是共現**（完全連鎖→併為單一 lineage）| 同上:102 |
| **n_independent（四配子）** | 一對 sSNV 同時見到 RR/RA/AR/AA 全 4 種 = perfect-phylogeny 違反;`n_independent_clean`=min(RA,AR)>1（非單雜訊）| :108-124（C2）|
| **cn** | 該區 copy-number 狀態（gain/loss/loh/neutral/unknown）;`SM_CNBED=""`→全 unknown | sm_region_integration.py:36-42 |

### B. 每區指標定義 + 計算法
| 指標 | 定義（計算式） | 源碼 |
|---|---|---|
| **n_clusters (c)** | `len([g for g in populations if "A" in g])` = 含≥1 個 ALT 的 genotype 向量數（germline 全-REF 不計）。**理論上界 c ≤ k+1**（perfect-phylogeny） | topology_analysis.py:157-158 |
| **topology_type** | 建樹（父=最大真子集）後依樹形分：`single`=1 節點 / `star(全姊妹)`=≥2 根且無分支 / `linear(全直系)`=maxd==節點數且無分支（單鏈）/ **`branched(直系+姊妹)`=有任一節點≥2 子（has_branch）** / mixed | topology_analysis.py:106-112 |
| **tree_shape** | 依 pairwise 邊分（較保守,需 read 確認）：`inconsistent`=有環 / `full_tree`=n_nested≥1 **且** n_sibling≥1 / `linear_nested`=n_nested≥1（深≥2）/ `sibling_only`=n_sibling≥1 / `co_linked_lineage`=n_colinked≥1 / `no_confirmed_structure`=皆無 | sm_region_integration.py:147-158 |
| **determinacy** | `incompatible`=有環∨topology incompatible∨(乾淨四配子違反 且 cn≠gain)〔C2〕 / `A_noisy`=去噪丟>10%〔C3〕 / `A_determined`=≥2 ALT 向量 且 總 read≥6 / `B_pairwise`=tree_shape∈{full_tree,linear_nested,sibling_only} / `C_underdetermined`=no_confirmed_structure / `A_ambiguous_order`=缺中間群 | topology_analysis.py:169-176 |
| **canonical shape** | 樹同構標準式（AHU）：`rec(n)='('+ ''.join(sorted(rec(子)))+')'`,從 ROOT 遞迴,**只看形狀不看標籤**;有環回 `(X)`。11 種 = 全 7 樣本合計不同值數 | 本分析 `canonical_shapes.json` |

### C. 統計指標定義 + 計算法（皆 c=2 子集）
| 指標 | 計算式 | 意義 |
|---|---|---|
| **branched%** | `count(topology_type=='branched' ∧ c=2) / count(c=2)` × 100 | 2 條 ALT 向量非嵌套的比例 = **幾何上界**（不需 read 共現） |
| **confirmed%** | `count(tree_shape∈{full_tree,linear_nested,sibling_only,co_linked_lineage} ∧ c=2) / count(c=2)` × 100 | 有 read 驗證結構的比例 = **下界** |
| **一致性:成環** | `count(canonical shape 含 '(X)')` per 樣本 → 全 0 = 無無效樹 | 建出的樹是否 acyclic |
| **一致性:c>k+1 違反** | `count(n_clusters > n_sSNV+1)` per 樣本 → 全 0 = 符 perfect-phylogeny | 群數是否超過理論上界 |
| **CramérV** | `sqrt(chi2 / (N·(min(rows,cols)−1)))`,對「c=2 branched-vs-linear × 7 樣本」列聯表（scipy chi2_contingency）| 跨樣本效應量（0.227=small-med）;**p 值不用**（pseudoreplication,真 n=7） |

> **計算流程（可重現）**：① 7 樣本各跑 `sm_region_integration.py`(C2)→`topology_analysis.py`(C3) 產 topology_per_region.json（`SM_CNBED=""` 保 cn=unknown）→ ② main-loop python 讀 7 檔算上述指標落 `three_axis.json`/`canonical_shapes.json` → ③ Read 回驗證 → ④ 寫本報告 + HTML 面板（§13.0 分批,數字全 grep-able）。

## 來源
- 資料:7 樣本 topology_per_region.json（C2/C3 修正後,今日重生）;三軸中間 `_assets/20260702_topology_three_axis/three_axis.json`（cross-tab + 卡方 + confirmed% + CN 分層,皆重算 grep-able）;canonical 中間 `_assets/20260702_topology_three_axis/canonical_shapes.json`（11 形狀 + per-sample + 一致性）。
- 對抗驗證:workflow `wf_ebf15c81-927`（3 agents）。
- 關聯:`InterSubMod/docs/methodology/20260701_topology_algorithm_audit_findings_01.md`（C2/C3/D4 + branched=多根 lineage）、memory `project_topology_workstation_features_state`。
