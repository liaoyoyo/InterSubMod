---
title: 拓撲結構三軸觀察（群數 c × 樹拓撲 × 樣本）— 統計分析 + 對抗驗證
date: 2026-07-02
status: analysis
build_branch: research/subclonal-reconstruction-202606
evidence_grade: L3 / ⭐3（descriptive characterization, single-bulk, partly-artifact）
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json, /big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{sample}/topology_per_region.json
method: 7 樣本 C2/C3 修正後同版本;main-loop cross-tab + scipy 卡方;3-agent workflow(獨立重算 + confound skeptic + 收斂)獨立驗證
---

# 拓撲三軸觀察 — c × 拓撲 × 樣本

## 🔴 一句話裁決
數字**逐位可重現**（獨立 scipy），跨樣本差異**真實存在但只是部分生物訊號（`partly_artifact`）**；證據級 **L3/⭐3 = 描述性 characterization，非 subclone 判別**。兩個「branched-heavy」離群樣本（COLO829、HCC1954）**成因不同、不可合併敘述**。

## TL;DR
- **c=1 → single 100%（機械）**;c=2 是唯一有 branched-vs-linear 自由度的格。
- **c=2 branched% 隨樣本顯著變**：COLO829 87.0 / HCC1954 78.6 vs 其他 5 樣本緊聚 53–63%。卡方 chi2=525.9、CramérV=**0.227**(small–medium)。
- 🔴 **但**：branched% 是**read-未驗證的幾何上界**;CN 是決定性未控 confound;p 值因 pseudoreplication 無意義（真 n=7）;COLO829 的高 branched 主要是**低 coread artifact**。

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
- 資料:`scratchpad/canonical_shapes.json`（per-sample shape counts + descriptors + consistency,皆重算 grep-able）。
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

## 來源
- 資料:7 樣本 topology_per_region.json（C2/C3 修正後,今日重生）;三軸中間 `scratchpad/three_axis.json`（cross-tab + 卡方 + confirmed% + CN 分層,皆重算 grep-able）。
- 對抗驗證:workflow `wf_ebf15c81-927`（3 agents）。
- 關聯:`InterSubMod/docs/methodology/20260701_topology_algorithm_audit_findings_01.md`（C2/C3/D4 + branched=多根 lineage）、memory `project_topology_workstation_features_state`。
