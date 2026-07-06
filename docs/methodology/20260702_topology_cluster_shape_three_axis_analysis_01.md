---
title: 拓撲結構三軸觀察（群數 c × 樹拓撲 × 樣本）— 統計分析 + 對抗驗證
date: 2026-07-02
updated: 2026-07-06（重算到 07-05/06 資料版）
status: analysis
build_branch: research/subclonal-reconstruction-202606
evidence_grade: L3 / ⭐3（descriptive characterization, single-bulk, partly-artifact）
data_version: 2026-07-05/06 重生（gap#1 subcube 救回 + incompatible 重分類）;總 38655 區
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json, /big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{sample}/topology_per_region.json, docs/methodology/_assets/20260702_topology_three_axis/three_axis_v2_6288.json, docs/methodology/_assets/20260702_topology_three_axis/canonical_v2_6288.json
method: 7 樣本同版本;main-loop cross-tab + scipy 卡方;3-agent workflow(獨立重算 + confound skeptic + 收斂)R1-R4 DRY 驗證(方法學)
---

# 拓撲三軸觀察 — c × 拓撲 × 樣本

> ## 🔄 2026-07-06 資料刷新（重算到 07-05/06 版）
> 本報告原以 **07-02 資料(總 22676 區)**做 + R1-R4 對抗驗證。並行 session **07-04~06** 把資料重生（**gap#1 subcube 救回** + **incompatible 重分類** 118→18）：總區 **22676→38655**、新增 **~16k 個 c=0「no_genotype_vectors」區**（gap#1 救回·無 ALT 向量·非拓撲分析對象）、新 determinacy 類別（`E_subcube_recovered(gap#1)`、`recurrence_*`）、incompatible 全降。
> 🟢 **核心裁決 ROBUST（重算確認）**：**c≥1 區的 c=1/2/≥3 幾乎不變**（HCC1395 c2=1744 兩版相同）;**c=2 branched-vs-linear 跨樣本 CramérV 0.227→0.222**（近乎不變）;COLO829 87.1%/conf 39.3%、HCC1954 79.0%/conf 52.6% 皆與舊版一致;**CN gain 67.6→68.2 vs LOH 27.5→29.8**（robust）。→ **partly-artifact/L3 主裁決不變**。
> 🔴 **改變的**：總區數、incompatible%（因分母增 +重分類：合計 5.9%→**0.9%**、逐樣本 0.4–19.1%→**0.0–2.8%**）、canonical 形狀種類（11→**32**，H2009 複雜區增）。以下表格已更新到新版;方法學/caveat 不變。

## 🔴 一句話裁決
數字**逐位可重現**（獨立 scipy），跨樣本差異**真實存在但只是部分生物訊號（`partly_artifact`）**；證據級 **L3/⭐3 = 描述性 characterization，非 subclone 判別**。兩個「branched-heavy」離群樣本（COLO829、HCC1954）**成因不同、不可合併敘述**。

## TL;DR（數字＝07-05/06 資料版）
- **c=1 → single 100%（機械）**;c=2 是唯一有**乾淨二元對比**（branched-vs-linear）的格,故用 c=2 控 c↔topology 機械相依。**新增 ~40% 區為 c=0（無 ALT 向量·gap#1 救回）＝非拓撲對象**。
- **c=2 branched% 隨樣本顯著變**：COLO829 **87.1** / HCC1954 **79.0** vs 其他 5 樣本緊聚 55–63%。CramérV=**0.222**(small–medium;07-02 版 0.227,近乎不變)。
- 🔴 **但**：branched% 是**read-未驗證的幾何上界**;CN 是決定性未控 confound;p 值因 pseudoreplication 無意義（真 n=7）;COLO829 的高 branched 主要是**低 coread artifact**。

## 📊 數據總表（07-05/06 版；群數 c 分布 + 有效性 + read 驗證；依 c=2 branched% 排序）
| 樣本 | 總區 | c=0(無向量)% | c=1% | c=2% | c≥3% | 形狀種類 | incompatible% | **c=2 branched%**(幾何上界) | c=2 confirmed%(read) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| COLO829 | 6613 | 43 | 39 | 19 | 0 | 5 | 0.0 | **87.1** | 🔴 39.3 |
| HCC1954 | 2837 | 30 | 29 | 39 | 2 | 10 | 0.0 | **79.0** | 52.6 |
| H1437 | 8169 | 42 | 31 | 25 | 2 | 18 | 0.8 | **63.0** | 58.7 |
| HCC1937 | 1886 | 13 | 32 | 52 | 3 | 9 | 0.0 | **58.0** | 67.7 |
| HCC1395 | 6288 | 38 | 32 | 28 | 2 | 13 | 0.3 | **57.7** | 57.6 |
| HCC1395_DORADO | 3418 | 30 | 34 | 35 | 1 | 9 | 0.1 | **57.0** | 60.1 |
| H2009 | 9444 | 55 | 19 | 22 | 3 | 26 | 🔴 2.8 | **55.7** | 66.7 |
| **合計/均** | **38655** | — | — | — | — | **32 種** | **0.9** | — | — |

**欄位**：**c=0(無向量)**=no_genotype_vectors（gap#1 subcube 救回·無 ALT 向量·非拓撲對象,新版才有）;c=1/2/≥3=有 ALT 向量的群數分布;形狀種類=此樣本不同 canonical 樹形數(全 7 樣本合計 **32 種**);**incompatible%**=determinacy=四配子/perfect-phylogeny 違反率(誠實有效性訊號;07-04 重分類後 0.0–2.8%);**c=2 branched%**=兩 ALT 向量非嵌套(幾何上界·不需 read);confirmed%=c=2 有 read 驗證結構的比例(對照上界)。

**三組比較（c=2 branched% 主軸，穩健於資料版）**：
- **偏平行(branched% 高·conf 低)**：COLO829(87.1→conf 39.3🔴 缺口最大=低 coread artifact)、HCC1954(79.0→conf 52.6·undetermined)。
- **均衡**：H1437/HCC1937/HCC1395/DORADO/H2009 branched 55–63% ≈ confirmed。
- **結構複雜度**：COLO829 最簡(僅 5 種形狀·c≥3≈0);**H2009 最複雜**(26 種形狀·incompatible 2.8% 最高·c=0 佔 55%=gap#1 救回最多)。
- **跨樣本效應**：c=2 branched-vs-linear CramérV=**0.222**(small-med;07-02 版 0.227);**不報 p 值**(pseudoreplication·真 n=7)。
> 🔴 **枚舉完整**（32 種 canonical 形狀·0 未分類）;**有效性**：incompatible 合計 **0.9%**（逐樣本 0.0–2.8%,H2009 最差;07-04 重分類把多數舊 incompatible 轉 recurrence/隱藏祖先）。**c=2 branched% 跨樣本差異仍受 branched=幾何上界 / CN未控 / coread 三 confound**（§②③④）— 統計正確 ≠ 生物差異證據。

## 前置：7 樣本已同版本（C2/C3 修正後）
**資料版本(07-05/06)**：7 樣本經 gap#1 subcube 救回（新增 c=0 no_genotype_vectors 區）+ incompatible 重分類（把多數舊 incompatible 轉 recurrence/隱藏祖先樹）。c≥1 的核心區保留，故 c=2 分析穩健（見上方刷新 banner）。

## ① 三軸主結果（每樣本 c=2：branched% 與 read-confirmed% 並列；07-05/06 版）
| 樣本 | c=2 n | **c=2 branched%**(幾何上界) | **c=2 confirmed%**(tree_shape read-驗證) | 缺口 |
|---|--:|--:|--:|--:|
| COLO829 | 1225 | **87.1** | **39.3** | 🔴 47.8 |
| HCC1954 | 1095 | **79.0** | 52.6 | 26.4 |
| H1437 | 2053 | 63.0 | 58.7 | 4.3 |
| HCC1937 | 985 | 58.0 | 67.7 | −9.7 |
| HCC1395 | 1744 | 57.7 | 57.6 | 0.1 |
| DORADO | 1184 | 57.0 | 60.1 | −3.1 |
| H2009 | 2093 | 55.7 | 66.7 | −11.0 |

→ **COLO829 branched% 87% 但 read-confirmed 僅 39.3%（缺口 47.8,全樣本最大）** = 幾何斷言但 read 未驗證。（與 07-02 版幾乎相同 → c=2 分析穩健。）

## ①b Canonical 樹形 + 有效性（2026-07-02;🔴 R1 對抗審計修正 — 原把「枚舉完整」誤陳為「100% 有效」）
用**樹同構標準式（AHU canonical form，只看形狀不看標籤）**統計所有觀察樹。**兩件事必須分開講**（R1 審計 M1-M3）:

**(a) 枚舉完整性 ✅（正確、可放心）**
- 全 7 樣本合計 **32 種**不同 canonical 樹形（07-02 版 11 種;新版因 gap#1 救回更多複雜區—尤 H2009 26 種—而增）、**0 未分類**（每棵有 edge 的樹恰對應一種形狀）。仍以淺樹為主。
- 空間仍小的真正原因 = **經驗淺薄（多數 c≤2、深度 1-2）** + perfect-phylogeny/laminar 約束。⚠ **不宜寫「非 n^n」**：n^n 是 *labeled* 基準,32 是 *unlabeled*(AHU) 形狀,類別不同不可直接比;真驅動是淺薄。
- top 形狀（合計數·樹同構式）：`(())` 單群 11490 / `(()())` 2平行 6542 / `((()))` 2直系 3688 / `((())())` 2平行+1直系 513 / `(()()())` 3平行 112 / `((()()))` 92 / `(((())))` 3直系 31 / 尾端 c=4+ 各 1-19。

**(b) perfect-phylogeny 有效性 ⚠（非 100%，這才是誠實訊號）**
- 🔴 **「成環=0 / 100% 有效樹」是結構性恆真（tautology），不是驗證**：stored edges 是**去環後的 acyclic 近似**,AHU 的環標記結構上不可能出現 → 恆 0。同理 **「c>k+1 違反=0」是 off-by-one 的 vacuous check**（c 界應 ≤k）。
- **真正的有效性訊號 = determinacy `incompatible`（四配子/perfect-phylogeny 違反）**,07-05/06 版重算:
  - 合計 **355/38655 = 0.9%** incompatible（07-02 版 5.9%;下降因分母增 +gap#1 + **07-04 重分類把多數舊 incompatible 轉 recurrence/隱藏祖先樹**);`has_cycle` 631。
  - 🔴 **逐樣本 incompatible%**：COLO829 **0.0%** / HCC1954 0.0% / HCC1937 0.0% / DORADO 0.1% / HCC1395 0.3% / H1437 0.8% / **H2009 2.8%**（仍最差,cyclic 473）。
  - ⚠ **新 determinacy 類別（07-04~06）**：`E_subcube_recovered(gap#1;pairwise樹/欠定/含衝突)` ≈16k、`recurrence_candidate/artifact/LOH_unresolved/required`（Model A recurrence m-通道拆分）;舊「incompatible」大量移入這些桶。
- 資料:`_assets/20260702_topology_three_axis/canonical_v2_6288.json`（per-sample n/c分布/inc/cyc/nshapes + CN 分層,皆重算 grep-able）。
> 結論:**枚舉可放心當定論**（32 形狀·0 未分類）;**有效性**：incompatible 合計 **0.9%**（逐樣本 0.0–2.8%,H2009 最差）,**不可**寫「100% 有效樹」。跨樣本 c=2 差異另受 branched=幾何上界/CN/coread confound（§②③④）。

## ② topology_type vs tree_shape — 兩套不同計算，該用哪個
- **topology_type='branched'（topology_analysis.py:107 判定 has_branch + 111 命名,從 genotype 向量樹）= 確定性幾何陳述**：2 條被 call 的 ALT 向量彼此**非嵌套**→ 必各自掛 ROOT 分岔（非 a→b 線性）。含義 = 「≥2 條非嵌套 / 平行 / 多根相容譜系」,**不需 read 共現證據**、完整可重現。
- **tree_shape（sm_region_integration,從 pairwise 巢狀/姊妹邊）= 較保守**：需 read 共現確認 nested/sibling 邊,證據弱→`no_confirmed_structure`。
- **兩軸回答不同問題**：branched% = 「2 向量非嵌套」的幾何比例（不需 read）;confirmed% = 「有無被 read 確認結構」的比例（含線性/嵌套/姊妹確認）。⚠ **兩者非夾住單一真值的上/下界**（confirmed% 不是平行性下界,故 4/7 樣本 confirmed>branched → 缺口為負）;branched%−confirmed% 是**兩個非嵌套量的描述性差**。發散是 by construction 非 bug。
- 🔴 **報告/顯示用法**：以 branched% 為描述性主軸（完整可重現）**但務必永遠與 confirmed% 並列**;措辭「≥2 條非嵌套 ALT 向量（平行/多根相容）」,**禁**「confirmed parallel subclones / branched biology」。
- ⚠ **branched% 與 incompatible% 正交但部分重疊**（R2 揭露；07-04 重分類後 incompatible 大降故重疊也降）：branched 與 incompatible 各為獨立 single-value 分類、各區各計一次,分母**不互扣**;H2009 仍有部分 c=2-branched 落 cyclic → 那些不可詮釋為乾淨平行 lineage。headline 離群 COLO829/HCC1954 的重疊極小,不影響結論。

## ③ Confound 裁決（3-agent 對抗驗證）
| Confound | 裁決 | 關鍵證據 |
|---|---|---|
| 數字重現性 | 非 confound | 獨立 scipy 全小數點級吻合 |
| sSNV 密度 | **排除**（僅低密度帶） | c=2 mean n_sSNV：COLO829 2.35 / HCC1954 2.13 ≈ 低密度帶(2.09–2.41,**5/7 樣本**);⚠ H1437 4.25 / H2009 17.41 遠高但 branched% 普通 → 密度不追蹤 branched%,故排除 |
| Basecaller | **排除** | HCC1395 56.7% vs DORADO 56.5%（同樣本兩 basecaller,近乎相同）|
| **Coverage/coread** | 🔴 **COLO829 主導(多為 artifact)/ HCC1954 否決** | COLO829 within-sample:medCoread=16(最低)、51% c=2 群<6 reads、branched 中僅 4.7% 成對確認(sib-edge)、confirmed 僅 39.3%。⚠ 跨樣本 Pearson(medCoread, br-sib%)=+0.867 **僅 n=7 且 leverage-driven**（去 COLO829→0.76 p=0.079、再去 HCC1937→0.54 非顯著）→ 佐證但**裁決不依賴它**,靠 COLO829 within-sample 證據。HCC1954 medCoread=40、38.8% sib-confirmed 健康 |
| **CN/LOH** | 🔴 **material 且未控(決定性)** | HCC1395(唯一有 CN)c=2 branched%(07-05/06 版):**穩健對比 gain(n=1189) 68.2 vs LOH(n=483) 29.8 = 差 38 點**;neutral(n=65) 75.4 / loss(n=7) 42.9 為小 n 參考。**6/7 樣本 cn=unknown → 無法 CN-adjust**（注:6 樣本部分已補 SAVANA cna-only CN,但本表仍以 HCC1395 SEQC2 為準）|
| Pseudoreplication | 確認 | chi2 把 10230 區當獨立,實則同樣本共享 purity/coverage/CN → **p=2.2e-110 不可解讀**,真 n=7 |
| **Mapping-bias / MAX_SNV=8 截斷** ⭐R1補 | **未控,待查** | genotype_cap=8:H2009 17.7%(752/4243) 區真 sSNV 遠>8 被截（部分截為 c=1）、H1437 4.6%、HCC1395 1.1%。高突變密度區可能重疊低 mappability → 建議做 mappability sensitivity 再談生物意義 |
| c↔topology 機械相依 | 真實但已由 within-c=2 控制 | c=1 強制 single;c=2 內乾淨二元對比=branched-vs-linear,回流 coread confound（c≥3 佔 3.6% 有多向自由度） |
| TP/FP 標籤 | 弱 | 僅 HCC1395 SEQC2 真值;HCC1954 fp=0、COLO829 fp=255,6 樣本 census-derived |
| H2009 資料品質 | 旗標 | **incompatible 2.8%(仍全樣本最高) + cyclic 473(最多) + canonical 26 種形狀(最複雜) + c=0 佔 55%(gap#1 救回最多)** → noise/FP 疑慮,生物結論前單獨審查 |

## ④ 兩離群的分解（務必分開講）
- **COLO829 87%** = **主要低 coread artifact**（medCoread=16、51% 無法成 link、僅 4.7% 成對確認、confirmed 僅 39.3%）。
- **HCC1954 79.0%** = **非 coread artifact,但 UNDETERMINED**（coread 充足、confirmed 52.6%;可能真分岔、也可能 CN-gain multiplicity,因 cn=unknown 無法排除）。

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
3. 呈現 per-sample 表 + CramérV=0.222（07-05/06 版）,**不**呈現 chi2 p 值（或明標 pseudoreplication·真 n=7）。
4. CN caveat 醒目 + 附 HCC1395 CN 分層（穩健錨 gain 67.6% vs LOH 27.5% = 差 40 點;neutral n=65 小樣本參考）作 confound 存在證明。
5. 兩離群分開（COLO829=低 coread artifact / HCC1954=undetermined）。
6. H2009 資料品質旗標。
7. 頂部 ribbon：characterization / single-bulk / weak census labels / no CN control / L3-⭐3。
8. 已排除項（sSNV 密度、basecaller）也寫,示已對抗檢驗。

## 附錄:定義與計算方法（每指標精確定義 + 公式 + 源碼位置，可重現）

### A. 資料單位與上游欄位
| 名詞 | 定義 | 來源 |
|---|---|---|
| **region（區）** | 最大單分子連鎖區（≤50kb 的 somatic-sSNV 鏈）;分析只取 `n_sSNV≥2` 的區 | topology_analysis.py:141（`n_sSNV>=2` filter）+ sm_region_integration.py:19,200 |
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
| **n_clusters (c)** | `len([g for g in populations if "A" in g])` = 含≥1 個 ALT 的 genotype 向量數（germline 全-REF 不計）。🔴 **正確上界 c ≤ k**（觀察到 germline root 時;原寫 k+1 誤把未計的 root 也算進 → off-by-one） | topology_analysis.py:157-158 |
| **topology_type** | 建樹（父=最大真子集）後依樹形分：`single`=1 節點（108）/ `star(全姊妹)`=≥2 根且無分支（109;本資料 0 例）/ `linear(全直系)`=maxd==節點數且無分支·單鏈（110）/ **`branched(直系+姊妹)`=有任一節點≥2 子·has_branch（107 判定+111 命名）** / `mixed`（112;0 例）/ **`incompatible`（有環→topology_analysis.py:83;本資料 267×）** | topology_analysis.py:106-112 |
| **tree_shape** | 依 pairwise 邊分（較保守,需 read 確認）：`inconsistent`=has_cycle / `full_tree`=n_nested≥1 **且** n_sibling≥1 / `linear_nested`=n_nested≥1 且深≥2 / `single_nested`=n_nested≥1 但深<2（結構上不可達→0 例）/ `sibling_only`=n_sibling≥1 / `co_linked_lineage`=n_colinked≥1 / `no_confirmed_structure`=皆無 | sm_region_integration.py:147-158 |
| **determinacy**（07-05/06 版,全 7 樣本占比） | `E_subcube_recovered(gap#1;pairwise樹/欠定/含衝突)` ≈41%（07-04~06 gap#1 救回,最大桶）/ `A_determined` 26.2% / `B_pairwise` 11.6% / `other` 10.8% / `C_underdetermined` 6.3% / `recurrence_*`（candidate/artifact/LOH_unresolved/required·Model A m-通道）2.1% / **`incompatible` 0.9%** / `A_ambiguous_order` 0.7%。⚠ 07-02 版的 6 桶（A_determined 44.2%…incompatible 5.9%）已被 gap#1 + recurrence 重分類取代 | topology_analysis.py:169-176 |
| **canonical shape** | 樹同構標準式（AHU）：`rec(n)='('+ ''.join(sorted(rec(子)))+')'`,從 ROOT 遞迴,**只看形狀不看標籤**;stored edges 已去環故 `(X)` 恆不出現。**32 種** = 全 7 樣本合計不同值數（07-02 版 11;枚舉完整,非有效性保證） | 本分析 `canonical_v2_6288.json` |

### C. 統計指標定義 + 計算法（皆 c=2 子集）
| 指標 | 計算式 | 意義 |
|---|---|---|
| **branched%** | `count(topology_type=='branched' ∧ c=2) / count(c=2)` × 100 | 2 條 ALT 向量非嵌套的比例 = **幾何上界**（不需 read 共現） |
| **confirmed%** | `count(tree_shape∈{full_tree,linear_nested,sibling_only,co_linked_lineage} ∧ c=2) / count(c=2)` × 100 | 有 read 驗證結構的比例 = 「有無被 read 確認結構」的比例。⚠ **非平行性下界**（含線性/嵌套確認,故 4/7 樣本 confirmed>branched → 缺口為負）;branched%−confirmed% 是兩個非嵌套量的描述性差,非夾住單一真值 |
| 🔴 **incompatible%**（有效性訊號） | `count(determinacy=='incompatible') / n` × 100 | 四配子/perfect-phylogeny 違反率。**07-05/06 版 合計 0.9%,逐樣本 0.0–2.8%**（07-02 版 5.9%/0.4–19.1%,07-04 重分類轉 recurrence/隱藏祖先後降）。**取代原「成環=0」「c>k+1=0」** tautology |
| **枚舉完整性** | **32 種** canonical 形狀·0 未分類 per 全樣本（07-02 版 11） | 每棵樹都對應一種形狀;與「有效性」正交 |
| **CramérV** | `sqrt(chi2 / (N·(min(rows,cols)−1)))`,對「c=2 branched-vs-linear × 7 樣本」列聯表（scipy chi2_contingency）| 跨樣本效應量（**07-05/06 版 0.222**;07-02 版 0.227,近乎不變=穩健）;**p 值不用**（pseudoreplication,真 n=7） |

> **計算流程（可重現）**：① 7 樣本各跑 `sm_region_integration.py`(C2)→`topology_analysis.py`(C3) 產 topology_per_region.json（`SM_CNBED=""` 保 cn=unknown）→ ② main-loop python 讀 7 檔算上述指標落 `three_axis.json`/`canonical_shapes.json` → ③ Read 回驗證 → ④ 寫本報告 + HTML 面板（§13.0 分批,數字全 grep-able）。

## 來源
- 資料:7 樣本 topology_per_region.json（**07-05/06 重生版**,gap#1 subcube 救回 + incompatible 重分類）;三軸中間 `_assets/20260702_topology_three_axis/three_axis_v2_6288.json`（per-sample c分布 + branched/confirmed% + CramérV）+ `canonical_v2_6288.json`（n/c分布/inc/cyc/nshapes + CN 分層）,皆重算 grep-able。**07-02 版中間檔（three_axis.json/canonical_shapes.json,11 形狀/5.9%）保留為歷史。**
- 對抗驗證:workflow `wf_ebf15c81-927`（3 agents）。
- 關聯:`InterSubMod/docs/methodology/20260701_topology_algorithm_audit_findings_01.md`（C2/C3/D4 + branched=多根 lineage）、memory `project_topology_workstation_features_state`。
