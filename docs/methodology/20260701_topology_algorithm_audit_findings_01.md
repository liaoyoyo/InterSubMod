---
title: 克隆樹建樹演算法方法學稽核（Q1-Q4）+ 待辦 ②④ 記錄
date: 2026-07-01
status: audit_findings
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/{topology_analysis.py,sm_multilocus_combinations.py,backbone_resolution.py}
method: dynamic-workflow(4 research + 4 adversarial-verify agents)+ 主回合資料量化;所有數字重算自 topology_per_region.json
---

# 克隆樹建樹演算法方法學稽核

> ⚠ **2026-07-01 同日更新（並行 session fe9d8003 的 C2/C3/D4 修正,commit 979c54d/44ed19c/e8be1bb）**：本稽核跑在 pre-fix 資料;修正後 determinacy = **A_determined 1741（原 1812）/ incompatible 118（原 12）/ A_ambiguous 62（原 76）/ B_pairwise 943（原 958）/ C_underdetermined 544（原 550）**;situation_dist = 衝突 118 / 已確定 2178 / 跨HP 91 / pairwise 892 / 欠定 544 / 順序 62。🔴 **本稽核的待辦 ② 四配子 hoisting 已由 C2 獨立落地**（把非 CN-gain 乾淨四配子違反從靜默丟改判 incompatible → 這是 12→118 暴增的主因）;C3 修「丟 60.9% reads 仍標 A_determined」（本稽核 Q1/Q5 flag 的 truncated-A 失真）。以下正文的具體計數是稽核當時 pre-fix 快照,結論（verdict/盲點）不變,現行計數以本 banner 為準。單一真值:`InterSubMod/docs/methodology/20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md`。

## TL;DR
對建樹演算法四個核心設計做了「查真實文獻定理 + 從資料量化 + 對抗驗證」稽核。**結論：核心在群/向量層穩健，在樹拓樸推論有系統性盲點。** 最該顯著標註的一點：所有 perfect-phylogeny/IDPP 建構的有效性都依賴 infinite-sites，癌症 LOH/CNV（本專案 82-91%）系統違反。已落地「透明度」與「正名」到工作站；動到上游建樹演算法的兩項（②④）記錄於此待日後排程。

## 摘要（3 行）
Q1 缺值 complete-case = 正確保守（別 impute）但太激進 + 未記錄損失；Q2 laminar-on-rows 與 4-gamete-on-columns 不可比、且畫不出深分支（最大盲點）；Q3 父指派最大真子集正確、99% ambiguity 是 CN artifact；Q4 A 門檻不是瓶頸、read-span 才是（>50kb 區 A 率僅 11%）。

---

## Q1 — 部分覆蓋 read 靜默丟棄（complete-case，不 impute）
**判定：sound_with_caveats。** 對「單分子共現」估計目標，丟棄未跨全部位點的 read 是正確保守選擇——無法觀察分子沒覆蓋位點的 linkage，impute 會捏造 phased linkage（chr17 bug 已證：未覆蓋→REF impute 造假中間群 RAA + 假線性樹）。
**理論**：① 根已知（germline 全-REF）→ 缺值補全屬 IDPP（Pe'er-Pupko-Shamir-Sharan 2004）多項式可解；無根版 NP-complete（Gramm et al. 2009）。② complete-case 在「完整案例機率 ⊥ outcome given covariates」時 unbiased（White & Carlin 2010）；此處 missingness 由 read span 驅動非 allele 驅動 → 對共現結構 unbiased。③ haplotype assembly（MEC/WhatsHap，Majidian 2020）標準是用 `-` 標未覆蓋、read 仍貢獻覆蓋子位點——不整條丟。
**盲點/量化**（HCC1395，源 topology_per_region.json 7143 區）：40.4%(2887/7143) n_sSNV≥2 區向量全空被丟；42 密集區截斷丟 1132 sSNV；co-read 隨距離衰減（88@<1kb→11@20-50kb）=span-informative missingness；逐區丟棄比例未記錄。🔴 甲基插補救援已測失敗（nearest-centroid LOOCV balanced acc 0.566 < 多數 baseline，僅 29.5% 勝）。
**verify 更正**：row nested-or-disjoint 對有根情況是 perfect phylogeny 的**精確**充要（非近似）；但 IDPP 有效性依賴 infinite-sites，LOH 破壞 → 不可 oversell 成便宜填補。

## Q2 — laminar（row ALT-set）vs 4-gamete（column）
**判定：suboptimal。** 不相容 screen 用錯物件：`_laminar`（topology_analysis.py:58-61）查每 read ALT 集合（行），Gusfield 四配子要求 character 對（列）。🔴 verify：兩者**不可比**（既漏檢真四配子違反 30/32，又可能錯誤丟棄合法深分支 read）。
**最大盲點（對主軸最關鍵）**：row-laminar 建樹**結構上畫不出 subclone-of-subclone**——1112/1113 個「branched」只在 germline 根分支。理論：非空節點兩 child 都含 parent ALT 集→交集→laminarity 逼 nested→不能 sibling。→「branched」應叫「多根 lineage」（已落地正名）。
**量化**（源 topology_per_region.json 3885 detail 區）：determinacy = A_determined 1812 / B_pairwise 958 / C_underdetermined 550 / other 477 / A_ambiguous 76 / incompatible 12。12 incompatible 全 has_cycle（9 CN-gain+3 loh），與重算的 32 個真四配子違反**不重疊**（不同機制、兩者都需要）。四配子 test 已在 backbone_resolution.py:22-30 但被 gate 只事後跑。

## Q3 — 父指派（最大真子集 + 字串 tiebreak）
**判定：sound_with_caveats。** 在 laminar 前提下最大真子集 = 可證明正確的最大簡約 parent，非 hack、非瓶頸：只 1.98%(77/3885) 受影響，**~99% 是 CN-gain multiplicity artifact**（CN-neutral 只 1 區）。唯一明確缺陷：把 jump>1（parent 仍唯一、缺中間 Steiner 節點）與真多父 tie（38 區）混進同一 ambig flag。
**更好做法**：報告層分開 jump>1 vs tie（下游 backbone 已 prototype）；tie 用 CN-corrected CCF + pigeonhole（Tarabichi 2021；pos_vaf 全有）但 tie 本身是非樹結構、CCF 是 band-aid；LOH 刪突變（24/77）用 Dollo parsimony。🔴 絕不用甲基排序（循環）。多樹 posterior（Pairtree）單 bulk 不可行 → 維持等機率候選集。

## Q4 — A_determined 門檻 + read-span（主回合實測補完）
**≥6 read 門檻不是瓶頸**：有 ≥2 ALT 向量的 1890 區，0 個卡在 <6 read。
**真正限制是 read-span**：51%(1995/3885) 有向量區只有單一 ALT 向量（結構上達不到 A 需要的 ≥2）；A 率隨 span 崩（<2kb 47% / 2-10kb 70% / 10-50kb 41% / >50kb **11%**）；sSNV 三桶 linked 61.0% / underpowered 15.4% / isolated 23.5%（38.9% 因距離>read-span 連不起來）。→ 確認「可建單分子樹的區受 read 長度限制」，這是 ⭐3「regional partition ≤read-span 非 genome-wide tree」的根本原因，屬**物理極限非門檻設太嚴**。

---

## 已落地（工作站，commit 見 git log 2026-07-01）
- **① 透明度面板**：整體觀察區加「這張工作站沒顯示什麼/已知損失」摺疊（40.4% 空向量、read-span 衰減、深分支上限、四配子事後、infinite-sites caveat）。
- **③ 正名**：topology_type「branched(直系+姊妹)」顯示改「多根 lineage」+ STAT_DICT 標明 row-laminar 畫不出 subclone-of-subclone（1112/1113 只在根分支）。
- **LOH 標籤/篩選**：逐區 LOH tag（list + detail KV）+ f_loh 改「僅 LOH 區」用 longphase-TO LOH.bed overlap（補 6 樣本 cn=unknown 的洞）。

## 待辦（②④ 動上游建樹演算法，較大改動，記錄待排程）
| 項 | 內容 | 成本 | scope 限制 |
|---|---|---|---|
| **② 四配子 screen** ✅**已由 C2 落地(07-01)** | C2(commit 44ed19c)把 `independent`(4-gamete)從靜默丟 → 計數 + 非 CN-gain 乾淨違反判 incompatible(12→118) | 已完成 | ✅ 非 CN-gain 乾淨違反已改判;CN-gain multiplicity 假違反仍待 CN-aware genotyping(未做) |
| **② CN-aware genotyping** | 對 CN-gain loci 建模 allele multiplicity（SPRUCE/PhyloWGS/Canopy 標準） | 中等（需 allele-specific CN，部分可從 SM_CNBED/SAVANA） | 可合法清 22 CN-gain 假違反 |
| **④ 記錄丟棄** | 每區輸出 n_reads_spanning_any / n_full_cov_reads / dropped_partial_fraction | 極低 | 純透明；把靜默丟棄變可稽核 |
| **④ partial fragment 拓樸** | 讓 multi-locus 拓樸吃 partial fragment（已在 pairwise 層做）救回部分空向量區 | 中等 | 先用 pairs_eps2 ∩ 空區算真正可救比例（<40.4%，未算），不可把 40.4% 當可救目標 |
| Dollo parsimony | normal-polarized 救 LOH 違反（Q2 9 + Q3 24 LOH），當 annotation 非靜默改樹 | 中等 | 少數區 |

## 相關研究問題（2026-07-01 附）
- **SAVANA 能判一區有幾群 CNV 嗎？** ❌ 不能。SAVANA（Nature Methods 2025, s41592-025-02708-0）做 purity/ploidy（BAF@LOH）+ allele-specific CN，但**假設單一 clonal 解**（best purity-ploidy fit 全基因組套用），**不做 subclonal 建模**；bulk 有亞克隆時它**平均** CN 而非解出各 clone。要「幾群 CNV」需 subclonal CN caller（Battenberg/ASCAT-subclonal）輸出 fractional CN，但單 bulk 仍受限（與本專案「定不出來」一致）。
- **甲基能幫分群/算差異嗎？** 分群 ❌（cis-ASM 循環：甲基群對齊你已用 genotype 定的軸、非獨立；LOOCV 甲基連 genetic label 都救不回，balanced acc 0.566<baseline）。算差異 ⚠️：能算 Δβ 但屬 cis-ASM（對齊 genotype/HP 軸），非獨立判別器 = bounded-auxiliary（CROSS-HP ~35% 可控、SAME-HP 結構性無解）。06-28 cis-control 定案。Reopen：C1 多樣本 / C3 single-cell methylation ground truth。

## 來源
- 動態 workflow（4 research + 4 adversarial verify）transcript：session subagents/workflows/wf_39723c87-867。
- 真實定理：Pe'er et al. SIAM J.Comput. 2004(IDPP)；Gusfield 1991(4-gamete)；White & Carlin Stat.Med. 2010(CCA vs MI)；Majidian et al. PLOS ONE 2020(MEC)；Tarabichi 2021(CCF ordering)；SAVANA Nature Methods 2025。
- 數字全重算自 topology_per_region.json（3885 detail / 7143 n_sSNV≥2 / 42 truncated）。
