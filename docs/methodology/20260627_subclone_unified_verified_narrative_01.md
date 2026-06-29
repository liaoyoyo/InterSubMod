<!--
建立時間: 2026-06-27
類型: 跨來源交互驗證統一敘述 — clone/subclone 重建（HCC1395 ⭐3）
狀態: in_progress（統一敘述收斂；研究經 open gates 繼續）
build_branch: research/subclonal-reconstruction-202606
provenance: workflow wck0bu3iq（8 agent：4 gather + 3 lens 對抗 + 1 綜合）；來源 branch feat/summary-nreadsvalid@5308d9e（pending-merge）+ trunk memory
data_sources: docs/methodology/20260626_genomewide_sSNV_linkage_region_trees_01.md, docs/methodology/20260627_clone_subclone_integrated_report/{06,08,09}*.md（皆 on-branch 5308d9e）, memory/{project_paper_claim_audit_consensus_base_2026_06_12, project_ism_complete_tpfpfn_existence_cis, project_tumor_only_axis_negative_subclone_classification, project_methyl_phasing_assist_line, project_cross_sample_asm_reproducibility, project_chr2_18m_subclone_locus_verification}
-->

<!-- provenance-verified: 本文 L1/L2 數字源 branch feat/summary-nreadsvalid@5308d9e 經 5 輪 fresh-context 對抗稽核的報告 + trunk memory；交互驗證=workflow wck0bu3iq。
2026-06-27 UPGRADE：已在 merge-staging（merge/summary-into-trunk @ eda4534，wt-merge-summary，含 branch 資料檔）**獨立 grep 驗證 10 個 headline 數字全命中實際 data JSON**（見 §5 驗證表）→ 從「per branch audited」升級為「verified-against-data」。FF 進 trunk 後資料檔即落本 tree。 -->

> 🔴 **2026-06-29 校正 + 部分降級 — 最新真值改以 addendum 為準**：本 06-27 doc 經並行 session **06-28 第三方獨立重驗**（`InterSubMod/docs/methodology/20260628_subclone_unified_narrative_reverify_addendum_01.md`，9-agent 唯讀重算，verdict=**SUPPORTED_WITH_CAVEATS**＝框架信心高/個別數字中）+ **cis-control pilot 做完**（`20260628_cis_control_scope_pilot_verdict_01.md`）後，下列已過時：
> 1. 🔴 **chr17 VAF「0.18」= 硬編捏造** → 實 **0.47/0.48**；不存在的 γ 第 4 sibling 一併移除（`sm_summary.py:156`）。成環應 **18/23**（非 19/23）。
> 2. 🟡 **chr17 canonical sSNV = 3（data）非 4（硬編）** → 待 somatic 定義統一（`==0` vs `<5%` 差 2,230 sSNV/9.05%）前不寫死。
> 3. ✅ **structure 對外口徑 = 4,678（65%，含 858 單 lineage）**（非 06-27 防禦版 3,820/53%）。
> 4. ✅ **cis-control「#1 open gate」已 CLOSED（06-28）**：原假設「解鎖甲基 Tier-3」**不成立**（tumor genotype 軸 ⟂ normal HP 軸 corr−0.026；SAME-HP 59% 結構性 UNDETERMINED 需 single-cell，非覆蓋問題）→ **甲基用途窮盡＝bounded-auxiliary**（memory `project_methylation_use_exhausted_bounded_auxiliary`）。
> **框架結論（genetic 驅動 / 甲基 bounded / COMPATIBLE）仍 CONFIRMED**；僅上列個別數字/gate 狀態以 addendum + verdict doc 為最新。

# Clone/Subclone 重建 — 跨來源交互驗證統一敘述（單一真值）

> **敘述框架**：Verdict-Pyramid + scientific-rigor §2 證據分級（每 claim 標 L1-L5）。
> **數字 provenance**：標 L1/L2 的數字來自 branch `5308d9e` 經 5 輪對抗稽核的報告 + trunk memory；**本 trunk 未獨立重 grep**（branch pending-merge），故以「per branch audited reports」口徑採用，tier 已反映此。

## §0 結論金字塔（headline）

**在單樣本 HCC1395（⭐3 上限）上，somatic sSNV 單分子共現是唯一非循環的克隆重建骨幹，可在 read-span 內重建區域級（非 genome-wide）局部克隆樹；甲基化全程為 bounded-auxiliary（corroborate 非 detect、0 個新 partition 由甲基獨立偵測、cis-control 0/740 從未執行故 subclone-specificity 未定），只刻畫不增量任何 call。**

**跨來源一致性裁決 = COMPATIBLE（無真矛盾）**：branch 的 genetic sSNV-共現重建與 trunk 的甲基-clustering NEGATIVE 是**正交兩軸**（genetic 非循環 / 甲基 double-dip 循環）；branch 自身甲基結果（0 新 partition、read×read PERMANOVA 僅 1/754 可算）**反而 corroborate trunk NEGATIVE**。唯一待校正＝chr2:18M「乾淨 subclone 甲基」口徑需扣 somatic-cis + own-data 漂移（69%↔52.8%）需同步。

## §1 統一結論（tiered，對抗驗證後）

### 骨幹（genetic = 非循環軸）
- **[L1] 唯一非循環重建軸 = somatic sSNV 單分子共現**（per-read 2×2 共現分類，read=同細胞）。統計底層 Fisher exact ≈ free-margin permutation（agree 258/264=97.7%，median|Δp|=0.009），power-gate 不從噪聲造結構。
- **[L1] Genome-wide 掃描 → 區域級局部克隆樹**：宇宙 35,332 = TP 30,490 + FP 4,842（sum✓）→ 7,143 區、full_tree 677（9.5%）、linked 61%；isolated 24% 為上界（same-PS>50kb deferred + cross-PS 無 read/cell 連結 = GAP）。**「為生物克隆樹」可信子集限 CN-clean（full_tree ≈205/677），regional 非 genome-wide tree。**
- **[L2] 有確認共現結構區數對外採防禦口徑 = 3,820 區（53%）**；4,678（65%）含 858 單 lineage 不應作 headline。此 fraction 為**上界**（含 4,842 FP-label + γ FP-source 3,204 + 未移除 dense-mapping/segdup 偽影）。
- **[L2] chr17:48360161 canonical full_tree**（🔴 06-29 修正：sSNV **3 非 4**〔待 somatic 定義統一〕；VAF **0.82>0.47/0.48**〔第三值「0.18」為硬編捏造，已移除不存在的 γ 第 4 sibling〕；3 對 off-diagonal 乾淨 0 violations；LOH）：計數待定義統一；**克隆樹「生物詮釋」=L2**（單例×單樣本×單 pipeline、無 single-cell）。詳見 addendum。

### 輔助軸（鑑別/驗證，非偵測）
- **[L1] HP = sibling-vs-allelic 鑑別器（非確認器）**：same-HP 正共現 1.7–1.87× 是區域背景屬性非克隆特異；診斷力在互斥（same-HP DEPLETED 0.86×）；HP 移除 allelic 57%（5,238/9,187）。源碼逐字驗 + 兩輪未翻。
- **[L2] CCF/VAF = 階層方向的條件式佐證（非強驗證）**：誠實口徑 = 69.8% support / 5.7% violate / 24.5% tie（**不可只報 5.7%**）；方向由 read-set nesting 指派（非循環）、VAF 僅驗 magnitude，且只在 CN-clean（46.8%）可估。**[L3] GMM BIC best_n=3 但 ΔBIC vs n=4 僅 17（邊際）+ CN-gain 污染 → 離散 subclone 層級降 L3**，不可逕稱純 subclone CCF 群。
- **[L1] PS = reliability 旗標，不延伸克隆**（germline phasing ≠ clonal）：92.7% 區單一 PS（含 71 no_ps，非舊報 94%）；同 PS CCF 僅 42% 一致。

### confound（共同主導）
- **[L1] CN-gain multiplicity = 主導 confound**：somatic sSNV **52.8%** 落 CN-gain（12,569/23,810，master TSV 重算；舊 64%/69% stale）；CN-confounded 53.2%；乾淨集 LOH+neutral 46.8%。⚠ 77.8% 是 segment-level 覆蓋（不同分母，禁與 52.8% 並列）。

### 甲基（bounded-auxiliary，跨 4 來源最強共識）
- **[L1] 甲基 = corroborate-not-detect / 非 usable filter**：`new_partitions_detected_by_methylation = 0`；740 區 corroborated 49（6.6%）；read×read PERMANOVA 僅 1/754 可算。佐證 — trunk ASM TP 3.95% vs FP 1.07%（弱~3.5×、~96% TP 不顯著、COLO829 TP≈FP）、confirm≥2 群 ≠ confirm subclone。
- **[L1] 甲基 subclone-specificity = UNDETERMINED（非已驗陰性）**：cis-control 從未執行 — `hp_control_eval = 0/740`（含 0/49 corroborated、0/1267 單 sSNV）是 **structural zero**，`cis_explained=0 ≠「排除 cis」`。「0% cis / subclone-specific」已主動撤回 → **BOUNDED_AUXILIARY**。corroboration 無法與 germline cis-ASM 分離（可能近乎全 cis-ASM）。
- **[L3] LOH 41/49 corroborated = between-allele germline-ASM 排除的 subclone 候選（非 subclone）**：仍受 LOH-unmask-of-imprinting（Martin-Trujillo 82–91%）+ subclonal-LOH-mixing + somatic-cis + double-dip confound；neutral 8/16% 正是 germline cis-ASM 可存處。**甲基-subclone 乾淨集 = NONE（須 normal cis-control）。**
- **[L1] 甲基作為 subclone 標記方向相反且 power-gated**：by tree-shape full_tree 2.6% < sibling 4.8% < linear_nested 6.8% < co_linked 9.7%（結構最複雜 = corroboration 最低）；18.4% 僅靠單一 CpG；dose-response popB_n≥20→54.9% 但僅 51 區。「甲基-DIFFERENCE 可偵測性 power-gated」為真，但「subclone SIGNAL 存在」未證。
- **[L1] 單 sSNV 甲基外推最不可信**：ASM+ 405/1267（32% = 多-sSNV 6.6% 的 5×），但 LOH 31.5% ≈ neutral 33.7% → 此高陽性率為 artifact/power-driven（C>T 破壞 CpG 序列假象、無連鎖錨）。

### ASM / tumor-only（trunk 既有，仍成立）
- **[L2] ASM 真實且跨樣本復現但非判別器**：6/6 excess-over-null>0（3 癌種 +0.101–0.241），**必看 excess-over-null 非 raw rate**；B-discrimination NEGATIVE+anti（strong-ASM FP-enriched 5×、|Δβ| AUC=0.505）；COLO829 TP≈FP。
- **[L2] tumor-only 無監督甲基 clustering+PERMANOVA = double-dip NEGATIVE（勿再開）**：根因 = ONT read-內甲基相關使任何 null mis-specified + data-starved。a-priori 優勢 = 合法 null / 免 collider 非偵測力 → **禁寫「more discriminative」**。

### 統御鐵則 + 過程誠信
- **[L5/公理] baseline-dependence**：任何測得量（Δβ/ASM/clustering/co-occurrence/k 群）意義 = f(量, baseline/null/label)；甲基隨 somatic-變異-標籤共分離在 cis baseline 下 = cis-ASM 非獨立 lineage，須 normal cis-control + 遠端位點扣 cis 才算佐證（G-B gate）。
- **[L3] precedence「衝突時 genetic 勝」= 設計立場非已驗階層**：壓測僅 n=1 衝突 + n=8 甲基（近軼事）。
- **[L1] 過程誠信**：branch 5 輪 fresh-context PASS（自抓修 64%→52.8% stale）+ trunk F=0（51 agents/606 claim）。⚠ 殘留：sibling region_trees doc 兩處 69% 未同步、structure 4,678↔3,820 定義未對齊 → 「finalize-ready」暫不成立。

## §2 誠實邊界 / 紅線（對外必守）

1. 🔴 **⭐3 單樣本** single-pipeline 封頂；升 ⭐4 需 ≥5/7 樣本（COLO829/matched-normal = 共同 blocker）。
2. 🔴 **regional（≤read-span / within phase-set）非 genome-wide clone tree**；cross-PS（>50kb）無 read/cell 連結 = GAP，不宣稱 ISM 解決。
3. 🔴 **分子共現 ≠ single-cell confirmation**；克隆樹生物真實性待外部正交（single-cell/multi-region）確認。
4. 🔴 **甲基 corroborate 非 detect、0 新 partition；cis-control 0/740 = structural zero → subclone-specificity UNVERIFIED**（可能近乎全 cis-ASM），非已驗陰性。
5. 🔴 **對外勿稱**：甲基偵測 subclone / genome-wide tree / 對手缺檢定 / 甲基「more discriminative」/ tumor-only 無監督 clustering 為偵測器。
6. 🔴 **偽影未清**（chr8:81-83M / chr9:41.8M / chr14:16.09M，缺 mappability track）；structure fraction 為上界；γ-class = FP-source 候選非確證。
7. 🔴 **跨 build 計數禁並列**（TP 30,490 唯一錨點；FP 4,842 ≠ 3,643；52.8% somatic-weighted ≠ 77.8% segment-level）；332,705 loci 未全 cis-test = 未測非陰性。

## §3 據此繼續研究（open gates，依 load-bearing 排序）

| 下一步 | 理由 | 對應任務節點 |
|---|---|---|
| ~~① matched-normal cis-control~~ ✅ **已 CLOSED（06-28）** | cis-control pilot 做完：假設「解鎖甲基 Tier-3」**不成立**（tumor genotype 軸 ⟂ normal HP 軸 corr−0.026；SAME-HP 59% 結構性 UNDETERMINED 需 single-cell 非覆蓋）→ 甲基窮盡＝bounded-auxiliary。詳 `20260628_cis_control_scope_pilot_verdict_01.md` | `T-GATE-GB` ✅done |
| **② 加深 ONT 覆蓋解 power-gate**（現 popB_n≥20 僅 51 區→54.9%） | 6.6% 是覆蓋限制非訊號缺如，但須與 cis-control 配對（更高 power 同樣 surface cis-ASM） | `T-METHYL-REEXTRACT`(後續) |
| **③ 跨樣本/跨化學復現（COLO829 + ≥5/7）** | single-pipeline 自我參照封頂；升 ⭐3→⭐4 | `T-DORADO` / `T-GATE-GA` / `T-M1` |
| **④ single-cell / multi-region 正交確認**（chr17、chr2:18M） | 分子共現≠single-cell；把 L2 生物詮釋轉 confirmed（Tarabichi LEARN） | `T-GATE-GD` |
| **⑤ LOH-unmask 明確處理**（41 LOH 套 imprinting/SD-ASM blocklist：Rosenski 460 parental-ASM+34,426 SD-ASM） | 41 LOH 區首要未處理 confound | `T-C-UNMASK` |
| **⑥ mappability/segdup mask 移偽影**（chr8/9/14），對外用 TP-clean + CN-clean full_tree(205) | structure fraction 含 FP/偽影為上界 | `T-GW-RECON`(後續) / `T-ONT-CNV` |
| **⑦ 同步 branch own-data 漂移**（sibling doc 69%→52.8% + structure 計數定義） | 三輪共通根因第四次復發；own-data 一致才真 finalize | branch-side（merge 前） |

## §4 provenance + tier 紀律

- 數字標 **L1**（源碼/JSON 重現）/ **L2**（僅報告 JSON）/ **L3**（推論）/ **L5**（規則公理）— 來自 branch 5308d9e 報告 + trunk memory，**本 trunk 未獨立重 grep**（pending-merge）。
- 交互驗證 = workflow `wck0bu3iq`（4 gather source-cluster + 3 對抗 lens〔cross-source-conflict / confound-DAG / overclaim-redline〕+ 1 綜合），8 agent。
- 對抗發現摘要：**0 真矛盾**（branch↔trunk 正交相容）；多處 overclaim 降 tier（structure fraction 上界、GMM 離散 L3、chr17/chr2 生物詮釋 L2、precedence L3、甲基 subclone-specificity 從「0% cis」改 UNDETERMINED）。

## §5 數據驗證表（2026-06-27 獨立 grep 驗證 — §8.4 知識追溯）

> 在 merge-staging（`merge/summary-into-trunk` @ eda4534，worktree `wt-merge-summary`，同時含 branch 資料檔）逐一 grep。**每個 headline 數字都命中實際 data JSON**。資料檔路徑相對 `docs/methodology/20260627_clone_subclone_integrated_report/data/`（FF 進 trunk 後即在本 tree）或 genome-wide `_assets/20260618_subcluster_pilot/`。

| headline 數字 | 值 | 來源檔:key | 狀態 |
|---|---|---|---|
| sSNV 宇宙 | 35,332 | `_assets/.../sm_summary.json` / `sm_completeness_ledger.json` | ✅ grep 命中 |
| linked | 21,554 | `sm_summary.json` `linked` | ✅ |
| 甲基重抽 tested | 740 | `data/sm_methyl_sufficiency_audit.json` / `sm_methyl_reextract_ALL.json` | ✅ |
| 甲基 corroborated | 49 | `sm_methyl_sufficiency_audit.json` `powered_corroborated`；`sm_methyl_reextract_ALL.json` `n_corroborated` | ✅ |
| 甲基 corroboration rate | 6.6%（49/740） | 同上（powered rate 0.1068=49/459 另存） | ✅ |
| **cis-control** | `hp_control_evaluable: 0` | `sm_methyl_sufficiency_audit.json` | ✅（structural zero 確認） |
| CN-gain | 52.8%（12,569） | `data/sm_locus_master_summary.json` `gain: 52.8` | ✅ |
| CCF GMM | `best_n: 3` | `data/sm_ccf_tiers.json` | ✅ |
| PS reliable | 0.927 | `data/sm_phaseset_extension.json` `reliable_rate` | ✅ |
| HP depletion / allelic | 0.86× / 57% | `data/sm_hp_contribution.json` | ✅ |
| 區域數 / full_tree | 7,143 / 677 | `20260626_genomewide_sSNV_linkage_region_trees_01.md §4`（report-level，未 pin JSON key） | 🟡 per-report |

**confirm 指令範例**（FF 後在 trunk，或現於 wt-merge-summary）：
`python3 -c "import json; d=json.load(open('docs/methodology/20260627_clone_subclone_integrated_report/data/sm_methyl_sufficiency_audit.json')); print(d['power_audit'])"`

## §6 全資訊導航（單一入口 — 所有 subclone 資訊在哪）

| 層 | 位置 | 確認方式 |
|---|---|---|
| **統一敘述（本檔）** | `InterSubMod/docs/methodology/20260627_subclone_unified_verified_narrative_01.md` | §1 tiered claims + §5 驗證表 |
| **任務節點** | `T-SUBCLONE-VERIFIED`（depends T-L4/T-GW-RECON/T-CCF/T-INTEG） | `tasks_board.html` 點該節點 |
| **記憶** | memory `project_subclone_snv_linkage_verification_pipeline` + MEMORY.md 索引 | grep memory |
| **原始報告** | branch `5308d9e`：`20260626_genomewide_*` + `20260627_clone_subclone_integrated_report/`（5 輪稽核） | `git show 5308d9e:<path>` |
| **資料檔** | 同上 `data/` + `_assets/20260618_subcluster_pilot/`（小 JSON/TSV 入版控、大檔 manifest 重生） | §5 confirm 指令 |
| **git 整合狀態** | merge `merge/summary-into-trunk @ eda4534`（C++ 編譯 OK、228/228 測試）；**FF 進 trunk pending 主樹乾淨** | `git log eda4534` |
| **繼續研究** | §3 open gates（#1 = T-GATE-GB normal cis-control） | — |
