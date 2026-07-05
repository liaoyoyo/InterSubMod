<!--
建立時間: 2026-07-01
類型: 研究主軸 + 流程 + 檔案地圖 handoff（給其他 session AI 快速 onboarding）
狀態: authoritative snapshot（HCC1395 ⭐3；反映到 2026-07-01 全部修正）
build_branch: feat/summary-nreadsvalid
data_sources: docs/methodology/20260627_clone_subclone_integrated_report/data/*.json
-->

# 論文研究 Handoff — 給其他 session AI 的快速 onboarding

> **一句話讀懂**：碩論《**Subclonal reconstruction using somatic haplotagging and methylation profiles**》（HCC1395 單樣本、ONT 長讀、**⭐3 exploratory**）。**用 somatic sSNV 單分子共現重建局部克隆結構（真引擎）；haplotag 鑑別 allelic vs clonal；甲基是有界的表觀 characterize（非偵測器）**。新 session 讀本檔即可掌握主軸/流程/檔案；細節進各層文件。

---

## §0 主軸與敘述架構（A-framing，最重要，勿違反）

**論文分工（摘要首段就要點明）**：
| 資訊源 | 角色 | 一句話 |
|---|---|---|
| **somatic sSNV 單分子共現** | 🟢 **重建引擎**（reconstruction 歸此）| 唯一非循環的克隆共現證據 |
| **haplotag（HP tag）** | 🟢 **鑑別器**（allelic vs clonal，非確認器）| 移除等位假 subclone |
| **甲基（methylation）** | 🟡 **有界 characterize**（非 driver、非 detector）| 對已確認結構的表觀註解 |

🔴 **紅線（對外/論文絕不可越）**：regional partition（≤read-span）非 genome-wide tree；甲基 **corroborate 非 detect**；分子證據非 single-cell confirmation；勿稱「甲基偵測 subclone / genome-wide tree / 對手缺檢定」。

---

## §1 研究流程（三段，帶 verified 數字 + 檔案路徑）

> 所有數字可在 `InterSubMod/docs/methodology/20260627_clone_subclone_integrated_report/data/*.json` grep（§13-A 注入）。報告夾＝ `InterSubMod/docs/methodology/20260627_clone_subclone_integrated_report/`（`00_INDEX.md` 為導覽 SoT）。
> ⚠ **用詞範圍**：本段「確認結構 / 局部克隆樹」指**讀段內 sSNV 連鎖可解析出巢狀/分支**（遺傳連鎖層事實），**非生物學 subclone 確認**——單樣本分子證據的 confirm 天花板見 §5（需 single-cell/multi-region）。

### 第一段：somatic 共現骨幹 → 重建 subclone 結構（區域 partition + nested/sibling）
- **35,332 sSNV**（30,490 TP + 4,842 FP）→ **7,143 區域**：**4,678 有確認結構**〔含 **677 full_tree** + linear_nested + sibling + 858 co_linked 單 lineage〕、**2,465 sparse**。其中 **3,820（53%）具分支/階層結構**（`20260626_genomewide_sSNV_linkage_region_trees` 口徑，排除 858 單 lineage）。同一 ONT read 上多 sSNV 的 2×2 共現 → 互斥/巢狀/共連 → 局部克隆樹。
- 檔：`01_locus_master.md`（per-locus 主表 35,332）、`06_integrated_narrative.md`；資料 `data/sm_locus_master.tsv`、`_assets/20260618_subcluster_pilot/sm_linkage_genomewide.json`。
- 腳本：`sm_linkage_genomewide.py` → `sm_evolution_build.py` → `sm_region_integration.py`。

### 第二段：haplotag 鑑別 allelic vs clonal
- **HP = 鑑別器非確認器**。不分 HP 時誤判 **9,187** subclone → HP gate 後 same-HP **保留 3,949（候選，非確認）**、**移除 5,238（57%）等位（diff-HP）**。🔴 same-HP 高多為「區域背景」非克隆證據（故 3,949 是「過等位 gate 的候選」非已確認 subclone）；HP 的真判別力在**互斥**（mutual_excl DEPLETED 0.86×）。
- 檔：`02_hp_contribution.md`、`00b_methods_grounding_HP_ISM.md`；資料 `data/sm_hp_contribution.json`。
- 🔑 **HP tag 語意（L1 核對 `src/core/ReadParser.cpp:121-154`；genome-wide V1 交叉表 2026-07-04 佐證）**：`1/2`=germline（乾淨、非循環；V1: 0% 碰 somatic ALT）；`1-1/2-1`=germline HP 上帶 somatic 的 haplotype（V1: ~97% confident ALT）；`3`=HP3 = **經過任一 somatic ALT，但 germline 定相未確認/HP 衝突**（V1: ~97.5% confident ALT）；`0`=**unphase = REF-only 定相失敗**（沒碰 germline 或 germline 衝突；V1: 99.2% 不碰 ALT；ReadParser 輸出字串 `"0"` 非 `"unphase"`）。判別軸 = **是否經過 confident somatic ALT**（germline 0% vs somatic 家族 ~97%）；HP3≠unphase = 兩個不同 read 族群的定相失敗，勿混談。LongPhase-S **完整 9 態**（含 `4` 與 `1-2/2-2` multi-subclone）。🔴 `1-1/2-1/3` 當獨立 subclone 證據=循環。**V2 分樹決策**：HP3 併入帶該突變的 lineage 樹；unphase 不進 somatic 樹（背景）。

### 第三段：甲基輔助 characterize（有界，非偵測）
- 從 tumor BAM MM/ML 重抽 genome-wide **740 區**（驗證 corr=1.000 vs ISM）：僅 **49（6.6%）corroborate**（power-gated：高功率區 54.9%，僅 51 區達）；**0 個新 partition 由甲基獨立發現**。
- **cis-control（修正後）**：14/49 可評估，**11 解析為 germline 等位特異甲基化（非 subclone）**、3 候選 → where testable，甲基多為 allele-specific 非 subclone。
- 檔：`05_methylation_corroboration.md`、`08_methylation_sufficiency_audit.md`；資料 `data/sm_methyl_reextract_ALL.json`。

---

## §2 甲基能/不能幫什麼（L9 + L10，最新大規模驗證）

**兩層鐵律**（`10_two_layer_validation.md`，754 區三軸驗證）：
| 層 | 訊號 | 大規模數據 | 有用？ |
|---|---|---|---|
| **Haplotype（germline ASM）** | 強、普遍、非循環 | neutral **62%**、Δβ **0.97** | ✅ 對 **phasing/read 招募**有用 |
| **Subclone（somatic-specific）** | 稀少 | 全基因組僅 **7–9 區**（chr21 neutral 最乾淨）| 🟡 需 matched-normal |

**read 招募力**（`sm_read_recruitment.py`）：甲基錨點把 read 招募到 **HP 準確度 96.3%**、可救 **334** 個 unphased/HP3 read（擴充 phasing）；但 **within-HP subclone 招募僅 10.5%** → **甲基串接是 phasing 層擴充非 subclone 層**。

**somatic HP 子相位**（`10_two_layer_validation.md §9`，`sm_somatic_hp_methyl.py`）：測 `1` vs `1-1`（同 germline HP 內 somatic 軸，非 germline 循環）→ **680 區可評估、15% 甲基對齊 somatic 相位、Δβ 0.96**。🔴 但 94% LOH、與突變共變 = **somatic-cis + subclone 混合**，非已證 subclone，需 normal 分離。

**單 sSNV / 無法建樹位點**：單位點 ASM 32% 但**最不可解釋**（無錨+C>T 破 CpG 假象）；read×read 甲基分群 **1/754 data-starved**。→ 補無法建樹位點靠**遺傳（更多 sSNV/更深 read）非甲基**。

---

## §3 ISM 是什麼（方法 + 目標 + 創新點口徑）

- **方法**：per region 產 read×CpG 甲基矩陣 → read-read 距離矩陣 → clustering → **PERMANOVA 結構存在性檢定**；甲基 corroborate 用 **genotype-anchored**（sSNV 定群，非無監督 clustering，杜 double-dip）；cis-control 用 HP 軸 + normal-baseline cis-test。
- **目標**：不是「用甲基找 subclone」；是**遺傳重建 + 甲基/HP 佐證/刻畫**。
- **論文創新點口徑（投稿必守）**：「**無監督 read×read 距離矩陣結構 PERMANOVA + normal-baseline cis-test + somatic-subclone 目標**」。🔴 **禁**用「對手二代定序 / 對手缺顯著性檢定」當差異（cvlr/ASMS/MethylBERT 都 ONT-capable 且有 randomization 檢定）。

---

## §4 檔案地圖（其他 session 從這裡找）

| 類別 | 路徑 |
|---|---|
| **整合報告（主）** | `InterSubMod/docs/methodology/20260627_clone_subclone_integrated_report/`（`00_INDEX.md` 導覽 + `00b`基礎 + `01`-`10` 各層 + `index.standalone.html` 互動版）|
| **canonical 數據** | 同上 `data/*.json`/`*.tsv`（每數字 grep-able）|
| **骨幹腳本（26 個 sm_*.py）** | `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/scripts/sm_*.py` |
| **重生 manifest** | `.../_assets/20260618_subcluster_pilot/README_sm_linkage_pipeline.md` |
| **論文 framing** | `InterSubMod/docs/paper_focus/00_PAPER_MASTER_subclonal_reconstruction_haplotagging_methylation.md` |
| **HP/工具權威** | `src/core/ReadParser.cpp:121-154`（HP tag 語意）；KB `05_tools/longphase-s.md`、`longphase-to.md`（MCP knowledge_get_doc）|
| **外部文獻庫** | `/big7_disk/liaoyoyo2001/external_validation/`（74 源，repo 外；入口 `REGISTRY.md`）。⚠ 橋接索引 `docs/method_comparison/20260613_external_validation_library_index_01.md` 目前只在 main repo branch，本 worktree 分支 `feat/summary-nreadsvalid` 未同步（跨 session 查此索引請切 main 或用 external_validation/REGISTRY.md）|
| **postmortem（本 session 修的 bug）** | `InterSubMod/docs/postmortems/20260629_hp_tag_str_int_comparison_bug_postmortem.md` |

**relevant memory**（`.claude/.../memory/`）：`project_clone_subclone_integrated_report_finalized`、`project_subclone_snv_linkage_verification_pipeline`、`project_thesis_writing_architecture`、`project_subclonal_reconstruction_paper_focus`、`feedback_baseline_dependence_not_result`、`feedback_feature_name_vs_definition_rule`。

---

## §5 關鍵 caveat + 教訓 + 下一步

- 🔴 **str/int bug 教訓（2026-06-29）**：pysam `get_tag("HP")` 對 `HP:Z:` 回 **str**；早期腳本用整數 `in (1,11)` 比對恆 False → `hp_control_eval=0` 假象曾誤寫「cis-control 不可評估」。**任何 tag 比對必 `str()` 正規化；任何「N=0/從不發生」先疑 parse bug 非生物學**（見 postmortem + `feedback_feature_name_vs_definition_rule`）。
- 🔴 **baseline-dependence 陷阱**：甲基隨 somatic 變異共分離、扣 cis 前 = cis-ASM 非獨立佐證；「招募到 HP」≠「招募到 subclone」。
- **⭐3 天花板**：單樣本分子證據只能 characterize 不能 confirm subclone（需 single-cell/multi-region）。
- **唯一推進路徑（B 軸，deferred）**：**matched-normal cis-control** 把 germline cis-ASM 殘差化 → 存活的 tumor-specific 甲基才可能是真 subclone（把 7–9 候選 + 15% somatic-HP 對齊從「somatic-cis 疑似」升級）。normal 端 5/6 樣本 ready（COLO829 缺）。

---

## §6 三十秒總結（給趕時間的 session）

> **遺傳（sSNV 共現）建骨幹 → HP 分 allelic/clonal → 甲基有界 characterize**。甲基主訊號是 **germline 等位特異甲基化（haplotype 層，62% neutral、Δβ 0.97）→ 對 phasing/招募有用（HP 招募 96%、救 334 read）**；**subclone-specific 甲基全基因組僅 7–9 區 + somatic-HP 對齊 15%（多 somatic-cis，需 normal）→ 甲基不是 subclone 偵測器**。單樣本 ⭐3 = characterize 非 confirm。
