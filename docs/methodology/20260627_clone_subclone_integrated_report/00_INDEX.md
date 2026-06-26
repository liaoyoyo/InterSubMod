<!--
建立時間: 2026-06-27
類型: 索引 + 定義 + 格式規範 — HCC1395 clone/subclone 整合分析報告
狀態: in_progress（HCC1395 ⭐3, Tier-R）
build_branch: feat/summary-nreadsvalid
-->

# 00 — 索引 / 定義 / 格式（HCC1395 clone/subclone 整合分析）

> **報告目的**：用 **sSNV + ONT read** 建構整體 clone/subclone，逐層（HP→CCF→PS→甲基）**單變量驗證**各資訊源的 data-supported 貢獻，統整成完整可驗證的分析。HCC1395 單樣本 ⭐3，Tier-R。
> **可信度先後鐵則**（見 `00b`）：**sSNV 連鎖 > HP > 甲基**；衝突 genetic 勝。

## §1 檔案索引

| 檔 | 層 | 內容 | 圖 | 數據 |
|---|---|---|---|---|
| `00_INDEX.md` | — | 本檔：索引/定義/格式 | — | — |
| `00b_methods_grounding_HP_ISM.md` | 基礎 | HP/ISM 定義 + precedence + 錯誤模式 | — | (源碼) |
| `01_locus_master.md` | L0 | per-locus 主表（35,332 sSNV 全狀態 join）| `01_*.png` | `sm_locus_master.{tsv,json}` |
| `02_hp_contribution.md` | L1 | HP ablation（移除 57% allelic）| `02_*.png` | `sm_hp_contribution.json` |
| `03_ccf_tiers.md` | L2 | 祖先≥後代 VAF 96% + CCF 峰 | `03_*.png` | `sm_ccf_tiers.json` |
| `04_phaseset.md` | L3 | PS-reliability 旗標 + Tier-PS 不延伸 | `04_*.png` | `sm_phaseset_extension.json`, `region_ps_flag.json` |
| `05_methylation_corroboration.md` | L4 | 甲基既有輸出覆蓋不足 + 小樣本 | `05_*.png` | `sm_methyl_corroboration.json` |
| `06_integrated_narrative.md` | L5 | 整合敘述（對抗稽核後）| — | — |
| `07_cross_sample_capability.md` | L6 | 跨樣本能力（只 HCC1395 ready）| — | — |
| `index.standalone.html` | L7 | **整合 HTML**（嵌圖 + 全層 §13-A 注入）| — | — |

## §2 定義詞典（meaning）

| 詞 | 定義 |
|---|---|
| **sSNV 連鎖** | 同一 ONT long read 上多 somatic SNV 的物理共現（非循環、ground truth）|
| **region（區域）** | 最大可關聯區域 = ≤50kb 內相連的 somatic sSNV chain（Tier-R）|
| **tree_role** | ancestor / descendant / intermediate / sibling / region_member(co_linked併或區內未連) / isolated |
| **HP / HP1-1** | longphase-S 單倍型 tag；HP1-1 = somatic phased 到 HP1 germline 背景（HP1∪HP1-1=同單倍型）|
| **same-HP / diff-HP** | 兩 sSNV ALT reads 同/異 germline 單倍型 → subclone / allelic 鑑別 |
| **CCF tier** | CN-corrected VAF 估的 subclone 細胞分率層（clonal~0.9 / major~0.45 / rare~0.05；僅 CN-clean）|
| **PS-reliability** | region 跨幾個 phase-set；單一=HP 可信，多 PS=phase-switch 風險不可信 |
| **甲基 corroborate** | 用 sSNV 定好的群（genotype-anchored）看甲基是否區分；**非 detect**（單獨循環）|
| **cis-ASM / double-dip** | genotype-同質群內甲基子結構 = 循環假象（無 a-priori 錨）|

## §3 格式規範（format）
- **每層 .md**：frontmatter `data_sources:` 指向 `data/` JSON（§13 provenance）；含「方法 / 數據觀察(verified) / data-supported 結論 / 限制」。
- **數字**：一律可在 `data/*.json` 或 `*.tsv` grep 到（§13 反捏造）；HTML/報告由 JSON 注入（§13-A）。
- **單變量紀律**：每層只加一變量 + ablation/對照 + 落檔，未驗證不堆疊。
- **證據層級**：claim 標 L1（源碼/數據重現）/ L2（僅 JSON）/ L3（推論）。
- **紅線**：⭐3 單樣本；regional 非 genome-wide tree；甲基 corroborate 非 detect；對外勿稱「甲基偵測 subclone / genome-wide tree / 對手缺檢定」。

## §4 如何驗證（之後研究可重用）
1. 任一數字 → 查對應 `data/*.json`（frontmatter data_sources 指明）→ grep。
2. 逐位點查 `data/sm_locus_master.tsv`（35,332 列）；逐區域旗標 `data/region_ps_flag.json`。
3. 重現指令見骨幹 manifest `_assets/20260618_subcluster_pilot/README_sm_linkage_pipeline.md` + 本報告各層腳本 `_assets/.../scripts/sm_*.py`。
4. 對抗稽核紀錄見 `06_integrated_narrative.md`（L5 fresh-context evaluator verdict）。

## §5 層狀態
L0 ✅ · L1 ✅(HP=鑑別器非確認器) · L2 ✅(GMM BIC n=3) · L3 ✅(PS-reliable 92.7%) · L4 ✅(**已從 BAM 重抽 genome-wide 740 區，6.6% 弱 corroborate**) · L5 ✅(**對抗稽核 5 輪 → PASS**，修正全套用) · L6 ✅ · L7 ✅(整合 HTML)。
全部 data-supported；**每個數字皆可從 `data/*.json`/`*.tsv` grep（§13-A，5 輪 fresh-context 稽核三方驗證 prose=JSON=HTML 一致）**。詳見 `06_integrated_narrative.md` §4 五輪對抗稽核紀錄。
