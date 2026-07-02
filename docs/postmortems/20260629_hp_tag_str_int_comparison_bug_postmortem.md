<!--
建立時間: 2026-06-29
類型: postmortem — 分析腳本 str/int 型別 bug 導致假陰性結論
影響: L8/L9 cis-control 結論（HCC1395 clone/subclone 整合報告）
狀態: RESOLVED（4 腳本已修 + 全基因組重跑 + 文件全修正）
build_branch: feat/summary-nreadsvalid
-->

# Postmortem — HP tag str/int 比對 bug 造成「cis-control 不可評估」假象

## 摘要

4 個 Python 甲基分析腳本以**整數** `hp in (1, 11)` 比對 read 的 HP tag，但 pysam 對此 HCC1395 BAM 的 **`HP:Z:` 字串** tag 回傳 **str**（`'1'`/`'1-1'`），`'1' in (1, 11)` 恆 `False` → **HP 從未被讀到** → `hp_control_eval` 普遍為 0。這使 L8/L9 得出**錯誤結論**：「cis-control 0/740 不可評估」「HP 排解 cis-ASM 不可行」。實際 BAM 充滿 HP reads（chr17 抽樣 `HP:Z:1`×70,709 + `HP:Z:2`×112,070）。

## 根因

| 層 | 事實 |
|---|---|
| 工具 | 此 tumor BAM 由 **longphase-s** tag（字串 `HP:Z:1/2/1-1/2-1/3`），非 longphase-to（整數 11/21/33）|
| pysam | `get_tag("HP")` 依 tag 儲存型別回傳：`Z`→`str`、`i`→`int`。此 BAM 是 `Z`→回 **str** |
| 腳本 | `reads[rn][2] in (1, 11)`（int 集合）；`str("1") in (1,11)` = `False` → 全歸 "other" |
| 結果 | `hp1`/`hp2` 永遠空 → `hp_control_eval` 恆 0；非生物學「無 HP」，是型別不匹配 |

受影響腳本：`sm_methyl_reextract.py`、`sm_single_locus_methyl.py`、`sm_methyl_genetic_concordance.py`、`sm_hp_disambiguation.py`（皆 `hp in (1,11)` / `hpclass()`）。

## 偵測

使用者要求「確認 HP-tag 定義 + 甲基輔助是否合理」→ 對 **ReadParser.cpp 源碼**核對 HP 語意（正確）→ 對 **KB longphase-s 文件**核對（9 態詞彙）→ **抽真實 BAM** 發現 HP tag 是 `HP:Z:` 字串 → pysam 實測 `get_tag` 回 str 且 `in (1,11)` 全 False → 定位 bug。**不是程式報錯，是靜默假陰性**（`0` 看起來像合法結果）。

## 影響（blast radius）

**受影響（假象，已修正）**：
- `hp_control_eval = 0/740`（reextract）→ L8 §4「structural-zero / subclone-specific 撤回」框架。
- `hp_eval = 0/1267`（single_locus）→ L9 §5「HP 排解 cis-ASM 不可行」。
- 所有「cis-control 不可評估 / HP 無法排解 cis-ASM」結論。

**不受影響（HP-independent，結論不變）**：
- 6.6% genotype-anchored corroboration（49/740）— 用 REF/ALT 等位非 HP。
- 單位點 ASM 32%、read×read data-starved 1/754、LOH 41/49（CN bed）、power dose-response。

## 修正後真值（全基因組重跑，corroboration 完美重現 754→740→49）

| 指標 | bug 版 | 修正後 |
|---|---|---|
| hp_control_eval 可評估 | 0/49 | **14/49**（neutral 8/8 + LOH 6/41）|
| cis_explained＝germline-cis（甲基非 subclone）| 0 | **11** |
| subclone/somatic-cis 候選 | — | **3** |

→ 結論**反轉且更強**：HP 排解 cis-ASM **可行**；where testable，甲基 corroboration **多為 germline 等位特異甲基化而非 subclone** —— 比原「structural-zero 推託」更坐實「甲基非 subclone 偵測器」。BOUNDED_AUXILIARY 裁決方向不變。

## 修正

- 4 腳本改 `str(reads[rn][2]) in ("1","1-1","11")`（HP1）/ `("2","2-1","21")`（HP2）—— `str()` 同時兼容 longphase-s 字串與 longphase-to 整數，對齊 `ReadParser.cpp:135-139` 的正規化意圖。
- 全基因組重跑 → canonical 數據（`sm_methyl_reextract_ALL.json/_perregion.tsv`）替換。
- L8/L9 + 05/06/00_INDEX + index.standalone.html 全數修正（prose=JSON=HTML 三方一致）。

## 教訓（reusable）

1. **🔴 pysam tag 型別陷阱**：`get_tag` 依儲存型別回 str（`Z`/`H`）或 int（`i` 等）。比對 HP/任何 tag **必先正規化**（`str()` 或統一型別），尤其同一 pipeline 可能混 longphase-s（字串）與 longphase-to（整數）兩種 tag 編碼。
2. **🔴 靜默 universal-zero = 可疑 bug 非生物學**：任何「N/總數 = 0 / 從未發生一次」的結果，先質疑「輸入是否真的被讀到」（加一行 assert「at least one non-default」），別直接當成生物學陰性寫進結論。本案 `hp_control_eval=0` 被當「LOH 無 HP」寫了兩份報告，實際 BAM 滿是 HP。
3. **feature 定義必查源碼 + 真實資料**：`feedback_feature_name_vs_definition_rule` 的延伸 —— 不只查 C++ 定義，還要**抽真實 BAM 驗證該欄位實際長相**（型別/取值），定義對不代表 parse 對。

關聯：`feedback_feature_name_vs_definition_rule`、`feedback_researcher_claim_needs_empirical_verification`、`project_clone_subclone_integrated_report_finalized`。
