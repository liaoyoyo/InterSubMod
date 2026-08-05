# V3 — G-B / cis-control：甲基節點級標記能否升 confirmer

**軌**: V3 (G-B / T-GATE-GB matched-normal cis-control)
**樣本**: HCC1395 單樣本 ⭐3 · tumor.bam + normal.bam (big8) · hg38 ONT 5mCG
**日期**: 2026-07-04 · **狀態**: COMPLETED（新跑 correct-null PoC，43 區→39 有 CpG）

---

## 一句話裁決

**甲基節點級 ASM 標記無法從 annotation 升 confirmer。** 我在 43 個乾淨 CROSS-HP 區
（`cleanest_applicable_set`）上**首次實際執行**了過去只分類、從未跑過的 **HP-oriented
matched-normal correct null**。結果：正確 null **對混淆無槓桿**——corr(tumor Δβ, normal Δβ)
= **+0.008**（逐區中位 −0.002），扣除不但沒縮小訊號反而**放大變異**（median |Δβ| 0.076 →
|residual| 0.175），84% 的 tumor-significant CpG 在扣除後**存活**。因此 T-GATE-GB 的驗證閘
**不開**：normal-baseline 甲基**未能**分離出「獨立於突變 cis」的成分，只是注入 germline-ASM 噪音。
甲基維持 **Tier-3 bounded-auxiliary / L3 候選旗標**，與 framework §4 紅線一致。

---

## Q1｜G-B 正確 null 跑了沒？用對 null 沒？

### 裁決：**改前 = NO（正確 null 從未在可用區跑過）；改後 = 我補跑了，結論見下。**

先前 verdict（`20260628_cis_control_scope_pilot_verdict_01.md`）的所有數字是**scope 分類**，
不是正確 null 的執行結果：

| 事實（原始碼/JSON 溯源） | 值 | 來源 |
|---|---|---|
| pilot 用的 normal 對照 = `nH1 − nH2`（germline-het HP 差）| — | `scripts/pilot_cis_control.py`（`nH1=per_cpg(nmeth,set(h1)); nd=nH1[c]-nH2[c]`）|
| pilot 6 個有 CpG 的區 **全 LOH** | 6/6 loh | `cis_control_scope_summary.json → pilot_regions_with_cpg` |
| LOH → 兩 subclone 群同一 HP（SAME-HP）→ `nH1−nH2` **正交（by construction）** | corr −0.026 | verdict §4；`joint_distribution.corr` |
| 43 個乾淨 CROSS-HP 區僅有**區名清單**，**無 residual 計算** | 43（僅 examples 8 個列名）| `cis_control_scope_summary.json → cleanest_applicable_set`（無 residual 欄）；`clean_subset.json`（只有 shape 計數）|

**判別（回答 PI 原問「正確 null 還是 germline-het 錯對照」）**：pilot 用的是 germline-het
`nH1−nH2` 對照，但**套在 SAME-HP（LOH）區**——那裡兩 subclone 群共用單一 germline HP，
germline-ASM **本來就自動抵消**，`nH1−nH2` 對 tumor 的 A−B 差**貢獻結構性為零**（verdict
§10.5 角度1 代數已證：`(A−N)−(B−N)=A−B`，matched-normal 在 SAME-HP 唯一非冗餘角色 =
per-arm 極性錨 `A−N`/`B−N`，非差異對照）。所以 corr≈0 是**恆等式，不是「無混淆」的好消息**。
→ 既有結論 = **germline-het 對照套在它結構上無效的區**；真正該跑正確 null 的 CROSS-HP 區
**從未跑過**。

---

## Q2｜正確 null 未跑 → 我實際跑了一個小 PoC（43 乾淨 CROSS-HP 區）

### 設計（`scripts/gb_correct_null_poc.py`，複用 pilot/classify 的讀取函式，無漂移）

- 對每個乾淨 CROSS-HP 區：tumor 依 genotype 取 top-2 ALT 群（pop A/B），確認**讀層 dominant HP
  不同**（hpA≠hpB，真 read-level CROSS）；
- normal 取**同兩條 HP** 的 reads，**方向對齊**（A→hpA, B→hpB）算 germline-ASM 基線 `nHa−nHb`
  = **正確 null**；
- 每 CpG：`residual = (tA−tB) − (nHa−nHb)`。此為 CROSS-HP 唯一結構上有效的 normal 扣除
  （popA/popB 分居兩 HP → 其 Δβ 真的含 HPa-vs-HPb germline 成分）。

### 結果（每個數字溯源 `gb_correct_null_poc.json`）

| 量 | 值 | 意義 |
|---|--:|---|
| 輸入區 / 有 CpG 區 | 43 / **39** | 4 區無跨四組共有 CpG（`funnel`）|
| CpG 對總數 | **8,494** | `summary.n_cpg_pairs_total` |
| **corr(tumor Δβ, normal Δβ)** | **+0.008** | ≈0；逐區中位 **−0.002**（範圍 −0.177~+0.147）|
| median \|tumor Δβ\| | 0.076 | 扣除前 |
| median \|normal Δβ\|（germline-ASM）| 0.093 | 基線本身有可觀 spread（p90=0.719）|
| **median \|residual\|** | **0.175** | 扣除後 **放大**（>兩輸入）|
| tumor-significant CpG (\|Δβ\|≥0.2) | 2,006（23.6%）| 分母 |
| **↳ residual 存活 \|≥0.2\|** | **1,686（84.0%）** | 正確 null **移除不掉**的 somatic 訊號 |
| ↳ residual 崩解 \|<0.1\| | 139（6.9%）| germline-ASM 能解釋的極少數 |
| 逐區「多數 sig 崩解」的區數（≥5 sig CpG）| **0 / 39** | 無任何區被 germline-ASM 主導 |
| 逐區「多數 sig 存活」的區數 | **39 / 39** | 全數殘差主導 |

### 判讀（為何 corr≈0 + 高存活 = confirmer 失敗，而非好消息）

1. **正確 null 是 inert 的，不是「無混淆」**：CROSS-HP 區的 tumor A−B 差與 normal Ha−Hb 差
   **統計獨立**（corr≈0）。若 germline-ASM 是 tumor 訊號裡的真實共享成分，應為**正相關**；
   實測 ≈0 → normal 基線**沒有共享混淆可扣**。
2. **扣除放大變異（Var(A−B)=Var A+Var B−2cov；cov≈0 → 加總）**：median |Δβ| 0.076→|res|0.175
   → 拿正確 null 當「清洗」步驟會**讓甲基訊號更糟**，不是更乾淨。
3. **存活的 84% 不是「被驗證的 subclone marker」**：它 ≈ 未變的 tumor A−B 訊號，而該訊號本身
   = **somatic-cis footprint**（chr17 α_cis=23 ≫ lineage=6，framework §1）；單一 bulk
   **無法識別** footprint vs lineage（verdict §10.5：normal 無 somatic 軸）。存活只代表「訊號
   是真的、是 somatic」，**不代表「是獨立的譜系證據」**。

---

## Q3｜甲基節點級標記（Thread-2 ③④）能否升 confirmer？

### 裁決：**不能。維持 annotation / Tier-3 bounded-auxiliary（L3 候選旗標）。**

T-GATE-GB 驗證閘（framework §4 紅線#4）有**兩個**條件：

1. ✅ **cis-control 跑通** — 我跑了（39/43 區、8,494 CpG）。
2. ❌ **證 normal-baseline 甲基獨立於突變 cis** — **失敗**。扣除是 inert（corr≈0），殘差 ≈ tumor
   訊號 = 突變-cis 主導；normal 基線**沒有**分離出獨立於 cis 的成分，只注入 germline-ASM 噪音。

→ 兩條件必須同時成立才升級；條件 2 不成立 → **閘不開**。

confirmer 需要 normal-same-HP 基線能**獨立確認**某甲基差異是 subclone-specific（非 germline、
非 somatic-cis 後果）。PoC 直接證明正確 null **既不能移混淆（無共享混淆可移）也不能確認
（殘差留在 cis-footprint vs lineage 不可識別的同一狀態）**。因此：

- **audit 文件的 C′**（`within_hp_subclone_permanova`，a-priori germline-tag/carrier-tag 標籤）
  維持 **"mark for confirmation" 非 confirmer** — 本 PoC 證實所缺的 cis-control 步驟
  **不能**把它從標記升成確認。
- 甲基節點級 ASM 標記在拓樸 DAG 上維持 **auxiliary annotation（Tier-3, w↓, cis-control-pending
  永遠標 pending）**，永不 override genetic/HP（framework §4 紅線#1）。

---

## 誠實邊界 + BLOCK 缺口

- **⭐3 單樣本**：本 PoC 是 HCC1395 單一 bulk 的 exploration，非跨樣本驗證。
- **結構性不可識別（非覆蓋問題）**：即使正確 null 跑通，單一 bulk 仍**無法**區分殘差是
  subclone lineage marker vs somatic-mutation-cis footprint（normal 無 somatic 軸重建
  within-HP subclone）。
- **BLOCK — 要真正升 confirmer 需要（皆現不可及）**：
  1. **single-cell 甲基真值**（錨定 lineage vs footprint）；
  2. **COLO829 / 第二樣本 multi-region**（≥5/7 樣本，framework §4 紅線#6）；
  3. **per-read CN-aware 分群**解 CN-gain 610 區的 multiplicity（解鎖更多 CROSS-HP，但仍不解
     單-bulk 不可識別）。
- 這些是 framework §9 / scientific-rigor §8.3.1 的 **Reopen 條件 C1/C3**，非本輪可補。

---

## 產物（皆存 `docs/methodology/_assets/20260704_V_verification/V3_GB_cis_control/`）

- `gb_correct_null_poc.py` — 正確 null PoC 腳本（複用 pilot/classify 讀取函式）
- `gb_correct_null_poc.json` — 全輸出（summary + 39 區 + 8,494 CpG 對，每數字 grep-able）
- `_v3_clean_cross_hp_regions.json` — 43 區清單（從 `hp_alignment_fullscan.json` 重生）
- `gb_correct_null_poc.png` — 3-panel 診斷圖（見下）

**溯源**：所有輸入凍結於 `hp_alignment_fullscan.json`（1874 區分類）+ tumor/normal BAM
（big8）+ 同 germline phased VCF（md5 一致，verdict §10.5 角度2 已證 tumor/normal HP 指同一
實體染色體，98.85% allele 一致）。
