<!--
建立時間: 2026-06-23
類型: 方法說明 (methodology spec)
狀態: in_progress (exploratory)
tier: ⭐2-3 (single-sample HCC1395, 無外部真值)
影響 track: tumor-only subclone 向 (merged 當 cis-control 對照)
data_sources: _assets/20260618_subcluster_pilot/phylo_weak_structure_counts.json, _assets/20260618_subcluster_pilot/geom_summary.json, _assets/20260618_subcluster_pilot/geom_records.json, _assets/20260618_subcluster_pilot/phylo_cpp_wg_records.json
adversarial_verification: workflow wmgsefb9d (4-lens, SHIP_WITH_CORRECTIONS)
-->

# 弱結構保留分類 + geometry-vs-null 分歧 — 方法說明

> ⚠ **PARTIAL / 探索級**：單樣本 HCC1395、tumor-only。全基因組 5 態分類涵蓋全 34,736 位點；geometry 層為 **~400 位點抽樣**（非全量）。無外部真值，結論為 characterization 不是 subclone 宣稱。tier ⭐2-3。

## 0. TL;DR

phylo-v4.1 的嚴 verdict（`coarse_ng`）把 TP 的 **71.01%** 標成「無結構」，但這層把**弱結構與真乾淨混在一起**。本方法把每位點拆成 **互斥 5 態決策流**，顯式分出弱結構候選（不丟進黑洞）；再加一個正交的 **`geometry_ng`**（scipy 0.7×max 純幾何平切）量化「**肉眼/幾何看得到、但 null 閘擋掉**」的分歧。抽樣 400 位點重跑得 **124/400（31%）分歧**。

🔴 **對抗驗證（workflow wmgsefb9d）關鍵更正**：弱結構候選（S4）**TP/FP=1.15 = chance-level，與噪音統計不可分** → 是「待驗證候選」不是「已保留的真結構」；`geometry_ng` 是「null 閘保守度」的診斷，**本就會過切單一 clone，非 subclone**。

---

## 1. 動機：兩種裁切準則 + 不漏掉弱結構

同一棵 **UPGMA 樹**（同一 BERNOULLI 距離），可用兩種準則裁切歸群：

| 準則 | 機制 | 看資料？ | 噪音行為 |
|---|---|---|---|
| **C++ phylo null 閘** | 遞迴 between/within `r≥SEP_MIN` 且 `r>` 欄內置換 null 第 N 百分位 | 看甲基欄結構 + 顯著性 + 平衡 + 穩定 | 擋不過 → 報 1 群（保守） |
| **scipy 幾何平切** | `color_threshold = 0.7×max(Z[:,2])` 單一高度切 | 只看樹幾何（合併高度） | 高度有變異就切 → 一定 ≥2（會過切） |

需求：使用者要「**不漏掉弱結構**」（可能是有用資訊）+「**可觀察確認**」。→ 本方法保留弱訊號欄位 + 顯式分桶 + 加幾何欄 + 建逐位點觀察工作站。

---

## 2. 變數清單（有哪些變數）

### 2.1 每位點輸入變數（`phylo_cpp_wg_records.json`，C++ native，34,736 位點）

| 變數 | 定義 | 型別 |
|---|---|---|
| `chrom`, `pos` | 染色體 + region 起點（= SNV − 5000，window 起點） | str/int |
| `set` | `TP` / `FP`（pileup track；FP=germline/artifact-enriched 負類） | enum |
| `n` | 通過 peel（filter-for-complete-matrix）的 read 數 | int |
| `coarse_ng` | **null95 嚴閘**群數，K=10 modal 多數決（主 verdict） | int |
| `fine_ng` | **null90 寬閘**群數，單一 base seed（候選低信心） | int |
| `n_other` | descend-quarantine 殘離群數（≥MINSZ 才標 other） | int |
| `unstable` | `modal_frac<0.7`（K=10 seed 群數分歧） | bool |
| `hidden_het` | `n_other/n > 0.30`（隱藏異質） | bool |
| `aligned` | `max(V_hp,V_allele) ≥ 0.3`（cis-ASM 指紋，**非 subclone 檢定**） | bool |
| `V_hp`, `V_allele` | CramerV(cluster, hp) / CramerV(cluster, alt_support) | float |

### 2.2 geometry 子集新增變數（`geom_records.json`，~400 位點）

| 變數 | 定義 |
|---|---|
| `geometry_ng` | `fcluster(Z, 0.7×max, 'distance')` 後 size≥MINSZ 的群數（= scipy 預設樹配色等價，RNG-free） |
| `geom_flat_raw` | 平切原始群數（含 singleton） |
| `geom_cut` | 0.7×max 切高 |
| `divergence` | `geometry_ng≥2 AND coarse_ng<2`（幾何看到但 null 擋掉） |
| `geom_gt_coarse` | `geometry_ng > coarse_ng` |

### 2.3 算法常數（`scripts/phylo_v4.py`）

`MINSZ=3` · `SEP_MIN=1.3` · `RNULL=40` · null95(coarse)/null90(fine) · `K=10` modal · `MODAL_CONF=0.7` · `HIDDEN_HET=0.30` · 距離=BERNOULLI(min_common=3) · 樹=UPGMA(average linkage) · 幾何=scipy `color_threshold=0.7×max(Z[:,2])`。

---

## 3. 分類決策流（有哪些分類）

**互斥 5 態**（ordered first-match；partition 已獨立驗證：0 重疊 / 0 缺口 / sum-check ✓）：

```
每位點:
  if coarse_ng≥2 and not unstable and aligned      → S1 確認多群·對齊
  elif coarse_ng≥2 and not unstable and not aligned → S2 確認多群·不對齊
  elif coarse_ng≥2 and unstable                     → S3 不穩定多群
  elif coarse_ng<2 and fine_ng≥2                    → S4 次閾值候選
  elif coarse_ng<2 and fine_ng<2 and (n_other>0 or hidden_het) → (S5 殘留)
  else                                              → S6 乾淨單群
```

| 態 | 意義 | 性質 |
|---|---|---|
| **S1** 確認多群·對齊 | 顯著+穩定+對齊 hp/allele | `aligned`=cis-ASM 指紋，**非 subclone** |
| **S2** 確認多群·不對齊 | 顯著+穩定，無標籤對應 | 結構但未對齊 |
| **S3** 不穩定多群 | coarse 報多群但 seed 分歧 | borderline，留候選 |
| **S4** 次閾值候選 | 嚴閘併成 1、寬閘切 ≥2 | 🔴 **chance-level（見 §5），待驗證候選非真結構** |
| ~~**S5** 殘留~~ | coarse<2 且有殘離群 | **結構性恆 0**（`n_other`/`hidden_het` 只在 coarse≥2 出現，與 coarse<2 矛盾 = unsatisfiable；親驗 1896/1896 `n_other>0` 皆 coarse≥2）→ 實為 5 態 |
| **S6** 乾淨單群 | 兩解析度都 1 群、無殘留 | 唯一真「無弱訊號」桶 |

---

## 4. 結果數量與比例

### 4.1 全基因組 5 態（N=34,736 = TP 30,077 + FP 4,659）

| 態 | TP n | TP % | FP n | FP % |
|---|---:|---:|---:|---:|
| S1 確認多群·對齊 | 5,538 | 18.41 | 1,349 | 28.95 |
| S2 確認多群·不對齊 | 1,878 | 6.24 | 312 | 6.70 |
| S3 不穩定多群 | 1,302 | 4.33 | 177 | 3.80 |
| S4 次閾值候選 | 1,796 | 5.97 | 242 | 5.19 |
| ~~S5 殘留~~ | 0 | 0.00 | 0 | 0.00 |
| S6 乾淨單群 | 19,563 | 65.04 | 2,579 | 55.36 |

- confident-multi（S1+S2）= TP **24.66%** / FP **35.65%** → **FP 反而更多（反判別）**。
- 弱結構候選（S3+S4）= TP **10.30%** / FP 8.99%。
- 真正乾淨（S6）= TP **65.04%**。

### 4.2 舊二分 vs 新分類（弱結構不再被黑洞吃掉）

舊「structure vs no_structure」把 coarse<2 一律當「無結構」= TP **21,359（71.01%）**。新分類拆出其中 **1,796（5.97%）= S4 有弱訊號**，真正乾淨只剩 S6 65.04%。

### 4.3 geometry-vs-null 分歧（子集抽樣 400 位點）

抽樣組成：S6 n≥12 共 200（TP150/FP50）+ S4 120（TP90/FP30）+ S1/S2/S3 對照 80。重跑 binary 重生矩陣 → 算 `geometry_ng`。

- **分歧（`geometry_ng≥2` 但 `coarse_ng<2`）= 124 / 400 = 31.0%**；分布 S6 64 + S4 60。
- `geometry_ng≥2` 總 185；`geometry_ng > coarse_ng` 總 134；coarse==geom 一致 233/400；反向（coarse≥2 但 geom<2）19。

**coarse(null95) × geometry 交叉表**（橘 = 分歧區）：

| null95 ＼ geom | geom1 | geom2 | geom3 | geom4 | geom5+ |
|---|---:|---:|---:|---:|---:|
| **coarse1** | 196 | **51** | **39** | **11** | **23** |
| **coarse2** | 15 | 27 | 4 | 3 | 2 |
| **coarse3+** | 4 | 10 | 11 | 4 | — |

→ `coarse1` 列右側（geom≥2）合計 **124 = 分歧數**。代表案例：`chr11:61576011`（S6, n76, null95=1 fine=1 **幾何=8**，同質單群被幾何切 8）、`chr7:154546475`（S4, n153, null95=1 fine=4 幾何=8）。

---

## 5. 🔴 對抗驗證修正 + 必掛 caveat

workflow `wmgsefb9d`（4-lens）裁決 **SHIP_WITH_CORRECTIONS**，以下為**引用本數據時必須同時陳述**的修正：

1. **S4 ≠ 已保留的真弱結構**：只過 null90（α=0.10）、被 null95（α=0.05）擋；**TP 5.97% vs FP 5.19% = ratio 1.15（chance-level）**，與噪音統計不可分。正名為 **sub-threshold candidates pending validation**。
2. **S4/S6 的 `aligned=False` 非實證**：`aligned` 只在 coarse≥2 計算（0/24,180 coarse=1 被算過），對 coarse<2 應標 **NA/未測**，非 False。
3. **confident-multi 反判別**：FP 35.65% > TP 24.66%（**預期**：germline het 帶 allele-specific 甲基 → 真實 cis-ASM）。`aligned`=cis-ASM 指紋，**不是 subclone 檢定**。本分類是**觀察/刻畫層，不可當 TP/FP filter 或自動 subclone caller**。
4. **`geometry_ng` 是診斷非宣稱**：純 0.7×max 平切無 null，**本就會過切單一 clone**（gap/silhouette 已證 1 clone 被判 k=3-4）。它量化 null 閘的保守度，**非 subclone 結構**。
5. **無多重檢定校正**：34,736 位點 + 遞迴 within-locus 分裂皆無 FDR；`RNULL=40` 的百分位 Monte-Carlo SE 大 → null90 vs null95 區分本身粗糙。
6. **partition 與計數可信**：獨立重算 0 重疊/0 缺口/sum-check 到個位數吻合；**修正在「詮釋與標籤」非算術**。
7. **scope**：單樣本 HCC1395 探索 ⭐2-3，無外部真值；裁定 cis-ASM vs subclone 需 normal cis-control / multi-sSNV CCF / single-cell（本方法皆未做）。baseline-dependence：甲基分群結構在扣 cis/allele baseline 前非獨立結果。

---

## 6. 格式選擇（使用哪種格式適合）

本工作的三層產物各有適配格式，**分工而非擇一**：

| 格式 | 檔案 | 適合什麼 | 為什麼 |
|---|---|---|---|
| **Markdown 方法說明**（本檔） | `20260623_weak_structure_classification_and_geometry_divergence_01.md` | 規格 + 數據表 + 推理 + caveat 的**留底 SoT** | grep-able、git diff、INDEX 登錄、跨 session 引用、純文字校對 |
| **互動 standalone HTML 工作站** | `_assets/.../20260623_geometry_divergence_observation_dashboard.standalone.html` | **逐位點肉眼觀察 + 人工裁決** | path-based 400 圖 lazy-load、篩選（set/態/★分歧）、judgment localStorage + JSON/CSV 匯出 |
| **JSON 資料層** | `phylo_weak_structure_counts.json` / `geom_records.json` | **機讀真值 / 防捏造注入源** | 數字 SoT、可程式查詢、§13-A 由構造注入（缺值 refuse） |

> **建議**：方法/比例/分類 → **本 .md**（留底、引用、論文素材）；要「查詢與觀察分類」→ **HTML 工作站**；要重算/再分析 → **JSON**。此即本專案 `html-report-build` + `verify-workstation` + §13-A 反捏造的標準三件式。給 PI 終版時，可由本 .md 經 `/html-report-build --mode standalone` 生 companion HTML（sticky TOC + tier badge）。

---

## 7. 可重現性

```bash
# 全基因組 5 態分類（純讀 records，零 compute）
python3 scripts/phylo_weak_structure_classify.py
#   → phylo_weak_structure_counts.json + phylo_weak_structure_classified.json

# geometry 子集重跑（抽樣 400 → mini-VCF(snv=rec_pos+5000) → binary → geometry_ng）
python3 scripts/geom_divergence_run.py      # ~42s, 背景安全
#   → geom_records.json + geom_summary.json（矩陣保留於 output/_geomdiv_{TP,FP}）

# 渲染 400 觀察圖（Pool 平行，雙色欄 C++ null95 vs 幾何 + 6 側欄）
python3 scripts/geom_render.py               # → figs_geomdiv/

# 生成觀察工作站（JSON 注入，§13-A 不手打數字）
python3 scripts/build_geom_dashboard.py
```

binary = `build/bin/inter_sub_mod`（feat/summary-nreadsvalid，phylo c186092）；BAM/REF/VCF 見各腳本頂部路徑。

## 8. Provenance 溯源表（metric → 來源）

| metric | 值 | 來源檔案 |
|---|---|---|
| WG N / TP / FP | 34736 / 30077 / 4659 | `phylo_weak_structure_counts.json` `n_total`,`TP.n`,`FP.n` |
| 5 態計數 TP/FP | §4.1 表 | `phylo_weak_structure_counts.json` `TP/FP.states.*` |
| confident-multi TP/FP | 24.66 / 35.65 | `..counts.json` `*.multi_confident_pct` |
| 弱結構候選 / 乾淨 TP | 10.30 / 65.04 | `..counts.json` `TP.weak_candidate_pct`,`truly_clean_pct` |
| 舊無結構 / 弱藏其中 | 21359(71.01) / 1796(5.97) | `..counts.json` `TP.old_no_structure_*`,`weak_hidden_*` |
| 分歧 124/400(31%) + 分布 | 124 / S6 64 S4 60 | `geom_summary.json` `n_divergence`,`divergence_by_wstate` |
| coarse×geometry 交叉表 | §4.3 表 | `geom_summary.json` `coarse_vs_geometry_xtab` |
| 反向 19 / 一致 233 | 19 / 233 | `geom_records.json`（重算） |
| S4 TP/FP 1.15 / aligned tautology / FDR 缺 | — | workflow `wmgsefb9d` Lens1（4-lens 對抗驗證） |

## 9. 限制與下一步

- **限制**：單樣本、geometry 層僅 400 抽樣（非全 ~21,852 分歧候選池）、無 FDR、無 normal cis-control。
- **下一步（可選）**：(a) 擴跑全分歧池；(b) 加 FDR/null-CI 校正；(c) normal cis-control 區分 cis-ASM vs subclone；(d) 經 `/html-report-build` 生 PI standalone companion。

---
> 關聯：memory `[[project_phylo_subcluster_labeling_doubledip_fix]]`、`[[project_tumor_only_axis_negative_subclone_classification]]`；上游 spec `20260622_phylo_v3_method_spec_and_statistics_01.md`（主repo）。
