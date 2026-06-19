<!--
建立時間: 2026-06-18
類型: cpp-change spec（待核准 → 進 /cpp-change 6 步 PDD）
狀態: pending_approval
目標 branch: feat/summary-nreadsvalid（within-HP + Bootstrap class 在此；研究 branch research/subclonal-reconstruction-202606@b761336 無此碼）
依據: docs/methodology/20260616_structure_label_verdict_false_negative_audit_01.md §7（B+D）
data_sources: docs/experiments/in_progress/2026/06/20260614_u1_maxdist_vs_skip/method_audit/within_hp_clean_pilot.json（@4ec339b）
-->

# within-HP bootstrap-stability + cluster×label ARI emit — 最小 cpp-change spec

> 修 §7 盤點出的 ③② gap。**最小**：只動 2 處、復用既有 `Bootstrap` class，不啟用昂貴的 main-clustering bootstrap、不把完整 typology 移進 C++。

## 0. 目的與範圍

| 變更 | 修哪個 gap | 風險 | 順序 |
|---|---|---|---|
| **C1** within-HP 子分群加 bootstrap-stability | ③：within-HP 用 raw silhouette → 14%↔41% 閾值漂移（循環） | 中（改分群認列行為） | 後 |
| **C2** 算 + emit `ARI(cluster, label)` | ②：對齊強度從沒被量化 | 低（純加觀測欄） | 先 |

**明確不做（留待後續）**：① 啟用 main-clustering `enable_bootstrap`（昂貴 + cluster-first 不進 verdict，無收益）；② 把 `StructureNoLabel/MultiGroupNoLabel` 三 case typology 移進 C++（先留 offline Python，靠 C2 的 ARI 欄即可離線分型）。

---

## 1. C2 — `ARI(cluster, label)` 計算 + emit（先做，低風險純加）

### 1.1 問題
「幾何群是否對齊標籤」目前只有離散 FFH p / CramersV，沒有對齊指標。`adjusted_rand_index()` 已實作但只用於 bootstrap 比較兩次分群（`Bootstrap.cpp:154`），**從沒拿來比「幾何群 vs 標籤」**。

### 1.2 設計（2 行計算 + 2 欄輸出）
- 位置：`SignificanceAnalyzer::analyze`，`cluster_labels`（silhouette 幾何群）與 `hp_merged_labels`/`allele_labels`（已於 `:284-304` 建好）皆在 scope。
- 計算：
  - `result.ari_cluster_hp     = Bootstrap::adjusted_rand_index(cluster_labels, hp_merged_labels);`
  - `result.ari_cluster_allele = Bootstrap::adjusted_rand_index(cluster_labels, allele_labels);`
  - 缺一邊標籤（如全 -1）時回 NaN，不阻斷。
- struct：`SignificanceResult` 加 `double ari_cluster_hp = NAN; double ari_cluster_allele = NAN;`
- emit：CSV header（`RegionProcessor.cpp:1338` 區）+ 值（`:1477` 區）加 `ARI_Cluster_HP,ARI_Cluster_Allele`。

### 1.3 三 case 可離線分型（不需 C++ typology）
有了 ARI 欄，offline 即可分：
- **對齊（1↔1）**：`ARI_Cluster_HP` 或 `ARI_Cluster_Allele` 高（如 ≥0.5）。
- **沒對齊（多標籤→1 群 / 結構非標籤）**：`optimal_k≥2` 且 ARI≈0。
- **標籤內次結構（1 標籤→多群）**：`WithinHP_CleanMultigroup=true`（已有欄）。

### 1.4 Step→Verify（C2）
- `cd build && make -j$(nproc)`；`ctest --output-on-failure` 全綠。
- ARI ∈ [-1,1]；抽 legacy `Strong` 位點 ARI 應偏高、`StructureNoLabel`/`Subclone` 應≈0（與 reclassify_v2 typology 對照）。
- regression 雙守護（SKIP+MAX_DIST golden）：**新增欄不得改既有欄位值**（純加）。

---

## 2. C1 — within-HP 子分群加 bootstrap-stability（後做，改行為）

### 2.1 問題（鐵證）
`compute_within_hp_substructure`（`RegionProcessor.cpp:2079`，呼叫 `:1180-1188`）用 raw silhouette + heuristic 閘（sil≥0.5+balanced+gap>0.15）判「裂幾群」。實測 loose(sil>0.3)=**98.25%** / clean=**14.2%**（隨閾值 14→27→41% 漂移）→ silhouette 過切，「裂幾群」無原則答案。

### 2.2 設計（復用 `Bootstrap`，限 within-HP）
- 在 `compute_within_hp_substructure` 內、得到候選 `ngroups≥2` 後，加**穩定度驗證**：
  - 對該 HP-family 的 reads bootstrap 重抽 B 次（重抽 reads，B 預設 **50**），每次以**同分群函式**重新分群，算 co-clustering 一致性（pairwise 同群比例平均）或對 original 的 mean ARI。
  - 只當 `stability ≥ τ`（預設 τ：mean co-cluster ≥ **0.8** 或 ARI ≥ **0.5**，二選一，spec 取 co-cluster≥0.8）才認列 `within_hpX_ngroups`；否則 collapse 回 1。
  - 復用既有 `Bootstrap`（其 `compute` 已接受 `cluster_func`）；within-HP 專用輕量 config（B、seed、τ），**與主分群 bootstrap 各自獨立**。
- 新 config：`within_hp_bootstrap_enable`（預設 **true** — 這是修正核心，非預設關）、`within_hp_stability_tau`、`within_hp_bootstrap_B`。
- struct：加 `double within_hp1_stability; double within_hp2_stability;`
- emit：CSV 加 `WithinHP1_Stability,WithinHP2_Stability`。
- `within_hp_clean_multigroup`（`:1188`）的定義改為：PATTERN 多群**且通過穩定度** OR level bimodal（level 軸本就靠 mean-β 雙峰，較穩，維持）。

### 2.3 對應檔
- `src/core/RegionProcessor.cpp:2079`（`compute_within_hp_substructure` 內加 bootstrap）、`:1180-1189`（呼叫處 + clean_multigroup 定義）、`:1338`/`:1477`（header+值加 stability 欄）。
- `include/core/RegionProcessor.hpp:252-257`（struct 欄）、Config 加 3 個 within-HP bootstrap 參數。
- 復用 `src/core/Bootstrap.cpp`（不改）。

### 2.4 Step→Verify（C1）
- `make` + `ctest` 全綠。
- **預期方向**：within-HP `clean_multigroup` 比例**從 14–27% 下降**（過切假群被穩定度擋掉）；對 **chr2:18M HCC1395 已知 subclone 位點**（`project_chr2_18m_subclone_locus_verification`），within-HP stability 應**偏高、仍保留**（真 subclone 撐得住重抽）。→ 真降假、留真。
- regression 雙守護 PASS；C1 會改 `WithinHP*` 欄與依賴它的 verdict（`:1265-1270` OR 條件）→ golden 需**有意更新**並記錄差異（非 bit-identical）。

---

## 3. 風險 / 平衡

- **compute 成本**：within-HP bootstrap 只在「該 region 有 ≥2 HP-family 且每 family reads ≥ min」時跑，B=50 限制；對 30k region 屬可接受（主 compute 仍是距離矩陣）。先在單樣本量 wall-clock，再決定 B。
- **τ 閾值仍是參數**：但「重抽穩定度」比「單次 silhouette 高低」原則得多（雜訊假群在重抽下散掉）；τ 可用 chr2:18M 真 subclone + permutation-null 兩端校準，不憑空設。
- **不碰既有 verdict 主邏輯**（D1/D2 已落地）：本 spec 只「把 within-HP 認列加穩定度把關」+「加 ARI 觀測欄」。

## 4. 驗收標準（總）

- [ ] `ctest` 221/221（或當前基準）全綠
- [ ] regression：C2 bit-identical（純加欄）；C1 golden 有意更新 + 差異記錄
- [ ] within-HP `clean_multigroup` 比例下降（過切被擋），chr2:18M subclone 位點仍保留（stability 高）
- [ ] `ARI_Cluster_HP/Allele` ∈[-1,1]，Strong 高 / StructureNoLabel 低
- [ ] 全基因組 F1 / 分類分布 no-regression（對照修前）
- [ ] 單樣本 wall-clock 增幅可接受（B=50）

## 5. cpp-change 6-step 對應

1. **C2 先**（低風險純加觀測欄）：struct 欄 → 計算 2 行 → CSV emit → ctest → 一 commit。
2. **C1 後**（改行為）：Config 參數 → `compute_within_hp_substructure` 加 bootstrap → clean_multigroup 重定義 → CSV stability 欄 → ctest + chr2:18M 驗證 → golden 更新 → 一 commit。
- 每 commit 過 pre-commit 編譯 Hard Gate；在 `feat/summary-nreadsvalid` 上做（或新開 `feat/within-hp-bootstrap-ari` 從該 branch 切，避免並行污染）。

## 6. 用戶核准

**核准進 /cpp-change？** [x] **C2 done**（committed `64e68d9` on `feat/cluster-label-ari`，2026-06-18；ctest **224/224** 含 3 ARI 測試；observation-only 無 F1 影響；worktree `/big7_disk/liaoyoyo2001/wt-cluster-label-ari` 保留）/ [ ] **C1（within-HP bootstrap）待核准** / [ ] 調整 τ・B
**理由 / 調整**：C2 純加觀測欄先落地；待 genome-wide run 看 ARI 三 case 分布 + 用 C2 欄離線驗證後，再決定 C1。
**未完成（交付/收尾）**：C2 commit 未 push / 未 merge（feature branch 隔離）；evidence_ledger 未記（observation-only 無 delta_f1，且 shared-state 待用戶）；worktree 待 C1 或驗證後清理。
