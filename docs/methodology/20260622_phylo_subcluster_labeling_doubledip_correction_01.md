<!--
建立時間: 2026-06-22
問題類型: 統計方法（演化群 phylo labeling 的 null 設計）
影響 track: subclonal reconstruction pilot（tumor-only 切群層）
狀態: validated_correction（v2 已修並驗證；C++ 內建決策 deferred）
data_sources: _assets/20260618_subcluster_pilot/phylo_v2_final.json,_assets/20260618_subcluster_pilot/phylo_v1v2_compare.json,_assets/20260618_subcluster_pilot/phylo_groups.json
build_branch: feat/summary-nreadsvalid
-->

# 演化群 phylo labeling — double-dip 缺陷更正（v1→v2）

> **一句話**：演化群階層 labeling（`phylo_groups.py`）的 split 顯著性檢定有 **double-dipping**（null 沿用真實樹「已最大化分離」的分割、且遞迴 null 只算一次全域）。修正（per-subgroup 重分群 null + RNULL=40）後，**pilot 的兩個「subclone 候選」消失，6/30 多群位點全部對齊 germline 軸（cis-ASM），零 subclone 候選** — 與全基因組「subclone 候選非 somatic-特異」結論完全一致。先前 v1 的 subclone 訊號是統計假象。

## 0. 來龍去脈（三方裁決互相修正）

針對「演化群 phylo 是否健全 / 成果是否正確 / 能否內建進 C++」三問，出現三個不同裁決，**親測（純噪音校準 + v1/v2 對照）才是真相**：

| 議題 | 主 agent ① 裁決 | workflow（3 驗證 agent） | 🟢 親測真相（本文件） |
|---|---|---|---|
| 方法健全？ | 「SOUND」**過強**（漏測噪音校準） | 「UNSOUND」**過強** | **機制有缺陷，但影響被 CpG 數閘住** |
| 純噪音 FP (n=20) | 沒測 | 「41.5%」**未重現** | root 1.0% / 完整遞迴僅 **C≤10 才 16–24%** |
| double-dip 真假 | 漏看 | 真 ✓ | **真**（重分群對照 FP→0%）✓ |
| pilot 成果 | 「清楚正確」 | 「只信 n≥64」過嚴 | **root 多群全可信；深層標籤是假象** |
| C++ 內建 | 「GO」**過早** | 「PROBE / 先修方法」✓ | **先修（v2 已修）再 decide C++** ✓ |

> 教訓：主 agent 不該在沒跑純噪音校準前宣稱 SOUND；但也不該全信 workflow 的 41.5%/UNSOUND 而不獨立驗證（§13.7）。**親測 reconcile 兩者**才得正解。

## 1. 缺陷（mechanism）

`phylo_groups.py::phylo_label` 的 `split_real(c1,c2)`：
- 真實分離比 `r = between(S1,S2)/within`，其中 **S1,S2 是真實 UPGMA 樹「為了最大化此分割」選出的子節點**。
- null `ns=[_bw(Dn,S1,S2) for Dn in nulls]` **沿用同一組 S1,S2** 套在欄內置換的 `Dn` 上 → 但 S1,S2 是為原排列最佳化的，置換後不再匹配 → null 比降到 ~1.0 → `r`（被最佳化膨脹）> null95 = **偽陽性（double-dip）**。
- 額外：`nulls` 只在**全 locus 算一次**，深層子分裂被拿去跟「全域打散」比，而非「子群內打散」→ depth≥1 用錯虛無假設。
- `SEP_MIN=1.3` 是未校準 OK 繃（純噪音曾 r=1.233，僅因 <1.3 沒誤判）；`RNULL=12` 反保守（≈12 抽樣 max）。
- 專案先前在 `cluster_redesign.py::stab_excess` 已修掉同型 double-dip，phylo 又踩回。

## 2. CpG-gated 偽陽性校準（`phylo_noise_calibration.py`，seed=20260622）

純噪音 = 每 read 對每 CpG 獨立抽 Bernoulli(per-CpG 率)，有邊際結構但零 read×read 結構。回傳 ≥2 群 = 偽陽性。

**測試 1 — root split FP（C=40, 25% 缺值, TRIALS=200）**：

| n | 現法(v1) FP | 重分群對照(v2) FP |
|---|---|---|
| 20 | 1.0% | 0.0% |
| 40 | 0.5% | 0.5% |
| 80 | 0.0% | 0.0% |

**測試 2 — 完整遞迴 phylo_label FP × CpG 數（TRIALS=200）**：

| n | C=10 | C=20 | C=40 | C=76(中位) |
|---|---|---|---|---|
| 20 | 16.5% | 11.0% | 0.5% | 0.0% |
| 40 | 24.0% | 3.0% | 0.0% | 0.0% |
| 80 | 22.5% | 1.5% | 0.0% | 0.0% |

→ workflow 的 41.5% **未重現**（高估）；但 double-dip **方向為真**，在 **C≤20（少 CpG）才嚴重**。
→ **pilot 全 30 位點 CpG min=38 / median=93 / max=398，C≤20 的位點數 = 0**（C≤40 僅 1）→ 在此密度 FP≈0%，故 root-level 多群呼叫可信、深層遞迴標籤才是缺陷高發區。

## 3. 修法（v2）+ 前後對照（`phylo_compare_v1v2.py` / `phylo_v2_final.py`）

**v2 修法**（真實分割不變 — UPGMA 節點的子節點本就是該子樹最佳 2-way 分割；只改 null）：
1. **per-subgroup 重分群 null**：對該子群 read 欄內置換 → **重新 UPGMA 分群** → 取其 root 分裂比（消 double-dip）。
2. **遞迴 null 改子群內**（depth≥1 用對的 null）。
3. **RNULL 12→40**（對齊 `kprofile_stability.py`）。

**前後對照（30 位點）**：v1 平均 1.30 群 → v2 平均 1.20 群；6 多群 → 6 多群，但 **3 個位點本質改變**：

| 位點 | n | C | v1 | v2 | 變化 |
|---|---|---|---|---|---|
| chr21_10353822 | 255 | 72 | **4 群**（含 2 個 hp1-1/ALT subclone 候選）| **1 群** | 🔴🔴 完全塌掉（全假象）|
| chr20_30274614 | 52 | 72 | **3 群**（含 within-ALT/hp1-1 群 `1-2-1`）| **2 群** | 🔴 同-allele 第二群消失 |
| chr20_21855867 | 26 | 70 | 2 群 `1-1-1`/`1-2`（depth-3）| 2 群 `1-1`/`1-2` | ⚠️→✅ 深度修淺 |
| chr21_40743336 | 24 | 38 | 1 群 | **2 群** REF/ALT | v2 多救一個 clean cis-ASM（C 最低需留意）|

**v2 最終 6 個多群位點 — 全部 germline 軸分裂（cis-ASM）**：

| 位點 | groups（dominant hp / allele）| 型態 |
|---|---|---|
| chr21_40743336 | hp1-1/REF · hp2-1/ALT | REF/ALT |
| chr22_44981696 | hp2/REF · hp2-1/ALT | REF/ALT |
| chr20_56191455 | hp1/REF · hp2-1/REF(n6) | hp 軸（皆 REF）|
| chr20_21855867 | hp2/REF · hp2-1/ALT | REF/ALT |
| chr20_30274614 | hp1-1/ALT · hp2/REF | REF/ALT |
| chr20_24575172 | hp1/REF · hp2/REF | hp 軸（皆 REF）|

> **🔴 決定性更正**：v1 唯二的「same-germline 內多群」subclone 候選（chr20_30274614 `1-2-1`、chr21_10353822 `2-1-2-1`+`2-2`）在 v2 **全部消失**。修正後 **pilot 零 subclone 候選，每個多群分裂都對齊 germline（REF/ALT 或 hp1/hp2）= cis-ASM**。這與全基因組（REAL_NOVEL subclone 候選非 somatic-特異、單樣本無法當 subclone）+ memory `project_tumor_only_axis_negative_subclone_classification` 完全一致。

## 4. C++ 內建決策（deferred — 先修方法是對的）

- 主 agent ① 評估 C++ 機械成熟度高（UPGMA/BERNOULLI/PERMANOVA-置換範本/Bootstrap 槽/Tree+cycle-guard/methylation 在記憶體 ≈70% 已存在），新增僅 per-CpG null + 遞迴分裂 + 階層標籤。**此評估不變、隨時可用**。
- **但本更正證明：若把 v1 硬化進 Hard-Gate C++ 層，會把假 subclone 候選烤進 binary。** 順序必須是 **Python 把 v2 定稿 → 重驗 → 才 decide C++**（同 workflow PROBE 建議）。
- 輸出格式建議（內建後）：獨立 `phylo_groups.tsv`（read_id/phylo_label/is_outlier/group_n/seed/sep_min/rnull）+ `phylo_groups_summary.json`（per-group dominant hp/allele）；**不塞進 reads.tsv**（穩定 schema vs 衍生結果分層）；Python 端只讀檔畫圖。
- C++ 現有 cluster 輸出用 `TreeCutter::find_optimal_clusters`(silhouette best_k) = pilot 證實不可靠的無監督 k 選法；phylo v2 是更有原則的替代（若內建）。

## 5. 驗證表（每 headline 數字 → 來源 + 重算 + L 級）

| # | 數字 | 來源 | L 級 | 狀態 |
|---|---|---|---|---|
| 1 | 純噪音 FP C=10 n=40 = 24.0% | `phylo_noise_calibration.py` 測試2（seed 20260622, TRIALS=200）| L2（模擬，可重跑）| ✓ |
| 2 | 純噪音 FP C=76 = 0.0% | 同上 | L2 | ✓ |
| 3 | pilot CpG min=38/median=93/C≤20=0 | `phylo_v2_final.json` 各位點 C 欄 + 快取 methylation.csv header | L1（從矩陣重算）| ✓ |
| 4 | v1 1.30→v2 1.20 平均群 | `phylo_v1v2_compare.json`（30 位點）| L1 | ✓ |
| 5 | chr21_10353822 v1 4 群→v2 1 群 | `phylo_v1v2_compare.json` + `phylo_v2_final.json` | L1 | ✓ |
| 6 | v2 6/6 多群皆 germline 軸（0 subclone 候選）| `phylo_v2_final.json` align 欄 + 對齊判讀 | L1 | ✓ |

**限制**：單樣本 HCC1395 tumor-only、30 位點抽樣（chr20-22）、⭐2-3 characterization 非 subclone 確認。v2 修法已驗證消 double-dip，但「群=cis-ASM 非 subclone」仍待 normal cis-control。C=38 的 chr21_40743336（v2 新增多群）在 FP 略高區，clean 100%REF/100%ALT 對齊支持為真但標 borderline。

## 5.5 四方法對照（`phylo_methods_compare.py` → `methods_compare_30loci.json`）

同 30 pilot 位點 + 純噪音 FP，比 silhouette best_k(C++現用) / 4-gate / phylo-v1 / phylo-v2：

| 方法 | null gate | 「無結構」verdict | 輸出 | 純噪音 FP(C=76,n=20/40/80) | 30 位點多群 |
|---|---|---|---|---|---|
| silhouette best_k（C++ raw）| ❌ | ❌ 強制 k≥2 | 單 k + 平標籤 | 100%（構造）| 30/30 |
| **silhouette + PERMANOVA**（C++ 實際）| ⚠️ double-dip | ❌ | 同上 | 🔴 **74.7 / 80.7 / 90.7%** | **30/30 全顯著** |
| 4-gate（cluster_redesign）| ✅ clusterboot-excess | ✅ NO_CLEAR | 5 類 + coarse/fine | ✅ 1.3 / 0.7 / 0.7% | 24/30 |
| phylo-v1（double-dip）| ⚠️ 沿用樹分割 | ✅ | 階層標籤+對齊 | ⚠️ CpG-gated(C≤10:16-24%) | 6/30 |
| **phylo-v2**（修法）| ✅ per-subgroup 重分群 | ✅ | 階層標籤+對齊 | ✅ **~0%** | 6/30 |

> **🔴 關鍵**：C++ 現用 silhouette+PERMANOVA **純噪音 FP 74.7-90.7%（隨 n 升）** — silhouette 挑最大化分離的 k、PERMANOVA 再測那個被最佳化的分割 = double-dip 極端版 → 30/30 位點全判「有顯著結構」=系統性過度宣稱。**這是 pilot 全程證實 silhouette k-選不可靠的量化根因**。
> **specificity 排序**：phylo-v2(~0%) ≈ 4-gate(~1%) ⋙ silhouette+PERMANOVA(75-91%) > silhouette-raw(100%)。
> **4-gate vs phylo-v2 互補非誰贏**：4-gate=「有無可重現結構」(含 diffuse,寬)；phylo-v2=「有無清楚離散群」(要 gap+子群 null,嚴) → 對 subclone 離散 lineage 問題 phylo-v2 判準更對 + 多階層標籤/對齊。silhouette 應淘汰。

## 5.7 第二個缺陷：FM1 復發（單離群吃群位）+ v3.1 修正（用戶肉眼複核抓出）

用戶在 v2 工作站肉眼複核，旗標 4 個「v2 判 1 群但眼睛看到 2-4 群」位點 → 親查 trace 發現 **phylo 遞迴第二個 bug**：

- **bug**：phylo 只測**樹根的兩個子節點**。當根分裂是 `(tiny_outlier, n-1)`（`min<MINSZ`）→ `split_real` 直接 return False → **整棵樹判 1 群，遞迴從不往下**。但平衡切割就在下一層（`phylo_recursion_trace.py` 證實）：

| 位點 | root 分裂 | 平衡 k 最大兩群 r | Jaccard | V_allele | v2 | v3.1 |
|---|---|---|---|---|---|---|
| chr22_26939195 | (1,153) 拒 | **1.944** | 0.999 | 0.96 | 1 群 | **2 群** ✅ |
| chr22_30454004 | (1,94) 拒 | k4 | 0.96 | 0.63 | 1 群 | **4 群** ✅ |
| chr21_40426852 | (1,36) 拒 | **2.095** | 0.992 | 1.0 | 1 群 | **2 群** ✅ |
| chr20_42981498 | (2,78) 拒 | 1.306 (<null95 1.34) | 0.85 | 0.99 | 1 群 | 1 群（真 borderline）|

> 這是用戶最初提的 **FM1「單一離群吃群位」在 phylo 遞迴層復發**（cluster_redesign/4-gate 用 quarantine 修過，phylo 沒做）。**非 v1→v2 回歸**（v1/v2 都犯，與 double-dip 修正正交）。

**v3.1 修法**（`phylo_v3_validate.py`）：`rec` 進 node 時 **descend 剝離 (tiny,big) caterpillar 單離群至「兩子群皆≥MINSZ」的平衡節點**再測 split；且 **quarantine 物有所值** — 只在分裂成立時才把剝離的標離群，分裂不成立則整群保留（避免侵蝕 caterpillar 子樹）。

**驗證**（全 L1，`phylo_v3_validate.json`）：
- (A) 用戶 4 旗標：3 個救回（chr22_26939195→2、chr22_30454004→4、chr21_40426852→2，全 germline 對齊 cis-ASM），chr20_42981498 真 borderline（r=1.306 < 噪音 null95 1.34）仍 1 群 = 可辯護的嚴格拒切。
- (B) 全 30 位點：**v2 多群 6 → v3 多群 10，救回 6，減少 0**（無侵蝕），v3 平均 1.47 群。
- (C) 純噪音 FP（C=76）**v3 = 0.0%（n=20/40/80）** → 修法不過切、specificity 完整保留。

**驗證方法學（回答「如何驗證方法與判別」）**：一個候選切割是否真實 = **clusterboot Jaccard 穩定（Hennig >0.75）+ 對齊 CramérV（hp/allele）+ 平衡 between/within r ≥ 噪音 null95 + 全域純噪音 FP**，四者一致才採信。v3.1 救回的 6 個 Jaccard 0.85-0.999 + V 0.63-1.0 + r 1.3-2.1 + 噪音 0% → 確認真實。

> **subclone 結論不變**：v3.1 救回的多是 germline 軸（cis-ASM）；chr22_30454004 的 3 個 hp2-1/ALT 子群是 within-allele 多群（subclone 候選 vs cis-ASM 子結構待肉眼+normal cis-control 判）。單樣本仍 ⭐2-3 characterization。

## 5.8 門檻數字的數據解釋（敏感度掃描 `phylo_sensitivity.py` → `phylo_sensitivity.json`）

回答「各門檻能否用數據解釋」：跑掃描看「換門檻會怎樣」，區分 data-derived（資料推導）vs convention（約定，但落安全平台）。

**掃描結果（全 L1，C=76, noise TRIALS=80）**：

| 掃描 | 變動 | 噪音 FP | 30 位點多群 | 4 旗標 |
|---|---|---|---|---|
| **SEP_MIN** | 1.0 / 1.1 / 1.2 / 1.3 | **全 0%** | **全 10** | [2,4,1,2] 一致 |
| | 1.4 / 1.5 | 0% | 9 / 6 | 開始砍真群 |
| **null 百分位** | 90 / 95 / 99 | 全 0% | 12 / 10 / 8 | 控敏感度 |
| **RNULL** | 12→80（chr22_26939195）| — | **5 次全 2 群（穩）** | — |
| | 12→80（chr20_42981498）| — | **1/2/3 翻動（RNULL=80 仍翻）** | — |
| **MINSZ** | 3 / 4 / 5 | 全 0% | 全 10 | — |

**逐門檻數據判定**：

| 門檻 | 值 | 類型 | 數據解釋 |
|---|---|---|---|
| **null95（per-subgroup 重分群）** | 95pct | ★ **data-derived（方法核心）** | 真正控 FP 的閘 — SEP_MIN 1.0-1.5 全 0% 證明 FP 來自此 null 非 SEP_MIN。百分位控敏感度（90→12 / 99→8 多群，都 0% FP）|
| **SEP_MIN** | 1.3 | convention（**安全平台**）| 1.0-1.3 全給 10 多群+0% FP（冗餘於 null95）；1.4↑ 才砍真群。可降到 1.0 無害，1.3 留安全邊際 |
| **MINSZ** | 3 | convention（穩健）| 3/4/5 結果不變（都 10 多群）→ 非關鍵；3=最小可估群統計量 |
| **RNULL** | 40 | data-informed | 穩定位點 RNULL=12 已收斂；**borderline 位點任何 RNULL 都翻** → 40 對穩定位點足夠，borderline 需另以「多 seed 一致性」標 ambiguous（非加 RNULL）|
| Jaccard >0.75 | — | literature（Hennig 2007 clusterboot）| 業界 cluster 穩定性標準 |
| CramérV≥0.3 / Cochran e≥5 | — | literature（效應量小門檻 / χ² 可靠性）| 標準統計約定 |
| BERNOULLI weight 2\|p-0.5\| | — | data（資訊量加權）| CpG 越偏離 0.5 越具判別力 |

**🔴 兩個由數據導出的方法洞察**：
1. **null95 是 specificity 的唯一來源**（SEP_MIN/MINSZ 是非關鍵約定，落安全平台）→ 對外解釋方法時應強調 per-subgroup null，而非 SEP_MIN=1.3。
2. **borderline 位點（chr20_42981498）verdict 跨 seed 不穩（RNULL=80 仍 1/2/3 翻）= 資料在決策邊界，非參數問題**。正解 = **加 instability flag（跑 k seeds，verdict 翻動 → 標 ambiguous）**，誠實標「資料無法決定」而非強判。此為 v3 後的下一個改進點。

## 6. 檔案

- 缺陷校準：`scripts/phylo_noise_calibration.py` → 純噪音 FP 表
- 前後對照：`scripts/phylo_compare_v1v2.py` → `phylo_v1v2_compare.json`
- v2 最終：`scripts/phylo_v2_final.py` → `phylo_v2_final.json`（含對齊）
- v1（缺陷版，保留供對照）：`scripts/phylo_groups.py` → `phylo_groups.json` + `figs_phylo/`
- v1│v2 並排圖：`scripts/phylo_render_v1v2.py` → `figs_phylo_compare/cmp_*.png`（4 變動位點）
- 四方法對照：`scripts/phylo_methods_compare.py` → `methods_compare_30loci.json`
- FM1 復發診斷+trace：`scripts/phylo_undersplit_diagnostic.py` / `phylo_recursion_trace.py`
- v3.1 修法+驗證：`scripts/phylo_v3_validate.py` → `phylo_v3_validate.json`；render `scripts/phylo_v3_render_all.py` → `figs_phylo_v3/`
- v3 肉眼工作站：`phylo_v3_verify_workstation.standalone.html`（repo 外 `/big7_disk/liaoyoyo2001/ism_html_review/20260622_phylo_v3_verify_workstation.html`）
- 門檻敏感度：`scripts/phylo_sensitivity.py` → `phylo_sensitivity.json`（SEP_MIN/null-pct/RNULL/MINSZ 掃描）
