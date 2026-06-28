<!--
建立時間: 2026-06-28
類型: methodology decision record — sSNV 單分子共現 2×2 分類門檻定案(ε=2%) + somatic/HP3 原則
狀態: in_progress(定義定案;pipeline code 套用待 branch 重跑)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data
provenance: 來源 branch feat/summary-nreadsvalid@5308d9e(pending-merge)凍結資料;三路驗證 workflow wf_f2b070ea-64c
-->

# sSNV 單分子共現分類門檻定案：ε=2%（決策紀錄）

> 框架：ADR（決策紀錄）。每數字 grep-able（凍結 `data/` JSON/TSV）。本決策以**可驗證方法**為主、統計推論為輔（小樣本 + 低 read，統計不可靠）。

## §0 決策（一句話）

**一個 2×2 cell 算「真」當且僅當 `count > coread × 0.02`（ONT 噪聲底線，隨定序深度縮放）。保留最低 1 條（低 coread 時單讀仍算）；高 coread 的單讀（1 ≤ coread×2%，即 coread≥50）判 noise。** 結構數字一律報「現行(≥1) ↔ ε2%」band。

## §1 為何是 ε=2%（方向由 3 路印證 + 值錨 ONT；**非 3 路收斂於同一數**）

> 🔴 對抗稽核校正：原寫「三路皆指 ε=2%」是 overclaim。實況 = **Path A 定方向（ε 對、VAF 錯）、Path B 封頂（別太高）、值 2% 錨在 ONT 錯誤率**。三路不獨立收斂於「2」這個數，2% 是 A（往高推）與 B（往低封）之間、錨在 ONT 噪聲的折衷。

### 路徑 A — FP 比例裁判（外部 SEQC2，僅觀察不進前處理）→ 定**方向**非值
弱判定（單讀）both-FP **57.9%** vs 強判定 30.1%（1.9×）。分離度「丟FP% − 留FP%」：ε=1% +14.8 / ε=2% +17.9 / ε=3% +24.6。⚠ **此序列單調遞增**（留-FP% 56.2→49.0→41.2）→ A 只說「ε 越高丟越多 FP」即**方向對**（vs VAF 調整 **−10.3 backwards** 出局），**不指定 2 這個值**。

### 路徑 B — 最後 region 結構穩定性（封**上限**，不靠外部標籤）
門檻傳播到 region 層重推樹形：full_tree **677→616**（ε2%，−9%）、相鄰 ε 換 shape 僅 **2.5–2.8%/step**（穩定非 knife-edge）。**region 級零單讀檢定（off≥2 flat）full_tree 仍存 525（77%）** → headline full_tree **非單讀假象**（即使全砍單讀）。⚠ 此為**對凍結 edge 的 post-hoc 重推**（非 pipeline 重跑）、113/3627 edge 查不到保守當存活。B 的作用 = 不讓 ε 推太高（保住結構），對 A 的「越高越好」封頂。

### 路徑 C — 塌陷集中度（偽影特徵，印證打對地方）
現行→ε2% 塌陷的 234 個 structured 區 vs 維持的 2,351 區：edge 中位 coread **74 vs 27**、CN-gain **74% vs 60%**。→ 門檻移除「高 coread + CN-gain」偽影嫌疑（dense-cluster/multiplicity），不碰乾淨低-coread 訊號。

**值 2% 的錨**：ONT R10/Dorado substitution 錯誤 ~1–2%，且 pileup `min_base_quality=0`（無品質過濾）推高至 ~2–3% → 取 **2% 為中央、1–3% 為 band**。A 定方向、B 封頂、C 印證、**ONT 錯誤率定值**。

> ⚠ **scope**：此 ε=2% 是 **HCC1395 單樣本 + SEQC2 標籤校準**的值（⭐3 上限）；跨樣本/跨平台需重新校準 ε（列為 TODO，不可逕用於多樣本宣稱）。

## §2 被否決的替代方案

| 方案 | 否決理由（可驗證） |
|---|---|
| 純絕對 flat ≥2 | 太鈍：丟全部單讀(10,270)含低-coread 合理者；生 1,502 sparse |
| VAF 調整（coread×VAF）| 🔴 backwards：FP 分離度 −10.3，丟 TP(高VAF/高覆蓋)留 FP；組合規則把 co_linked 推到 19,299 比 flat≥2 還嚴 |
| binomial / Fisher 檢定 | 小樣本（決定格中位僅 1 條、coread 中位 36）下 underpowered + 需假設 ε，**不可驗證**；改用 count>coread×ε 的算數門檻 |

## §3 最終 band（ε=2% 套用 38,049 pairs）

| config | 現行(≥1) | ε2% |
|---|--:|--:|
| co_linked | 11,750 | **16,048** |
| nested（合） | 13,113 | 11,352 |
| independent | 6,281 | 3,425 |
| mutual_excl | 6,905 | 7,223 |

重分類 5,966（15.7%）；單讀保留(coread<50) 5,159、單讀→noise(coread≥50) 5,111。
逐筆可驗：`data/pairs_eps2_annotated.tsv`（每對附 RR/RA/AR/AA + coread + floor_2pct + 強/弱/noise + orig→eps2）。

## §4 連帶原則定案

- **somatic 定義**：build = ClairS→longphase-S 的 **TP∪FP union**；somatic 由 **normal 比對**定（`is_somatic` normal-VAF）；**SEQC2 TP/FP 只觀察、不進前處理/定義**（程式碼已遵守）。統一 chr17 builder 的 `==0` → 用 normal-based 一致定義 → **chr17 = 3 個 somatic sSNV**。normal-VAF 仍輸出供下游標「疑似 germline 滲漏」。
- **源頭降噪**：pileup `min_base_quality=0` → 改 **≥10**，從源頭降 ε（比事後門檻乾淨、可驗證）。
- **HP3**：HP3('3') 暫不強改；**以 sSNV 共現資訊優先**，僅當 sSNV 證據與 HP3 分群衝突才回頭修。標記 sibling_sameHP **~6.3%（1317/20815 HP3 linked-somatic）** 污染待 sSNV 驗證。

## §5 實作狀態 + 待辦

- **現況**：ε=2% 為**對凍結 per-pair 的 post-hoc 重分類**（驗證/定義用），尚未進 pipeline。
- **待 code-fix（branch 重跑）**：① `classify()` 改 `count > coread×0.02`（對稱、取代不對稱 aa<2/==0）② pileup `min_base_quality≥10` ③ somatic 統一 ④ 輸出加 `coread×VAF` 描述欄 + confidence flag。重跑後所有下游 region/structure 數字以 ε2% band 重報。
- VAF 只當**描述欄**（報 `coread×VAF`「該有幾條」），**不當 cutoff**。

## §6 可驗證產物

`data/`：`pairs_eps2_annotated.tsv`(逐對) · `eps2_final_band.json`(band) · `threshold_comparison.json`(FP 裁判) · `region_threshold_impact.json`(結構穩定/集中) · `weak_pair_observation.json` · `fp_and_vaf_threshold.json` · `sSNV_combination_enumeration.json`。
`scripts/`：`apply_eps2_canonical.py` · `compare_thresholds.py` · `region_threshold_impact.py` · `observe_weak_pairs.py` · `fp_and_vaf_threshold.py` · `enumerate_sSNV_combinations.py`（皆可重跑複算）。
