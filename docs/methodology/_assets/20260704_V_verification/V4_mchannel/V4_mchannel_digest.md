# V4 驗證數字報告 — Model A 群組斯坦納樹 recurrence-required 的 m-通道拆分

**軌**: V4 (m-channel / recurrence-required 拆分) · **樣本**: HCC1395 (⭐3 單樣本上限) · **狀態**: **completed**
**日期**: 2026-07-04 · **sub-agent 驗證**

---

## 裁決一句話

Model A 的獨立多重度通道 **m_f(x) 已存在、已落地 (LANDED)**,且其 CN 來源是**外部正交** (SEQC2 ground-truth + SAVANA WGS),**非自共現 B 導出** — 通過 formal_statement §4/§8 的非循環判準。舊 118 incompatible 已按 Model A 完整重判;本軌額外用 **SAVANA allele-specific CN 落地了 50 個 LOH-unresolved 的拆分** (⭐3 候選升級),使全 118 區三分乾淨。

---

## 1. 舊 118 incompatible 的 Model A 重判 (已驗證,已落地)

來源: `topology_per_region.json` (07-04 重生) + `20260704_incompatible_reclassification_hidden_node_finding_01.md` (狀態 LANDED, 7 樣本重跑, 回歸乾淨)。

| 分流 | 數 | 判準 | 處置 |
|---|---|---|---|
| **救回 (隱藏祖先深分支樹)** | **30** | nic==0 & no cycle | `build_hidden_node_tree` → A_determined **23** + A_ambiguous **7** |
| **recurrence-required** | **70** | nic>0 & cn≠gain & no cycle | 送 m-通道拆 (見 §2) |
| **維持 incompatible** | **18** | has_cycle 或 (nic>0 & cn==gain) | 真 artifact,不建樹 |

- 30 救回已在資料中確認: `topology_type.branched(需隱藏祖先)=23` + `branched(順序未定,需隱藏祖先)=7`。
- 18 incompatible cross-tab (實測): has_cycle=True **12** (cycle_cause: CN-gain-multiplicity **9** + other-pairwise-cycle **3**); cn 分布 **15 gain + 3 loh**。
- ⚠ 對抗驗證 (wf_520ddd16) 揪出並修復 1 個 builder bug (H1437 chr20 homoplasy 假樹);HCC1395 的 30 救回全乾淨,不受影響。

## 2. m-通道拆分 70 recurrence-required (已落地,SEQC2 整數 CN)

`m_channel_split()` @ `topology_analysis.py:76-91`。**CN 來源 = SEQC2 ground-truth 整數拷貝數 BED** (`/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_{gain,loss}_cn.bed`)。判準:m_f = major-allele CN ceiling。

| verdict | 數 | 判準 | Model A 動作 |
|---|---|---|---|
| `recurrence_artifact(m>1;CN-amp)` | **11** | 整數 CN≥3 (gain) | 多重度 artifact **(棄)** — 與舊 Model B 結論一致 |
| `recurrence_candidate(m=1)` | **9** | neutral CN2 / loss CN≥1 | **候選真 recurrence (留)** — Model B 漏掉的新結構 |
| `recurrence_LOH_unresolved` | **50** | copy-neutral LOH total=2, allele-specific CN 無 | m∈{1,2} 未定 (VAF L3 軟旗標) |

- 9 個 m=1 候選: 8 neutral (cn_int=2) + 1 loss (chr4:18469698, cn_int=1)。範圍 n_sSNV=2~39。
- 11 個 m>1 artifact: 整數 CN 3~5 (chr17:82593283 CN=5 最高)。
- 50 LOH 的 VAF 軟旗標: 35 likely_artifact(高VAF>0.7) / 15 likely_recurrence(低VAF)。

## 3. ⭐3 候選升級 — SAVANA allele-specific CN 落地 50 個 LOH-unresolved

**問題**: 50 個 copy-neutral LOH 因「allele-specific CN 無」只能給 VAF 軟旗標。**但 SAVANA WGS 早已產出 allele-specific CN,只是未接線** (`HCC1395_smcnbed.bed` 缺;整合腳本 `rerun_cn_integration_from_savana.sh` 未對此 recurrence 拆分跑)。

**本軌落地**: 讀 `HCC1395_segmented_absolute_copy_number.tsv` (SAVANA, purity ρ=0.96, ploidy=2.79),對 50 區用標準多重度估計 `m̂ = VAF·(C_tot + 2(1-ρ)/ρ)`,cap 至 major-allele CN。

| 結果 | 數 |
|---|---|
| SAVANA 段對映覆蓋 | **50/50** |
| 確認 copy-neutral (total=2, minor=0) | **49/50** (1 個 total=3) — 證實 m-通道對 LOH 假設的 total=2 正確 |
| allele-specific verdict: artifact(m>1) | **35** |
| allele-specific verdict: candidate(m=1) | **15** |
| 與 VAF-only 軟旗標一致率 | **100%** (35→35, 15→15) |

**解讀 (誠實)**: SAVANA 升級**不翻轉**計數 — 因 49/50 為 total=2,m̂ 對 VAF 單調,硬門檻 (m≥2 ⟺ VAF≥~0.72) 與軟門檻 (0.7) 幾乎重合。SAVANA 的**真實貢獻 = 把 total=2 backbone 從「假設」變成「正交量測證實」**,並提供有原則的多重度估計 → 這 15 個低 VAF LOH 可從「未解 (軟旗標)」**升級為 allele-specific-CN-confirmed 的候選真 recurrence (m=1)**,退役 `LOH_unresolved` 桶。

## 4. V4 綜合結果 (SAVANA 升級後,118 三分)

| 類別 | 數 | 組成 |
|---|---|---|
| **候選真 recurrence (m=1,Model B 漏掉的新結構)** | **24** | 9 (SEQC2 neutral/loss) + 15 (SAVANA LOH copy-neutral) |
| **多重度 artifact (m>1,棄,與 Model B 一致)** | **64** | 11 (SEQC2 CN-amp) + 35 (SAVANA LOH 高VAF) + 18 (incompatible: gain-mult/cycle) |
| **救回隱藏祖先樹 (A_determined/ambiguous)** | **30** | 23 + 7 |
| **合計** | **118** | ✓ |

## 5. 非循環性 (formal_statement §4/§8 硬需求) — 通過

- **SEQC2 CN**: 外部 benchmark truth BED,與共現 B 完全正交。✓
- **SAVANA CN**: 由 read-depth + germline het-SNP BAF 呼叫,非由 somatic sSNV 共現 B。正交於 B。✓
- **多重度內的 VAF**: per-site 邊際等位基因分數,**非**共現 B (B = 跨位點 AA/RA/AR/RR 連鎖)。在正交 CN backbone 內定位多重度合法 → ∂(tree set)/∂φ=∂/∂d=0 保證不破。✓

**結論**: Model A **可安全落地** — 獨立 m 通道存在且非循環。formal_statement §12.4 的開放問題「118 中有多少翻案成真 recurrence」得答: **24 個候選真 recurrence** (弱主張,⭐3 單樣本,需 orthogonal 佐證且非甲基確認)。

## 6. group-Steiner 形式化措辭升級建議

1. **§4 的「無獨立 m → over-build」風險已解除**,但應在 spec 標明 m 通道的**兩級證據**: (L-hard) 整數 CN≥3 → artifact 硬判; (L-soft) copy-neutral LOH → allele-specific CN 確認 backbone + VAF 定多重度 (弱,需標 ⭐3 候選)。
2. §7 解法 1「directed Steiner tree / arborescence」應升級為 **group-Steiner**: recurrence-required 區的每個「翻兩次的位元」對應一組候選 Steiner 端點 (m=1 真 recurrence 節點 vs m>1 artifact 座標),m 通道即 group-membership 的外部裁決器。建議措辭: 「terminals 分組,每組至多選一個代表 (m≤1);m>1 組整組視為 multiplicity artifact 剔除」。
3. §10 輸出格式的 `m_check_status` 應增列 `allele_specific_cn_confirmed` 值,區分 SEQC2-categorical (loh 只知類別) vs SAVANA-allele-specific (知 total/minor 整數)。

---

### 資產盤點 (m-通道可用性)

| 資產 | 路徑 | 狀態 |
|---|---|---|
| SEQC2 整數 CN (gain/loss) | `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_{gain,loss}_cn.bed` | ✓ 已接線 (m_channel_split) |
| SEQC2 gain/loss/loh categorical | `.../ngs_benchmark_cnvs_gain_loss_loh.bed` | ✓ 已接線 (sm_region_integration cn 欄) |
| SAVANA WGS allele-specific CN | `/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv` | ✓ 存在, ⚠ **未接線** → 本軌落地於 50 LOH |
| SAVANA fitted purity/ploidy | 同上目錄 `_fitted_purity_ploidy.tsv` | ρ=0.96, ploidy=2.79 |
| 自導出 cn_proxy (BAF/depth) | `cn_proxy_annotation.py` | 輔助代理,**不進 m 通道** (誠實隔離) |

### 必讀清單交叉核對 — `backbone_resolution.json`
任務指定先讀之 `backbone_resolution.json` 為 **pre-C2-fix 快照** (06-29/07-01 era):`A1_incompatible.n=12` (舊 12→118 修正前值)、`cause_dist: CN-gain 9 + artifact/其他 3`、全 12 標 `fixable_artifact`「標 artifact 不建樹」。**無 recurrence / m_channel / cycle / Steiner 概念** — 早於 Model A。已被 `topology_per_region.json` (07-04 post-fix) 取代 (本軌所用)。**不改任何 V4 數字**,但**獨立佐證 CN-gain 主軸**:舊 12 區 9/12=75% CN-gain,與現行 18 incompatible 的 cycle_cause 9/12=75% CN-gain-multiplicity 一致。

### 交付檔案 (host: `docs/methodology/_assets/20260704_V_verification/V4_mchannel/`)
- `v4_msplit_summary.json` — 全 118 三分 + 非循環判準 + 資產盤點
- `v4_118_regions.json` — 88 區 (70 recurrence + 18 incompatible) 逐區 CN/determinacy/m_channel
- `v4_loh_savana_msplit.json` — 50 LOH 的 SAVANA allele-specific m-拆分逐區明細
