<!--
建立時間: 2026-07-05
類型: 觀察分析結果 (跨樣本延伸)
狀態: in_progress
主軸: Subclonal reconstruction — 甲基 read 分群 vs CN 關聯（是否只偵測 CNV / 獨立）
data_sources: /big7_disk/liaoyoyo2001/big7_disk_output/canonical/{H1437,H2009,HCC1954}/paired_*/*_complete_matrix/intersubmod_tp/significance_summary.csv, /big7_disk/liaoyoyo2001/cnv_sv_work/{H1437,H2009,HCC1954}/savana_wgs/*_smcnbed.bed, docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/20260621_kism_vs_cn_perread_result_01.md
build_branch: research/subclonal-reconstruction-202606
關聯 memory: project_ont_cnv_sv_subclone_verification_feasibility, project_savana_cna_only_pipeline_and_colo829_purity
-->

# 甲基 read 分群 vs CN — 是否只偵測 CNV / 獨立？（跨樣本延伸）

> **版本** v1.0 (2026-07-05) · 證據 tier：**HCC1395 = L1（⭐3）**；跨樣本 proxy = **L2-L3**（粗 metric）
> 腳本 scratchpad `xsample_methyl_structure_vs_cn.py`（結果 json 同層）；權威前結果見 §data_sources

## 0. 一句話結論

**甲基 read 分群（非監督）≠ 只偵測 CNV** — 大致**獨立於 CN**（HCC1395 ρ≈0；跨樣本 2/3 |r|<0.1）。但需區分兩種「甲基」信號：**① 非監督 read 分群 = CN-獨立**（有獨立內容、來源未定）；**② allele-ASM = CN/LOH 依賴**（需兩 allele，LOH 區測不了 → r=−0.57~−0.80）。→ 甲基有獨立於 CN 的內容，但**bulk 下無法證明是 subclone**（來源未定，⭐3 bounded-auxiliary）。

## 1. §7.1 Pre-registration（confirmatory 延伸）

- **H**：甲基 cluster 結構率獨立於 CN 狀態（LOH≈neutral≈gain），複現 HCC1395。
- **否證**：≥2/3 樣本結構率強依賴 CN（|r|>0.2）。
- **結果**：**否證未觸發**（僅 H1437 1/3 樣本 |r|~0.22>0.2）→ 獨立性未被推翻，但 H1437 有弱 CN 依賴。

## 2. 權威前結果（HCC1395, L1, ⭐3 — 不重造）

（`docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/`）
- k_ISM vs k_CN **無相關**（Spearman ρ=−0.038，效應≈0）。
- 只 **9.5%** somatic 區有非監督甲基 cluster 結構（PERMANOVA-gated；TP 9.5% vs FP 3.8%=2.5×）。
- 結構區僅 **0.64% 對齊 germline-HP（CN 可解釋）**、het 區 **85% unaligned**、LOH 區 72% candidate_beyond_CN。
- 結構「存在與否」與 coverage/CN **全 |ρ|<0.05**（既非 coverage 也非 CN 假象）。
- → **甲基分群非 CN 鏡子（有獨立內容），但結構大多不對齊任何遺傳軸 = 來源未定 epigenetic 異質**。

## 3. 跨樣本延伸（H1437/H2009/HCC1954，L2-L3）

> 🔴 **caveat**：新樣本 ISM run（20260315）**早於非監督 cluster-PERMANOVA 落地（6月）** → `ClusterPermanovaValid` 全 false → **無法用 HCC1395 同 metric**。改用舊 run 有的 **`PassedGating`（CramersV+GlobalP heuristic gate，over-permissive）** 當粗 proxy。COLO829/HCC1937 排除（purity mis-fit）。H1437/HCC1954=paired_pileup、H2009=paired_full（mode 不一）。

### 3a. 非監督甲基異質（gated）率 by CN state + 相關

| 樣本 | gated% | LOH% | neutral% | gain% | r(LOH) | r(gain) | r(cov) |
|---|---|---|---|---|---|---|---|
| H1437 | 35.3 | 11.9 | 42.2 | 40.2 | **−0.225** | +0.219 | +0.109 |
| H2009 | 12.2 | 7.6 | 13.5 | 14.3 | −0.081 | +0.051 | +0.041 |
| HCC1954 | 29.4 | 34.8 | 31.0 | 28.8 | +0.032 | −0.031 | −0.066 |

→ **2/3 樣本（H2009/HCC1954）甲基異質率近乎獨立於 CN（|r|<0.1）**，複現 HCC1395。**H1437 有弱 CN 依賴**（gain/neutral 42% > LOH 12%，r~0.22）— 但未達「≥2 樣本」否證門檻。coverage 相關全 <0.11（非 coverage 假象）。

### 3b. label-level allele-PERMANOVA（ASM）率 by CN state

| 樣本 | asm% | LOH% | neutral% | gain% | r_asm(LOH) |
|---|---|---|---|---|---|
| H1437 | 86.2 | 26.2 | 97.0 | 99.2 | **−0.803** |
| H2009 | 84.0 | 47.0 | 99.7 | 98.2 | **−0.574** |
| HCC1954 | 99.2 | 99.8 | 100.0 | 99.6 | +0.017 |

→ **ASM（allele 層）強烈 CN/LOH 依賴**（H1437/H2009 LOH 區 ASM 率暴跌、r=−0.57~−0.80）**機制清楚**：ASM 測兩 allele 甲基差，**LOH 區一 allele 丟失 → 測不了 → ASM-valid 率崩**。→ allele-ASM 信號**內在被 CN/allele-availability 綁定**，非純 epigenetic。（HCC1954 purity 0.66 殘留 normal allele 故 LOH 區仍可測 → r≈0，反印證機制。）

## 4. 逐一回答四向問題

| 問 | 答 | 依據 |
|---|---|---|
| **甲基有偵測到 read 分群差異？** | ✅ 有（但少）：HCC1395 9.5% 區有結構（TP 2.5×FP）；跨樣本 gated 12-35% | §2, §3a |
| **是否只偵測 CNV？** | ❌ **不是**：非監督 read 分群大致獨立於 CN（HCC1395 ρ≈0；2/3 樣本 |r|<0.1）| §2, §3a |
| **其實不只？** | ✅ 有 CN-獨立內容，**但大多不對齊任何遺傳軸 = 來源未定** | §2 |
| **甲基單獨、關聯不高？** | ✅ **非監督分群 = 大致單獨（低 CN 關聯）**；🔴 **但 allele-ASM 例外 = 強 CN/LOH 依賴** | §3a vs §3b |

## 5. §4 DAG + 質疑（confound）

```
CN/LOH ──強(r=−0.8)──> allele-ASM-valid   [機制:LOH 丟 allele → ASM 測不了]
CN     ──弱/無(ρ≈0)──> 非監督甲基分群       [HCC1395 L1;跨樣本 2/3 確認]
非監督甲基分群 ──?──> subclone            [🔴 來源未定:cis-ASM/技術/真 subclone 不可分]
```
- 🔴 **「獨立於 CN」≠「是 subclone」**：CN-獨立的結構仍可能是 cis-ASM（單分子順式）、技術、或真 subclone — **bulk 不可分**（⭐3 天花板）。
- 跨樣本是 **L2-L3 粗 proxy**（heuristic gate 非 PERMANOVA cluster；mode 不一；n=3）— 弱於 HCC1395 L1。
- H1437 弱 CN 依賴（r~0.22）需留意（gain 區異質較高 — 可能擴增區 read 多、gate 較易過）。

## 6. 發想（next）

1. **證明 CN-獨立結構是 subclone**：normal-cis control（同結構是否在 matched normal 也有？有=非 tumor subclone）+ multi-sSNV CCF + single-cell（⭐4 需求）。
2. **正確跨樣本複現**：對新樣本**重跑 ISM 含非監督 cluster-PERMANOVA**（20260315 run 缺）→ 才能用 HCC1395 同 metric 得 9.5%/independence。
3. **ASM=CN-driven 是乾淨負對照**：r=−0.8 跨樣本確認「ASM 内在被 CN 綁」→ 支持「ASM ≠ subclone」（allele-ASM 不可當 subclone 證據）。

## 7. Verdict（§2 證據分級）

- **L1（HCC1395 ⭐3）**：非監督甲基分群獨立於 CN、來源未定 bounded-auxiliary。
- **L2-L3（跨樣本粗 proxy）**：2/3 複現獨立性 + 乾淨確認 ASM 是 CN/LOH-driven。
- **總裁決**：**甲基 read 分群不是只偵測 CNV（大致獨立於 CN，有獨立內容）；但獨立內容來源未定、bulk 不能證明 subclone；而 allele-ASM 信號則確實被 CN/LOH 混淆。** 甲基定位 = **有界輔助（characterize 非 detect subclone）**，與既有主軸裁決一致。
