---
title: per-CpG 多軸歸屬（定位差異甲基位點 + 軸歸屬）— BRCA2 pilot 驗證 + ISM 輸出方向
date: 2026-06-20
status: in_progress
tier: 3
sample: HCC1395（單樣本；BRCA2 已知 ASM 驗證案；⭐3 characterization 非 subclone）
scope: pilot（單位點 chr13:32315128 BRCA2；零重跑 Python + 1-region binary）
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/brca2_idea1_out.json,docs/methodology/_assets/20260618_subcluster_pilot/brca2_multiaxis_out.json,research/tsg_promoter_asm_reviewer/output/step4_ism_results.json
build_commit: e57a7ab
observation_standard: true
---

# per-CpG 多軸歸屬 — BRCA2 pilot 驗證

> **想法（Idea 1）**：用「以甲基距離切群 + 標籤對齊」的分類結果，回頭做**附近 CpG 的篩選與定位 + 軸歸屬** — 找出哪些 CpG 與標籤高度相關、與**哪個軸**（HP / carrier / allele）相關、強度與數量，輸出 {CpG 座標 + 軸 + Δβ + FDR}，輔助 ASM 分析與下游研究。
> **驗證標的**：BRCA2 chr13:32315128（已 verified ASM，HP1 germline vs HP1-1 somatic Δβ=−0.1219, wilcoxon p=6.1e-11, somatic HYPOMETHYLATED）。

## 1. 雙重交叉驗證（兩條獨立路都指向 carrier 軸）

| 路徑 | 結果 |
|---|---|
| 我的 decisionflow（tumor） | 此位點 = **⑤確認真結構、aligned=carrier**；perm carrier F=15.4 / hp F=4.0 / allele F=8.7（p=0.01） |
| 既有真值（step4 stats HP1_vs_HP1-1） | region Δβ=**−0.1219**, wilcoxon p=6.1e-11, "somatic HYPOMETHYLATED" |

→ 無監督切群 + per-CpG 定位 + 既有 ASM 真值，三者一致指向 carrier（germline vs somatic-haplotype）軸。

## 2. 想法1 定位結果（兩種讀法，結論一致）

**(a) 從 cpg_records（tsg-reviewer，germline HP1 vs somatic HP1-1，每組≥5 read）**

| 指標 | 值 |
|---|---|
| 可測 CpG | 197 / 412 |
| 顯著 CpG（BH-FDR q<0.05） | **25** |
| 方向符 −0.122（somatic hypo, Δβ<0） | **25/25 = 100%** |
| 顯著 \|Δβ\| 中位 | 0.622（範圍 0.36–0.95） |
| 定位 dist_to_var 中位 | **−416**，範圍 **[−600, −169]**（變異上游 ~430bp 焦點低甲基塊） |
| pool 全體 Δβ 均值 | −0.102（對照 region −0.122） |

**(b) 從 my-binary read 矩陣（tumor reads；somatic 在此 BAM 標 HP2-1）**

per-CpG 三軸各測 BH-FDR：**hp 4 / carrier 18 / allele 5** 顯著（carrier 主導）。

> 兩讀法 carrier 顯著數 25 vs 18 差異來源：不同 tagged BAM（somatic = HP1-1〔tsg〕vs HP2-1〔ClairS-pileup binary〕）+ 不同 read 集 + 聚合 count vs per-read binarized；同量級、同結論。

## 3. 多軸歸屬 — 一個 CpG 可同時對齊多軸（你的直覺成立）

**軸組合（21 顯著 CpG 聯集，tumor）**

| 組合 | CpG 數 | 意義 |
|---|---:|---|
| carrier only | 14 | somatic 單倍型整體低甲基（不分 allele）= 單倍型層 ASM |
| carrier + allele | 2 | 多軸：carrier 與 allele 都顯著 |
| hp + carrier + allele | 2 | 三軸全中（最 confounded） |
| hp only | 2 | 物理單倍型差 |
| allele only | 1 | 僅 REF/ALT |

**🔑 carrier × allele 共線診斷（tumor）**：germline = **全 REF（30）**、somatic = ALT(13) + **REF(13)** 混。
- germline⟹REF（完全共線）→ 但 somatic **同時有 ALT 和 REF**，故 carrier 與 allele **可部分分離**。
- **14 個 carrier-only CpG** = somatic-REF 與 somatic-ALT read **都**低甲基 → 訊號來自 **somatic 單倍型本身**（haplotype-level ASM），**非** ALT allele 特異。
- 4 個多軸（carrier+allele / 三軸）= allele 也參與，最 read-identity-confounded（adversarial review 已標 allele 軸最 confounded）。

## 4. 圖（figs_brca2/，gitignored 可重生）

- `track_carrier_dbeta.png`：per-CpG carrier Δβ vs dist — 全 ±5000 窗只有 −600~−169 段有顯著（紅，全負），其餘 Δβ≈0 → 訊號**焦點化非全域**。紫菱=兼 allele 多軸（落在最強效 CpG）。
- `dualpanel_brca2.png`：read×CpG 甲基（RdBu_r）+ read×read 距離（magma），排序 carrier→allele→β，側欄 carrier/HP/ALT/T-N/strand — somatic（carrier 紅）reads 在焦點窗清楚低甲基塊。

## 5. 驗證表

| 數字 | 值 | 來源 | 重算 | L |
|---|---|---|---|---|
| 可測 CpG（cpg_records） | 197 | brca2_idea1_out.json.n_tested | 412 中兩組各≥5 | L1 |
| 顯著 CpG | 25 | .n_sig | BH-FDR q<0.05 | L1 |
| 方向符 −0.122 | 25/25 | .n_neg | Δβ<0 計數 | L1 |
| carrier 顯著（my-binary） | 18 | brca2_multiaxis_out.json.sig.carrier | per-read Fisher+BH | L1 |
| 多軸聯集 | 21 | sig 三軸 union | hp4+carrier18+allele5 去重 | L1 |
| carrier×allele cross-tab | germ-REF30/som-ALT13/som-REF13 | brca2_multiaxis.py stdout | reads.tsv 計數 | L1 |
| region truth Δβ | −0.1219 | step4_ism_results.json[0].stats | 既有 verified | L1 |

## 6. 問題方向 + ISM 能否輸出這些座標/分析/觀察

**方向**：從「位點層分類（decision flow 5 態）」→「**CpG 層歸屬/註解**」。每位點輸出一張 CpG 表 {座標, dist_to_var, 各軸 Δβ+q, 多軸 flag, primary 軸}，= 一條 **ASM 註解 track**，餵下游 + 論文「甲基 characterize」。

**ISM 現況**：
- 區域層 `significance_summary.csv` 已輸出（但只 region 級 3 數）。
- **`PerCpgAsm.cpp:274-319` 已做逐 CpG Fisher + BH-FDR，但只測 HP-family 一軸 + per-CpG 向量未落檔**（只 aggregate 成 region 3 數）→ **機制已存在，只是丟掉了 per-CpG 向量**。

**兩條輸出路**：
1. **Python 後處理層（建議先做）**：對既有 `methylation.csv` + `reads.tsv` 矩陣（29,754 region 全有，零 C++ 風險）跑本 pilot 的多軸歸屬 → 輸出 per-CpG 表。本 pilot 即此法的單點實證。
2. **ISM C++ 原生輸出（若廣泛有用再做）**：擴 `PerCpgAsm` → (a) 三軸（hp/carrier/allele）、(b) emit per-CpG 向量檔 {coord, axis, Δβ, p, q, n_axes_sig}、(c) 多軸 flag。屬 cpp-change（gated）。

**建議**：先 Python 層全基因組跑（證明價值 + 數量是否夠多），再決定是否 productionize 進 ISM。輸出須內建誠實邊界 caveat（見 §7）。

## 7. 誠實邊界

- **單樣本 cis-ASM characterization，非 subclone**：carrier 軸 = germline vs somatic-haplotype；BRCA2 單 TSG 無 linkage、normal unphased → 只能 ASM 定位。
- **allele 軸最 confounded**（read-identity）；多軸 CpG 的 allele 歸屬須 normal cis-control 才能確認。
- **carrier×allele 部分共線**：能分離靠 somatic-REF read（此例 13 個），數量少時不可分。
- 顯著 CpG 數依「每組≥N read」門檻；本 pilot N=5。

## 8. 重生

```bash
cd InterSubMod/.claude/worktrees/ism-review-infra
# 1-region binary（取 read 矩陣）：mini-vcf chr13:32315128 → output/_brca2_region/
python3 docs/methodology/_assets/20260618_subcluster_pilot/scripts/brca2_idea1_percpg.py       # cpg_records 版
python3 docs/methodology/_assets/20260618_subcluster_pilot/scripts/brca2_multiaxis_attribution.py # 三軸版
python3 docs/methodology/_assets/20260618_subcluster_pilot/scripts/brca2_plots.py              # track + dual-panel
```
