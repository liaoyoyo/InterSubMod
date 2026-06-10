<!--
建立時間: 2026-06-11
狀態: in_progress (大規模方法學驗證 — 816 loci; Phase B 延伸)
報告類型: largescale_methodology_verification
受眾: PI · 方法學決策(Fisher 敏感度 / normal-HP 比對是否加)
framework: Verdict-Pyramid + 3 問題(Q1 率vs結構 / Q2 normal-HP 價值 / Q3 Fisher 敏感度)
data_sources:
  - runs/largescale_summary.json + runs/largescale_perlocus.tsv (816 loci 實跑, 已讀回+sanity)
  - research/tsg_promoter_asm_reviewer/pipeline/cache/level1/*.tsv.gz (816 Level-1 cache)
  - cis_scan_full.json (tier)
provenance_note: 數字自 largescale_summary.json + perlocus.tsv 實跑讀回+sanity 核對; struct_delta 為 Python NHD d_between-d_within proxy(非 ISM binary PERMANOVA, 99-perm p下限0.01); tumor/normal absdbeta 為不同軸(HP1vsHP1-1 somatic / HP1vsHP2 germline)的 magnitude 比較, matched-CpG 為更乾淨指標。單樣本 HCC1395。本檔 Write 與分析不同 batch(§13.0)。
-->
<!-- provenance-verified: 數字引 runs/largescale_summary.json + perlocus.tsv 實跑(816 loci)+ sanity 重算核對; 撰寫與分析分批。 -->

# 11 — 大規模方法學驗證（816 loci）：Fisher 敏感度 + normal-HP 比對價值

> **這份回答用戶 3 個方法學疑問**（大規模、真資料）：① Fisher（率）vs PERMANOVA（結構）—— 「肉眼看得到但率測提取不出」在全基因組成立嗎；② 加 **normal haplotype-tag 比對**（HP1-vs-HP2 germline 基線）值不值得；③ Fisher 敏感度漏多少。

---

## L0 — 三句話結論

1. **「肉眼看得到結構、率測提取不出」是真的，但少數**：293 個結構顯著 loci 中 **25 個（8.5%）率測弱（Fisher frac<5%）**，集中在 **T0/T2 低 tier**。Fisher 抓到 83%。→ 結構（PERMANOVA）是補抓那 8.5% 的必要互補；**不能說 Fisher 沒事，但也沒「大規模盲掉」**。
2. **🔴 normal-HP-tag 比對 = 值得加（強烈建議）**：normal 的 germline 等位基線 **|Δβ(HP1 vs HP2)| median = 0.237，比 tumor somatic 軸 0.181 還大**，且目前**沒被當基線控制** → 現行分析可能把 germline 等位甲基化誤算成 somatic ASM。matched-CpG 顯示 **73% loci 有真 tumor-specific CpG**（tumor 差但 normal 平），normal 參考能**直接提取**它們。
3. **Fisher 敏感度 gap 真實但有界**：17% 結構 loci 率測不強、8.5% 率測弱 → effect-size 門檻 / beta-binomial-with-delta 仍值得試（抓那 8.5–17%）。

---

## L1 — 三問題結果表

| 問題 | 量化（816 loci，442 testable，299 有 normal）| 判定 |
|------|---------------------------------------------|------|
| **Q1 率 vs 結構** | 結構顯著(d_between>d_within, perm p≤0.05) = **293**；其中率測弱(Fisher frac<5%) = **25 (8.5%)**；率測強(≥10%) = 243 (83%) | 🟡「肉眼看得到提取不出」**真但少數**，集中 T0/T2 |
| **Q2 normal-HP 價值** | tumor\|Δβ\|(HP1vsHP1-1) median **0.181**；**normal\|Δβ\|(HP1vsHP2) median 0.237（更大）**；tumor>normal+0.05 = 68(23%)；matched-CpG **209/285(73%) loci 有 tumor-specific CpG**(中位 2/locus) | 🔴 **值得加** — germline 基線大且未控；能提取真差異 |
| **Q3 Fisher 敏感度** | 結構 loci 中率測**不強** 50/293(17%)、率測**弱** 25/293(8.5%) | 🟡 gap 真實有界 → effect-size 門檻可補 |

---

## L2 — 逐問題判讀

### Q1 — 「肉眼看得到結構、率測提取不出」🟡 真但少數
- 293 個 loci 在**讀層結構**（NHD d_between−d_within > 0，99-perm p≤0.05）顯著；其中 **25 個（8.5%）per-CpG Fisher 卻測弱**（顯著 CpG <5%）。
- 這 25 個的 tier = **T2(12) + T0(13)** —— 正是原本 cis-scan 用率/cis 判準**降權或漏掉**的低 tier 位點。
- **意義**：用戶的觀察成立 —— 確有「結構看得到、率測提取不出」的位點（這正是 487-latent 現象的 generalization），但**是少數（~8.5%）**，且集中在低 tier。Fisher 對 83% 的結構 loci 仍抓得到。
- → **結構（PERMANOVA）是補抓這 8.5% 的必要互補**；對這些位點，**換敏感度導向標準**（effect-size 門檻 / beta-binomial-with-delta / 讀層結構）才提取得出。

### Q2 — normal haplotype-tag 比對 🔴 值得加（核心發現）
- **germline 等位基線很大**：normal 的 HP1-vs-HP2（兩條親本 haplotype）median **|Δβ|=0.237**，**比 tumor somatic 軸（HP1 vs HP1-1）的 0.181 還大**。
- **目前沒被當基線控制**：現行 cis-test 只用 normal **HP1**（d_drift），**沒用 normal HP1-vs-HP2 的等位基線** → 在 germline ASM/imprinting 強的位點，tumor 的「等位差異」可能**部分來自 germline 等位甲基化**而非 somatic。
- **matched-CpG（更乾淨）**：同 CpG 上「tumor |Δ|≥0.3 但 normal |Δ|<0.2」= 真 tumor-specific → **209/285 (73%) loci 至少有 1 個**（中位 2/locus）。
- **意義（直接回應用戶目標）**：
  1. **驗證**：肉眼「block 在 somatic tag 甲基化」→ 用 normal HP1-vs-HP2 確認「normal 是否本來就這樣」→ normal 平 = 真 tumor-specific。
  2. **提取**：normal-flat 過濾出每位點 ~2 個真 tumor-specific CpG，把「率測漏的讀層差異」surface 出來。
  3. **濾假**：germline 基線(0.237)≥tumor(0.181) 的位點（77%）提醒：沒有 normal 等位基線，會把 germline ASM 誤算成 somatic。
- → **強烈建議納入 normal-HP1-vs-HP2 基線**（doc 10 #13）。

### Q3 — Fisher 敏感度 gap 🟡 真實有界
- 結構顯著 loci 中：率測**不強**（Fisher frac<10%）= 50/293（17%）；率測**弱**（<5%）= 25/293（8.5%）。
- → 純看 p 的 Fisher 會漏這 8.5–17%（多在覆蓋邊際 / 讀層 block 差異）。**effect-size 門檻（|Δβ|≥閾值並列）或 beta-binomial-with-delta** 值得試來補敏感度（doc 10 #8b）。**這推翻了我先前「Fisher 不需改」的過早結論的一半** —— 假陽性端 OK，敏感度端有 8.5–17% gap。

---

## L3 — 對 doc 10 改進計畫的影響

| doc 10 項 | 大規模驗證後 |
|-----------|------------|
| #13 normal-HP 比對 | 🔴 **升為強建議**（germline 基線 0.237 > tumor 0.181 未控；73% loci 有可提取的 tumor-specific CpG）|
| #8b Fisher 配 delta（敏感度）| 🟡 **保留 OPEN**（8.5–17% 結構 loci 率測漏；effect-size 門檻可補）|
| 🔴1 多重檢定校正 | 不變（正交）|
| ❌ 「Fisher 不需改」| **半推翻**：FP 端不需改，**敏感度端需要**（換標準/加 effect-size）|

## 誠實邊界（caveat）
- `struct_delta` = **Python NHD d_between−d_within proxy + 99-perm**（非 ISM binary 真 PERMANOVA；p 下限 0.01）→ 結構顯著數是近似，production 應用 binary 重跑確認。
- tumor vs normal **absdbeta 是不同軸**（HP1-vs-HP1-1 somatic / HP1-vs-HP2 germline）的 magnitude 比較；**matched-CpG（n_tumor_specific_cpg）是更乾淨的指標**，兩者一致指向「normal 基線大且有可提取的 tumor-specific 訊號」。
- **單樣本 HCC1395**；somatic 軸依賴 HP1-1 tag 品質；442/816 testable（其餘覆蓋/CpG 不足）。
- 這是**方法學方向驗證**，非已證 production uplift；要落地仍須 /methodology-audit + A/B。

## Provenance
- `runs/largescale_summary.json` + `runs/largescale_perlocus.tsv`（816 loci 實跑 `largescale_verify.py`）+ sanity 重算核對。
- tier：`cis_scan_full.json`。輸入：816 個 `pipeline/cache/level1/*.tsv.gz`（含 tumor HP1/HP1-1 + normal HP1/HP2）。
