<!--
建立時間: 2026-04-06 23:45
目標: 每張 slide 的佈局規格 + 口頭講稿
處理範圍: 2026-03-31 ~ 2026-04-06
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/04/20260406_研究週報_LOH雙定義與特徵探索全面關閉/pptx_config.json
-->

# Slide 佈局規格 + 口頭講稿

---

## Slide 1 — 封面（Dark BG）

**佈局**：`add_title_slide` — 全頁深色背景 #1E2A44，右側 accent panel
- 左：kicker badge + 標題(Georgia 28pt) + 副標題 + summary line
- 右：導讀 bullets（3 條）+ 報告日期/作者

**講稿**：
> 各位好，這是 03/31 到 04/06 的研究週報。本週主題是 LOH 雙定義與特徵探索全面關閉。我們完成了 166 張圖表、16 項判定，最終得出一個決定性的戰略結論。

---

## Slide 2 — 本週重點結論

**佈局**：`add_key_findings_slide` — 4 cards (2×2 grid)
- 左上(綠)：SEQC2 Jaccard=0.928
- 右上(紅)：TO 特徵全面關閉
- 左下(紫)：甲基化 NEGATIVE + TO NO-GO
- 右下(藍)：戰略轉向 characterization

**講稿**：
> 先看四個重點。第一，我們找到 FDA 金標準 SEQC2，驗證 LOH.bed Jaccard=0.928，定義可信。第二，LOH 區域 10/10 filter 策略全失敗，Non-LOH 也沒有特徵超過 0.58 門檻。第三，甲基化三維度全被 confound 否決，TO 60 多個 VCF 特徵也全部小於 0.64。第四，我們的結論是 ISM 的價值在 characterization 而非 filtering。

---

## Slide 3 — LOH 雙定義與四象限

**佈局**：`add_framework_slide` — 左 3-layer stack + 右 2×2 quadrant
- 左：三層分工（LOH 定義驗證 → ISM 影響 → 全面關閉）
- 右：四象限定義（Q1-Q4 各佔比）

**講稿**：
> 本週的分析框架是圍繞 LOH 的兩種定義建立的。ISM 有兩種 LOH 定義——LOH.bed 來自 LongPhase-TO，HP_Ratio 來自 ISM 自己的計算。我們建立了四象限：Q1 兩者都判 LOH 佔 26.7%，Q2 只有 HP 判的佔 15.2%，Q3 只有 LOH.bed 判的近乎為零，Q4 都不判的佔 58.1%。

---

## Slide 4 — 目標與完成度

**佈局**：`add_goal_status_slide` — 表格 + 3 summary cards
- 表格：6 目標（goal/subgoal/completion/conclusion）
- 底部 3 cards：LOH.bed 可信 / TO 全面關閉 / 轉向 characterization

**講稿**：
> 本週有 6 個主要目標全部完成。LOH 雙定義三波分析、甲基化三維度、TO VCF 特徵、Self-phasing 因果鏈、QS 修正、ASM 驗證。重點看底部三張卡：定義可信、特徵全面關閉、方向轉向。

---

## Slide 5 — 時間軸

**佈局**：`create_timeline_slide` — 4 時間卡片 + 連結線
- 03-31：甲基化 NEGATIVE
- 04-01~02：TO NO-GO + Self-phasing
- 04-03~05：LOH Wave 1+2 + Deep Study
- 04-06：Wave 3 + 全面關閉

**講稿**：
> 時間線。3/31 完成甲基化三個假說否決。4/1-2 完成 TO 特徵 NO-GO 和 self-phasing 因果鏈驗證。4/3-5 建立 LOH 四象限和 SEQC2 驗證。4/6 完成 Wave 3 Non-LOH 分析，正式全面關閉。

---

## Slide 6 — SEQC2 外部驗證

**佈局**：`add_seqc2_validation_slide` — 3-column + bottom callout
- 左：SEQC2 描述（FDA, 320 entries, 多平台）
- 中：LongPhase-TO 描述（ONT, 方法學零重疊）
- 右：驗證指標（Jaccard, Sensitivity, Precision, F1, HP_Ratio AUC）
- 底：結論 callout（深色背景）

**講稿**：
> 這是本週最重要的驗證。左邊是 SEQC2——FDA 主導的多平台多中心共識。右邊是 LongPhase-TO——我們用 ONT 的 phased genotype。兩者方法學完全獨立。結果 Jaccard=0.928，Sensitivity 和 Precision 都在 0.96 以上。這代表 LOH.bed 作為分層基礎是可靠的。

---

## Slide 7 — ISM 在 LOH 區域全面失效

**佈局**：`add_ism_loh_impact_slide` — 4 stat cards (2×2) + root cause + fix
- 左上：PERMANOVA Valid Rate 5-6%
- 右上：LOH 甲基化 AUC ~0.50
- 左下：Filter 10/10 FAIL
- 右下：LOH FP rate 0.239
- 底左：根因說明
- 底右：已修正（QS mode-aware, 綠色 badge）

**講稿**：
> 但 LOH 可信不代表 ISM 在 LOH 裡有用。四個關鍵數字：PERMANOVA valid rate 只有 5-6%，因為 94% 的位點只有一條 haplotype 有 reads。甲基化 AUC 約 0.50 跟隨機一樣。10 種 filter 全失敗。而且 LOH 裡 FP rate 反而比 Non-LOH 低——它是 TP-enriched 的。我們已經實作了 QS mode-aware 修正。

---

## Slide 8 — Non-LOH + 多特徵組合

**佈局**：`add_non_loh_closure_slide` — table + voting + Simpson's Paradox + callout
- 左：Non-LOH 特徵 AUC 表
- 右上：Voting AUC=0.577
- 右下：cnLOH Simpson's Paradox
- 底：結論 callout

**講稿**：
> 那 Non-LOH 呢？佔了 58.1%，如果可以區分就能集中資源。但看左邊表格，corrected AUC 沒有超過 0.58 的。三特徵 Voting 最高只到 0.577。右下角特別注意 Simpson's Paradox：cnLOH 整體看似 0.587 超過門檻，但 per-sample mean 是 0.50，被 H2009 主導。問題是全域性的，不是 LOH 特異的。

---

## Slide 9 — 甲基化 NEGATIVE + TO NO-GO + ASM

**佈局**：`add_negative_results_slide` — methylation rows + TO card + ASM dark card + summary
- 左：O11-O13 三行（before→after, cause）
- 右上：TO NO-GO 統計
- 右下(dark)：ASM 唯一亮點（32-66%）
- 底：結論

**講稿**：
> 同時進行的其他方向。甲基化三個假說——O11 heterogeneity 表面 0.845 是 n_reads confound，修正後 0.530。O12 是 AF confound。O13 是 shared read confound。TO VCF 特徵 60 多個全部低於 0.64。但 ASM 是亮點——5 種方法驗證 32-66% 的 SNV 位點有 allele-specific methylation。ISM 的 PERMANOVA 是唯一能量化這個現象的工具。

---

## Slide 10 — 14 結論總表

**佈局**：`add_conclusion_table_slide` — full-width table
- 14 行（C1-C14）× 5 欄（#, 結論, 判決, 穩定度, 本週變動）
- NEW 標記的用粗體

**講稿**：
> 總結來看，14 個結論。本週新增 C5 LOH.bed 可信、C6 不可作 filter。三個新 NEGATIVE：甲基化、TO 特徵、Non-LOH。C12 多特徵組合不可行。C14 確認 LOH 的價值在 stratification。6 個 CONFIRMED、5 個 NEGATIVE、1 個 POSITIVE、1 個 CONDITIONAL。

---

## Slide 11 — 下一步（Dark BG）

**佈局**：`add_next_steps_slide` — 3 columns
- 左(淺色)：目前尚缺
- 中(白色)：下一步行動
- 右(accent)：不建議做的

**講稿**：
> 下一步。左邊是尚缺的——最重要的是 PON-only phasing 後的真實分佈。中間是行動——Phase 2A 的 PON-only 全量重跑和 normal reference baseline。右邊是明確不建議的——不要再嘗試 LOH filter、不要用甲基化過濾 TO、不要混合 paired/TO 結論。

---

## Slide 12 — 複查路徑

**佈局**：`add_paths_slide` — 3 columns
- 左：核心文件路徑
- 中：數據來源路徑
- 右：重生指令

**講稿**：
> 最後是複查路徑。來源週報、LOH 分析報告、研究全景索引都在這裡。如果需要重新生成這份簡報，右邊有完整指令。
