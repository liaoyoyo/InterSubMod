<!--
建立時間: 2026-05-30
作者: Claude (Opus 4.8, 1M) — 綜合 Agent (盤點 + 可行性 + 三鏡頭對抗整合)
報告類型: strategic-feasibility-assessment
任務類型: D handoff / 決策支援 (給主 agent 帶回用戶討論)
框架: Pyramid (BLUF → verdict → 支撐) + Cynefin (域分類 gate)
樣本 scope: ⚠ partial — 證據鏈核心為單樣本 HCC1395 paired_full；discrimination 軸跨 TO/paired/V6 多樣本
驗證原則: feedback_existing_artifacts_must_verify — 每個 verdict 數字本 session grep source + cross-check ledger
證據分級: 每 claim 標 L1 (親自讀碼/重算) / L2 (從碼推論) / L3 (推測)
read-only: 本評估無改 C++、無跑 >1min compute、未碰 genome_survey_v2/
-->

# ISM「甲基 × haplotype × somatic」全鏈策略可行性評估

> **一句話 BLUF**：InterSubMod 能達成「meth × HP × somatic 全鏈 observe→discover→verify→analyze→output」這個高價值目標 ——
> **但僅在 A-characterization 框架下成立（CAN-with-known-fixes），在 B-discrimination 框架下是橫跨 3 pipeline、6+ cycle 收斂的 concluded NEGATIVE，達不到且不值得再投入。**
> 區分這兩個框架是整份評估的軸心；混淆它們是過去最大的策略風險。

---

## 1. 最終目標精確定義 — 把模糊的「全鏈」收斂成兩個正交框架

使用者高價值目標的字面陳述：「用 InterSubMod 完成甲基(methylation) × haplotype × somatic 三者關聯的全鏈 observe → discover → verify → analyze → output」。這個陳述**語意上同時涵蓋兩個科學目標，但它們的可行性結論相反**，必須先拆開：

| 框架 | 精確目標陳述 | 對應五目標 | 成功標準 |
|------|------------|-----------|---------|
| **A — characterization (表徵)** | 「哪些 somatic 位點周邊有 allele-specific methylation (ASM)？方向 (hypo/hyper)、空間分佈、是否抗多重檢定？」 | **目標 1**（per-CpG 甲基多標籤關聯評分）+ 部分目標 2/3 | 可解釋性成功：能說明甲基模式與 HP/SNV 關聯，case study 可重現 + 跨樣本方向一致 |
| **B — discrimination (判別)** | 「能否用 methylation/ASM 訊號把 somatic TP 與 FP 分開、過濾 FP / 救回 TP，提升 caller F1？」 | **目標 5**（整合 evidence panel 提升 F1） | F1 提升成功：precision/recall 可量化改善 + 跨樣本泛化 |

> **依據**（L1, `project_research_vision_five_goals.md:13-17`）：專案五目標白紙黑字把「目標 1 = per-CpG 甲基位點多標籤關聯性評分」與「目標 5 = 整合 evidence panel 提升 F1」列為**兩個獨立目標**，且明文「兩者都算成功，不能混用」。A 對應目標 1，B 對應目標 5。

**全鏈五段定義**（兩框架共用語彙）：
- **observe** — 同時讀到 per-read 的 5mC/5hmC call + HP:Z haplotype tag + somatic sub-haplotype (HP1-1/HP2-1) + allele type (ref/alt)
- **discover** — 計算 per-CpG β、跑 dual-axis (HP-axis non-LOH + ALLELE-axis 含 LOH) 對比、找位點聚集
- **verify** — Wilcoxon paired + Cohen's d + negative control + 雙獨立實作 cross-check + 多重檢定校正
- **analyze** — 對單位點 / 全基因組下科學結論（A: ASM 表徵；B: TP/FP 判別力）
- **output** — TSV / IGV session / HTML / 圖 / ledger entry

---

## 2. 現有結論誠實盤點 — POSITIVE / NEGATIVE table

> 全部 cross-check 自 `research/autoresearch/evidence_ledger.jsonl`（53 entries, append-only, L1 本 session 親讀）+ MEMORY「Concluded」區。

### 2.1 A-characterization 軸 — POSITIVE 但 modest + 待修

| 結論 | verdict | tier | 證據 | 級別 |
|------|---------|------|------|------|
| BRCA2/ZAR1L chr13:32,315,128 ASM 真實可重現抗 Bonferroni | `positive_but_modest` | ⭐3 | ledger `ASM_TSG_PROMOTER_ZAR1L_BRCA2`；HP1 vs HP1-1 Wilcoxon p=1.8e-10 (Bonferroni-robust) | L1 |
| 全基因組 ASM 廣泛存在 | POSITIVE | — | 22 autosomes / 18,149 records / Bonferroni 90 / 45 spatial clusters | L1 (ledger + report §5) |
| ASM 32-66% somatic 位點存在 | POSITIVE | — | MEMORY `project_snv_methylation_association` | L1 |

**A 軸的四個誠實 caveat**（ledger entry `caveat` 欄逐字）：(1) magnitude 方法依賴 −0.05 (MSA REF-anchored) ~ −0.12 (pysam read-anchored 高估)；(2) **TSG promoter enrichment = 0.00×**（582 sig 中 0 個落 TSG promoter，但 n=3 underpowered，非「證明無關」）；(3) 方向 **hyper(11) > hypo(8)**，BRCA2 是少數 hypo；(4) **single-sample HCC1395，跨樣本一致性 NOT verified**。

### 2.2 B-discrimination 軸 — 🧱 concluded NEGATIVE 硬牆（不可繞、不該 reopen）

| 結論 | verdict | 證據 | 級別 |
|------|---------|------|------|
| TO germline FP 過濾 G1-G7 全 AUC<0.64 | NO-GO | MEMORY `germline_fp_nogo`；ledger H001/H002/H003/H004/H005/H007 全 `negative` | L1 |
| methylation 對 F1 是 vestigial，全靠 caller_af proxy | NEGATIVE | ledger `H016+H_M1a+H_A1: ism_vestigial_caller_af_dominant` | L1 |
| FP filter 無 sample-level 泛化 (LOSO ΔF1≈−0.00004) | NEGATIVE | ledger `LOSO_sample_level_validation: negative_filter_no_sample_level_generalization` | L1 |
| Phase2 cross-sample 反向 | NEGATIVE | ledger `H_PHASE2_C1_GLOBAL_FP_FILTER: negative_cross_sample_with_caveat` | L1 |
| ASM 廣泛但 FP(65.6%) > TP(35.9%) 重疊大、AUC<0.60 無法區分 | NEGATIVE | MEMORY `project_snv_methylation_association` + `beyond_auc_exhaustion` | L1 |

> **硬牆性質**：B-NEGATIVE 不是單一 cycle 的失敗，而是**橫跨 TO / paired / V6 三條 pipeline、≥6 個獨立 cycle 全部收斂**到 methylation 無 TP/FP 判別力。符合 MEMORY `feedback_productive_failure_reopen_threshold` 的 already-concluded 情境。

### 2.3 看似給 B 留活路、實則站不住的兩個誘惑（誠實打掉）

- **誘惑 1**：B4-S4 entry `POSITIVE_CALLER_DRIVEN`、LR AUC 0.717。**打掉**：verdict 名稱本身寫明 `CALLER_DRIVEN` — 主導特徵是 AF/AlleleDelta，**不是 methylation**（L1）。
- **誘惑 2**：HCC1395 global FP filter ΔF1=+0.022 `positive_strong_single_sample`。**打掉**：同 hypothesis_id 的後續 entry `negative_cross_sample` + `ism_vestigial_caller_af_dominant` 證明那 +0.022 是 sample-level circularity、methylation 貢獻 vestigial（L1）。

---

## 3. 研究流程 + 驗證方法健全度 — sound vs 需修

> 依據三份本 session L1 覆核的 audit：`agentA_method_validity_audit.md`（method validity）/ `msa_audit_synthesis.md`（tool fitness）/ `20260529_..._verification_01.md`（既有 artifact 驗證，自掛 🟠 PROVISIONAL banner）。

### 3.1 Sound（可信、方向結論成立）

- **OBSERVE 鏈兩端齊備（L1 本 session 親讀碼）**：上游 tag 生產 = longphase-to-mod `HaplotagProcess.cpp` getVote/HAPLOTYPE1_1（MEMORY `reference_longphase_getvote_source` 行號）；下游 tag 消費 = MSA `MethylHaploExtractor.cpp:339-373` —— 我直接 grep 確認該 .cpp **明確 strcmp 比對 `"1-1"` / `"2-1"` 字串**（行 348-351）+ `origHP` fallback（行 368-373），HP:Z 五類 sub-haplotype 不 collapse。header `:46/62/71/82` 有 extractMethylation / extractHaplotypeTag / determineAlleleType / classifyMethylationState。
- **方向性結論 robust**：BRCA2 somatic hypomethylated 在兩條獨立 pipeline（pysam script03 + MSA C++）× 三種 collapse 口徑（raw / MAX / 退化 MIN 除外）全部同號（agentA §2b, L1）。
- **多重檢定校正做對**：未校正 p<0.05=3,630 → Bonferroni 90 → Bonferroni+|Δβ|≥0.1=19，誠實縮水（L1, report §5.1）。
- **cross-check 雙實作方向一致**：Python Δβ−0.122 vs MSA −0.055，方向一致（L1, synthesis §2）。

### 3.2 需修（影響 magnitude / 顯著性 / 泛化，但**不影響方向**）

| # | 問題 | 嚴重度 | 影響什麼 | 修法 | 是否改 C++ |
|---|------|--------|---------|------|-----------|
| (a) | **Level1 5mC/5hmC 重複發射 bug** — 同 read 同 CpG 同 strand 把 p 與 1−p 各記一列，2,395/2,395 互補對；n 灌水 2.4× (10,632→4,362)、p=6e-11 被人為膨脹、絕對 β 隨 collapse 漂 0.0~0.34「不可信」 | **高 (L1)** | n / p / 絕對 β（**不影響 Δβ 符號**） | dedup-by-MAX collapse per-read-per-CpG（script18 已含 workaround；C++ export 端未修） | C++ 修屬 Hard Gate；Python workaround 免改 |
| (b) | **script03 CpG-anchoring 瑕疵** — `C+m?` unspecified mode 把 27.5% 非 CpG 'm' call 當 CpG，膨脹 magnitude | 中 (L1) | magnitude（script03 −0.122 = 2.3× 高估） | 棄用 script03 為定量，僅留 cross-check；報告數字一律用 MSA REF-anchored | 免改 |
| (c) | **negative control 統計力不足** — n=5 中 2 個 n_paired=0，有效 null 僅 3 點；且用 script03 瑕疵口徑與 test 口徑不一致 | 中 (L1) | empirical null 上界不可信 → 「BRCA2 是 outlier」近乎無統計意義 | n≥50-100 matched（TVAF + 局部 CpG density + coverage matching），修正口徑重做 | 免改 |
| (d) | **ALLELE-axis confound 未控** — collider/selection bias + ALT-read mapping bias；18_dual_axis_pivot.py 未過濾 MAPQ、未控 read length；ALLELE-axis 與 HP-axis 非等價對比群 | 中 (L1+L2) | ALLELE-axis 只能 exploratory，不可單獨下定量 ASM | 加 MAPQ/read-length 過濾；non-LOH 先驗兩軸群體重疊 + 符號一致 | 免改 (script) |
| (e) | **MSA RAM/disk dead param** — `--max-ram-gb` 不生效，全基因組 Level1 估 50-200 GB、累積 vector 可爆 100-300 GB | infra | 全基因組 batch 可行性 | per-chr 刪檔 + --max-read-depth 500 + --window 1000（script19 已 workaround） | C++ chunked streaming 未實作（可選） |
| (f) | **跨樣本泛化未做** — 全 A 軸建在單樣本 HCC1395 paired_full；canonical/ 下他樣本 *tagged.bam 不在同 layout，MEMORY 載 V3F/V5 4 樣本 BAM 不存在 | scope | A 軸不可升 ⭐4 | fan-out 既有 multi-sample-consistency + parallel-benchmark 到 4+ 樣本（known-cost 工程，非科學未知） | 免改 |

> **健全度總評**：A 軸的方法**方向上 sound、定量上待修**。所有「需修」項中，(a)(b)(c)(d) 影響的是 magnitude / 顯著性強度 / null 嚴謹度，**沒有一個動搖「有沒有 ASM、朝哪個方向」這個 characterization 核心問句**。(f) 是真正的科學未知（跨樣本一致性）。**沒有任何 hard MISSING 阻斷全鏈。**

---

## 4. 可行性裁決 — per sub-goal × framing (CAN / PARTIAL / CANNOT)

> 下表為三鏡頭對抗（optimist steelman + skeptic steelman）整合**後**的修正版。修正點在表下說明。

| sub-goal | framing | verdict | 理由（壓縮） | 主要 blocker | 證據級別 |
|----------|---------|---------|------------|-------------|---------|
| observe | A | **CAN** | tag 生產(HaplotagProcess) + 消費(MethylHaploExtractor:339-373) 兩端 L1 齊備 | — | L1 |
| observe | B | **CAN** | 同上；B 不是卡在 observe | — | L1 |
| discover | A | **CAN-for-HP-axis / PARTIAL-for-ALLELE** | HP-axis(non-LOH) 設計正確 agentA §3 PASS；genome survey 22 autosomes 已實跑 = CAN 強證；ALLELE-axis 是 exploratory extension | ALLELE-axis confound (d) | L1 |
| discover | B | **CANNOT** | FP>TP 重疊 confound-driven | concluded confound | L1 |
| verify | A | **PARTIAL (direction-CAN, magnitude/null-PARTIAL)** | 方向 robust = CAN；絕對 β/p/null 受 (a)(c) 拖累 = PARTIAL | Level1 dup (a) + null n=3 (c) | L1 |
| verify | B | **CANNOT** | LOSO circularity + cross-sample 反向 | concluded NEGATIVE | L1 |
| analyze | A | **CAN-with-known-fixes** | BRCA2 ⭐3 end-to-end existence proof + genome survey；天花板 ⭐3 因 single-sample + 小效應 | tier3 ceiling = 待跨樣本 (f) | L1 |
| analyze | B | **CANNOT** | ≥6 cycle 跨 3 pipeline AUC<0.60 | concluded NEGATIVE | L1 |
| output | A | **CAN (載體) / PARTIAL (可信定量內容)** | TSV/IGV/HTML/ledger 載體 EXISTS 且已產出；但繼承 (a) 的 β/p 數字須 dedup 重出正式版 | Level1 dup (a) 繼承 | L1 |
| output | B | **CANNOT** | 內容 NEGATIVE，只剩 postmortem 價值 | no positive artifact | L1 |

### 4.1 對抗驗證後的關鍵修正（相對於初版 verdict）

1. **analyze/A 的 evidence「7-sample POSITIVE」是 conflation，已移除**（skeptic major catch，verdict_survives=false on that cell）。報告 §16 + ledger caveat(4) 白紙黑字 ASM 全程**單樣本 HCC1395**。專案裡的 7-sample 工作是 (i) V6 phasing/LR-filter（與 per-CpG ASM 無關）或 (ii) `project_snv_methylation_association` 那條「FP>TP、7/7 一致但 AUC<0.60 無法區分」—— **後者反而是 B-NEGATIVE 的證據，不是 A 跨樣本 POSITIVE 的證據**。修正：analyze/A 維持 CAN-with-known-fixes，但 rationale 改為「BRCA2 單樣本 end-to-end existence proof」，**不引用 7-sample**。
2. **A 軸整體從「PARTIAL leaning CAN」改寫為「CAN-with-known-fixes」**（optimist + skeptic 都接受）：科學上可行、僅工程/口徑未竟；兩個可修項（dedup-by-MAX + ALLELE confound 控制）皆 no-code 或 ≤50 行，方向結論不受影響。
3. **verify/A、output/A 明標 blocker 為 magnitude/口徑層級**，避免讀者誤以為 characterization 本身不成立。direction = CAN，quantitative = PARTIAL until dedup。
4. **discover/A 區分 HP-axis(non-LOH)=CAN vs ALLELE-axis=exploratory-extension**：把 ALLELE-axis confound 從「核心 discover blocker」降為「延伸雄心的未竟」。
5. **B 軸全列 CANNOT 維持，且證據比初版更強**：skeptic 鏡頭也攻不破，所有 optimist 翻案嘗試（S4 / +0.022 / 低純度潛力）被 ledger L1 全數擊倒。

---

## 5. 核心 Verdict — InterSubMod 能否達成此高價值目標？

**能 —— 但只在 characterization 框架，且需補一個科學未知 + 修兩個工程口徑。對誰高價值要說清楚。**

InterSubMod 已經把「甲基 × haplotype × somatic」全鏈在 A-characterization 框架下**走通過一次完整的 end-to-end**：observe (C++ tag 兩端 L1 齊備) → discover (genome survey 22 autosomes 實跑、HP-axis 設計正確) → verify (方向 robust、抗 Bonferroni、cross-check 一致) → analyze (BRCA2 ⭐3 existence proof) → output (TSV/IGV/HTML/ledger 全產出)。剩下的不是「可行性未知」，而是兩類**known-cost** 工作：(1) **科學未知**僅一項 —— 跨樣本一致性（A 軸全在單樣本 HCC1395，這是升 ⭐3→⭐4 的唯一實質門檻，但用既有 multi-sample-consistency + parallel-benchmark 基礎設施 fan-out 即可，是工程不是科學賭注）；(2) **工程口徑**兩項 —— Level1 dedup-by-MAX + ALLELE-axis confound 控制，皆 no-code 或 ≤50 行，且**只動 magnitude/顯著性，不動方向**。

反之，B-discrimination 框架是**確定的 CANNOT**：用 methylation/ASM 區分 TP/FP、過濾 FP、提升 caller F1，已被橫跨 TO/paired/V6 三條 pipeline、≥6 個獨立 cycle 反覆釘死（germline_fp_nogo AUC<0.64 / LOSO ΔF1≈−0.00004 / ism_vestigial caller_af-dominant / cross-sample 反向）。這不是「還沒做好」，是**做過很多次、機制上 FP 與 TP 在 methylation 空間重疊**的 concluded NEGATIVE，符合 reopen-3-條件不該 reopen 的情境。

**對誰高價值（Cynefin 域分類）**：
- **A-characterization → 高價值，且是 Complicated 域**（因果可分析、有 expert-knowable 答案）：對「InterSubMod 作為可解釋 epigenetic evidence 整合框架」這個論文定位高價值 —— 它能回答「哪些 somatic 位點有 ASM、方向、空間分佈」這類**機制描述**問題，是論文的 mechanism/case-study 章節素材。**對 PI / 論文 reviewer 高價值**。
- **B-discrimination → 零增量價值，是 Clear 域**（答案已知 = NO）：對「提升 caller F1」這個工程目標無價值，不應再投入算力。**對任何想用甲基化做 somatic caller 的人 = 明確勸退**。

---

## 6. 建議路徑 — 若可行怎麼走 / 硬牆怎麼繞 / reopen 條件

### 6.1 A-characterization 的 known-cost 推進路徑（建議走）

```
[P0 修口徑]  (1) Level1 dedup-by-MAX 落地到正式 pivot (Python 免改 C++；或 /cpp-change 修 export = Hard Gate)
             (2) script03 正式棄用為定量，報告數字一律 MSA REF-anchored
             (3) negative control 重做 n≥50-100 matched (TVAF+CpG density+coverage)
   │  ← 全部修「magnitude/null」，方向結論預期不變
   ▼
[P1 跨樣本]  用既有 /multi-sample-consistency + parallel-benchmark fan-out BRCA2 + genome-survey
             流程到 4+ 樣本 (COLO829 / H2009 / HCC1954 ...)
   │  ← 唯一真實科學未知；方向一致性 ≥4/N 樣本 → 升 ⭐3→⭐4
   │  ⚠ 前置: 確認他樣本 *tagged.bam (帶 HP1-1 + 5mC) 是否存在於同 layout (MEMORY 載 V3F/V5 4 樣本 BAM 不存在 — 需先盤點)
   ▼
[P2 analyze] ALLELE-axis 控 confound (MAPQ/read-length filter) → non-LOH 先驗兩軸符號一致再採信
             LOH 區 ALLELE-axis 結論需獨立 phasing 驗證
   ▼
[P3 output]  正式 magnitude HTML/TSV/ledger (移除 🟠 PROVISIONAL banner)
             + spatial clustering 從 ad-hoc ~50 行固化為 pipeline component
```

### 6.2 B-discrimination 硬牆 — 不繞，明確封閉

不建議任何「繞過」嘗試。methylation 對 TP/FP 的判別力是機制性 NEGATIVE（FP 與 TP 在 methylation 空間重疊），不是參數調校問題。任何宣稱「用甲基化提升 F1」的新提案，**必須先通過 reopen gate**。

### 6.3 reopen NEGATIVE 需要的 C1/C2/C3 條件（MEMORY `productive_failure_reopen_threshold`）

要重開 B-discrimination，至少需滿足一項：
- **C1 新數據** — 例如更高純度/深度的新樣本，或新平台（非 ONT）的 orthogonal methylation；現有 HCC1395 + 6 樣本不算。
- **C2 新方法** — 例如非「平均 β 比較」的新表徵（如 read-level methylation haplotype block entropy），且須先過 `/auc-confound-guard`（within-group OLS + AF-bin + permutation）。
- **C3 新前置** — 例如 caller 升級後 FP 結構改變、或限定到某個特定 zone（如低純度 subclone）而非全域。MEMORY `project_read_level_germline_fp` 留了「低純度樣本有潛力」一線，但須 C1 新樣本支撐才成立。

---

## 7. 要與使用者確認的關鍵決策點（3-5 個）

> 以下為帶回主 agent → 用戶討論的 gate；每個都標影響度/信心度供暫停判斷。

1. **【框架確認 — 最關鍵】** 使用者要的「全鏈」是 **A-characterization（表徵哪些位點有 ASM，論文 mechanism 素材）** 還是 **B-discrimination（用甲基提升 caller F1）**？
   → 若 A：往 §6.1 走（可行）。若 B：這是 concluded NEGATIVE，需先談 §6.3 reopen 條件，否則建議不投入。（影響: 高 / 信心: 高 — 兩框架結論相反，走錯方向浪費整個 cycle）

2. **【跨樣本 BAM 盤點 gate】** A 軸升 ⭐4 的唯一科學門檻是跨樣本，但 MEMORY 載 V3F/V5 4 樣本 tagged BAM 不存在、canonical/ 下他樣本 *tagged.bam 不在同 layout。**是否先花 <30min 盤點哪些樣本有帶 HP1-1+5mC 的 tagged BAM**，再決定 fan-out 範圍？（影響: 高 / 信心: 中 — 若 BAM 不存在，跨樣本路徑需先補 tag 生產 pipeline，成本量級不同）

3. **【Level1 dedup 修法選擇】** 修 5mC/5hmC 重複發射 bug 有兩條路：(a) Python pivot 端 dedup-by-MAX（免改 C++，快，但下游 C++ Level2/3 輸出仍錯）；(b) `/cpp-change` 修 C++ export（根治，但 **C++ commit = Hard Gate** 需編譯驗證 + 用戶確認）。**走哪條？**（影響: 中 / 信心: 高 — 影響正式 magnitude 數字可信度與 output 範圍）

4. **【negative control 重做的算力授權】** null 從 n=5(有效3) 擴到 n≥50-100 matched 需要跑額外 MSA。**是否授權這次（中等）長計算**，還是先用現有 3 點 null 出 PROVISIONAL 版？（影響: 中 / 信心: 高 — 不重做則「BRCA2 是 outlier」這個 verify 結論不可信，但背景有 genome_survey_v2 在跑，需排程避撞）

5. **【B 軸是否正式封箱】** 是否同意把 B-discrimination 在 ledger / CURRENT_FOCUS 正式標為 **concluded NEGATIVE（reopen 需 C1/C2/C3）**，停止任何「甲基提升 F1」的探索，把算力全投 A-characterization？（影響: 高 / 信心: 高 — 這是策略聚焦決定，避免重複踩 6+ cycle 的同一面牆）

---

## 8. Provenance（本 session L1 親自覆核）

| claim | source | 級別 |
|-------|--------|------|
| HP:Z 1-1/2-1 sub-haplotype 解析 | `/big8_disk/.../src/msa/core/MethylHaploExtractor.cpp:339-373`（本 session grep strcmp "1-1"/"2-1"+origHP）| **L1 親讀碼** |
| extractMethylation/AlleleType/HaploTag | `MethylHaploExtractor.h:46,62,71,82` | L1 親讀 header |
| BRCA2 ASM verdict + 4 caveat | ledger `ASM_TSG_PROMOTER_ZAR1L_BRCA2`（本 session json.tool dump）| L1 |
| B-discrimination NEGATIVE 硬牆 | ledger H001-H007 + LOSO + ism_vestigial + cross-sample（本 session grep）| L1 |
| Level1 dup bug / script03 anchoring / null n=3 / ALLELE confound | `agentA_method_validity_audit.md §1-§4`（本 session 全讀）| L1 (agentA 親測) |
| MSA tool fitness / 3 調整免改 C++ | `msa_audit_synthesis.md §0-§5`（本 session 全讀）| L1 (synthesis) |
| 五目標 A/B 對應 | `project_research_vision_five_goals.md:13-17`（⚠ memory 48 天，目標定義穩定不受時效影響）| L1 |
| 報告自掛 PROVISIONAL | `20260529_..._verification_01.md:14`（本 session 全讀）| L1 |
| longphase HaplotagProcess 上游 tag | MEMORY `reference_longphase_getvote_source`（行號未本 session 重驗，grep 該路徑下檔案未即時命中）| L2 |

**未驗證 / TBD**：跨樣本 tagged BAM 實際存在性（決策點 2 待盤點）；longphase-to-mod HaplotagProcess.cpp 確切路徑（grep 未即時命中，依 MEMORY L2）；dedup 修正後正式 magnitude（待 P0）。

---

*業界框架對照：Pyramid Principle (Minto, BLUF 先給結論) + Cynefin domain classification (Snowden — A=Complicated/B=Clear) + Productive Failure reopen gate (Kapur, C1/C2/C3)。*
