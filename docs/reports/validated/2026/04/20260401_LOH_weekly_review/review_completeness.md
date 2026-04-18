<!--
建立時間: 2026-04-01 18:00
目標: 審查 LOH 週報 10 份文件的敘述完整性（五環敘述鏈），產出逐文件評分與缺失項目
處理範圍: 00_background.md ~ 09_conclusions_and_actions.md
關聯檔案:
  - docs/reports/validated/2026/04/20260401_LOH_weekly_review/00_background.md
  - docs/reports/validated/2026/04/20260401_LOH_weekly_review/09_conclusions_and_actions.md
-->

# LOH 週報敘述完整性審查報告

**審查日期**：2026-04-01
**審查範圍**：`docs/reports/validated/2026/04/20260401_LOH_weekly_review/` 目錄下 10 份 `.md` 文件
**審查標準**：五環敘述鏈完整性

---

## 五環敘述鏈定義

| 環 | 名稱 | 審查要點 |
|:---:|------|---------|
| 1 | **動機** | 為什麼做這個分析？研究問題是否明確陳述？ |
| 2 | **方法** | 怎麼做？參數定義、公式、操作步驟是否完整？ |
| 3 | **觀察** | 觀察到什麼？數據表格、圖表引用是否充分？ |
| 4 | **推論** | 為什麼這樣解讀？邏輯推論是否清晰？替代解釋是否被考慮/排除？ |
| 5 | **結論** | 結論是否有明確 statement？影響與後續行動是否說明？ |

**評分標準**：
- ✅ 完整：該環節有充分、清晰的敘述
- ⚠️ 部分缺失：該環節存在但有明顯遺漏或不夠深入
- ❌ 缺失：該環節完全未涉及或極度不足

---

## 總覽評分表

| 文件 | 動機 | 方法 | 觀察 | 推論 | 結論 | 整體評語 |
|:-----|:---:|:---:|:---:|:---:|:---:|:---------|
| 00_background.md | ✅ | ✅ | ✅ | ⚠️ | ⚠️ | 背景鋪陳極佳，但作為定位文件缺乏獨立結論 |
| 01_hp_integer_tag_fix.md | ✅ | ✅ | ✅ | ✅ | ✅ | 五環俱全，堪稱範本 |
| 02_loh_evidence_panel_rounds1_4.md | ✅ | ✅ | ✅ | ✅ | ✅ | 四輪邏輯鏈清晰，方法與推論皆詳盡 |
| 03_post_hp_fix_loh_enrichment.md | ✅ | ✅ | ✅ | ✅ | ✅ | 數據完整，推論嚴謹 |
| 04_mechanism_investigation.md | ✅ | ✅ | ✅ | ✅ | ✅ | 逐步推導機制假說，替代解釋處理優良 |
| 05_systematic_observation_O1_O10.md | ✅ | ✅ | ✅ | ⚠️ | ✅ | 觀察極為豐富，但部分觀察缺乏替代解釋討論 |
| 06_methylation_hypothesis_negative.md | ✅ | ✅ | ✅ | ✅ | ✅ | negative result 報告的最佳實踐 |
| 07_qs_mode_aware_change.md | ✅ | ✅ | ⚠️ | ⚠️ | ⚠️ | 動機與方法清晰，但缺乏實際驗證數據 |
| 08_literature_survey.md | ✅ | ⚠️ | ✅ | ✅ | ⚠️ | 文獻覆蓋佳，但缺乏系統性方法論與總結性結論 |
| 09_conclusions_and_actions.md | ✅ | ❌ | ✅ | ⚠️ | ✅ | 整合文件不需方法環，但推論環偏薄 |

---

## 逐文件詳細審查

### 00_background.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | 清楚說明 ISM 的定位、LOH 核心問題、Paired vs TO 差異的重要性。Section 2 明確定義本週兩個核心研究問題。 |
| 方法 | ✅ | 完整描述 ISM 的三層資訊整合方法（甲基化、somatic SNV、haplotype）、距離度量、統計檢驗，以及 LOH-like 操作型定義（HP_Ratio 0.1/0.9 閾值）。 |
| 觀察 | ✅ | 提供 7 個樣本詳細表格、TP/FP 比例圖表引用、數據規模（748,391 rows x 116 columns）、關鍵術語表含具體數據。 |
| 推論 | ⚠️ | Section 4.3 列出 Paired/TO 區分的五個理由，但這些理由主要是**陳述事實**而非**推導為什麼會這樣**。例如「LOH 方向相反」被當作事實呈現，但本文件中未解釋機制。此外，Section 7「待驗證問題」中的五個問題暗示了推論空間，但並未在文中展開推導。 |
| 結論 | ⚠️ | 作為背景文件，未設置獨立的「結論」章節。「待驗證問題」和「認知門檻補充建議」部分地承擔了結論的角色，但缺乏一個明確的 summary statement（如「本文件建立了以下共識基線：...」）。 |

**具體缺失項**：
1. 缺乏 Section 4.3 中「為什麼 FP 組成完全不同」的機制推導（僅陳述 germline vs 技術性誤報）
2. 缺乏總結性結論段落，應明確聲明本文件建立的共識基線

---

### 01_hp_integer_tag_fix.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | Section 1 詳細說明了 bug 的根本原因（整數 vs 字串 HP tag 格式差異），直接後果（30-99% reads 被錯誤歸為 HP0），並明確劃分受影響/不受影響的範圍。 |
| 方法 | ✅ | Section 2 提供完整的修正程式碼（switch-case mapping）、每個 mapping 的語義解釋表格、防護機制說明。修正方式可完全重現。 |
| 觀察 | ✅ | Section 3 提供了 per-sample LOH eligibility 變化表、Tier 分佈變化表、hp_assign_rate 重新量化表、欄位翻轉比例表。7 個樣本均有數據，包含修正前後對比圖。 |
| 推論 | ✅ | Section 6「修正後的關鍵口徑變化」明確推導了 TO LOH 從「FP filter 候選」到「TP enrichment 結構訊號」的定位轉變。Section 7 仍然成立 vs 需改寫的結論分離清晰。 |
| 結論 | ✅ | Section 8「最終判決」提供了邊界清楚的結論：保留什麼、重算什麼、改寫什麼。Section 5 的文件修正分級（A/B/C 級）是額外的結論層次。 |

**無明顯缺失**。此文件是五環敘述鏈的範本。

---

### 02_loh_evidence_panel_rounds1_4.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | 文件開頭的核心問題清晰（「LOH 相關特徵是否能作為 FP 鑑別或過濾指標？」），四輪邏輯鏈圖（R1->R2->R3->R4）說明了每輪的問題演進。每輪均有獨立的「目的」段落。 |
| 方法 | ✅ | 每輪均有完整的方法描述：R1 的 LOH-like 定義與 HP_Ratio 計算公式、R2 的 Tier 分層定義與統計力依據、R3 的 HP0 filter 假設與 2x2 矩陣設計、R4 的 F1 基線修正方法。分析腳本路徑與 workspace 路徑均已記錄。 |
| 觀察 | ✅ | 四輪共引用 10+ 張圖表、20+ 個數據表格，涵蓋全域統計與 per-sample 分解。數據來源均標注至具體 TSV 檔案。 |
| 推論 | ✅ | R1 的 HP_Ratio=0.5 陷阱推論、R2 的 Tier B 反轉 Simpson's paradox 解釋、R3 的 HP0 filter 否定原因分析、R4 的 Tier A vs A+ 方向相反的機制推導（中等深度 vs 古老克隆性 LOH），均有充分邏輯鏈。 |
| 結論 | ✅ | 「四輪研究總結論」提供了 9 項結論一覽表，且有「LOH 的最終定位」提供具體的 LOH_FP_risk_score 公式與 evidence panel 建議。方法論注意事項（sample composition artifact 等）亦有涵蓋。 |

**微小建議**：
1. Round 3 的 LOH+HPMergedSig 7.4x 發現的「生物學解釋」（Section Round 3 第 2 點）可增加對替代解釋的排除討論（已在 R4 部分補足，但 R3 敘述中略顯跳躍）

---

### 03_post_hp_fix_loh_enrichment.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | Section 1 清楚說明 HP bug 影響範圍與為什麼需要重跑，並解釋 Round 1 舊數據為何「偶爾看似接近正確」但不可信。 |
| 方法 | ✅ | Section 2 描述完整的修正與重跑步驟（5 步），數據確認日期和分析腳本路徑齊全。LOH eligible 定義（eff_hp >= 30）清晰。 |
| 觀察 | ✅ | Table 1-4 提供 per-sample LOH eligibility、per-tier enrichment、per-sample enrichment、Tier A 統計顯著性等完整數據。圖表引用充分（enrichment heatmap、concordance 圖）。Section 5 單獨處理 COLO829 離群值。 |
| 推論 | ✅ | Section 6 的 Paired vs TO 方向翻轉分析包含 TP/FP 基礎差異的解釋、per-sample 對比表、方向對比總結。Section 7 推導 LOH penalty 對 TO QS 的反向影響。 |
| 結論 | ✅ | Section 8 以結構化表格呈現 6 項核心結論，每項含強度評級和證據來源。COLO829 標記為「低信心」的建議也屬於結論。 |

**無明顯缺失**。

---

### 04_mechanism_investigation.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | Section 1 以單一核心問題開場（「為什麼相同基因組位點在不同 calling mode 下 LOH enrichment 方向翻轉？」），動機極為聚焦。 |
| 方法 | ✅ | Section 2 描述同位點 concordance 分析方法（288,609 valid loci），Section 4-5 提供 HP balance 四象限分析方法和 AF 分層方法。 |
| 觀察 | ✅ | 同位點 concordance 數據（7 樣本一致 16-52x 偏向）、HP balance 四象限數據、AF 中位數差異數據、VerificationClass 分層數據、低 AF TP 數據，均以表格和圖表呈現。 |
| 推論 | ✅ | **此文件的推論環是 10 份文件中最出色的**。Section 3 以四步邏輯鏈逐步推導機制假說（TO phasing 差異 -> somatic allele HP imbalance -> FP 相對 HP balance 好 -> Paired FP 不同機制），每一步都標示前提、證據、推導，並在 Section 7 明確列出「不能宣稱的」三點。 |
| 結論 | ✅ | Section 7 以結構化表格區分「確認的結論」（5 項，含強度評級）和「不能宣稱的」（3 項）。TP rescue 潛力被適當標記為「初步」而非確定。 |

**無明顯缺失**。此文件的推論環可作為範本推廣。

---

### 05_systematic_observation_O1_O10.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | 開篇段落清楚說明為什麼需要系統性觀察（748K rows x 116 columns 的複雜度）、分拆為 9 個觀察單元的策略，以及觀察等級定義（Level A/B）。每個 O 觀察均有獨立「目的」。 |
| 方法 | ✅ | 每個觀察列出方法摘要（如 O1 對 20 個核心欄位繪製分佈圖、O2 的 Leave-One-Component-Out sensitivity analysis、O3 的 concordance matrix）。交叉驗證方法（零矛盾準則）有專門段落。 |
| 觀察 | ✅ | 9 個觀察共引用 18+ 張代表圖表，每個 O 觀察均有關鍵發現清單和代表圖表。Top 10 Cross-Observation Level A 發現表格是極好的整合呈現。 |
| 推論 | ⚠️ | **多數 O 觀察的推論僅為一句話解讀**，缺乏替代解釋的討論。例如：O4 發現「AF 在 TO 方向反轉 AUC=0.418」，僅附「高 AF = 更多 FP」的解讀，但未討論為什麼高 AF 在 TO 對應 FP（是 germline homozygous 的特徵？還是 caller bias？）。O8 發現「H2009 佔 36.2%」，但未深入推論 H2009 高佔比對其他統計結論的影響程度。O10 的「FP 傾向高甲基化區域」缺乏生物學機制推導。 |
| 結論 | ✅ | Top 10 發現表格、行動建議表格（P0/P1/P2 分級）、閱讀順序建議，結論層次分明。 |

**具體缺失項**：
1. O4 的 AF 反轉現象缺乏機制推論（高 AF TO FP 的生物學來源）
2. O8 的 H2009 樣本偏大問題對全域統計影響的定量推論不足
3. O10 的 FP 高甲基化傾向缺乏替代解釋討論（是 genomic context 效應？還是 caller 偏好？）
4. 部分觀察缺乏「排除替代解釋」的段落（如 O6 VerificationClass 無效——是 VC 設計問題？還是甲基化本身對 TP/FP 無區分力？）

---

### 06_methylation_hypothesis_negative.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | 開頭即聲明 negative result 的價值（三點），O11/O12 各有清晰假說表格，理論基礎說明充分（mQTL 穩定性 vs cancer stochasticity；三場景生物學差異）。 |
| 方法 | ✅ | O11 的 6 個 heterogeneity 特徵定義表格（含物理意義）、residualization 流程說明、O12 的四層 confound 控制（L0/L1/L2/L3）設計，方法論極為嚴謹。L3 優於 L2 的原因有專門解釋。13 個 novel features 定義完整。 |
| 觀察 | ✅ | O11 的 raw vs residualized AUC 對比表、confound 證據鏈（3 組數字）。O12 的 per-sample AUC 表、collider bias 案例表（8 例）、novel feature AUC 排名、整合穩定性排名。圖片清單含 15+ 張圖。 |
| 推論 | ✅ | O11 的 confound 推論鏈完整（n_reads -> epipolymorphism 機械性相關 -> AUC artifact）。O12 的 L2 collider bias 推導（Section 2.4）包含數學直覺、具體數字說明、教訓總結，是方法學貢獻的深度推論。AlleleDelta = AF confound 的推論也有充分支持。 |
| 結論 | ✅ | 兩個方向均有明確的「正式關閉」聲明。Section 三的總結論含正式排除表格、方法學收穫 3 點、唯一剩餘路徑、一致性確認 4 點。結論層次極豐富。 |

**無明顯缺失**。此文件是 negative result 報告的範本，L2 collider bias 發現尤其有方法學貢獻價值。

---

### 07_qs_mode_aware_change.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | Section 1 清楚呈現 TO QS AUC=0.497 的問題，根因分析（LOH penalty 和 verify bonus 方向反轉），佐證數據 3 點。 |
| 方法 | ✅ | Section 2 提供完整的 QualityScoreWeights struct 設計、Paired vs TO 權重對比表（含每項修改理由）、mode detection 邏輯程式碼、改動範圍（1 file, +48/-21）。可完全重現。 |
| 觀察 | ⚠️ | **缺乏修改後的實際驗證數據**。Section 3 的 AUC 改善僅為「預估」（「~0.546」帶波浪號），承認「需實際跑 benchmark 驗證」。圖表引用（Section 4）全部來自 O2 的修改前分析，沒有任何修改後的 before/after 圖。對於一個已 commit 的程式碼修改，缺少 post-implementation 驗證數據是顯著缺失。 |
| 推論 | ⚠️ | Section 3.2「效果有限的原因」有推論（其他 QS component 在 TO 也弱），但 Section 3.3「實際意義」的三點偏向**價值判斷**而非推論。缺乏對「為什麼移除兩個 component 只能提升 ~0.049」的定量推導（例如：兩個 component 的 weight 在總 QS variance 中佔多少比例？）。 |
| 結論 | ⚠️ | 缺乏獨立的「結論」段落。Section 3.3 部分承擔結論角色，但沒有明確的 summary statement。Section 5 的 commit 資訊是元數據而非結論。「待驗證問題」中的第 1 點（「修改後的 TO QS AUC 是否確實達到 ~0.546？」）實際上是本文件應該回答的核心問題，卻被放在「待驗證」裡。 |

**具體缺失項**：
1. **無修改後的 benchmark 驗證數據**（最嚴重缺失）——已 commit 的修改應有 before/after 對比
2. 缺乏兩個被移除 component 在 QS 總 variance 中佔比的定量分析
3. 缺乏獨立結論段落（帶有明確 summary statement）
4. effect_bonus_strong/moderate 在 TO 保持不變的決定缺乏數據支持（僅說「仍有參考價值」）

---

### 08_literature_survey.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | 隱含動機清晰：為本週觀察提供理論基礎和文獻對照。Section 1-2 的結構（LOH 生物學 + 甲基化區分文獻）與本週觀察（LOH enrichment + O11/O12 negative）直接對應。 |
| 方法 | ⚠️ | **缺乏文獻調查方法論**——未說明文獻檢索策略（搜尋關鍵字、資料庫、時間範圍、納入/排除標準）。讀者無法判斷文獻覆蓋的系統性和完整性。這在文獻調查文件中是標準要求。此外，Section 2.2 的工具比較表缺乏評選標準說明。 |
| 觀察 | ✅ | Section 1 的文獻共識表格（4 項）、Section 2.1 的 mQTL 效應表格（5 項）、Section 2.2 的工具比較表（5 工具）、Section 4 的 10 篇關鍵文獻摘要表。數據和引用充分。三場景理論的文獻對應描述詳細。 |
| 推論 | ✅ | Section 3「文獻如何支持/不支持本週觀察」是此文件的核心推論環，做得很好。四個子節（TO LOH 過判、甲基化區分力弱、O12 三場景不可區分、O11 heterogeneity 否決）逐一對照文獻與實驗觀察，並給出文獻支持度評級（強支持/部分支持/不矛盾但也不支持/間接支持 negative result）。Direction E 的推導有外部驗證支持。 |
| 結論 | ⚠️ | **缺乏獨立的總結論段落**。Section 3 的四個子節各有小結論，但沒有一個統合的 summary（如「文獻調查支持以下 N 項結論：...；質疑以下 M 項方向：...」）。Section 2.3 的「文獻的主要結論」是針對文獻本身的結論，而非本文件的結論。「待驗證問題」4 項和「認知門檻補充建議」承擔了部分結論角色，但不夠聚焦。 |

**具體缺失項**：
1. **缺乏文獻檢索方法論**（搜尋策略、資料庫、時間範圍、篩選標準）
2. 缺乏文獻覆蓋完整性的自我評估（如：「本調查涵蓋 2010-2026 年的 N 篇核心文獻，可能遺漏的方向包括...」）
3. 缺乏統合性總結論段落
4. Section 2.2 工具比較未說明選擇這 5 個工具的理由（是否有其他競爭工具被排除？）

---

### 09_conclusions_and_actions.md

| 環 | 評分 | 說明 |
|:---:|:---:|------|
| 動機 | ✅ | 作為整合文件，動機是彙整全週結論並提出行動建議，文件開頭的結構（確認 10 項 / 否決 4 項 / 待定 3 項）即是動機的體現。Section 7 的研究方向定位與 2026-Q2 策略對接。 |
| 方法 | ❌ | **作為整合結論文件，不需要實驗方法環，此項不適用但按框架標記為缺失**。然而，Section 2-4 的結論分類缺乏「如何從 Section 01-08 彙整出這 15+3+3 項結論」的方法論說明——例如：結論的強度如何判定？「確認」vs「否決」vs「待定」的分類標準是什麼？ |
| 觀察 | ✅ | Section 2 的 10 項確認結論（含數據摘要和來源）、Section 3 的 4 項否決假說（含否決依據和教訓）、Section 4 的 3 項待定問題（含所需數據和預估影響），數據引用完整。 |
| 推論 | ⚠️ | **行動建議的推論邏輯不夠清晰**。Section 5 的 P0/P1/P2 優先級排序缺乏判定依據（為什麼「建立 Paired/TO 分離策略」是 P0 而非 P1？）。Section 6 的下週工作排程缺乏資源估算和時間可行性論證。Section 5.2 的依賴關係圖有助推論，但缺乏對「如果 P0 延遲會如何影響 P1/P2」的討論。 |
| 結論 | ✅ | 本文件本身即是結論文件，Section 2-4 的三類結論清晰，Section 5 的行動建議結構化，Section 7 的研究方向定位與策略對接良好。 |

**具體缺失項**：
1. 結論分類標準（確認/否決/待定）的判定方法未明示
2. P0/P1/P2 優先級的判定邏輯缺乏說明
3. 下週工作排程缺乏工時估算與可行性論證
4. 缺乏「本週研究整體效率與方向正確性」的回顧性評估

---

## 整體發現

### 優勢

1. **數據引用極為豐富**：10 份文件共引用 40+ 張圖表、50+ 個數據表格，每個關鍵數字均可追溯至具體的 TSV/workspace。這在研究報告中罕見。
2. **推論鏈的典範**：04_mechanism_investigation.md 的四步推導和 06_methylation_hypothesis_negative.md 的 L2 collider bias 推導，是生物資訊報告中推論環的優秀範本。
3. **Negative result 的嚴謹處理**：O11/O12 的否決報告不僅說明了結果，還提供了方法學教訓和研究方向的正式關閉聲明。
4. **認知門檻補充建議**：每份文件末尾的「認知門檻補充建議」段落是非常有用的讀者輔助，降低了跨領域閱讀門檻。
5. **待驗證問題的誠實列舉**：每份文件均列出待驗證問題，標誌著作者對不確定性的坦誠態度。

### 共通弱點

1. **07_qs_mode_aware_change.md 缺乏 post-implementation 驗證**：已 commit 的程式碼修改應附帶 before/after benchmark 對比，而非僅有「預估」數字。
2. **08_literature_survey.md 缺乏檢索方法論**：文獻調查未說明系統性檢索策略，影響結論的可重現性。
3. **05_systematic_observation.md 的推論環偏薄**：9 個觀察中多數缺乏替代解釋的排除討論，部分發現（如 AF 反轉、H2009 偏大）停留在描述層面。
4. **部分文件缺乏獨立結論段落**：00_background.md、07_qs_mode_aware_change.md、08_literature_survey.md 沒有明確的「結論」或「Summary」段落。

### 建議改善順序

| 優先級 | 文件 | 改善項目 |
|:---:|------|---------|
| P0 | 07_qs_mode_aware_change.md | 補充修改後的 benchmark 驗證數據（before/after AUC 對比） |
| P1 | 08_literature_survey.md | 補充文獻檢索方法論段落、統合性結論 |
| P1 | 05_systematic_observation.md | 為 O4/O8/O10 補充替代解釋討論 |
| P2 | 00_background.md | 補充總結性結論段落 |
| P2 | 09_conclusions_and_actions.md | 補充結論分類標準與優先級判定邏輯 |

---

*本審查報告僅評估敘述完整性，不涉及數據正確性或科學結論的合理性判定。*
