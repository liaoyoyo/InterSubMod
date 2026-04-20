# 03 — P2 Precision Decisions（5 題）

> **建立日期**: 2026-04-19

> **優先等級**：P2（用詞精確化 / 方法論補強）
>
> **依賴**：多數獨立可執行；P2-C 需等 Phase 2 A+D

---

## 問題 P2-A: LOH 層次混用精確化

**問題描述**：
結論 9「62% TO LOH 消失」在不同文件引用時未精確區分兩層概念：
- **LOH.bed**（VCF coordinate，AF/VAF-based）：不受 self-phasing 影響，Jaccard=1.0
- **ISM HP_Ratio LOH**（BAM HP tag-based）：受 self-phasing 影響，62% 消失

`CURRENT_FOCUS.md:53` 與 `07_LOH_CN_AF.md:102` 未補註層次差異。

**影響範圍**：
- 受影響結論：**C09**（Self-Phasing）, **C13**（LOH.bed SEQC2）, **C21**（LOH.bed clean）
- 受影響功能：F1 Self-Phasing（主要論文貢獻，需嚴謹敘述）
- 若不處理的風險：論文 Results 章節語義歧義；reviewer 質疑「哪個 LOH？」

**現況證據**：
- Audit card: `cards/C09_Self_Phasing.md` D1「⚠️ LOH 層次」
- 修正位置：`docs/CURRENT_FOCUS.md:53`, `docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md:102`
- Memory: `project_self_phasing_causal_chain_confirmed.md`

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 兩個檔案補註「（此指 ISM HP_Ratio LOH，非 LOH.bed）」+ 全文 grep「62%」加註 | 20 分鐘 | 低 |
| B | 新建術語表文件（`docs/glossary/LOH_terminology.md`）供未來引用；既有文件加連結 | 1 小時 | 低 — 根本解決但需持續維護 |
| C | 不修現有敘述，只在 C09 audit card 加註（已做） | 0 分鐘 | 中 — 其他文件仍混用 |

**驗證標準（無論選 A/B/C）**：
- **必達**：`grep -rn "62%" docs/` 所有結果附近 ±3 行包含「ISM HP_Ratio」或「LOH.bed」任一明確標註
- **驗收指令**：`grep -B1 -A1 "62%" docs/ | grep -v "HP_Ratio\|LOH.bed"` 預期 0 行

**依賴關係**：無（獨立可執行）

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P2-B: AUC 0.58 門檻定義統一

**問題描述**：
多個結論使用 AUC 0.58 作為 POSITIVE/NEGATIVE 判定門檻，但未統一定義。有些結論用 0.55, 0.6, 0.62 作為 cutoff。跨結論可比性差。

**影響範圍**：
- 受影響結論：**C03**（TO AUC<0.58）, **C07**（G1-G7 <0.64）, **C11**（Phase 1A）, **C16**（HPFineNGroups）
- 受影響功能：F5 Variant Confidence + 所有 AUC 聲明
- 若不處理的風險：reviewer 質疑「為何此處用 0.58 彼處用 0.62？」；內部決策依據模糊

**現況證據**：
- 00_INDEX.md D5 統計列顯示「⚠️ 14 項」
- Feedback：未成文，但多次在研究報告中出現不同 cutoff

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 建立 `docs/standards/auc_threshold_definition.md`：定義 POSITIVE=AUC≥0.58, STRONG=AUC≥0.70, WEAK=0.55-0.58；grep 既有文件修正 | 2 小時 | 低 |
| B | 只在新產出強制統一，既有文件保留 | 30 分鐘 | 中 — 既有不一致保留 |
| C | 每個結論保留原 cutoff 並附「本結論採用 X 門檻」註腳 | 1 小時 | 高 — 仍不統一 |

**驗證標準（無論選 A/B/C）**：
- **必達**：`docs/standards/auc_threshold_definition.md` 存在
- **必達**：22 audit card 全部引用此標準（或註腳說明例外）
- **驗收指令**：`grep -rn "AUC.*0\.[4-9]" docs/reports/audit/cards/ | head -20` 數值合理解讀

**依賴關係**：無（獨立）

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P2-C: ASM germline vs somatic 區分

**問題描述**：
C12 報告 ASM 32-66%（跨樣本），但未區分：
- **Germline ASM**（imprinting + 其他）：文獻 15-30%
- **Somatic ASM**（tumor-specific）：文獻未明確量化

C12 數值高於 germline 文獻值，暗示有 somatic 成分，但目前無法驗證（需 Normal BAM）。

**影響範圍**：
- 受影響結論：**C12**（ASM POSITIVE）
- 受影響功能：F4 ASM Profiling（主力 Phase 2 方向）
- 若不處理的風險：ASM 機制解釋過於籠統；論文 Results 章節深度不足

**現況證據**：
- Audit card: `cards/C12_ASM.md` D6「⚠️ germline/somatic」
- Phase 2 Role 2 審查 C-BIO-1
- 相關文獻（knowledge base 02_DNA_Methylation_Biology）

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 啟動 Phase 2 A+D Normal BAM 整合後區分；本審查只註明「pending Phase 2」並引用 germline 文獻 15-30% | 2-3 週（Phase 2 啟動） | 低 — 與研究路徑對齊 |
| B | 立即用 imprinting 區域先驗做部分區分（chr15 / chr11 / 其他 known imprinted regions） | 3-5 天 | 中 — 只能部分區分，不完整 |
| C | 不做區分，論文敘述改為「ASM 整體分佈」 | 0 天 | 中 — 論文深度不足 |

**驗證標準（無論選 A/B/C）**：
- **必達**：C12 audit card 補上文獻引用 germline 15-30%
- **必達**：Phase 2 A+D 完成後，C12 更新版本包含 germline vs somatic split
- **驗收指令**：`grep "germline" cards/C12_ASM.md` 有實質內容；Phase 2 報告輸出分層 ASM TSV

**依賴關係**：
- **前置依賴**：Phase 2 A+D Normal BAM 整合（選 A）
- **被阻塞項**：論文 F4 章節完整性

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P2-D: NG=3 非單調生物機制

**問題描述**：
HPFineNGroups 的 NG=2 → NG=3 → NG=4 的 TP rate 呈**非單調**（NG=3 低於 NG=2 與 NG=4）。目前無生物機制解釋。可能是：
- 生物機制：NG=3 對應「2 次 mutation events 中的中間態」不穩
- Artifact：NR=80 樣本量不足導致中間態估計不準
- 特徵設計缺陷：NG=3 在 HPFineNGroups 定義下本就邊界模糊

**影響範圍**：
- 受影響結論：**C16**（HPFineNGroups）
- 受影響功能：F2 Subclone Marker
- 若不處理的風險：reviewer 質疑「非單調暗示 marker 不 biologically meaningful」

**現況證據**：
- Audit card: `cards/C16_HPFineNGroups.md` D2「⚠️ C-BIO-3」
- Phase 2 Role 2 生物學家 C-BIO-3 觀察

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 補充「生物假說 vs artifact」雙軌討論段落於 C16 + 實驗 pilot：NR=120 時 NG=3 是否仍非單調 | 1 週 | 低 — 小 pilot 可驗證 artifact 假說 |
| B | 僅加假說討論，不做 pilot | 2 小時 | 中 — reviewer 可能要求數據 |
| C | 不解釋，論文敘述只 focus NG=4 filter | 0 天 | 高 — 刻意隱藏，reviewer 會抓 |

**驗證標準（無論選 A/B/C）**：
- **必達**：C16 audit card 新增「NG=3 非單調」討論段落（含假說 + 對應實驗設計）
- **觀察（選 A）**：NR=120 pilot 後 NG=3 TP rate 是否穩定
- **驗收指令**：`cards/C16.md` 含完整「Biological hypothesis vs Artifact」段落

**依賴關係**：
- **前置依賴**：P0-A CovM fix（pilot 需乾淨 baseline）
- **被阻塞項**：F2 論文 Results 章節

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P2-E: Self-Phasing 演算法層面根因

**問題描述**：
C09 Self-Phasing 報告觀察層面的 62% TO LOH 消失與 17.3:1 phasing bias，但**未從 LongPhase-TO 演算法層面**解釋為何 tumor-only mode 會產生 self-phasing。目前僅有「使用 tumor-derived haplotype 做 phasing 是 circular dependency」的黑箱解釋。

**影響範圍**：
- 受影響結論：**C09**（Self-Phasing mechanism）
- 受影響功能：F1 Self-Phasing（最強原創貢獻）
- 若不處理的風險：論文 reviewer 質疑「你們展示了 what 但未解釋 how」；F1 深度不足

**現況證據**：
- Audit card: `cards/C09_Self_Phasing.md` D2「⚠️ 演算法根因」
- Phase 2 Role 2 C-BIO-2 觀察
- Knowledge base: 02_LongPhase_TO_document（已有）

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 閱讀 LongPhase-TO 源碼 + 演算法層面補充「為何 TO mode 會產生 17.3:1 bias」章節於 C09 | 3-5 天 | 低 — 知識庫已有文件基礎 |
| B | 引用 LongPhase-TO 官方 document 與文獻，不讀源碼 | 1 天 | 中 — 深度不足 |
| C | 黑箱觀察（現況），不補機制 | 0 天 | 高 — F1 論文貢獻深度受損 |

**驗證標準（無論選 A/B/C）**：
- **必達**：C09 audit card 新增「Algorithmic Mechanism」段落
- **觀察（選 A）**：段落內含具體源碼路徑或演算法流程圖
- **驗收指令**：`cards/C09.md` 含 Mermaid 或流程圖說明 phasing 流程

**依賴關係**：無（獨立研究任務）

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## P2 區段總結

**推薦總動作**（全選 A）：
1. **P2-A 立即可做**（20 分鐘）：LOH 層次註解
2. **P2-B 立即可做**（2 小時）：AUC 門檻文件 + 現有文件修正
3. **P2-C 綁 Phase 2 A+D**：現階段先補文獻註
4. **P2-D 綁 CovM fix 後**：1 週 NR=120 pilot
5. **P2-E 獨立研究任務**：3-5 天 LongPhase-TO 源碼

**推薦理由**：
- P2-A/B 純文字低成本可立即清理
- P2-C 對齊研究路徑（Phase 2 A+D 自然觸發）
- P2-D 與 CovM fix 後 pilot 共享環境
- P2-E 大幅提升 F1 論文深度

**若全選 A 後的輸出**：
- `docs/standards/auc_threshold_definition.md`（新）
- `CURRENT_FOCUS.md / 07_LOH_CN_AF.md` LOH 註解補齊
- `cards/C09.md` Algorithmic Mechanism 段落
- `cards/C12.md` germline 文獻註
- `cards/C16.md` NG=3 討論 + NR=120 pilot 結果

---

## 關聯文件

- EXECUTIVE_DECISION_BRIEF.md 第 4-5 節
- `cards/C09.md`, `C12.md`, `C16.md`
- `cross_cutting/Characterization_Functions.md`（F1/F2/F4 對應）
- CHECKLIST.md
