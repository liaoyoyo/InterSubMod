# Example: 教 ISM 概念給新生 — Feynman + ELI5 + CPP + Diátaxis

> 場景: 教學, Who=新生, Why=解釋, What=概念, When=2-5min 對話

```markdown
---
framework: Feynman Technique + ELI5 + CPP (Concept-Procedure-Principle)
audience: 新研究生（剛加入 lab，無生資背景）
when: 2-5min 1-on-1 explainer
applied_by: narrative-frame v0.1
---

# Inter-Subclonal Modification (ISM) — 5 歲版

## Step 1: 一句話講給 5 歲小孩

「我們有兩隊細胞分開長，但 DNA 修飾沒有完全同步 — 我們在找『不同步的痕跡』來區分『真癌症 mutation』vs『讀錯的雜訊』。」

## Step 2: 找漏洞（哪裡卡住）

- ⚠ 「DNA 修飾」是什麼？(methylation — 化學標籤)
- ⚠ 「兩隊細胞」是什麼？(subclones — 同 tumor 內不同遺傳子族)
- ⚠ 「讀錯的雜訊」是什麼？(false positive — sequencing error 或 normal contamination)

## Step 3: 簡化 + 比喻

### 比喻 1: 雙胞胎髮型

- 想像同卵雙胞胎本來髮型一樣
- 但今天兩個出門前各自換不同髮型（methylation 差異）
- 我們透過「髮型差異」確認他們不是同一個人 — 即使 DNA（基因）完全一樣
- 一樣的，**ISM** 透過「methylation 差異」確認兩個 subclone 不同 — 即使 mutation 一樣

### 比喻 2: 撲克牌

- 一手牌（read）= 一條 ONT sequencing read
- 牌上的字（base）= ATCG nucleotide
- 牌的顏色（methyl）= 5mC modification
- Caller 看「字」決定是不是 mutation
- ISM 看「字 + 顏色」更準確 — 因為兩堆撲克牌（subclones）顏色不同

## Step 4: 重述 + 測試

### 重述

ISM = 用 ONT 同時看 DNA 序列 + methylation 標籤，當兩個 subclone 的 methylation pattern 不同時，能更準分辨「真 somatic mutation」vs「sequencing error」— 因為真 mutation 會帶著「來源 subclone 的 methylation 簽名」。

### 自測

✅ 新生能用自己的話講回來？
✅ 能舉一個範例（不是上面我的）？

---

## CPP（Concept-Procedure-Principle）

### Concept（是什麼）

**ISM** = 跨 subclone modification 差異訊號 — 用 ONT epigenetic context 強化 cancer caller。

### Procedure（怎麼做）

1. ONT BAM 讀 5mC tag
2. 每 variant 看周圍 methylation pattern（per-CpG）
3. 跨 read 比對 methyl pattern 一致性
4. Pattern 不一致 → 可能跨 subclone → 真 somatic
5. Pattern 一致 → 可能 sequencing error → 過濾

### Principle（為什麼有效）

**因為 epigenetic state 在 cell division 時可遺傳但會 drift**：
- Single mutation 起源於 1 個 founder cell
- Founder cell 的 methylation 在 division 時帶到所有後代
- 不同 subclone 經獨立 division → methylation 自然 diverge
- ⇒ 「真 mutation 帶著一致 methyl signature；sequencing error 隨機分佈」

---

## Diátaxis 補充（用什麼文件學）

- **Tutorial**: `InterSubMod/docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md`（hands-on 跑第一次）
- **How-to**: `InterSubMod/.claude/rules/workflow-commands.md`（跑 pipeline 配方）
- **Reference**: `InterSubMod/include/core/*.hpp`（API 規格）
- **Explanation**: 本檔（concept why）

---

Framework: Feynman Technique (Gleick《Genius》1992) + ELI5 (Reddit) + CPP (Instructional Design) + Diátaxis (Procida)

業界引用「If you can't explain it to a 6-year-old, you don't understand it yourself」(Feynman)
```

---

## 為什麼選此 hybrid（N2 推薦理由）

**5W 識別**:
- Who: 新生（無背景）
- Why: 解釋概念
- What: 領域核心概念
- When: 2-5min 對話
- How: 口頭 + 簡單 .md companion

**主框架 Feynman** — 用 5 歲語言重述 + 找漏洞 + 簡化（自我深度測試）
**Sub: ELI5** — 比喻具體化（雙胞胎髮型 / 撲克牌）
**Sub: CPP** — 系統化（concept / procedure / principle 三層）
**Sub: Diátaxis** — 指引後續學習資源（4 type doc 對照）

**不選擇**:
- BLUF: 缺 buildup，新生不懂 why care
- SCQA: 過於 business / consulting 風格
- A3: 工程框架不適合概念教學
- Pixar Spine: 故事過於 narrative，缺概念精確度
