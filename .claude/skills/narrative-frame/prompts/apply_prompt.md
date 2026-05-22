# N5 — Apply Prompt

> AI 在 N5 套用 framework 時遵循的 prompt 指南。

## 任務

把 N3 萃取的素材 + N4 mapping 表，填入 `templates/<framework>_skeleton.md`，產出 structured narrative。

## 步驟

1. **讀 template** — `templates/<framework>_skeleton.md`
2. **填空** — 按 mapping 表 section ↔ 重點 對應填
3. **語氣校準**：
   - [F] Fact → 「已驗證」「確認為」「N=X 跨 Y 樣本」
   - [O] Observation → 「初步觀察」「N 不足需 Y 驗證」
   - [I] Inference → 「推測」「可能」「值得進一步觀察」
   - [U] Unconfirmed → 「待釐清」「需 X 才能確認」
4. **Source citation** — 每段標 `(source: InterSubMod/.../file:line)`
5. **Cohen ribbon + CI**（含數字 claim 必加）
6. **業界引用 footer**（Tier 3 必加）

## 輸出 mode 選擇

- `--mode=outline` — markdown outline（only headings + bullets）
- `--mode=script` — 口頭講稿（含 timing：[0:00-0:30] section X）
- `--mode=slide` — slide outline（每 section = 1 slide title + body bullets，含 speaker note）
- `--mode=companion` — .md 報告附 frontmatter `framework:` + section 對應 source citation
- `--mode=inline`（預設 Tier 2）— 對話 inline 回覆首行聲明 framework，後接結構化內容

## Tier 對應行為

- **Tier 1**: skip（無 framework）
- **Tier 2**: inline 模式 — 首行「用 <framework>：」+ 結構化內容（200-500 字）
- **Tier 3**: companion / slide / paper 模式 — 完整 section + source + cohen + provenance

## 自審強制

每 section 必有 ≥1 對應 source；無 → 標 ⚠ gap 進 N6。

## 輸出格式範例（Tier 2 inline）

```
用 SCQA（Barbara Minto, McKinsey 1960s）：

【Situation】上週 V6 修正 chr19 priority bug 後 cross-sample F1 ↑0.022 (HCC1395 / source: cycle 1 findings)
【Complication】但 4 樣本 caller-F1-ceiling 卡住（best τ=0.10 keep all / source: loso_hnew4 results）
【Question】要 pivot 嗎？還是 PROBE 找 low-F1 panel？
【Answer】Path B chr8-specific zone gate — 不依賴 LR threshold

⚠ Gap: §Situation 數據引用待補 file:line
```

## 輸出格式範例（Tier 3 companion .md）

```markdown
---
framework: SCQA + STAR per item
audience: PI
when: 5min PI 1-on-1
applied_by: narrative-frame v0.1
cycle_id: 20260522-XXXX
---

# <Title>

## §S — Situation

<上下文 paragraph>

(source: InterSubMod/research/.../findings.md:42)

## §C — Complication

<張力 paragraph>

(source: ...)

## §Q — Question

<核心問題>

## §A — Answer

<解答 / 建議>

(supporting evidence: ...)

### Provenance Footer

- Framework: SCQA (Minto 2020) + STAR per item
- Commit: <hash>
- Build: <timestamp>
- Tier: 3
```
