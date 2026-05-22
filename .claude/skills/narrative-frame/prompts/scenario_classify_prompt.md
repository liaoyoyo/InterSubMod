# N1 — Scenario Classification Prompt

> AI 在 N1 識別場景時遵循的 prompt 指南。

## 任務

對用戶當輪 raw request + 提供的 N 份 source docs，識別 5W 維度。

## 步驟

1. **Parse user prompt** — 找觸發詞、自然語句、明示場景
2. **Check 場景特例**（讀 `references/scenario_to_framework.md §3`）— 命中則直接路由 thin wrapper
3. **5W 維度填空**：

```yaml
Who: [PI / 教授 / 同儕 / 大眾 / 自己 / 學生 / 投資人 / 評審 / 下屬 / 客戶]
Why: [説服 / 解釋 / 報告進度 / 探索 / 紀錄 / 答辯 / 教學 / 比較 / 決策 / 回饋 / 行銷]
What: [新發現 / 進度 / 概念 / 比較 / 故事 / 決策 / 案例 / 數據 / 流程 / bug]
When: [30s / 2-5min / 18min / 45-60min / 紙本]
How: [口頭 / slide / .md / HTML / paper / email / Slack / 對話 inline]
```

4. **Cynefin domain 判斷**（front-gate）— 問「相同行動是否曾重複產生**可預測**結果？」
   - Yes → Clear / Complicated（best-practice OK）
   - Maybe → **Complex（強制 PROBE，禁套 deterministic framework）**
   - No → Chaotic（先穩定）

5. **Tier 判斷**：
   - 簡單問答 / single-line → **Tier 1 skip framework**
   - 200-500 字 / 跨 2-3 概念 → **Tier 2 inline 首行聲明**
   - ≥500 字 / 跨 ≥3 概念 / 多文件 → **Tier 3 完整 N1-N6**

## 輸出格式（給用戶 C1 確認）

```
[narrative-frame N1] 場景識別

Who: <選項> (理由: ...)
Why: <選項> (理由: ...)
What: <選項> (理由: ...)
When: <選項> (理由: ...)
How: <選項> (理由: ...)

Cynefin domain: <Clear / Complicated / Complex / Chaotic>
Tier: <1 / 2 / 3>

★ 是否正確？(y/n/edit <維度>=<值>)
```

## 互動模式

- **互動**: 暫停等用戶 ack / edit
- **全自動**: AI 預設標籤 → 30s 倒數通過

## 失敗模式

- 5W 無法判斷 → 問用戶 1-2 維度（不全問 5 維）
- Cynefin 模糊 → 預設 Complicated（中性 default）
- Tier 邊界（如 180 字）→ 偏 Tier 2（保守）
