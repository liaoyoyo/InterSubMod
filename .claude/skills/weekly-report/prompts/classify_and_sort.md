# classify_and_sort.md — W3+W4 → C2 prompt（合併版）

## 使用時機

W2 主線確認後，對每筆素材做：
- W3：4 層分類 [F]/[O]/[I]/[U]
- W4：5 維度評分 + 4 桶分流（PPT / 講稿 / 備註 / 暫存）

合併為一個 checkpoint C2，避免兩階段切割造成過多互動。

## 觸發前置條件

- C1 主線類型已鎖定
- W1 raw data 清單已備齊

## AI 自動處理（用戶 review）

對 W1 清單中每筆素材，AI 自動：

1. 套用 4 層分類決策樹（→ `references/LAYER_STRUCTURE.md` §C）
2. 計算 5 維度評分（→ `references/LAYER_STRUCTURE.md` §D）
3. 建議桶分流（PPT / 講稿 / 備註 / 暫存）
4. 列三維表格供用戶 review

## AskUserQuestion 模板

AI 先呈現分類+排序結果表格，再問：

```yaml
question: "本週 N 筆素材分類+分桶結果（見上方表格）。如何處置？"
header: "C2 分類+排序"
multiSelect: false
options:
  - label: "完全接受 AI 分類"
    description: "進入 W5+W6 邏輯紅旗檢查 + 教授追問預測。"
  - label: "個別調整（指出哪幾筆）"
    description: "用戶列「素材 X 改 [F] / 桶 Y」之類的調整指令。"
  - label: "重新評分（改評分權重）"
    description: "AI 重新計算 5 維度（如把『教授關心』權重提高）。"
  - label: "拆分某筆素材"
    description: "用戶覺得 AI 把多 sub-item 合一筆，要拆開分別評分。"
```

## 表格輸出格式

```markdown
### W3+W4 分類+排序結果

| # | 素材摘要 | [F]/[O]/[I]/[U] | 5 維分數 | 桶 | Tier |
|---|---------|----------------|---------|---|------|
| 1 | Phase 1A 7 樣本 ΔF1=+0.0112 | [F] | 5+5+5+5+4=24 | PPT | 1 |
| 2 | HPFineNGroups N=5 phasing | [O] | 4+3+4+5+3=19 | PPT | 1 |
| 3 | ClairS-TO 4.5.0 是否解決？ | [U] | 5+1+5+5+1=17 | 講稿 | 2 |
| 4 | LOH-constrained phasing 推測 | [I] | 3+2+3+3+2=13 | 講稿 | 2 |
| 5 | 自動化雜項 commit | (略) | <8 | 暫存 | (棄) |
| ... | | | | | |

**桶上限檢查**：
- PPT 桶：N 筆（上限 8）→ ✅ / ⚠ 超出
- 講稿桶：M 筆（上限 15）→ ✅ / ⚠ 超出
```

## 用戶處置流程

| 用戶選擇 | AI 行動 |
|---------|---------|
| 完全接受 | 寫入 internal state，進 W5+W6 |
| 個別調整 | 等用戶列調整指令 → 修改表格 → 再次顯示 → 確認 |
| 重新評分 | 詢問新權重 → 重算 → 顯示新表格 |
| 拆分 | 詢問如何拆 → 重分類+評分 → 顯示新表格 |

## fast-track 模式

C2 為**非必停** checkpoint（依 Q9 用戶決策）。fast-track 下：
- AI 完成自動分類+評分後，**不暫停 AskUserQuestion**
- 顯示表格 + 30 秒倒數，用戶可中斷修改，否則自動推進
- 結果存於 internal state，後續 C4 母稿確認時用戶仍可調整

## 過度宣稱紅旗自動掃描（W5 預先觸發）

C2 完成時，AI 同時掃描分類結果，標出潛在錯誤：
- 若 [F] 描述含「初步」「可能」→ 標 ⚠ 應降為 [O] 或 [I]
- 若 [O] 描述含「證實」「確認」→ 標 ⚠ 應改為觀察語氣
- 若 [F] 但 sample_count < threshold → 標 ⚠ 應降 [O]
- 若 PPT 桶 > 8 → 標 ⚠ 須降級

紅旗在 W5 階段正式處理（→ check_and_predict.md）。

---

## v2.1 升級：csv export option（修正點 3）

C2 表格大時，「個別調整」option 提供 csv 匯出：

```bash
# AI 在 internal state 寫出 csv
cat > /tmp/classify_sort.csv << CSV
id,summary,tag,score,bucket,tier
1,"Phase 1A 7 樣本 ΔF1=+0.0112",F,24,PPT,1
2,"HPFineNGroups N=5",O,19,PPT,1
...
CSV

# 用戶可用 sed/awk/csvkit 編輯
sed -i 's/^2,.*,O,/2,...,F,/' /tmp/classify_sort.csv

# AI ingest 編輯後 csv 重組分類
```

csv 比 free-text 結構化，避免「個別調整」需手寫 N 條指令。
