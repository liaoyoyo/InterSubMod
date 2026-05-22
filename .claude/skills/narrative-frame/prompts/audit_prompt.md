# N6 — Audit Prompt（自審 checklist）

> AI 在 N6 跑自審 + 5 秒測試 + 缺漏標 + 用戶確認。

## 8 項自審 checklist（必跑）

對 N5 產出的 structured narrative 逐項 ✅ / ⚠ / ❌：

1. **每 framework section 都有對應 source？** 無 → ⚠ gap（具體列哪一 section）
2. **重要素材沒漏？** N3 萃取 list vs N4 mapping 對照（comm-skip）
3. **5 秒測試**：給同事看 5 秒能否説出 takeaway？模擬 — 讀首段 + 主視覺 → 寫一句 takeaway
4. **3 秒法則**：標題 + 主視覺 3 秒內傳達主訊息？標題 ≤ 15 字（中）/ 80 字元（英）
5. **Assertion-Evidence**：每段標題 = 結論句（不只是 topic 標籤）？❌ 範例「結果」 → ✅「V6 修 100% chr19 priority bug victims」
6. **6-item pre-publish**（design_principles.md Rule 12）：data-ink / CRAP / 留白 / hierarchy / 色 / colorblind 全 ✓？
7. **過度宣稱 check**：[O]/[I] 用了「證實 / 確認 / 解決 / 全部 / 完全」→ 改謙詞
8. **路徑前綴**：所有 .md ref 用 `InterSubMod/...`？

## Gap 標格式

```markdown
## ⚠ Gap（N6 自審）

- **<framework section>** 缺 source: <具體哪一段> — 推測但無 evidence；補 file:line 或降為 [I]
- **5 秒測試 FAIL**: 首段太密（X 字過長）；建議拆 stat-grid 4 個焦點數字
- **過度宣稱**: §X 用「證實」但 tier=[O]；改「初步觀察到」
- **標題非結論句**: §Y 標題「結果」太空泛；改 Assertion-Evidence
- **路徑前綴**: line 42 `docs/X.md` 缺 InterSubMod/ → 改 `InterSubMod/docs/X.md`
```

## C3 用戶確認格式

```
[narrative-frame C3] 套用完成

Framework: <name>
產出: X 字 / Y sections / Z 個 source citations

自審 8 項:
1. Section source ✅ / ⚠ <N items>
2. 重要素材 ✅
3. 5 秒測試 ✅ / ⚠ <reason>
4. 3 秒法則 ✅
5. Assertion-Evidence ✅
6. 6-item pre-publish ✅ / ⚠ <which item>
7. 過度宣稱 ✅
8. 路徑前綴 ✅

⚠ Gap: <N> 項（見上）

★ 是否接受？(y/n/edit/換框架/補 gap)
```

## 互動模式

- **互動**: C3 必停（即使全自動）— 用戶可：
  - (y) 接受
  - (n) 重來
  - (edit <section>) 局部修
  - (換框架 <X>) 跳回 N5 重套
  - (補 gap <section>) 跳回 N3 補萃取

## 失敗模式

- 多 gap（≥5 項）→ 建議跳回 N3 補 source 而非硬補
- 5 秒測試 FAIL → 強制改寫首段 / 改 framework
- 過度宣稱 ≥3 處 → 建議降 Tier（Tier 3 → Tier 2）

## Tier 對應

- **Tier 1**: skip audit
- **Tier 2**: 跑 1, 3, 5, 7（4 項簡審）
- **Tier 3**: 跑全 8 項 + 業界引用 footer + provenance
