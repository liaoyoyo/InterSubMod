# N2 — Framework Pick Prompt

> AI 在 N2 推薦框架時遵循的 prompt 指南。

## 任務

依 N1 場景結果 + Cynefin domain，從 `framework_catalog.md` 挑 **主框架 + 1-2 備框架** + 解釋為什麼。

## 步驟

1. **讀 `scenario_to_framework.md` §1** — 5W → framework 對照矩陣
2. **讀 `scenario_to_framework.md` §2** — 自然語句直接路由（命中則直接用）
3. **跑 Cynefin gate**：
   - Complex 域 → **禁止挑 SCQA / Pyramid / STAR / 任何 deterministic framework**
   - 改用 §9 決策（WRAP / Pre-Mortem）+ §6 教學（Feynman 反問）
4. **多維度交叉查詢**：
   - 主框架 = 命中最多維度的（Who + Why + What）
   - 備 1 = 切換 1 維度後最佳（如 When 變 30s）
   - 備 2 = 切換 2 維度後最佳（如 Who 變大眾）
5. **檢查 §6 hybrid 組合** — 主 + 內嵌 sub-framework
6. **檢查 §7 反 pattern** — 避免不推薦組合

## 輸出格式（給用戶 C2 確認）

```
[narrative-frame N2] 框架推薦

主框架: <name>
- 結構: <one-line>
- 為什麼選: 對應 5W <Who> + <Why> + <What>；業界源 <author year>
- 業界引用: 「<one-line quote>」
- Template: `templates/<name>_skeleton.md`

備框架 1: <name>
- 差異: <when to switch>
- 例: 若 <Y> 變 <Z>，這個更適合

備框架 2: <name>
- 差異: <when to switch>

不選擇:
- <name>: 不適合，因為 <reason>（針對 5W 哪一維）
- <name>: 不適合，因為 <reason>

★ 確認用 <主框架> 嗎？(y/n/換 X/混合 X+Y)
```

## 互動模式

- **互動**: C2 必停（即使全自動模式也必停 — framework 選擇是核心 lock）
- **用戶 override**: 説「換 X」立刻跳 N5；説「混合 X+Y」用 §6 hybrid

## 失敗模式

- catalog 中無對應 framework（例外場景）→ 用最接近 + 標 ⚠ + 問用戶要否新增 entry
- 用戶不熟 framework 名 → 用 `scenario_to_framework.md §2` 自然語句解釋
- 多個 framework 同等適合 → 列 top 3 + 各自 1 句 pros/cons

## Cynefin Complex 域特殊處理

```
[narrative-frame N2] Cynefin = Complex

⚠ Complex 域因果未確立 — 禁套 deterministic framework

推薦 PROBE-first 路徑:
1. WRAP Widen options - 拓寬選項
2. Pre-Mortem - 反推失敗
3. Feynman 反問 - 確認因果未泛化

★ 接受 PROBE 路徑？或 override 改用 deterministic（不推薦）？
```
