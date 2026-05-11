# v2 PPT Rendering Rules

> 強制規則。任何修改 build_pptx.py 必須先讀完此檔。違反規則的 slide 必須在 screenshot_all.py 內被偵測為 issue。

## 1. Latin + CJK 雙字型 fallback per-character

每段文字的每個字元都必須同時宣告 latin typeface + eastAsia typeface。否則：
- 純 CJK 字型在英文/數字位置會渲染為方塊
- 純 Latin 字型在中文位置會渲染為方塊或缺字

實作方式（透過 lxml 操作 `<a:rPr>`）：

```python
rPr = run._r.get_or_add_rPr()
latin_el = etree.SubElement(rPr, qn("a:latin"))
latin_el.set("typeface", "Arial")
ea_el = etree.SubElement(rPr, qn("a:ea"))
ea_el.set("typeface", "Droid Sans Fallback")
```

且 segments 必須切割：依字元 codepoint 判斷 CJK / Latin 後 group 同類字元為一個 run。

## 2. 圖片等比 fit-within

`fit_image_within(slide, path, x, y, max_w, max_h)` 必須：

```python
ratio = min(max_w / iw_emu, max_h / ih_emu)
final_w = int(iw_emu * ratio)
final_h = int(ih_emu * ratio)
```

絕對禁止：
- `add_picture(path, x, y, width=W, height=H)` — 強制雙軸縮放會破壞縱橫比
- 取 `max(...)` — 會超出框

## 3. CJK matplotlib 字型注入

任何 matplotlib 自繪 schematic 必須在 import 後立即：

```python
import matplotlib.font_manager as fm
fm.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({
    "font.family": [<resolved_name>, "DejaVu Sans"],
    "axes.unicode_minus": False,
})
```

否則中文標題、軸標籤會渲染為方塊。

**禁止使用 emoji**（⭐🏅🎯🔴🟢🟡 等）— Droid Sans Fallback 不含 emoji glyph，會渲染為方塊。改用 ASCII / box symbols（≥ ≤ ★ □ ■ • × ✓）。

## 4. 每張 slide 必有 speaker note ≥ 350 字

`set_speaker_note(slide, text)` 內必檢查 `len(text) >= 350`，否則 raise。

speaker note 內容必須涵蓋：
- 開場 1 句（破題）
- 3 個重點（量化 / 機制 / 影響）
- 預備教授 Q&A（11 條中對應的）
- 結尾 transition（指向下一張）

speaker_script_v2.md 為人類講者用稿，每張 ≥ 460 字（含 Q&A 整合）；speaker note 為 PPT notes_text_frame 寫入版，可略簡但 ≥ 350 字硬性下限。

## 5. slide 文字最少（標題 ≤ 15 字；每頁 ≤ 5 視覺元素）

- 標題 ≤ 15 中文字（不含半形數字、半形字母）
- 每頁 visual element ≤ 5（box / image / table / arrow / text-block）
- bullet ≤ 5 條，每條 ≤ 15 字
- 細節入 speaker note，slide 不要塞滿

## 6. 強制口徑校準（reviewer 第二輪審核採納）

| 主題 | 必須這樣寫 |
|------|---------|
| Self-phasing 修補語法 | 「self-phasing 問題鏈被分層處理：V2b 解 phase scaffold；V3F/V5 解 tag 層 bug 與 HP33 directional reassignment。**V5 不修 self-phasing 本身**（V2b 已處理）」|
| 業界對照 | 「同實驗室相鄰工作」（不寫「業界共識」「標準替代」） |
| ISM 影響百分比 | 只列 count（14 / 共 85），來源報告誤寫 7% 須註明 |
| enum 行號 | `Util.h:20-25`（不是 21-25）|
| concordance | 「+8.3 pp（**clean PS blocks，全基因組**）」必須加 caveat |

## 7. S18 climax 防誤解 caveat

S18 V5max1 IGV 大圖必須在底部（不蓋圖）加小字 caveat：
「但 V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）」

## 8. S2 30 sec 限速 framing

S2 業界家族樹是 framing 不是比較細節，30 sec 帶過。細節留 S21 比較表。

## 9. screenshot_all.py 必檢項

- shapes within slide bounds（[0, 0, slide_w, slide_h]）
- speaker note ≥ 350 字
- 圖片 referenced 存在於 disk
- title 存在
- wireframe 渲染所有 shape + text（CJK + Latin 無方塊）
- 最終輸出 `[ISSUES] none detected`
