# PPT Rendering Rules · Self-Phasing 完整教授報告

> **適用範圍**：本目錄下所有 build_pptx.py / screenshot_all.py。修改前先讀完此檔。
> **教訓來源**：2026-04-23 PPTX 截圖教訓（feedback_pptx_screenshot_rendering_rules.md、feedback_matplotlib_cjk_font_rule.md）。

---

## 1. Latin + CJK 雙字型 fallback（per-character）

**問題**：
- 單一 CJK 字型（Droid Sans Fallback）會讓英文字母看起來變方塊或字距異常。
- 單一 Latin 字型（Arial / Calibri）會讓中文變方塊。

**解法**：在 `python-pptx` 寫文字時，對**每個字元**判斷是否為 CJK，並透過 `<a:rPr>` 設不同 typeface：
- Latin 字（U+0000–U+007F、U+0080–U+024F 等）→ `latin typeface="Arial"`
- CJK 字（U+4E00–U+9FFF、U+3000–U+30FF、U+FF00–U+FFEF 等）→ `eastAsia typeface="Droid Sans Fallback"`

實作見 `build_pptx.py::add_text_with_fallback`。**禁止**：對整段文字只設一種字型。

## 2. 圖片 fit-within（等比縮放，不擠壓）

**問題**：強制同時設定 `width=` 與 `height=` 會強制 PIL 拉伸圖片，失去原始比例。

**解法**：
```python
def fit_image_within(slide, path, x, y, max_w_emu, max_h_emu):
    img = PIL.Image.open(path)
    img_w_px, img_h_px = img.size
    # 將像素換算到 EMU（假設 96 DPI）
    img_w_emu = img_w_px * 9525
    img_h_emu = img_h_px * 9525
    ratio = min(max_w_emu / img_w_emu, max_h_emu / img_h_emu)
    final_w = int(img_w_emu * ratio)
    final_h = int(img_h_emu * ratio)
    # 居中
    cx = x + (max_w_emu - final_w) // 2
    cy = y + (max_h_emu - final_h) // 2
    slide.shapes.add_picture(path, cx, cy, width=final_w, height=final_h)
```

**禁止**：`add_picture(path, x, y, width=W, height=H)`（除非 W/H 已是等比例計算結果）。

## 3. matplotlib CJK 字型注入（如果在 build 內生新圖）

```python
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

CJK_FONT = '/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf'
fm.fontManager.addfont(CJK_FONT)
font_name = fm.FontProperties(fname=CJK_FONT).get_name()
plt.rcParams.update({
    'font.family': [font_name, 'DejaVu Sans'],
    'axes.unicode_minus': False,
})
```

**禁止字符**：⭐🏅⚠️🔴🟡🟢 等 emoji（Droid Sans Fallback 不含），改用：
- ⭐ → ★
- 🔴 → ●（紅色）
- 🟡 → ◐（黃色）
- 🟢 → ○（綠色）
- ⚠️ → ⚠

## 4. Speaker note 必填

每張 slide 必須呼叫：
```python
slide.notes_slide.notes_text_frame.text = "..."
```
**禁止空 notes**。screenshot_all.py 會檢查每頁 notes 長度 > 0。

## 5. Slide 文字密度限制

- 標題 ≤ 15 個字（中文）
- bullet ≤ 5 條
- 每條 ≤ 15 個字
- 長解釋一律放 speaker notes，**不上 slide**

違反時截圖驗證會 flag「文字溢出」。

## 6. 圖片可用清單

所有圖片在 `figures/`，禁止從外部複製。新圖片若需要，在 build_pptx.py 內用 matplotlib 自繪並暫存到 `output/auto_figs/`。

## 7. 截圖驗證標準

`screenshot_all.py` 會輸出每張 slide 的 PNG 並檢查：
1. **文字方塊比例** > 5%（OCR 不可用時改用顏色直方圖）
2. **白色佔比** > 95%（圖片缺失）
3. **shape 邊界溢出 slide 範圍**（PIL bounding box）
4. **notes 是否非空**（python-pptx 直接讀）

通過標準：`[ISSUES] none detected` 或所有 issues 有解釋為何不可解。

## 8. 不可變動的 ground truth

- `00_storyboard.md`
- `figures/*`
- `source_materials/*`

修改 build_pptx.py 時**不要修改以上檔案**。
