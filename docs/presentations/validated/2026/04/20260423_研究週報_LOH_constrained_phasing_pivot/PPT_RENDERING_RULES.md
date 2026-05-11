<!--
建立時間: 2026-04-23
源自: 2026-04-23 週報 0423 生成過程教訓 + 使用者明確指示
儲存目的: 本資料夾自包含；未來 AI 複製此資料夾時必讀，避免重蹈覆轍
對應 memory: feedback_pptx_screenshot_rendering_rules
-->

# PPT 渲染強制規則（未來每次生成必守）

本資料夾 2026-04-23 首次生成時發生兩項關鍵錯誤，使用者已明確指示「以後 PPT 產生不要發生這問題」。以下兩規則**必守**：

---

## 規則 0 · Python matplotlib 圖片產生時 CJK 字型必須預先注入（2026-04-23 新增）

### 現象

Python matplotlib 產生的 PNG 圖中，中文字顯示為方塊 / 豆腐字 / 缺 glyph 警告：
```
UserWarning: Glyph 32681 (\N{CJK UNIFIED IDEOGRAPH-7FA9}) missing from current font.
```

### 根因

Matplotlib 預設 `font.family="DejaVu Sans"`（Latin-only）。對 CJK 字元會 fallback 到 `.notdef` glyph（方塊）或 raise warning。與 wireframe screenshot 的 PIL 不同，matplotlib **不會自動走 font-manager 的 fallback chain**；必須主動注入 CJK 字型。

### 正確做法（每個 Python 生圖腳本都必須加）

```python
import matplotlib
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt

CJK_FONT_PATH = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
try:
    fm.fontManager.addfont(CJK_FONT_PATH)
    CJK_FONT_NAME = fm.FontProperties(fname=CJK_FONT_PATH).get_name()
except Exception:
    CJK_FONT_NAME = "DejaVu Sans"

plt.rcParams.update({
    "font.family":     [CJK_FONT_NAME, "DejaVu Sans"],
    "font.sans-serif": [CJK_FONT_NAME, "DejaVu Sans"],
    "axes.unicode_minus": False,
    # ... 其他風格
})
```

### Emoji 亦需避免

`⭐ 🏅 ❤️` 等 emoji 不在 Droid Sans Fallback 與 DejaVu Sans 中，會 warning。改用：
- ⭐ → ★（U+2605 BLACK STAR）
- 🏅 → ★ 或純文字 "(✓)"
- ⚠️ → ⚠（U+26A0，通常有）

### 自動偵測

```bash
python3 script.py 2>&1 | grep -i "Glyph.*missing" && echo "⚠ CJK/emoji issue"
```
無輸出才算過關。

### 驗證範例

本規則於 2026-04-23 S20/S23 schematic 圖製作時確立。參考實作：
`scripts/analysis/20260423_s_multi_figures.py` 頂部字型注入段落

---

## 規則 1 · Latin + CJK 混合字型 fallback

### 現象

wireframe screenshot 中**英文 / 數字 / 符號**變 ☐ 方塊（但中文正常）。

### 根因

PIL `ImageDraw.text(..., font=ImageFont.truetype("DroidSansFallbackFull.ttf"))` 單一字型 render 時：

- `DroidSansFallbackFull.ttf` 是 **CJK-only** fallback 字型，沒有 Latin alphabet glyph
- PIL 不支援 per-char font fallback（每次 `draw.text` 只用一個字型）
- 結果：Latin 字元（N、G、=、2、LOH 等）找不到 glyph → 顯示 .notdef（方塊）

### 正確做法

實作 per-char font fallback helper（見 `screenshot_all.py` 的 `_draw_mixed` / `_text_w` / `_segment`）：

| 字元 code point range | 用字型 |
|---------------------|-------|
| `0x4E00-0x9FFF` CJK Unified | `DroidSansFallbackFull.ttf` |
| `0x3400-0x4DBF` CJK Extension A | 同上 |
| `0x3000-0x303F`、`0x3040-0x30FF` | 同上 |
| `0xFF00-0xFFEF` Full-width | 同上 |
| `0x2E80-0x2EFF`、`0x2F00-0x2FDF` | 同上 |
| 其他（Latin / Symbol / 數字）| `DejaVuSans.ttf`（或 `-Bold.ttf`）|
| font name 含 `consol/courier/mono` | `DejaVuSansMono.ttf` |

**流程**：先 `_segment(text)` 切 `(is_cjk, substring)` runs → 對每 run 取對應字型 → 累加 x 繪製。

**測量寬度**也必須用 mixed-font 版本（否則 `_wrap_text` 斷行位置錯）。

### 偵測違反

若 wireframe 中英文變方塊但中文正常 → 規則 1 違反。

---

## 規則 2 · 圖片等比 fit-within（禁止擠壓）

### 現象

熱圖 / 條圖 / 架構圖被水平或垂直壓扁；正方形圖變長方形。

### 根因

```python
slide.shapes.add_picture(path, left, top, width=Inches(w), height=Inches(h))
```

同時指定 `width` + `height` 會強制拉伸到目標尺寸，忽略原始 aspect ratio。

### 正確做法

`add_image(slide, path, x, y, w, h)` helper 必做 fit-within 計算：

```python
from PIL import Image as PILImage

def add_image(slide, path, x, y, w, h, halign="center", valign="middle"):
    with PILImage.open(str(path)) as im:
        ow, oh = im.size
    box_ratio, img_ratio = w / h, ow / oh
    if img_ratio > box_ratio:           # 圖較扁 → 填寬留上下空白
        actual_w, actual_h = w, w / img_ratio
    else:                               # 圖較高 → 填高留左右空白
        actual_w, actual_h = h * img_ratio, h
    dx = 0.0 if halign == "left" else (w - actual_w) if halign == "right" else (w - actual_w) / 2
    dy = 0.0 if valign == "top"  else (h - actual_h) if valign == "bottom" else (h - actual_h) / 2
    slide.shapes.add_picture(str(path), Inches(x + dx), Inches(y + dy),
                             width=Inches(actual_w), height=Inches(actual_h))
```

允許 box 內有留白（使用者明示接受），**禁止強制填滿**。

### wireframe 同步要求

`screenshot_all.py` 中 `shape_type == PICTURE` 分支必須**實際 `Image.paste`** 圖片（已由 add_image 確保 PPTX 中 shape 尺寸已等比），不得用 `[IMG]` placeholder。

### 偵測違反

若截圖中可識別圖被明顯拉伸 / 壓扁 → 規則 2 違反。

---

## 每次 PPTX 生成 pre-flight checklist

在宣告 PPTX v1 完成前，**必須全部勾選**：

- [ ] `build_pptx.py` 的 `add_image` 是否走 fit-within 分支？
- [ ] `screenshot_all.py` 是否包含 `_draw_mixed` + `_text_w` + `_segment` 三個 helper？
- [ ] 抽查 **3 張** wireframe：英文段、中文段、圖片三者皆正常顯示？
- [ ] `python3 screenshot_all.py` 末尾輸出 `[ISSUES] none detected · all text fits · all images loaded ✓`？
- [ ] 若輸出 ISSUES：先全部修掉再宣告完成

---

## 系統相依性（本專案環境實測）

| 字型 | 路徑 | 用途 |
|------|------|------|
| DejaVu Sans | `/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf` | Latin 內文 |
| DejaVu Sans Bold | `/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf` | Latin 粗體 |
| DejaVu Sans Mono | `/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf` | 程式碼 / 等寬 |
| Droid Sans Fallback Full | `/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf` | CJK 唯一可用 |

**系統無 Noto Sans CJK / Source Han Sans**。若未來安裝，可統一用 Noto Sans CJK TC（同時支援 Latin + CJK）取代上述 fallback 機制。

Python 套件：`python-pptx==1.0.2` / `Pillow`（已裝）。PPTX→PDF→PNG 轉換用的 `libreoffice` / `soffice` **未安裝**，故採 PIL wireframe 預覽。真實 PPTX 在本機 PowerPoint / Keynote / LibreOffice 開啟時字型會走系統 fallback（Arial + 微軟正黑體 / PingFang / Noto Sans CJK），不受本 wireframe 限制影響。

---

**Last updated**: 2026-04-23
**Trigger**: 本規則於 0423 週報生成第一版時被使用者兩次訂正後定案。
