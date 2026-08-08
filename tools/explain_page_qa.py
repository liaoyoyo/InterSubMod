#!/usr/bin/env python3
"""解釋頁交付前 QA — 抓「渲染看起來沒錯但其實壞掉」的三類問題。

起因（2026-08-06）：11_system-map-overview 的 funnel SVG 內誤用了 HTML 的 <b> 標籤，
瀏覽器 parser 在 SVG 命名空間遇到 HTML 標籤會中止 SVG 解析，
導致該標籤之後的所有 SVG 元素降級成純文字掉出圖外 —— 但 curl 200、無 pageerror、
無 broken image，只有實際截圖用眼睛看才發現核心圖少了一半。

檢查項目：
  1. SVG 良構性 — 逐個 svg 區塊做 XML 解析（會抓到 HTML 標籤污染、未閉合標籤）
  2. SVG 內的 HTML 標籤 — <text>/<tspan>/<desc>/<title> 內不得出現 <b>/<i>/<br> 等
  3. 缺字風險字元 — 目標環境字型渲染不出的 emoji（會變成豆腐方框）
  4. 無障礙 — 每個 svg 應有 <title> 與 <desc>
  5. 外部依賴 — standalone 頁不得有 <link rel=stylesheet> / <script src> / 外部 img

用法: python3 tools/explain_page_qa.py <file.html> [file2.html ...]
退出碼: 0 全過 / 1 有 FAIL
"""
import re, sys, os

# 優先用 defusedxml（擋 XXE / billion-laughs）；未安裝時退回 stdlib。
# 本工具只解析本 repo 自產的 HTML，但它會隨交接包給外部使用者，故採防禦性解析。
try:
    from defusedxml import ElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

# 目標環境（Noto Sans CJK TC 等）渲染不出、會變豆腐方框的字元。
# ✔✘⚠→←·—①②③ 等已實測可正常渲染，不列入。
TOFU_RISK = {
    "\U0001F534": "🔴 紅圓 → 改用 ● 加 color 樣式",
    "\U0001F7E0": "🟠 橘圓 → 改用 ● 加 color 樣式",
    "\U0001F7E1": "🟡 黃圓 → 改用 ● 加 color 樣式",
    "\U0001F7E2": "🟢 綠圓 → 改用 ● 加 color 樣式",
    "⭐":     "⭐ 星 → 改用 ★ (U+2605)",
    "✅":     "✅ → 改用 ✔ (U+2714)",
    "❌":     "❌ → 改用 ✘ (U+2718)",
    "\U0001F50D": "🔍 → 改用文字",
    "\U0001F4CC": "📌 → 改用文字",
    "\U0001F9EA": "🧪 → 改用文字",
}
HTML_IN_SVG = re.compile(
    r"<(?:text|tspan|desc|title)\b[^>]*>[^<]*<(b|i|u|strong|em|br|code|span|div|p)\b")


def check(path):
    fails, warns = [], []
    s = open(path, encoding="utf-8").read()
    name = os.path.basename(path)

    # 1 + 4. SVG 良構 + 無障礙
    svgs = list(re.finditer(r"<svg\b.*?</svg>", s, re.S))
    if not svgs:
        warns.append("頁面沒有任何 <svg>")
    for i, m in enumerate(svgs, 1):
        blk = m.group(0)
        try:
            root = ET.fromstring(blk)
        except ET.ParseError as e:
            line = s[: m.start()].count("\n") + 1
            fails.append("SVG #%d (約第 %d 行) XML 解析失敗 — %s" % (i, line, e))
            continue
        ns = "{http://www.w3.org/2000/svg}"
        has_title = root.find(ns + "title") is not None or root.find("title") is not None
        has_desc = root.find(ns + "desc") is not None or root.find("desc") is not None
        if not has_title:
            warns.append("SVG #%d 缺 <title>（螢幕閱讀器無法描述）" % i)
        if not has_desc:
            warns.append("SVG #%d 缺 <desc>" % i)

    # 2. SVG 內的 HTML 標籤
    for m in HTML_IN_SVG.finditer(s):
        line = s[: m.start()].count("\n") + 1
        fails.append("第 %d 行：SVG 文字元素內含 HTML 標籤 <%s> — "
                     "會中止 SVG 解析、其後元素全部掉出圖外。改用 <tspan>" % (line, m.group(1)))

    # 3. 缺字風險
    for ch, hint in TOFU_RISK.items():
        n = s.count(ch)
        if n:
            fails.append("含 %d 個缺字風險字元 U+%04X — %s" % (n, ord(ch), hint))

    # 5. stroke path 缺 fill="none" — SVG path 的 fill 預設是黑色，
    #    折線（≥3 個頂點）會被填成實心色塊。直線看不出來，折線一定露餡。
    for m in re.finditer(r"<path\b[^>]*>", s):
        tag = m.group(0)
        if "stroke=" in tag and "fill=" not in tag:
            d = re.search(r'd="([^"]*)"', tag)
            # 數頂點：L/M/H/V 指令 ≥3 個即為折線
            if d and len(re.findall(r"[MLHVmlhv]", d.group(1))) >= 3:
                line = s[: m.start()].count("\n") + 1
                fails.append("第 %d 行：折線 <path> 缺 fill=\"none\" — "
                             "會被預設黑色填滿成實心色塊" % line)

    # 6. 外部依賴
    #    先剝除 <style> 內容與 HTML 註解再檢查 —— 內嵌的 CSS 註解裡可能提到 <link ...>
    #    字面字串（例如記錄「原本用 link 外部引用」），那不是真的依賴。
    scrubbed = re.sub(r"<style\b.*?</style>", "", s, flags=re.S | re.I)
    scrubbed = re.sub(r"<!--.*?-->", "", scrubbed, flags=re.S)
    for pat, desc in [
        (r'<link[^>]+rel=["\']?stylesheet', "外部 CSS"),
        (r"<script[^>]+src=", "外部 JS"),
        (r'<img[^>]+src=["\']https?://', "遠端圖片"),
    ]:
        if re.search(pat, scrubbed, re.I):
            fails.append("standalone 頁含%s依賴，離線開啟會壞" % desc)

    print("── %s" % name)
    print("   svg=%d  size=%.1f KB" % (len(svgs), len(s.encode("utf-8")) / 1024))
    for f in fails:
        print("   [FAIL] %s" % f)
    for w in warns:
        print("   [WARN] %s" % w)
    if not fails and not warns:
        print("   ✔ 全部通過")
    elif not fails:
        print("   ✔ 無 FAIL")
    return len(fails)


if __name__ == "__main__":
    total = sum(check(p) for p in sys.argv[1:])
    print("\n總計 %d 個 FAIL" % total)
    sys.exit(1 if total else 0)
