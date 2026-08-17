#!/usr/bin/env python3
"""把解釋頁內嵌的手刻 SVG 抽成獨立檔，供 GitHub README / Wiki 使用。

為什麼需要這支：
  - GitHub README 不吃內嵌 <style>，但吃 <img src="*.svg">；
    本專案的 SVG 全部用 presentation attributes（fill="#xxx" 直接寫在元素上），
    因此 GitHub 的 SVG sanitizer 不會破壞它們。
  - GitHub 有暗色模式，而這些 SVG 是淺色設計且原本沒有底色 —— 在暗色模式下
    深色文字會直接看不見。抽出時一律補一塊與頁面同色的背景 rect。
  - 檢視者的系統若缺中文字型，SVG 會顯示豆腐字；因此同時渲染一份 PNG
    （字型在本機已驗證正確）供 README 引用，SVG 留檔供編輯與論文使用。

用法: python3 tools/extract_svg_for_github.py
輸出: docs/images/<name>.svg + <name>.png
"""
import os
import re
import subprocess
import sys

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod"
SRC = os.path.join(ROOT, "docs", "explain")
OUT = os.path.join(ROOT, "docs", "images")
BG = "#FAF9F5"          # 與解釋頁 --c-bg 相同
SCRATCH = "/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/64a5bba6-e59f-4d9c-989b-9885a695e7bf/scratchpad"

# (來源檔, 該檔第幾個 svg, 輸出檔名, 用途說明)
PLAN = [
    ("11_system-map-overview.standalone.html", 1, "architecture-overview",
     "五層系統全景：資料 → 上游工具 → 兩支 C++ → Python → HTML"),
    ("11_system-map-overview.standalone.html", 2, "methylation-circularity",
     "為什麼甲基化不能用來重建譜系（cis-ASM 循環）"),
    ("11_system-map-overview.standalone.html", 3, "funnel-7samples",
     "全 7 樣本 funnel：469,849 個突變如何流失到 63,506 個單一拓撲"),
    ("12_intersubmod-io.standalone.html", 1, "ism-internal-pipeline",
     "InterSubMod 內部八階段處理鏈"),
    ("13_longlineage-io.standalone.html", 1, "longlineage-funnel",
     "LongLineage 位點流失漏斗（為什麼輸出 0 棵樹）"),
    ("14_upstream-data.standalone.html", 1, "upstream-toolchain",
     "上游前處理鏈與 sidecar 串流設計"),
    ("15_python-html-layer.standalone.html", 1, "workstation-refuse-design",
     "工作站生成器的拒絕渲染設計（必填指標缺失時 fail-closed；不驗證科學真實性）"),
    ("16_how-to-run.standalone.html", 1, "howto-six-steps",
     "操作六步驟與各自的驗收點"),
]

SVG_RE = re.compile(r"<svg\b.*?</svg>", re.S)


def add_background(svg: str) -> str:
    """在最底層插入背景色塊，避免 GitHub 暗色模式下深色文字看不見。"""
    m = re.search(r'viewBox="([\d.\s\-]+)"', svg)
    if not m:
        return svg
    parts = m.group(1).split()
    w, h = parts[2], parts[3]
    rect = f'<rect x="0" y="0" width="{w}" height="{h}" fill="{BG}"/>'
    # 插在 <desc>/<title> 之後、第一個繪圖元素之前
    for anchor in ("</desc>", "</title>"):
        if anchor in svg:
            i = svg.index(anchor) + len(anchor)
            return svg[:i] + "\n  " + rect + svg[i:]
    i = svg.index(">") + 1        # 沒有 title/desc 就緊接 <svg ...> 之後
    return svg[:i] + "\n  " + rect + svg[i:]


def render_png(svg_path: str, png_path: str, width: int) -> bool:
    """用 playwright 渲染，確保字型與本機驗證過的解釋頁一致。"""
    html = os.path.join(SCRATCH, "_svg_render.html")
    svg = open(svg_path, encoding="utf-8").read()
    # 去掉 XML 宣告，讓它能直接嵌進 HTML
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    open(html, "w", encoding="utf-8").write(
        '<meta charset="utf-8"><style>'
        'html,body{margin:0;padding:0;background:%s}'
        'svg{display:block;width:%dpx;height:auto}'
        '</style>%s' % (BG, width, svg))
    subprocess.run(
        [sys.executable, os.path.join(ROOT, "tools", "render_html_shot.py"),
         html, "-o", png_path, "--width", str(width), "--full"],
        capture_output=True, text=True, timeout=180)
    if not os.path.exists(png_path):
        return False
    trim(png_path)
    return True


def trim(png_path: str, pad: int = 12) -> None:
    """裁掉與背景同色的外框留白 —— --full 會抓整頁高度，SVG 下方常留一大片空白。"""
    try:
        from PIL import Image, ImageChops
    except ImportError:
        return
    im = Image.open(png_path).convert("RGB")
    bg = Image.new("RGB", im.size, tuple(int(BG[i:i + 2], 16) for i in (1, 3, 5)))
    box = ImageChops.difference(im, bg).getbbox()
    if not box:
        return
    l, t, r, b = box
    l, t = max(0, l - pad), max(0, t - pad)
    r, b = min(im.size[0], r + pad), min(im.size[1], b + pad)
    im.crop((l, t, r, b)).save(png_path)


def main():
    os.makedirs(OUT, exist_ok=True)
    manifest = []
    for src, idx, name, caption in PLAN:
        path = os.path.join(SRC, src)
        blocks = SVG_RE.findall(open(path, encoding="utf-8").read())
        if len(blocks) < idx:
            print("  [SKIP] %s 只有 %d 個 svg，取不到 #%d" % (src, len(blocks), idx))
            continue
        svg = add_background(blocks[idx - 1])
        if "xmlns=" not in svg.split(">")[0]:
            svg = svg.replace("<svg", '<svg xmlns="http://www.w3.org/2000/svg"', 1)
        svg_out = os.path.join(OUT, name + ".svg")
        open(svg_out, "w", encoding="utf-8").write(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            "<!-- %s\n     來源: docs/explain/%s (svg #%d)\n"
            "     由 tools/extract_svg_for_github.py 產生，勿手改；改上游解釋頁後重跑。 -->\n"
            % (caption, src, idx) + svg + "\n")

        vb = re.search(r'viewBox="[\d.\s\-]*?\s([\d.]+)\s([\d.]+)"', svg)
        w = int(float(vb.group(1))) if vb else 900
        png_out = os.path.join(OUT, name + ".png")
        ok = render_png(svg_out, png_out, min(max(w, 700), 1100))
        size_svg = os.path.getsize(svg_out)
        size_png = os.path.getsize(png_out) if ok else 0
        print("  %-28s svg %6.1f KB | png %s" %
              (name, size_svg / 1024, ("%6.1f KB" % (size_png / 1024)) if ok else "FAILED"))
        manifest.append((name, caption, size_svg, size_png))

    print("\n共產出 %d 組" % len(manifest))
    return manifest


if __name__ == "__main__":
    main()
