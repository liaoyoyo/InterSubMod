#!/usr/bin/env python3
"""lint_figure.py — 方法解釋圖客觀自檢（無 rasterizer 時的 pixel-free 設計稽核）.

對一個 rendered .svg（+ 可選同目錄 data.json/figure_spec.json）跑設計規範檢查，
產出 PASS / WARN / FAIL，讓 renderer 改→render→lint→fix 自我迭代到合規。
檢查（對齊 detail_levels.md + genomics_figure_conventions.md + explainer_checklist.md）：
  C1 結構：viewBox + <title> + <desc>
  C2 溢出：所有 x/y/cx/cy + rect(x+width)/(y+height) 在 canvas 內（用戶「文字超出範圍」防呆）
  C3 字型：title >=14；可讀文字 >=8（<8 FAIL，8-8.5 WARN micro）
  C4 認知負擔：<text> 元素數 >門檻 → WARN（busy）
  C5 label-flip：含 HP1-1/HP2-1 → 必含 'label-flip' 或 '≡'
  C6 色彩約定：haplotype 圖（含 HP1/HP2/germline/somatic/haplotag）→ 須有藍(2F5597/2E75B6)+橘(C55A11/ED7D31)
  C7 示意標註：data.json synthetic:true → svg 須含 '示意' 或 'synthetic'
退碼 0=無 FAIL / 1=有 FAIL。
用法: python3 lint_figure.py <svg> [--max-text 80]
"""
import argparse, json, os, re, sys

TEXT_BUDGET = 80   # <text> 元素數軟上限（busy 門檻）
FONT_MIN_READABLE = 8.0
TITLE_MIN = 14.0


def num(x):
    try:
        return float(x)
    except Exception:
        return None


def lint(svg_path, max_text):
    svg = open(svg_path, encoding="utf-8").read()
    ddir = os.path.dirname(os.path.abspath(svg_path))
    res = []  # (level, check, msg)

    def add(level, check, msg):
        res.append((level, check, msg))

    # C1 結構
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    if not m:
        add("FAIL", "C1", "缺 viewBox")
        W = H = None
    else:
        W, H = float(m.group(1)), float(m.group(2))
    add("PASS" if "<title>" in svg else "WARN", "C1", "有 <title>" if "<title>" in svg else "缺 <title>（accessibility）")
    add("PASS" if "<desc>" in svg else "WARN", "C1", "有 <desc>" if "<desc>" in svg else "缺 <desc>")

    # C2 溢出 canvas
    if W:
        over = []
        for attr, lim in (("x", W), ("cx", W), ("y", H), ("cy", H)):
            for v in re.findall(rf'\b{attr}="(-?[\d.]+)"', svg):
                f = num(v)
                if f is not None and (f < -2 or f > lim + 2):
                    over.append(f"{attr}={v}>{lim:.0f}")
        # rect x+width / y+height
        for rx, ry, rw, rh in re.findall(r'<rect x="(-?[\d.]+)" y="(-?[\d.]+)" width="([\d.]+)" height="([\d.]+)"', svg):
            if num(rx) + num(rw) > W + 2:
                over.append(f"rect x+w={num(rx)+num(rw):.0f}>{W:.0f}")
            if num(ry) + num(rh) > H + 2:
                over.append(f"rect y+h={num(ry)+num(rh):.0f}>{H:.0f}")
        add("FAIL" if over else "PASS", "C2", ("溢出 canvas: " + "; ".join(over[:6])) if over else "無元素溢出 canvas")

    # C3 字型
    fonts = sorted({float(x) for x in re.findall(r'font-size="([\d.]+)"', svg)})
    if fonts:
        too_small = [f for f in fonts if f < FONT_MIN_READABLE]
        add("PASS" if max(fonts) >= TITLE_MIN else "WARN", "C3", f"最大字級 {max(fonts)}pt (title 應 >={TITLE_MIN})")
        add("FAIL" if too_small else "PASS", "C3", (f"字級 <{FONT_MIN_READABLE} 不可讀: {too_small}") if too_small else f"字級皆 >={FONT_MIN_READABLE}（{fonts}）")
    else:
        add("WARN", "C3", "無 font-size（純圖形？）")

    # C4 認知負擔
    ntext = svg.count("<text")
    add("WARN" if ntext > max_text else "PASS", "C4", f"<text> 數={ntext} (軟上限 {max_text}{'，過於密集' if ntext>max_text else ''})")

    # C5 label-flip
    if re.search(r"HP[12]-1", svg):
        ok = ("label-flip" in svg) or ("≡" in svg)
        add("PASS" if ok else "FAIL", "C5", "含 somatic 子標且有 label-flip 註" if ok else "含 HP1-1/HP2-1 但缺 label-flip 註明")

    # C6 色彩約定（haplotype 圖）
    is_hap = bool(re.search(r"HP1|HP2|germline|somatic|haplotag|haplotype|Germline|Somatic", svg))
    if is_hap:
        blue = bool(re.search(r"2F5597|2E75B6", svg, re.I))
        orange = bool(re.search(r"C55A11|ED7D31", svg, re.I))
        add("PASS" if (blue and orange) else "WARN", "C6",
            f"haplotype 色約定 藍={blue} 橘={orange}" + ("" if (blue and orange) else " ← 應藍 germline/橘 somatic"))

    # C7 示意標註
    dj = os.path.join(ddir, "data.json")
    if os.path.exists(dj):
        try:
            d = json.load(open(dj, encoding="utf-8"))
            syn = d.get("schematic", {}).get("synthetic") or d.get("meta", {}).get("synthetic")
            if syn:
                ok = ("示意" in svg) or ("synthetic" in svg.lower())
                add("PASS" if ok else "FAIL", "C7", "synthetic 已標註" if ok else "data synthetic:true 但 svg 未標『示意/synthetic』")
        except Exception:
            pass

    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("svg")
    ap.add_argument("--max-text", type=int, default=TEXT_BUDGET)
    a = ap.parse_args()
    res = lint(a.svg, a.max_text)
    fails = [r for r in res if r[0] == "FAIL"]
    warns = [r for r in res if r[0] == "WARN"]
    icon = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}
    print(f"=== lint {os.path.basename(a.svg)} ===")
    for level, check, msg in res:
        print(f"  {icon[level]} [{check}] {msg}")
    verdict = "FAIL" if fails else ("WARN" if warns else "PASS")
    print(f"  → {verdict}  ({len(fails)} fail / {len(warns)} warn)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
