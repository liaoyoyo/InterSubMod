#!/usr/bin/env python3
"""build_gallery.py — 把 examples/*/method_explainer.svg 併成單頁 gallery 供一次眼驗.

讀每個 example 的 figure_spec.json(title/subtitle) + method_explainer.svg(inline 嵌入) + lint verdict，
產 examples/_gallery.html。用法: python3 build_gallery.py [examples_dir]
"""
import glob
import json
import os
import subprocess
import sys

ORDER = ["c1", "u1", "u2", "u3", "u4", "u5", "u6", "u7"]


def keyfn(p):
    name = os.path.basename(os.path.dirname(p))
    for i, pre in enumerate(ORDER):
        if name.startswith(pre + "_") or name == pre:
            return (i, name)
    return (99, name)


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "examples")
    base = os.path.abspath(base)
    svgs = sorted(glob.glob(os.path.join(base, "*", "method_explainer.svg")), key=keyfn)
    lintp = os.path.join(os.path.dirname(__file__), "lint_figure.py")
    cards = []
    for svg in svgs:
        d = os.path.dirname(svg)
        name = os.path.basename(d)
        spec_p = os.path.join(d, "figure_spec.json")
        title, sub = name, ""
        if os.path.exists(spec_p):
            sp = json.load(open(spec_p, encoding="utf-8"))
            title, sub = sp.get("title", name), sp.get("subtitle", "")
        verdict = "?"
        try:
            r = subprocess.run([sys.executable, lintp, svg], capture_output=True, text=True)
            verdict = "PASS" if r.returncode == 0 else "FAIL"
            for ln in r.stdout.splitlines():
                if "→" in ln:
                    verdict = ln.split("→")[1].strip()
        except Exception:
            pass
        svg_content = open(svg, encoding="utf-8").read()
        badge = "#2E7D52" if verdict.startswith("PASS") else ("#D98E04" if "WARN" in verdict else "#C0392B")
        cards.append(
            f'<section class="card"><div class="hd"><span class="tag">{name}</span>'
            f'<span class="badge" style="background:{badge}">lint {verdict}</span></div>'
            f'<h2>{title}</h2><p class="sub">{sub}</p><div class="fig">{svg_content}</div></section>')
    html = (
        '<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="UTF-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">'
        '<title>methods-example 圖例 gallery</title><style>'
        ':root{--bg:#FAF9F5;--ink:#141413;--soft:#6B6B66;--line:#E4E2DA;}'
        'body{margin:0;background:var(--bg);color:var(--ink);font-family:"Noto Sans CJK TC","Noto Sans TC",system-ui,sans-serif;line-height:1.6;}'
        '.wrap{max-width:980px;margin:0 auto;padding:26px 22px 80px;}'
        'h1{font-size:22px;margin:0 0 4px;}.lead{color:var(--soft);font-size:13.5px;margin:0 0 22px;}'
        '.card{background:#fff;border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:16px 0;box-shadow:0 1px 6px rgba(20,20,19,.05);}'
        '.hd{display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;}'
        '.tag{font-family:monospace;font-size:12px;color:var(--soft);}'
        '.badge{color:#fff;font-size:11px;font-weight:700;padding:2px 8px;border-radius:20px;}'
        'h2{font-size:16px;margin:2px 0 2px;}.sub{color:var(--soft);font-size:12.5px;margin:0 0 12px;}'
        '.fig{overflow-x:auto;}svg{max-width:100%;height:auto;display:block;}'
        '</style></head><body><div class="wrap">'
        '<h1>methods-example 圖例 gallery</h1>'
        f'<p class="lead">{len(cards)} 個方法解釋圖單元（C1 真值 pilot + U1-U7 模板）。全部 figure_spec.json + data.json → renderer 注入；lint 客觀自檢。請眼驗排版/可讀性/分級顯示。</p>'
        + "\n".join(cards) + '</div></body></html>')
    out = os.path.join(base, "_gallery.html")
    open(out, "w", encoding="utf-8").write(html)
    print(f"[OK] gallery → {out}  ({len(cards)} figures)")


if __name__ == "__main__":
    main()
