#!/usr/bin/env python3
"""Render presenter_v4, measure per-slide overflow (scrollHeight vs 648), screenshot key slides."""
import os
from playwright.sync_api import sync_playwright

D = "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2"
OUT = "/tmp/v4_shots"
os.makedirs(OUT, exist_ok=True)

errors = []
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1280, "height": 800})
    page.on("pageerror", lambda e: errors.append(f"[pageerror] {e}"))
    page.goto(f"file://{D}/presenter_v4.html", wait_until="networkidle")

    # per-slide overflow audit
    data = page.evaluate("""() => {
      return [...document.querySelectorAll('section.slide')].map((s,i) => {
        const tag = (s.querySelector('.tag')?.textContent || s.querySelector('.eyebrow')?.textContent || 'cover').trim();
        return {i, tag, scrollH: s.scrollHeight, clientH: s.clientHeight,
                over: s.scrollHeight - s.clientHeight};
      });
    }""")
    print("=== per-slide overflow (over>5px = content clipped in 16:9) ===")
    for d in data:
        flag = "  ⚠ OVERFLOW" if d["over"] > 5 else ""
        print(f"  [{d['i']:2d}] {d['tag'][:20]:20s} scroll={d['scrollH']:4d} client={d['clientH']:4d} over={d['over']:+4d}{flag}")

    # screenshot the actual verdict (idx 10) + fig1 (3) + fig2 (4)
    slides = page.query_selector_all("section.slide")
    for idx, name in {3: "fig1_dbeta", 4: "fig2_oldnew", 10: "verdict_G", 11: "s4_c3"}.items():
        slides[idx].scroll_into_view_if_needed()
        slides[idx].screenshot(path=f"{OUT}/pres_{idx:02d}_{name}.png")
    browser.close()

print("=== pageerrors ===", errors or "(none)")
print(f"OUT={OUT}")
