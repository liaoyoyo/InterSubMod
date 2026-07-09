#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""batch_render_qa.py — batch artifact-QA over many HTML deliverables (one browser).

Renders each HTML headlessly, captures pageerror + broken-image(naturalWidth==0) + blank-page
(near-empty body) + render time, writes a PNG per file, and prints a summary table so Claude
can then Read the flagged PNGs. This is the batch form of the render->Read closed loop
(verify-workstation W7 / HARNESS_FAILURE_LOG H-001).

Usage: python3 tools/batch_render_qa.py <out_dir> <html1> [html2 ...]  [--max-mb 40] [--width 1500]
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("html", nargs="+")
    ap.add_argument("--max-mb", type=float, default=40.0)
    ap.add_argument("--width", type=int, default=1500)
    ap.add_argument("--timeout", type=int, default=60000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    from playwright.sync_api import sync_playwright

    rows = []
    with sync_playwright() as pw:
        b = pw.chromium.launch(args=["--no-sandbox"])
        for h in a.html:
            # namespace by parent dir to avoid PNG collision (layered/H2009 vs topology/H2009)
            name = f"{Path(h).parent.name}__{Path(h).stem}"
            rec = {"file": h, "name": name}
            if not os.path.exists(h):
                rec["status"] = "MISSING"; rows.append(rec); continue
            mb = os.path.getsize(h) / 1048576
            rec["mb"] = round(mb, 1)
            if mb > a.max_mb:
                rec["status"] = "SKIP_TOO_BIG"; rows.append(rec); continue
            if os.path.getsize(h) < 200:
                rec["status"] = "EMPTY_FILE"; rows.append(rec); continue
            errs = []
            try:
                pg = b.new_page(viewport={"width": a.width, "height": 1000})
                pg.on("pageerror", lambda e: errs.append(str(e)))
                t0 = time.time()
                pg.goto(Path(h).resolve().as_uri(), wait_until="networkidle", timeout=a.timeout)
                pg.wait_for_timeout(1000)
                # scroll through page to trigger lazy-loaded images (else fold-below img
                # stay complete==false and broken-image undercounts; HARNESS_FAILURE_LOG H-010)
                pg.evaluate("""async () => { await new Promise(res => { let y = 0;
                    const t = setInterval(() => { window.scrollBy(0, window.innerHeight);
                    y += window.innerHeight;
                    if (y >= document.body.scrollHeight) { clearInterval(t); res(); } }, 60); }); }""")
                pg.wait_for_timeout(1500)
                pg.evaluate("window.scrollTo(0, 0)")
                broken = pg.evaluate("Array.from(document.images).filter(i=>i.complete && "
                                     "i.naturalWidth===0).map(i=>i.getAttribute('src'))")
                nimg = pg.evaluate("document.images.length")
                blen = pg.evaluate("document.body ? document.body.innerText.trim().length : 0")
                png = f"{a.out_dir}/{name}.png"
                pg.screenshot(path=png, full_page=False)  # viewport-only (big pages)
                pg.close()
                rec.update({"render_s": round(time.time() - t0, 1), "n_img": nimg,
                            "broken": broken[:8], "n_broken": len(broken),
                            "body_len": blen, "png": png})
                if blen < 40:
                    rec["status"] = "BLANK_PAGE"
                elif len(errs):
                    rec["status"] = "PAGEERROR"; rec["errors"] = errs[:3]
                elif len(broken):
                    rec["status"] = "BROKEN_IMG"
                else:
                    rec["status"] = "OK"
            except Exception as ex:
                rec["status"] = "RENDER_FAIL"; rec["exc"] = str(ex)[:160]
            rows.append(rec)
        b.close()

    json.dump(rows, open(f"{a.out_dir}/qa_summary.json", "w"), ensure_ascii=False, indent=1)
    print(f"{'status':14}{'mb':>6}{'img':>5}{'brk':>4}{'blen':>7}{'s':>6}  name")
    order = {"BLANK_PAGE": 0, "PAGEERROR": 1, "BROKEN_IMG": 2, "RENDER_FAIL": 3,
             "EMPTY_FILE": 4, "SKIP_TOO_BIG": 5, "MISSING": 6, "OK": 9}
    for r in sorted(rows, key=lambda x: order.get(x["status"], 8)):
        print(f"{r['status']:14}{r.get('mb',0):>6}{r.get('n_img','-'):>5}{r.get('n_broken','-'):>4}"
              f"{r.get('body_len','-'):>7}{r.get('render_s','-'):>6}  {r['name']}")
    bad = [r for r in rows if r["status"] in ("BLANK_PAGE", "PAGEERROR", "BROKEN_IMG", "RENDER_FAIL", "EMPTY_FILE")]
    print(f"\n{len(rows)} total | {sum(1 for r in rows if r['status']=='OK')} OK | "
          f"{len(bad)} flagged | {sum(1 for r in rows if r['status']=='SKIP_TOO_BIG')} skipped(too big)")
    if bad:
        print("FLAGGED (Read these PNGs):")
        for r in bad:
            print(f"  {r['status']}: {r.get('png', r['file'])}")


if __name__ == "__main__":
    main()
