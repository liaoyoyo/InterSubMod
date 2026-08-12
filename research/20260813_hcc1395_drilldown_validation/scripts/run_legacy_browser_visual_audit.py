#!/usr/bin/env python3
"""Compare the legacy standalone and a current drilldown page.

Two current modes are explicit:

* overlay: fulfill index.html in-memory after replacing its embedded UI assets;
  this is a PARTIAL fixture because the shell/payload was not regenerated.
* generated: load --current-url exactly as served, without altering CSS, JS, DOM,
  or responses.  Use this for a source-complete staging/generated page.

Neither mode edits source HTML or an external bundle.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from bs4 import BeautifulSoup
from playwright.sync_api import Page, sync_playwright


ASSET_ORDER = [
    "dashboard.js",
    "tree.js",
    "locus.js",
    "methylgraph.js",
    "panel.js",
    "methyl.js",
    "detail.js",
    "analysis.js",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--legacy", type=Path, required=True)
    p.add_argument(
        "--current-mode", choices=("overlay", "generated"), default="overlay",
        help=("overlay replaces embedded UI assets in-memory and is PARTIAL; "
              "generated tests --current-url exactly as served"),
    )
    p.add_argument(
        "--current-index", type=Path,
        help="immutable source index used only by --current-mode overlay",
    )
    p.add_argument("--current-url", required=True, help="URL to test")
    p.add_argument(
        "--assets", type=Path,
        help="current dashboard asset directory used only by --current-mode overlay",
    )
    p.add_argument("--figures", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    return p.parse_args()


def patched_current_html(index: Path, assets: Path) -> str:
    soup = BeautifulSoup(index.read_text(encoding="utf-8"), "html.parser")
    style = soup.find("style")
    if style is None:
        raise RuntimeError("current page has no embedded stylesheet")
    style.string = (assets / "dashboard.css").read_text(encoding="utf-8")

    scripts = soup.find_all("script")
    if len(scripts) < 11:
        raise RuntimeError(f"unexpected embedded script count: {len(scripts)}")
    for script, asset_name in zip(scripts[2:10], ASSET_ORDER):
        script.string = (assets / asset_name).read_text(encoding="utf-8")

    analysis = json.loads(soup.select_one("#analysisData").string)
    summary = analysis.get("selfcheck", {}).get("summary", {})
    n_pass = int(summary.get("pass", 0))
    n_fail = int(summary.get("fail", 0))
    n_skip = int(summary.get("skip", 0))
    if n_fail:
        state = "fail"
        badge = f"自檢 BLOCKED · {n_fail} FAIL · {n_skip} SKIP"
        title = f"自檢未通過：{n_fail} 項不成立"
        message = f"有 {n_fail} 項不成立；另有 {n_skip} 項無法檢查。"
        callout = "stop"
    elif n_skip:
        state = "incomplete"
        badge = f"自檢 INCOMPLETE · {n_skip} SKIP"
        title = f"自檢未完整：{n_skip} 項無法檢查"
        message = f"已通過 {n_pass} 項；SKIP 不是 PASS。"
        callout = "warn"
    else:
        state = "pass"
        badge = f"自檢 PASS · {n_pass}/{n_pass}"
        title = f"自檢通過：{n_pass} / {n_pass} 項成立"
        message = "所有可執行自檢均通過。"
        callout = ""

    authority = soup.select_one(".authority")
    validation = soup.new_tag("span")
    validation["class"] = ["validation", f"validation-{state}"]
    validation["role"] = "status"
    validation.string = badge
    authority.append(validation)

    selfcheck = soup.select_one("#view-selfcheck")
    selfcheck.select_one("h2").string = title
    status = soup.new_tag("div")
    status["class"] = ["callout", callout, "selfcheck-status"]
    status["role"] = "status"
    status.string = f"{badge}。{message}"
    selfcheck.select_one("header").insert_after(status)

    shell = soup.select_one(".shell")
    nojs = soup.new_tag("noscript")
    nojs.string = "JavaScript 未啟用；請改讀 SELFCHECK.md。"
    shell.insert(0, nojs)

    scripts[10].string = r"""
(function () {
  var authority = document.querySelector('.authority');
  function syncStickyOffset() {
    if (!authority) return;
    document.documentElement.style.setProperty(
      '--authority-height', Math.ceil(authority.getBoundingClientRect().height) + 'px');
  }
  syncStickyOffset();
  window.addEventListener('resize', syncStickyOffset);
  if (window.ResizeObserver && authority) new ResizeObserver(syncStickyOffset).observe(authority);
  var segs = document.querySelectorAll('.seg button[data-view]');
  segs.forEach(function (b) {
    b.onclick = function () {
      segs.forEach(function (x) { x.setAttribute('aria-pressed', String(x === b)); });
      ['genome', 'cooccur', 'selfcheck', 'capability'].forEach(function (v) {
        var el = document.getElementById('view-' + v);
        if (el) el.style.display = (v === b.dataset.view) ? '' : 'none';
      });
      if (b.dataset.view === 'genome' && window.__DD && window.__DD.refresh) window.__DD.refresh();
    };
  });
})();
"""
    return str(soup)


def page_metrics(page: Page, sticky_kind: str) -> dict:
    widths = page.evaluate("""() => ({
      viewport: innerWidth,
      clientWidth: document.documentElement.clientWidth,
      docScrollWidth: document.documentElement.scrollWidth,
      bodyScrollWidth: document.body.scrollWidth,
      bodyFontSize: getComputedStyle(document.body).fontSize,
      bodyLineHeight: getComputedStyle(document.body).lineHeight
    })""")
    widths["horizontalOverflow"] = max(widths["docScrollWidth"], widths["bodyScrollWidth"]) > widths["viewport"] + 1
    semantics = page.evaluate("""() => {
      const vis = e => !!(e.offsetWidth || e.offsetHeight || e.getClientRects().length);
      const roles = [...document.querySelectorAll('[role="img"]')].filter(vis);
      const imgs = [...document.querySelectorAll('img')].filter(vis);
      const svgs = [...document.querySelectorAll('svg')].filter(vis);
      const all = [...document.querySelectorAll('body *')].filter(vis);
      const small = all.filter(e => parseFloat(getComputedStyle(e).fontSize) < 12);
      return {
        roleImgVisible: roles.length,
        roleImgLabelled: roles.filter(e => e.getAttribute('aria-label') || e.getAttribute('aria-labelledby')).length,
        svgVisible: svgs.length,
        svgWithTitleOrLabel: svgs.filter(e => e.querySelector('title') || e.getAttribute('aria-label')).length,
        imgVisible: imgs.length,
        imgMissingAlt: imgs.filter(e => !e.getAttribute('alt')).length,
        figureCaptionsVisible: [...document.querySelectorAll('figcaption')].filter(vis).length,
        visibleTextBelow12px: small.length,
        visibleElementCount: all.length,
        hasDenominatorText: document.body.innerText.includes('分母'),
        hasMissingText: /缺|未檢定|未對齊/.test(document.body.innerText),
        hasLegendText: /圖例|allele 軸|HP 軸|未對齊/.test(document.body.innerText)
      };
    }""")
    page.evaluate("window.scrollTo(0, 700)")
    page.wait_for_timeout(120)
    if sticky_kind == "legacy":
        sticky = page.evaluate("""() => {
          const b = document.querySelector('.bar').getBoundingClientRect();
          return {barTop:b.top, barBottom:b.bottom, overlap:false};
        }""")
    else:
        sticky = page.evaluate("""() => {
          const a = document.querySelector('.authority').getBoundingClientRect();
          const l = document.querySelector('.levelbar').getBoundingClientRect();
          return {authorityTop:a.top, authorityBottom:a.bottom, levelbarTop:l.top,
                  levelbarBottom:l.bottom, overlap:l.top + 1 < a.bottom};
        }""")
    page.evaluate("window.scrollTo(0, 0)")
    return {"widths": widths, "semantics": semantics, "sticky": sticky}


def new_page(browser, viewport: dict):
    context = browser.new_context(viewport=viewport, color_scheme="light")
    page = context.new_page()
    errors = []
    page.on("console", lambda msg: errors.append(f"console:{msg.type}:{msg.text}") if msg.type == "error" else None)
    page.on("pageerror", lambda err: errors.append(f"pageerror:{err}"))
    page.on("requestfailed", lambda req: errors.append(f"requestfailed:{req.url}:{req.failure}"))
    return context, page, errors


def audit_legacy(browser, url: str, viewport: dict, screenshot: Path) -> dict:
    context, page, errors = new_page(browser, viewport)
    page.goto(url)
    page.wait_for_load_state("networkidle")
    initial = page.locator("#cnt").inner_text()
    page.select_option("#fa", "4")
    page.wait_for_timeout(80)
    filtered = page.locator("#cnt").inner_text()
    cards = page.locator("details.card:visible")
    visible_cards = cards.count()
    if visible_cards:
        cards.first.locator("summary").click()
    detail_open = bool(visible_cards and cards.first.get_attribute("open") is not None)
    metrics = page_metrics(page, "legacy")
    page.screenshot(path=str(screenshot), full_page=False)
    result = {
        "url": url,
        "viewport": viewport,
        "initialCount": initial,
        "filteredCount": filtered,
        "visibleCards": visible_cards,
        "detailOpen": detail_open,
        **metrics,
        "browserErrors": errors,
        "screenshot": str(screenshot),
    }
    context.close()
    return result


def audit_current(browser, url: str, body: str | None, mode: str, viewport: dict, screenshot: Path,
                  selfcheck_screenshot: Path | None = None,
                  cooccur_screenshot: Path | None = None,
                  detail_screenshot: Path | None = None,
                  methyl_screenshot: Path | None = None) -> dict:
    context, page, errors = new_page(browser, viewport)
    if body is not None:
        page.route("**/index.html", lambda route: route.fulfill(
            status=200, content_type="text/html", body=body))
    page.goto(url)
    page.wait_for_load_state("networkidle")
    page.wait_for_selector(".vrow")
    scope = page.locator("#scope").inner_text()
    list_count = page.locator("#listCount").inner_text()
    selected_locus = None
    if methyl_screenshot:
        selected_locus = page.evaluate("""() => {
          const boot = JSON.parse(document.getElementById('bootData').textContent);
          const axis = Array.from((boot.l1 && boot.l1.axis) || []);
          let i = axis.findIndex(code => (Number(code) & 7) !== 0);
          if (i < 0) i = axis.findIndex(code => ![0, 8, 16].includes(Number(code)));
          if (i < 0) throw new Error('No ISM-associated locus available for methyl detail QA');
          window.__DD.select(i);
          return {
            index: i,
            axisCode: Number(axis[i]),
            chrom: window.__DD.CHROMS[window.__DD.chromOf[i]],
            pos: Number(window.__DD.posOf[i])
          };
        }""")
    else:
        page.locator(".vrow").first.click()
    page.wait_for_timeout(350)
    page.wait_for_function("""() => {
      const detail = document.querySelector('#detail');
      return detail && !detail.innerText.startsWith('載入 ') &&
             detail.innerText.includes('甲基分群');
    }""")
    detail = page.locator("#detail").inner_text()[:500]
    if detail_screenshot:
        page.locator("#detail").scroll_into_view_if_needed()
        page.wait_for_timeout(120)
        page.screenshot(path=str(detail_screenshot), full_page=False)
        page.evaluate("window.scrollTo(0, 0)")
    methyl_section_text = None
    if methyl_screenshot:
        methyl_summary = page.locator("#detail summary").filter(
            has_text="甲基分群 — 沿哪個軸？").first
        methyl_summary.scroll_into_view_if_needed()
        page.wait_for_timeout(120)
        methyl_section_text = methyl_summary.locator("xpath=..").inner_text()[:1400]
        if "此位點沒有 ISM 資料" in methyl_section_text:
            raise RuntimeError("Selected methyl detail QA locus has no ISM record")
        page.screenshot(path=str(methyl_screenshot), full_page=False)
        page.evaluate("window.scrollTo(0, 0)")
    metrics = page_metrics(page, "current")
    page.screenshot(path=str(screenshot), full_page=False)

    page.locator("button[data-view='cooccur']").click()
    page.wait_for_timeout(120)
    cooccur = {
        "tables": page.locator("#view-cooccur table").count(),
        "denominatorLabels": page.locator("#view-cooccur .denom").count(),
        "fakeFilterCells": page.locator("#view-cooccur td[data-x]").count(),
        "text": page.locator("#view-cooccur").inner_text()[:600],
    }
    if cooccur_screenshot:
        page.screenshot(path=str(cooccur_screenshot), full_page=False)

    page.locator("button[data-view='selfcheck']").click()
    page.wait_for_timeout(120)
    selfcheck = {
        "heading": page.locator("#view-selfcheck h2").inner_text(),
        "status": page.locator("#view-selfcheck .selfcheck-status").inner_text(),
    }
    if selfcheck_screenshot:
        page.screenshot(path=str(selfcheck_screenshot), full_page=False)
    result = {
        "url": url,
        "mode": mode,
        "evidence_status": "PARTIAL" if mode == "overlay" else "DIRECT_GENERATED",
        "fixture": (
            "in-memory current-source overlay on immutable HCC1395_v1 payload; "
            "shell/payload not regenerated"
            if mode == "overlay" else
            "direct generated/staging URL; no response or DOM overlay"
        ),
        "viewport": viewport,
        "scope": scope,
        "listCount": list_count,
        "selectedLocus": selected_locus,
        "detailText": detail,
        "methylSectionText": methyl_section_text,
        "cooccur": cooccur,
        "selfcheck": selfcheck,
        **metrics,
        "browserErrors": errors,
        "screenshot": str(screenshot),
    }
    context.close()
    return result


def main() -> int:
    args = parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    args.result.parent.mkdir(parents=True, exist_ok=True)
    if args.current_mode == "overlay":
        if args.current_index is None or args.assets is None:
            raise SystemExit(
                "--current-mode overlay requires both --current-index and --assets")
        current_body = patched_current_html(args.current_index, args.assets)
        current_stem = "current_patched"
        current_start = 8
    else:
        if args.current_index is not None or args.assets is not None:
            raise SystemExit(
                "--current-mode generated loads --current-url directly; "
                "do not pass --current-index or --assets")
        current_body = None
        current_stem = "current_generated"
        current_start = 13
    legacy_url = args.legacy.resolve().as_uri()
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        results = {
            "schema_name": "intersubmod.legacy_browser_visual_audit",
            "schema_version": "1.0.0",
            "current_mode": args.current_mode,
            "evidence_status": "PARTIAL" if args.current_mode == "overlay" else "DIRECT_GENERATED",
            "legacy_desktop": audit_legacy(
                browser, legacy_url, {"width": 1440, "height": 1000},
                args.figures / "06_legacy_standalone_desktop.png"),
            "legacy_mobile": audit_legacy(
                browser, legacy_url, {"width": 390, "height": 844},
                args.figures / "07_legacy_standalone_mobile.png"),
            "current_desktop": audit_current(
                browser, args.current_url, current_body, args.current_mode,
                {"width": 1440, "height": 1000},
                args.figures / f"{current_start:02d}_{current_stem}_desktop.png",
                cooccur_screenshot=(
                    args.figures / f"{current_start + 2:02d}_{current_stem}_cooccur_desktop.png"),
                detail_screenshot=(
                    args.figures / f"{current_start + 4:02d}_{current_stem}_detail_desktop.png"),
                methyl_screenshot=(
                    args.figures / f"{current_start + 5:02d}_{current_stem}_methyl_detail_desktop.png")),
            "current_mobile": audit_current(
                browser, args.current_url, current_body, args.current_mode,
                {"width": 390, "height": 844},
                args.figures / f"{current_start + 1:02d}_{current_stem}_mobile.png",
                selfcheck_screenshot=(
                    args.figures / f"{current_start + 3:02d}_{current_stem}_selfcheck_mobile.png")),
        }
        browser.close()
    args.result.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        k: {
            "viewport": v["viewport"],
            "scrollWidth": v["widths"]["docScrollWidth"],
            "overflow": v["widths"]["horizontalOverflow"],
            "sticky": v["sticky"],
            "browserErrors": v["browserErrors"],
        }
        for k, v in results.items() if isinstance(v, dict) and "viewport" in v
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
