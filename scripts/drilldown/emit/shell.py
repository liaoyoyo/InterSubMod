#!/usr/bin/env python3
"""HTML shell 組裝。

CSS/JS 從 assets/ 讀進來 inline，不寫在 Python 字串裡 —— 那是
build_exact_ps_layered_workstation.py 的痛點（718 行 JS 埋在 f-string 中，
Python 3.9 又不准 f-string 內有反斜線）。分開放讓 JS 能被 linter 正常處理。
"""
from __future__ import annotations

import html
import json
from pathlib import Path

ASSETS = Path(__file__).resolve().parent.parent / "assets"


def esc(s) -> str:
    return html.escape("" if s is None else str(s), quote=True)


def _asset(name: str) -> str:
    return (ASSETS / name).read_text(encoding="utf-8")


def _cap_rows(caps: list) -> str:
    out = []
    for c in caps:
        link = c["linkage"]
        link_s = (f"{link['numerator']:,} / {link['denominator']:,}"
                  f" <span class='denom'>({link['rate'] * 100:.1f}%)</span>") if link else "—"
        probes = "".join(
            f"<li class='{'' if p['ok'] else 'bad'}'>{'✓' if p['ok'] else '✗'} "
            f"{esc(p['stage'])} {esc(p['detail'])}</li>" for p in c["probes"])
        sha = ""
        for p in c["paths"]:
            if isinstance(p, dict) and p.get("sha256"):
                sha = f"<code>{esc(p['sha256'][:12])}…</code>"
                break
        out.append(
            f"<tr><td><span class='dot dot-{esc(c['state'])}'></span>{esc(c['title'])}"
            f"<ul class='probe-list'>{probes}</ul></td>"
            f"<td>{esc(c['tier'])}</td>"
            f"<td>{esc(c['state_label'])}</td>"
            f"<td>{esc(c['failed_stage'] or '—')}</td>"
            f"<td class='num'>{link_s}</td>"
            f"<td>{sha or '—'}</td></tr>")
    return "".join(out)


def _selfcheck_presentation(analysis=None) -> dict:
    """Turn the three-state selfcheck summary into honest, user-facing status.

    FAIL is blocking.  With no FAIL, any SKIP still means the validation is
    incomplete.  Only an all-PASS result may use the word "passed".
    """
    summary = ((analysis or {}).get("selfcheck") or {}).get("summary") or {}
    n_pass = int(summary.get("pass", 0) or 0)
    n_fail = int(summary.get("fail", 0) or 0)
    n_skip = int(summary.get("skip", 0) or 0)
    n_total = int(summary.get("total", n_pass + n_fail + n_skip) or 0)
    if n_fail:
        return {
            "state": "fail",
            "badge": f"自檢 BLOCKED · {n_fail} FAIL · {n_skip} SKIP",
            "title": f"自檢未通過：{n_fail} 項不成立",
            "callout": "stop",
            "message": (f"有 {n_fail} 項守恆等式不成立；此產物不可宣稱驗證通過。"
                        f"另有 {n_skip} 項無法檢查。"),
        }
    if n_skip:
        return {
            "state": "incomplete",
            "badge": f"自檢 INCOMPLETE · {n_skip} SKIP",
            "title": f"自檢未完整：{n_skip} 項無法檢查",
            "callout": "warn",
            "message": (f"已通過 {n_pass} / {n_total} 項；另有 {n_skip} 項因缺資料層無法檢查。"
                        "SKIP 不是 PASS。"),
        }
    if n_total:
        return {
            "state": "pass",
            "badge": f"自檢 PASS · {n_pass}/{n_total}",
            "title": f"自檢通過：{n_pass} / {n_total} 項成立",
            "callout": "",
            "message": f"本次可執行的 {n_total} 項守恆等式全部通過。",
        }
    return {
        "state": "unavailable",
        "badge": "自檢 UNAVAILABLE",
        "title": "自檢狀態不可用",
        "callout": "warn",
        "message": "此頁未附自檢摘要，不能宣稱驗證通過。",
    }


def _ratio_text(metric: dict) -> tuple[str, str]:
    if not metric or not metric.get("available"):
        return "UNAVAILABLE", "缺必要資料層；不可解讀為 0"
    num = int(metric.get("num", 0) or 0)
    den = int(metric.get("den", 0) or 0)
    pct = f"{num / den * 100:.1f}%" if den else "分母為 0"
    return f"{num:,} / {den:,}", pct


def _observation_scope_html(boot: dict) -> str:
    """Render the observation universe before filters or locus details.

    The useful information architecture from the historical standalone browser is
    its explicit universe → selected → displayed accounting.  Here it is generalized
    to runtime payloads, without inheriting that page's sample-specific claims.
    """
    scope = boot.get("observationScope") or {}
    asset_policy = boot.get("assetPolicy") or {}
    universe = scope.get("universe") or {
        "available": True, "num": int((boot.get("l1") or {}).get("n", 0) or 0),
        "den": int((boot.get("l1") or {}).get("n", 0) or 0),
        "label": "topology 觀察母體",
    }
    metrics = [
        ("觀察母體", universe, "topology"),
        ("ISM summary 可得", scope.get("summary_join") or {}, "ISM"),
        ("非循環軸可檢定", scope.get("axis_testable") or {}, "derived"),
        ("任一全域軸 raw p ≤ 0.05", scope.get("axis_significant_raw") or {}, "derived"),
        ("位點細節可得", scope.get("detail_available") or {}, "ISM"),
    ]
    cards = []
    for title, metric, source in metrics:
        value, note = _ratio_text(metric)
        label = metric.get("label") or title
        cards.append(
            f"<div class='scope-card{' unavailable' if not metric.get('available') else ''}'>"
            f"<div class='scope-card-top'><span>{esc(title)}</span>"
            f"<span class='src src-{esc(source.lower())}'>{esc(source)}</span></div>"
            f"<div class='scope-value'>{esc(value)}</div>"
            f"<div class='scope-note'>{esc(note)} · {esc(label)}</div></div>"
        )

    definitions = scope.get("axis_definitions") or []
    if definitions:
        axis_rows = "".join(
            "<tr>"
            f"<td><code>{esc(d.get('id'))}</code><br>{esc(d.get('label'))}</td>"
            f"<td><code>{esc(d.get('field'))}</code></td>"
            f"<td>{esc(d.get('scope'))}</td>"
            f"<td>{esc(d.get('definition'))}</td>"
            f"<td><span class='axis-availability {'yes' if d.get('available_in_run') else 'no'}'>"
            f"{'run 有此軸' if d.get('available_in_run') else 'run 無此軸'}</span></td>"
            "</tr>" for d in definitions
        )
        axis_contract = (
            "<details class='axis-contract'><summary>甲基 axis 定義與全域編碼範圍</summary>"
            "<div class='details-body'><div class='table-wrap'><table class='data'>"
            "<thead><tr><th>軸</th><th>p 欄位</th><th>在本頁的範圍</th><th>定義 / 限制</th>"
            "<th>本 run</th></tr></thead><tbody>" + axis_rows + "</tbody></table></div>"
            f"<div class='denom' style='margin-top:.45rem'>{esc(scope.get('multiplicity') or '')}</div>"
            "</div></details>"
        )
    else:
        axis_contract = (
            "<div class='capability-off'><b>甲基 axis 定義不可用。</b>"
            "此產物沒有內嵌 observationScope metadata；不應從顏色猜測軸語意。</div>"
        )

    availability = "" if scope.get("available") else (
        "<div class='callout warn'><b>ISM 觀察層不可用。</b>" +
        esc(scope.get("reason") or "未提供原因") +
        "。上方 UNAVAILABLE 不代表 0 個位點。</div>"
    )
    skipped_assets = []
    if asset_policy.get("igv") == "skip":
        skipped_assets.append("逐 region IGV read 對齊圖")
    if str(asset_policy.get("methylPanels", "")) in ("", "0"):
        skipped_assets.append("甲基雙面板 PNG")
    asset_notice = "" if not skipped_assets else (
        "<div class='callout warn asset-policy'><b>此 build 為輕量觀察版。</b>"
        + esc("、".join(skipped_assets))
        + "未在本輸出中生成；summary/axis 數字仍可觀察，但缺圖不可解讀為沒有 read-level 訊號。"
        "完整影像查證須用 asset_policy=build/all 的獨立新輸出。</div>"
    )
    return (
        "<div class='interpretation-rail' aria-label='觀察到細節的閱讀流程'>"
        "<div><b>01 觀察母體</b><span>先看總數、可測數與缺失，不從精選個案外推。</span></div>"
        "<div><b>02 明示收旂</b><span>篩選後分母會常駐；預設仍是全顯示。</span></div>"
        "<div><b>03 位點查證</b><span>選 sSNV 後再看 region、read、甲基與 provenance 細節。</span></div>"
        "</div>"
        "<div class='callout warn claim-ceiling' role='note'><b>Claim ceiling：觀察性下鑽。</b>"
        "本頁未使用 truth set / HC BED 做 benchmark；甲基軸是事後 characterize。"
        "raw p ≤ 0.05 僅用於找候選位點，不是 cohort-wide 顯著結論。"
        "<button type='button' class='view-jump' data-view-target='selfcheck'>先看選擇偏誤前哨</button></div>"
        "<div class='scope-ledger' aria-label='觀察分母與缺失帳本'>" + "".join(cards) + "</div>"
        "<div class='denom scope-ledger-note'>前四張是 topology → summary 可得 → 可檢定 → "
        "raw-significant 探索漏斗，每步都重新標分母；「位點細節可得」是平行的產物 gate。"
        f"{esc(scope.get('missingness') or '')}</div>" + availability + asset_notice + axis_contract
        + f"<div class='denom legacy-crosswalk'>{esc(scope.get('legacy_crosswalk') or '')}</div>"
    )


def _analysis_mounts(analysis=None):
    analysis = analysis or {}
    co = analysis.get("cooccur")
    sc = analysis.get("selfcheck") or {}
    co_html = "" if isinstance(co, list) and co else (
        "<div class='capability-off'><b>共出現分析不可用。</b>"
        "此頁未內嵌 cooccur payload；不顯示空矩陣冒充 0。</div>"
    )
    sc_html = "" if sc.get("checks") else (
        "<div class='capability-off'><b>自檢細節不可用。</b>"
        "只有摘要或連摘要也沒有；請改讀同目錄 SELFCHECK.md。</div>"
    )
    return (f"<div id='cooccurContent' class='analysis-mount'>{co_html}</div>",
            f"<div id='selfcheckContent' class='analysis-mount'>{sc_html}</div>")


def render(boot: dict, caps: list, meta: dict, analysis: dict = None) -> str:
    validation = _selfcheck_presentation(analysis)
    scope_html = _observation_scope_html(boot)
    cooccur_mount, selfcheck_mount = _analysis_mounts(analysis)
    body = f"""<!DOCTYPE html>
<html lang="zh-Hant"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(boot['sample'])} · 多層下鑽觀察儀表板</title>
<meta name="intersubmod-ui-contract" content="drilldown-v1">
<meta name="intersubmod-sample" content="{esc(boot['sample'])}">
<meta name="intersubmod-ssnv-count" content="{boot['l1']['n']}">
<style>{_asset('dashboard.css')}</style>
</head><body>

<div class="authority">
  <b class="authority-sample">{esc(boot['sample'])}</b>
  <span class="sep">│</span><span class="authority-primary">sSNV {boot['l1']['n']:,}</span>
  <span class="sep">│</span><span class="authority-secondary">region {meta.get('regions', 0):,}</span>
  <span class="sep">│</span><span class="authority-secondary">{esc(meta.get('cap_summary', ''))}</span>
  <span class="sep">│</span><span class="authority-secondary">build {esc(meta.get('built_at', ''))}</span>
  <span class="validation validation-{esc(validation['state'])}" role="status">{esc(validation['badge'])}</span>
</div>

<div class="shell">

  <noscript><div class="callout stop" role="alert"><b>JavaScript 未啟用。</b>
    此互動儀表板無法載入篩選、分片與自檢細節；請改讀同目錄的 <code>SELFCHECK.md</code>，
    或啟用 JavaScript 後再判讀本頁。</div></noscript>

  <div class="levelbar">
    <div class="seg" role="group" aria-label="檢視層級">
      <button type="button" data-view="genome" aria-pressed="true">全基因組</button>
      <button type="button" data-view="cooccur" aria-pressed="false">共出現</button>
      <button type="button" data-view="selfcheck" aria-pressed="false">自檢</button>
      <button type="button" data-view="capability" aria-pressed="false">資料能力</button>
    </div>
    <nav class="crumb" id="crumb" aria-label="目前位置"></nav>
  </div>

  <div id="shardError" class="callout stop" style="display:none"></div>

  <section class="section" id="view-genome">
    <header>
      <div>
        <div class="eyebrow">L1–L2 · 全基因組 → 區域</div>
        <h2>{boot['l1']['n']:,} 個 sSNV 全部顯示；勾選才收斂</h2>
      </div>
      <span class="src src-topology">◆ topology</span>
    </header>

    {scope_html}

    <div class="scope" id="scope"></div>

    <details open><summary>篩選維度（預設全勾 = 不篩掉任何東西）</summary>
      <div class="details-body"><div class="dims" id="dims"></div></div>
    </details>

    <details open><summary>顯示圖層（只改畫法，不改篩選結果）</summary>
      <div class="details-body"><div class="layers" id="layers"></div></div>
    </details>

    <div class="genome-wrap" id="genomeWrap">
      <canvas id="genomeCanvas" role="img"
        aria-label="全基因組 sSNV 分布；下方清單提供鍵盤可操作的等價內容"></canvas>
      <div class="genome-hint" id="genomeHint"></div>
    </div>

    <div class="browser" style="margin-top:1rem">
      <div>
        <div class="denom" style="margin-bottom:.3rem">清單 <span id="listCount"></span>（虛擬捲動，不截斷）</div>
        <div class="vlist" id="vlist" role="listbox" tabindex="0"><div class="vlist-inner" id="vlistInner"></div></div>
      </div>
      <div class="detail" id="detail">
        <div class="empty">在上方圖上點一個 sSNV，或從左側清單選一筆。<br>
          同一像素若疊了多個位點，會列出全部讓你自己選。</div>
      </div>
    </div>
  </section>

  <section class="section" id="view-capability" style="display:none">
    <header>
      <div>
        <div class="eyebrow">資料能力</div>
        <h2>這份頁面用了哪些資料層，缺的缺在哪一段</h2>
      </div>
    </header>
    <p class="denom">硬核心只有 topology；其餘缺了會把對應面板降級並在此標明，不會整份拒繪。</p>
    <div class="table-wrap">
      <table class="cap-table">
        <thead><tr><th>能力 / 探測歷程</th><th>層級</th><th>狀態</th><th>停在</th>
          <th class="num">與硬核心連結</th><th>sha256</th></tr></thead>
        <tbody>{_cap_rows(caps)}</tbody>
      </table>
    </div>
  </section>

  <section class="section" id="view-cooccur" style="display:none">
    <header><div><div class="eyebrow">共出現</div><h2>維度之間的交集</h2></div></header>
    {cooccur_mount}
  </section>

  <section class="section" id="view-selfcheck" style="display:none">
    <header><div><div class="eyebrow">自檢</div><h2>{esc(validation['title'])}</h2></div></header>
    <div class="callout {esc(validation['callout'])} selfcheck-status" role="status">
      <b>{esc(validation['badge'])}。</b>{esc(validation['message'])}
    </div>
    {selfcheck_mount}
  </section>

</div>

<script type="application/json" id="bootData">__BOOT__</script>
<script type="application/json" id="analysisData">__ANALYSIS__</script>
<script>{_asset('dashboard.js')}</script>
<script>{_asset('tree.js')}</script>
<script>{_asset('locus.js')}</script>
<script>{_asset('methylgraph.js')}</script>
<script>{_asset('panel.js')}</script>
<script>{_asset('methyl.js')}</script>
<script>{_asset('detail.js')}</script>
<script>{_asset('analysis.js')}</script>
<script>
(function () {{
  var authority = document.querySelector('.authority');
  function syncStickyOffset() {{
    if (!authority) return;
    document.documentElement.style.setProperty(
      '--authority-height', Math.ceil(authority.getBoundingClientRect().height) + 'px');
  }}
  syncStickyOffset();
  window.addEventListener('resize', syncStickyOffset);
  if (window.ResizeObserver && authority) new ResizeObserver(syncStickyOffset).observe(authority);

  var segs = document.querySelectorAll('.seg button[data-view]');
  segs.forEach(function (b) {{
    b.onclick = function () {{
      segs.forEach(function (x) {{ x.setAttribute('aria-pressed', String(x === b)); }});
      ['genome', 'cooccur', 'selfcheck', 'capability'].forEach(function (v) {{
        var el = document.getElementById('view-' + v);
        if (el) el.style.display = (v === b.dataset.view) ? '' : 'none';
      }});
      if (b.dataset.view === 'genome' && window.__DD && window.__DD.refresh) window.__DD.refresh();
    }};
  }});
  document.querySelectorAll('[data-view-target]').forEach(function (b) {{
    b.onclick = function () {{
      var target = document.querySelector(".seg button[data-view='" + b.dataset.viewTarget + "']");
      if (target) {{ target.click(); target.focus(); window.scrollTo({{top: 0, behavior: 'smooth'}}); }}
    }};
  }});
}})();
</script>
</body></html>
"""
    # boot payload 最後才塞，避免 f-string 處理巨大 JSON
    from payload import json_for_html
    body = body.replace("__BOOT__", json_for_html(boot))
    return body.replace("__ANALYSIS__", json_for_html(analysis or {}))


def write(out_path: Path, boot: dict, caps: list, meta: dict, analysis: dict = None) -> int:
    html_text = render(boot, caps, meta, analysis)
    tmp = out_path.with_suffix(out_path.suffix + f".stage.{id(boot)}")
    tmp.write_text(html_text, encoding="utf-8")
    tmp.replace(out_path)
    return len(html_text.encode("utf-8"))


_ = json
