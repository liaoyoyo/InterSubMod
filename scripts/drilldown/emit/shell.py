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


def render(boot: dict, caps: list, meta: dict) -> str:
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
  <b>{esc(boot['sample'])}</b>
  <span class="sep">│</span><span>sSNV {boot['l1']['n']:,}</span>
  <span class="sep">│</span><span>region {meta.get('regions', 0):,}</span>
  <span class="sep">│</span><span>{esc(meta.get('cap_summary', ''))}</span>
  <span class="sep">│</span><span>build {esc(meta.get('built_at', ''))}</span>
</div>

<div class="shell">

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
    <div class="capability-off">待 P3 實作。將包含 k×唯一解率、resolution×coarse、
      ISM 可得性×k（選擇偏誤前哨）等已驗證有結構的組合。</div>
  </section>

  <section class="section" id="view-selfcheck" style="display:none">
    <header><div><div class="eyebrow">自檢</div><h2>從結果確認整體沒有問題</h2></div></header>
    <div class="capability-off">待 P3 實作。將包含 C1–C8 守恆等式、覆蓋率缺口成因、
      輸入 provenance 清單。</div>
  </section>

</div>

<script type="application/json" id="bootData">__BOOT__</script>
<script>{_asset('dashboard.js')}</script>
<script>{_asset('tree.js')}</script>
<script>{_asset('detail.js')}</script>
<script>
(function () {{
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
}})();
</script>
</body></html>
"""
    # boot payload 最後才塞，避免 f-string 處理巨大 JSON
    from payload import json_for_html
    return body.replace("__BOOT__", json_for_html(boot))


def write(out_path: Path, boot: dict, caps: list, meta: dict) -> int:
    html_text = render(boot, caps, meta)
    tmp = out_path.with_suffix(out_path.suffix + f".stage.{id(boot)}")
    tmp.write_text(html_text, encoding="utf-8")
    tmp.replace(out_path)
    return len(html_text.encode("utf-8"))


_ = json
