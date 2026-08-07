/* Locus 排列圖：一條線上的多個突變位點，看得到彼此的實際距離與各自的 ALT 比例。
 *
 * 表格答不出「這三個位點離多遠、哪個 ALT 比例高」，所以這裡畫成：
 *   上排  座標軸（等真實 bp 間距）+ 位點刻度 + 相鄰間距標註
 *   下排  每個位點一張卡：ALT/REF 計數 + ALT 比例條
 *
 * ⚠ ALT 比例的分母是該位點 exact-PS 投影後的 ALT+REF read 數，
 *   不是 VAF、不是 CCF。本樣本共現區 94% 為 CN-altered，read 豐度受拷貝數影響。
 */
(function () {
    "use strict";
    var DD = window.__DD;
    if (!DD) return;

    function human(bp) {
        if (bp >= 1e6) return (bp / 1e6).toFixed(2) + " Mb";
        if (bp >= 1e3) return (bp / 1e3).toFixed(1) + " kb";
        return bp + " bp";
    }

    DD.locusStrip = function (row, focusPos) {
        var pos = (row.ap || []).slice().sort(function (a, b) { return a - b; });
        if (!pos.length) return "<div class='denom'>此 region 沒有 active position。</div>";

        var af = {};
        (row.af || []).forEach(function (a) { af[a[0]] = { alt: a[1], ref: a[2] }; });

        var W = 660, H = 96, padL = 26, padR = 26;
        var lo = pos[0], hi = pos[pos.length - 1], span = Math.max(hi - lo, 1);
        var x = function (p) { return padL + (W - padL - padR) * ((p - lo) / span); };

        var svg = ["<svg viewBox='0 0 " + W + " " + H + "' width='100%' style='max-height:" + H +
                   "px' role='img' aria-label='locus 排列與間距'>"];
        // 主軸
        svg.push("<line x1='" + padL + "' y1='42' x2='" + (W - padR) + "' y2='42' " +
                 "stroke='#c9c8bd' stroke-width='2'/>");
        // 相鄰間距標註
        for (var i = 0; i + 1 < pos.length; i++) {
            var xa = x(pos[i]), xb = x(pos[i + 1]), mid = (xa + xb) / 2, gap = pos[i + 1] - pos[i];
            svg.push("<line x1='" + xa + "' y1='60' x2='" + xb + "' y2='60' stroke='#98a09a' " +
                     "stroke-width='1' stroke-dasharray='3 2'/>");
            svg.push("<text x='" + mid + "' y='74' text-anchor='middle' font-size='9' " +
                     "font-family='ui-monospace,monospace' fill='#667069'>" + human(gap) + "</text>");
        }
        // 位點刻度
        pos.forEach(function (p, i) {
            var xx = x(p), a = af[p] || { alt: 0, ref: 0 }, tot = a.alt + a.ref;
            var frac = tot ? a.alt / tot : 0;
            var focus = p === focusPos;
            var col = frac >= 0.9 ? "#a94336" : frac >= 0.5 ? "#b66e20" : "#285f8f";
            svg.push("<line x1='" + xx + "' y1='30' x2='" + xx + "' y2='54' stroke='" + col +
                     "' stroke-width='" + (focus ? 3.4 : 2) + "'/>");
            if (focus) {
                svg.push("<circle cx='" + xx + "' cy='42' r='7' fill='none' stroke='#17231e' stroke-width='1.6'/>");
            }
            svg.push("<text x='" + xx + "' y='22' text-anchor='middle' font-size='9' " +
                     "font-family='ui-monospace,monospace' fill='#17231e'>S" + (i + 1) + "</text>");
            svg.push("<text x='" + xx + "' y='90' text-anchor='middle' font-size='8.5' " +
                     "font-family='ui-monospace,monospace' fill='" + col + "'>" +
                     (tot ? (frac * 100).toFixed(0) + "%" : "—") + "</text>");
        });
        svg.push("</svg>");

        // 逐位點卡片
        var cards = pos.map(function (p, i) {
            var a = af[p] || { alt: 0, ref: 0 }, tot = a.alt + a.ref;
            var frac = tot ? a.alt / tot : 0;
            var focus = p === focusPos;
            return "<div class='locus-card" + (focus ? " focus" : "") + "'>" +
                "<div class='site'><b>S" + (i + 1) + "</b>" + (focus ? "<span>目前</span>" : "") + "</div>" +
                "<div class='position'>" + DD.fmt(p) + "</div>" +
                "<div class='af'>" + (tot ? (frac * 100).toFixed(1) + "%" : "無覆蓋") + "</div>" +
                "<div class='reads'>ALT " + DD.fmt(a.alt) + " / REF " + DD.fmt(a.ref) +
                "　n=" + DD.fmt(tot) + "</div>" +
                "<div class='af-track'><span style='width:" +
                Math.max(0, Math.min(100, frac * 100)).toFixed(1) + "%'></span></div></div>";
        }).join("");

        return "<div class='locus-axis'>" + svg.join("") + "</div>" +
            "<div class='locus-strip'>" + cards + "</div>" +
            "<div class='denom' style='margin-top:.35rem'>" +
            "軸為真實 bp 間距（跨度 " + human(span) + "）。刻度顏色依 ALT 比例：" +
            "<b style='color:#a94336'>≥90%</b> / <b style='color:#b66e20'>50–90%</b> / " +
            "<b style='color:#285f8f'>&lt;50%</b>。" +
            "<b>ALT 比例的分母是該位點 exact-PS 投影後的 ALT+REF read 數 —— 不是 VAF、不是 CCF。</b>" +
            "本樣本共現區 94% 為 CN-altered，read 豐度受拷貝數影響。</div>";
    };
})();
