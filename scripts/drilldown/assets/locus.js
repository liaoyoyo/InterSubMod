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

    /* 甲基軸：與 payload.py 的 AXIS_BITS 一致。
       刻意不含「甲基自身分群」—— 那軸是循環論證（自檢 C10），
       畫上去幾乎每個位點都亮，等於沒有資訊。 */
    var AX = [
        { id: "hp", label: "HP", color: "#285f8f" },
        { id: "allele", label: "ALT/REF", color: "#a94336" },
        { id: "lineage", label: "lineage", color: "#6b5592" }
    ];
    var NS_COLOR = "#c9c8bd";      // 有測但都不顯著
    var NT_COLOR = "#e7e5da";      // 未測 / 無 ISM 資料

    function axisOf(chrom, pos) {
        var rec = (DD.L4 && DD.L4[chrom]) ? DD.L4[chrom][String(pos)] : null;
        if (!rec) return { state: "none", sig: [], tested: 0 };
        var sig = [], tested = 0;
        AX.forEach(function (a) {
            var p = rec[a.id + "_p"], v = rec[a.id + "_valid"], n = rec[a.id + "_n"];
            if (v === false || (n !== null && n !== undefined && n < 2) ||
                p === null || p === undefined) return;
            tested++;
            if (p <= 0.05) sig.push(a);
        });
        if (!tested) return { state: "untested", sig: [], tested: 0, rec: rec };
        return { state: sig.length ? "sig" : "ns", sig: sig, tested: tested, rec: rec };
    }

    DD.locusStrip = function (row, focusPos) {
        var pos = (row.ap || []).slice().sort(function (a, b) { return a - b; });
        if (!pos.length) return "<div class='denom'>此 region 沒有 active position。</div>";

        var chrom = row.c;
        var af = {};
        (row.af || []).forEach(function (a) { af[a[0]] = { alt: a[1], ref: a[2] }; });
        var axes = pos.map(function (p) { return axisOf(chrom, p); });

        var W = 660, H = 132, padL = 30, padR = 30;
        var lo = pos[0], hi = pos[pos.length - 1], span = Math.max(hi - lo, 1);
        var x = function (p) { return padL + (W - padL - padR) * ((p - lo) / span); };

        var svg = ["<svg viewBox='0 0 " + W + " " + H + "' width='100%' style='max-height:" + H +
                   "px' role='img' aria-label='locus 排列、間距與甲基軸'>"];
        svg.push("<line x1='" + padL + "' y1='42' x2='" + (W - padR) + "' y2='42' " +
                 "stroke='#c9c8bd' stroke-width='2'/>");

        for (var i = 0; i + 1 < pos.length; i++) {
            var xa = x(pos[i]), xb = x(pos[i + 1]), mid = (xa + xb) / 2, gap = pos[i + 1] - pos[i];
            svg.push("<line x1='" + xa + "' y1='60' x2='" + xb + "' y2='60' stroke='#98a09a' " +
                     "stroke-width='1' stroke-dasharray='3 2'/>");
            svg.push("<text x='" + mid + "' y='74' text-anchor='middle' font-size='9' " +
                     "font-family='ui-monospace,monospace' fill='#667069'>" + human(gap) + "</text>");
        }

        // 甲基軸帶：每個位點一格，被切成 N 段（N = 顯著軸數），顏色即軸別
        svg.push("<text x='" + padL + "' y='" + (H - 2) + "' font-size='8' " +
                 "font-family='ui-monospace,monospace' fill='#667069'>" +
                 "↑ 甲基軸：該位點的甲基沿哪個軸分群</text>");
        pos.forEach(function (p, i) {
            var xx = x(p), a = axes[i], bw = 15;
            if (a.state === "none" || a.state === "untested") {
                svg.push("<rect x='" + (xx - bw / 2) + "' y='88' width='" + bw + "' height='11' " +
                    "fill='" + NT_COLOR + "' stroke='#c9c8bd' stroke-width='.5' " +
                    "stroke-dasharray='2 1.5'><title>" + DD.fmt(p) + "：" +
                    (a.state === "none" ? "無 ISM 資料" : "有 ISM 資料但無軸可檢定") +
                    "　—— 這不是「不顯著」，是沒測</title></rect>");
            } else if (a.state === "ns") {
                svg.push("<rect x='" + (xx - bw / 2) + "' y='88' width='" + bw + "' height='11' " +
                    "fill='" + NS_COLOR + "'><title>" + DD.fmt(p) + "：已檢定 " + a.tested +
                    " 個軸，沒有一個顯著（p ≤ 0.05）</title></rect>");
            } else {
                var seg = bw / a.sig.length;
                a.sig.forEach(function (ax, k) {
                    svg.push("<rect x='" + (xx - bw / 2 + k * seg) + "' y='88' width='" +
                        seg + "' height='11' fill='" + ax.color + "'><title>" + DD.fmt(p) +
                        "：甲基沿 " + ax.label + " 軸顯著（已檢定 " + a.tested + " 軸中的 " +
                        a.sig.length + " 個）</title></rect>");
                });
            }
        });

        // 位點刻度（可點）
        pos.forEach(function (p, i) {
            var xx = x(p), a = af[p] || { alt: 0, ref: 0 }, tot = a.alt + a.ref;
            var frac = tot ? a.alt / tot : 0;
            var focus = p === focusPos;
            var col = frac >= 0.9 ? "#a94336" : frac >= 0.5 ? "#b66e20" : "#285f8f";
            svg.push("<a href='#p=" + DD.esc(chrom) + ":" + p + "' class='locus-tick'>");
            svg.push("<rect x='" + (xx - 9) + "' y='12' width='18' height='46' fill='transparent'/>");
            svg.push("<line x1='" + xx + "' y1='30' x2='" + xx + "' y2='54' stroke='" + col +
                     "' stroke-width='" + (focus ? 3.4 : 2) + "'/>");
            if (focus) {
                svg.push("<circle cx='" + xx + "' cy='42' r='7' fill='none' stroke='#17231e' stroke-width='1.6'/>");
            }
            svg.push("<text x='" + xx + "' y='22' text-anchor='middle' font-size='9' " +
                     "font-family='ui-monospace,monospace' fill='" +
                     (focus ? "#17231e" : "#285f8f") + "' text-decoration='underline'>S" +
                     (i + 1) + "</text>");
            svg.push("<title>S" + (i + 1) + "　" + DD.esc(chrom) + ":" + DD.fmt(p) +
                     "　ALT " + a.alt + " / REF " + a.ref +
                     "　點擊跳到此位點</title></a>");
            svg.push("<text x='" + xx + "' y='84' text-anchor='middle' font-size='8.5' " +
                     "font-family='ui-monospace,monospace' fill='" + col + "'>" +
                     (tot ? (frac * 100).toFixed(0) + "%" : "—") + "</text>");
        });
        svg.push("</svg>");

        var cards = pos.map(function (p, i) {
            var a = af[p] || { alt: 0, ref: 0 }, tot = a.alt + a.ref;
            var frac = tot ? a.alt / tot : 0;
            var focus = p === focusPos, ax = axes[i];
            var axTag = ax.state === "sig"
                ? ax.sig.map(function (s) {
                      return "<span class='tag' style='color:" + s.color + "'>" +
                             DD.esc(s.label) + "</span>";
                  }).join("")
                : ax.state === "ns"
                    ? "<span class='tag' style='color:var(--muted)'>各軸皆不顯著</span>"
                    : "<span class='tag' style='color:var(--muted);border-style:dashed'>" +
                      (ax.state === "none" ? "無 ISM" : "未檢定") + "</span>";
            return "<a class='locus-card" + (focus ? " focus" : "") +
                "' href='#p=" + DD.esc(chrom) + ":" + p + "'>" +
                "<div class='site'><b>S" + (i + 1) + "</b>" + (focus ? "<span>目前</span>" : "") + "</div>" +
                "<div class='position'>" + DD.fmt(p) + "</div>" +
                "<div class='af'>" + (tot ? (frac * 100).toFixed(1) + "%" : "無覆蓋") + "</div>" +
                "<div class='reads'>ALT " + DD.fmt(a.alt) + " / REF " + DD.fmt(a.ref) +
                "　n=" + DD.fmt(tot) + "</div>" +
                "<div class='af-track'><span style='width:" +
                Math.max(0, Math.min(100, frac * 100)).toFixed(1) + "%'></span></div>" +
                "<div class='locus-ax'>" + axTag + "</div></a>";
        }).join("");

        var anySig = axes.some(function (a) { return a.state === "sig"; });
        return "<div class='locus-axis'>" + svg.join("") + "</div>" +
            "<div class='locus-strip'>" + cards + "</div>" +
            "<div class='legend' style='margin:.35rem 0'>" +
            AX.map(function (a) {
                return "<span class='bar-key'><i class='swatch' style='background:" +
                    a.color + "'></i>甲基沿 " + DD.esc(a.label) + " 軸顯著</span>";
            }).join("") +
            "<span class='bar-key'><i class='swatch' style='background:" + NS_COLOR +
            "'></i>已檢定但皆不顯著</span>" +
            "<span class='bar-key'><i class='swatch' style='background:" + NT_COLOR +
            ";border-style:dashed'></i>未檢定 / 無 ISM 資料</span></div>" +
            "<div class='denom' style='line-height:1.7'>" +
            "軸為真實 bp 間距（跨度 " + human(span) + "）。<b>S1/S2… 與卡片可點</b>，" +
            "直接跳到該位點的完整資訊。刻度顏色依 ALT 比例：" +
            "<b style='color:#a94336'>≥90%</b> / <b style='color:#b66e20'>50–90%</b> / " +
            "<b style='color:#285f8f'>&lt;50%</b>。<br>" +
            "<b>下方「甲基軸」帶</b>標出每個位點的甲基沿哪個軸分群 —— " +
            "一格被切成幾段就是有幾個軸同時顯著。" +
            (anySig ? "" : "本 region 沒有任何位點的甲基軸顯著。") +
            "虛線框 = <b>未檢定，不等於不顯著</b>。" +
            "<b>不含「甲基自身分群」軸</b>（循環論證，見自檢 C10）。<br>" +
            "<b>ALT 比例的分母是該位點 exact-PS 投影後的 ALT+REF read 數 —— 不是 VAF、不是 CCF。</b>" +
            "本樣本共現區 94% 為 CN-altered，read 豐度受拷貝數影響。</div>";
    };

})();
