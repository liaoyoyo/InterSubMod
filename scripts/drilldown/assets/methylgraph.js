/* 甲基群關係圖：節點 = 群（標群內甲基群數），邊 = 兩群之間甲基有明顯差異。
 *
 * 群的定義來自 RegionProcessor.cpp:2806-2832：
 *     樣本(N/T) × HP家族(HP1/HP2) × 等位(REF/ALT)   最多 8 個
 *     HP1 = {"1","HP1","1-1"}，HP2 = {"2","HP2","2-1"}
 *     ⚠ HP3（unphased somatic）不屬於任何群，完全不進這張圖
 *     每群需 ≥ min_group reads 才保留
 *
 * 邊 = ComboDbeta_SigPairs 的 "g1~g2=Δβ(q)"，999 次置換 + BH-FDR 校正。
 * 只有顯著的配對才會出現在該欄位，所以「沒有邊」意思是「沒有通過 FDR 的差異」，
 * 而不是「檢定過且沒差異」—— nTested 會告訴你到底測了幾對。
 */
(function () {
    "use strict";
    var DD = window.__DD;
    if (!DD) return;

    var GC = {                      // 群的顏色：HP1 藍系、HP2 紫系，ALT 深、REF 淺
        "T.HP1.REF": "#8fb4d6", "T.HP1.ALT": "#285f8f",
        "T.HP2.REF": "#c2aede", "T.HP2.ALT": "#6b5592",
        "N.HP1.REF": "#b8ccb8", "N.HP1.ALT": "#5f8f6b",
        "N.HP2.REF": "#d8c9a8", "N.HP2.ALT": "#b66e20"
    };

    function parsePairs(s) {
        // "T.HP1.REF~T.HP1.ALT=0.148870(0.001000);…"
        var out = [];
        String(s || "").split(";").forEach(function (chunk) {
            chunk = chunk.trim();
            if (!chunk) return;
            var m = chunk.match(/^([^~]+)~([^=]+)=(-?[\d.]+)\(([\d.eE+-]+)\)$/);
            if (m) out.push({ a: m[1], b: m[2], dbeta: parseFloat(m[3]), q: parseFloat(m[4]) });
        });
        return out;
    }

    function nodesFrom(rec, pairs) {
        // 節點只能從邊與 HP family read 數推 —— summary 沒有逐群 read 數，
        // 所以沒出現在任何顯著邊上的群不會有節點。這點必須講明。
        var names = {};
        pairs.forEach(function (p) { names[p.a] = 1; names[p.b] = 1; });
        var list = Object.keys(names).sort();
        var within = {
            HP1: rec.WithinHP1_NGroups, HP2: rec.WithinHP2_NGroups
        };
        return list.map(function (n) {
            var fam = n.indexOf(".HP1.") >= 0 ? "HP1" : (n.indexOf(".HP2.") >= 0 ? "HP2" : null);
            var w = fam ? parseInt(within[fam], 10) : null;
            return { name: n, fam: fam, withinK: (isNaN(w) ? null : w) };
        });
    }

    function layout(nodes, W, H) {
        // 兩欄：HP1 左、HP2 右；欄內 REF 上、ALT 下
        var cols = { HP1: [], HP2: [], other: [] };
        nodes.forEach(function (n) { (cols[n.fam] || cols.other).push(n); });
        var pos = {};
        [["HP1", W * 0.28], ["HP2", W * 0.72]].forEach(function (pair) {
            var key = pair[0], x = pair[1], arr = cols[key];
            arr.sort(function (a, b) { return a.name.localeCompare(b.name); });
            arr.forEach(function (n, i) {
                pos[n.name] = { x: x, y: arr.length === 1 ? H / 2 : 52 + i * ((H - 104) / (arr.length - 1)) };
            });
        });
        cols.other.forEach(function (n, i) { pos[n.name] = { x: W / 2, y: 24 + i * 30 }; });
        return pos;
    }

    DD.methylGraph = function (rec) {
        var nTested = parseInt(rec.ComboDbeta_NTested, 10);
        var nSig = parseInt(rec.ComboDbeta_NSig, 10);
        var pairs = parsePairs(rec.ComboDbeta_SigPairs);

        if (!nTested || isNaN(nTested)) {
            return "<div class='capability-off'>此位點沒有做群間 Δβ 檢定 —— " +
                "通常是沒有任何一組 <code>樣本×HP家族×等位</code> 達到最小 read 數。" +
                "<b>這不代表群間甲基沒有差異，是沒測。</b></div>";
        }
        if (!pairs.length) {
            return "<div class='callout'>測了 <b>" + nTested + "</b> 對群組，" +
                "<b>沒有一對通過 BH-FDR</b>（Δβ 差異不顯著）。" +
                "<div class='denom' style='margin-top:.3rem'>群 = 樣本(N/T) × HP家族(HP1/HP2) × 等位(REF/ALT)；" +
                "HP3（unphased somatic）不在任何群內。</div></div>";
        }

        var nodes = nodesFrom(rec, pairs);
        var W = 520, H = Math.max(210, 66 * Math.max(2, Math.ceil(nodes.length / 2)) + 60);
        var pos = layout(nodes, W, H);
        var maxAbs = Math.max.apply(null, pairs.map(function (p) { return Math.abs(p.dbeta); }));

        var svg = ["<svg viewBox='0 0 " + W + " " + H + "' width='100%' style='max-height:" + H +
                   "px' role='img' aria-label='甲基群關係圖'>"];

        pairs.forEach(function (p) {
            var a = pos[p.a], b = pos[p.b];
            if (!a || !b) return;
            var wdt = 1.4 + 5 * (Math.abs(p.dbeta) / (maxAbs || 1));
            var col = p.dbeta > 0 ? "#a94336" : "#285f8f";   // 正 Δβ 紅、負藍
            svg.push("<line x1='" + a.x + "' y1='" + a.y + "' x2='" + b.x + "' y2='" + b.y +
                "' stroke='" + col + "' stroke-width='" + wdt.toFixed(2) + "' opacity='.72'>" +
                "<title>" + DD.esc(p.a + " ↔ " + p.b) + "\nΔβ = " + p.dbeta.toFixed(4) +
                "\nBH 校正後 q = " + p.q.toExponential(2) + "</title></line>");
            var mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
            svg.push("<rect x='" + (mx - 30) + "' y='" + (my - 9) + "' width='60' height='15' rx='3' " +
                "fill='#fffdf7' opacity='.88'/>");
            svg.push("<text x='" + mx + "' y='" + (my + 2) + "' text-anchor='middle' font-size='9' " +
                "font-family='ui-monospace,monospace' fill='" + col + "'>Δβ " +
                (p.dbeta > 0 ? "+" : "") + p.dbeta.toFixed(3) + "</text>");
        });

        nodes.forEach(function (n) {
            var p = pos[n.name];
            var fill = GC[n.name] || "#c9c8bd";
            svg.push("<circle cx='" + p.x + "' cy='" + p.y + "' r='26' fill='" + fill +
                "' stroke='#17231e' stroke-width='1.2'><title>" + DD.esc(n.name) +
                (n.withinK !== null ? "\n群內甲基分群數 = " + n.withinK : "") + "</title></circle>");
            var parts = n.name.split(".");
            svg.push("<text x='" + p.x + "' y='" + (p.y - 3) + "' text-anchor='middle' font-size='9' " +
                "font-family='ui-monospace,monospace' fill='#fffdf7'>" + DD.esc(parts[1]) + "</text>");
            svg.push("<text x='" + p.x + "' y='" + (p.y + 8) + "' text-anchor='middle' font-size='9' " +
                "font-family='ui-monospace,monospace' fill='#fffdf7'>" + DD.esc(parts[0] + "·" + parts[2]) + "</text>");
            if (n.withinK !== null && n.withinK >= 2) {
                svg.push("<circle cx='" + (p.x + 21) + "' cy='" + (p.y - 20) + "' r='9' fill='#b66e20'/>");
                svg.push("<text x='" + (p.x + 21) + "' y='" + (p.y - 17) + "' text-anchor='middle' " +
                    "font-size='9' font-weight='700' fill='#fff'>" + n.withinK + "</text>");
            }
        });
        svg.push("</svg>");

        var w1 = rec.WithinHP1_NGroups, w2 = rec.WithinHP2_NGroups;
        return svg.join("") +
            "<div class='denom' style='line-height:1.75;margin-top:.4rem'>" +
            "<b>節點</b> = 群（樣本 N/T × HP家族 × 等位 REF/ALT）。" +
            "琥珀角標 = <b>該 HP 家族內部的甲基分群數 ≥2</b>" +
            "（HP1 內 " + DD.esc(w1 || "—") + " 群、HP2 內 " + DD.esc(w2 || "—") + " 群）。<br>" +
            "<b>邊</b> = 兩群甲基有顯著差異；粗細 ∝ |Δβ|，紅=前者較高、藍=後者較高，" +
            "已做 999 次置換 + BH-FDR。<b>共測 " + nTested + " 對，" + nSig + " 對顯著</b>。<br>" +
            "⚠ 沒畫出來的群不是「沒差異」，是<b>沒出現在任何顯著配對上</b>；" +
            "HP3（unphased somatic）<b>完全不進這張圖</b>。</div>";
    };
})();
