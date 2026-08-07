/* 演化樹繪製。視覺文法沿用 layered_workstation 的 treeSvg / candidateUnionSvg，
 * 讓兩份產物看起來同一家族（顏色常數由 tests/test_drilldown_palette.py 鎖住）。
 *
 * 這裡把「三態」拆成節點與邊兩套，因為它們的語意完全不同：
 *
 *   節點三態  實測 / 補入 / 根
 *     實測 = label 不以 H_ 起頭，read 直接觀察到這個 pattern
 *     補入 = solver 為使樹合法而加的中間祖先，reads 沒有直接看到（全樣本 37.8%）
 *     根   = 全 REF 參考態，不是觀察對象
 *
 *   邊三態  必然 / 可替換 / 代表樹限定
 *     ⚠ 所有邊都是推論 —— 沒有任何一條邊被直接觀察到。
 *     必然     = 出現在全部 minimum trees
 *     可替換   = 只出現在部分
 *     代表樹限定 = 沒有候選 sidecar，只知道代表樹長這樣，無從判斷前兩者
 */
(function () {
    "use strict";
    var DD = window.__DD || (window.__DD = {});

    var C = {
        root: "#17231e",
        observedFill: "#fffdf7", observedStroke: "#176b58",
        hiddenFill: "#f5ead5", hiddenStroke: "#b66e20",
        edge: "#52625a",
        best: "#285f8f",
        methyl: "#d5663a",
        selected: "#176b58"
    };

    function isHidden(label) { return String(label || "").indexOf("H_") === 0; }

    /* 版面：depth → x，同層均分 y。代表樹實測最多 7 vertex / 6 edge，
       所以不需要通用 DAG layout 的那套複雜度。 */
    function layout(verts, edges, W, H) {
        var parent = {}, byId = {};
        verts.forEach(function (v) { byId[v[0]] = v; });
        edges.forEach(function (e) { parent[e[1]] = e[0]; });

        function depth(id) {
            var d = 0, seen = {}, cur = id;
            while (parent[cur] !== undefined && !seen[cur]) {
                seen[cur] = 1; cur = parent[cur]; d++;
                if (d > 64) break;
            }
            return d;
        }
        var levels = {};
        verts.forEach(function (v) {
            var d = depth(v[0]);
            (levels[d] = levels[d] || []).push(v);
        });
        var maxD = Math.max.apply(null, Object.keys(levels).map(Number));
        var pos = {};
        Object.keys(levels).forEach(function (dk) {
            var d = +dk, row = levels[d].slice().sort(function (a, b) { return a[0] - b[0]; });
            var x = 58 + (maxD ? d * ((W - 120) / maxD) : 0);
            row.forEach(function (v, i) {
                pos[v[0]] = { x: x, y: row.length === 1 ? H / 2 : 42 + i * ((H - 84) / (row.length - 1)) };
            });
        });
        return { pos: pos, byId: byId, maxD: maxD };
    }

    /* edgeClass(e) 回傳 'forced' | 'variable' | 'repOnly'
       目前沒有候選 sidecar，一律 repOnly；接上 candidates 能力後改由那層決定。 */
    function drawTree(row, opts) {
        opts = opts || {};
        var verts = row.v || [], edges = row.e || [];
        if (!verts.length) {
            return "<div class='capability-off'>此 region 沒有代表樹。" +
                "unit_status = <code>" + DD.esc(row.us || "?") + "</code> —— " +
                "無 active ALT / 分母為零 / 資源放棄 的 region 不會建樹。</div>";
        }
        var W = 520, H = Math.max(220, 46 * Math.max(2, verts.length));
        var L = layout(verts, edges, W, H);
        var edgeMode = opts.edgeMode || "repOnly";
        var dash = { forced: "", variable: "5 4", repOnly: "2 3" }[edgeMode];

        var svg = ["<svg viewBox='0 0 " + W + " " + H + "' width='100%' " +
                   "style='max-height:" + H + "px' role='img' aria-label='代表樹'>"];
        svg.push("<defs><marker id='ar' viewBox='0 0 10 10' refX='9' refY='5' " +
                 "markerWidth='6' markerHeight='6' orient='auto'>" +
                 "<path d='M0,0 L10,5 L0,10 z' fill='" + C.edge + "'/></marker></defs>");

        edges.forEach(function (e) {
            var a = L.pos[e[0]], b = L.pos[e[1]];
            if (!a || !b) return;
            svg.push("<line x1='" + a.x + "' y1='" + a.y + "' x2='" + b.x + "' y2='" + b.y +
                "' stroke='" + C.edge + "' stroke-width='2'" +
                (dash ? " stroke-dasharray='" + dash + "'" : "") +
                " marker-end='url(#ar)'><title>獲得 " + DD.fmt(e[2]) +
                "　edge score " + DD.esc(e[3] || "—") +
                "　此邊為推論（" + ({ forced: "出現在全部 minimum trees", variable: "只出現在部分 minimum trees",
                                    repOnly: "無候選 sidecar，無從判斷必然或可替換" }[edgeMode]) +
                "）</title></line>");
            var mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
            svg.push("<text x='" + mx + "' y='" + (my - 7) + "' text-anchor='middle' " +
                "font-size='9' font-family='ui-monospace,monospace' fill='#667069'>" +
                DD.fmt(e[2]) + "</text>");
        });

        verts.forEach(function (v) {
            var p = L.pos[v[0]]; if (!p) return;
            var id = v[0], label = v[1] || "";
            var root = id === 0, hid = !root && isHidden(label);
            var fill = root ? C.root : hid ? C.hiddenFill : C.observedFill;
            var stroke = root ? C.root : hid ? C.hiddenStroke : C.observedStroke;
            var title = root ? "根：全 REF 參考態，不是觀察對象"
                : hid ? "補入：solver 為使樹合法而加的中間祖先，reads 未直接觀察到"
                      : "實測：有 read 直接觀察到這個 pattern";
            svg.push("<circle cx='" + p.x + "' cy='" + p.y + "' r='18' fill='" + fill +
                "' stroke='" + stroke + "' stroke-width='2'" +
                (hid ? " stroke-dasharray='4 3'" : "") + "><title>" +
                DD.esc(label || "ROOT") + " — " + title + "</title></circle>");
            var txt = root ? "ROOT" : String(label).replace(/^H_/, "");
            svg.push("<text x='" + p.x + "' y='" + (p.y + 3.5) + "' text-anchor='middle' " +
                "font-size='" + (txt.length > 4 ? 8 : 9.5) + "' font-family='ui-monospace,monospace' " +
                "fill='" + (root ? "#fffdf7" : "#17231e") + "'>" + DD.esc(txt) + "</text>");
        });
        svg.push("</svg>");
        return svg.join("") + legend(edgeMode);
    }

    /* 圖例固定顯示在樹旁，不是 hover 才出現 —— 三態的意義不該要人猜。 */
    function legend(edgeMode) {
        var edgeNote = {
            forced: "實線 = 必然邊（出現在全部 minimum trees）",
            variable: "虛線 = 可替換邊（只出現在部分 minimum trees）",
            repOnly: "點線 = 代表樹限定（無候選 sidecar，無從區分必然／可替換）"
        }[edgeMode];
        return "<div class='denom' style='margin-top:.5rem;line-height:1.7'>" +
            "<b>節點</b>　" +
            "<span style='display:inline-block;width:10px;height:10px;border-radius:50%;background:" +
            C.observedFill + ";border:2px solid " + C.observedStroke + ";vertical-align:-1px'></span> 實測（read 直接看到）　" +
            "<span style='display:inline-block;width:10px;height:10px;border-radius:50%;background:" +
            C.hiddenFill + ";border:2px dashed " + C.hiddenStroke + ";vertical-align:-1px'></span> 補入（solver 推的中間祖先）　" +
            "<span style='display:inline-block;width:10px;height:10px;border-radius:50%;background:" +
            C.root + ";vertical-align:-1px'></span> 根（全 REF）<br>" +
            "<b>邊</b>　<b style='color:var(--danger)'>所有邊都是推論</b>，沒有任何一條被直接觀察到。" +
            edgeNote + "。邊上數字 = 獲得該突變的座標。</div>";
    }

    DD.drawTree = drawTree;
    DD.treeColors = C;
})();
