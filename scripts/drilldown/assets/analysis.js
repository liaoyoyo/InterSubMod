/* 共出現與自檢的前端渲染。
 *
 * 共出現矩陣的每一格都可以點 —— 點下去等同「把 X 維度限縮到該欄、Y 限縮到該列」，
 * 這就是「看維度之間的交集」的操作化。Shift+點 = 疊加到現有篩選。
 *
 * 所有比例在 n < 30 時抑制（NCHS 標準），只顯示原始計數。
 */
(function () {
    "use strict";
    var DD = window.__DD;
    if (!DD) return;
    var A = JSON.parse(document.getElementById("analysisData").textContent);

    /* Pearson 殘差 → 顏色。正殘差（超出預期）用綠、負用紅，
       絕對值 <2 視為雜訊不上色 —— 免得每一格都有顏色反而看不出重點。 */
    function residColor(r) {
        if (r === null || Math.abs(r) < 2) return "transparent";
        var t = Math.min(Math.abs(r) / 8, 1);
        return r > 0 ? "rgba(23,107,88," + (0.10 + 0.55 * t).toFixed(2) + ")"
                     : "rgba(169,67,54," + (0.10 + 0.55 * t).toFixed(2) + ")";
    }

    function cellHtml(c) {
        var pct = c.small ? "<span class='denom'>n&lt;30</span>"
            : (c.colPct === null ? "" : "<span class='denom'>" + c.colPct.toFixed(1) + "%</span>");
        var res = (c.small || c.resid === null) ? "" :
            "<span class='denom' style='opacity:.8'>r " + (c.resid > 0 ? "+" : "") + c.resid + "</span>";
        return "<div style='font-family:var(--mono);font-weight:600'>" + DD.fmt(c.n) + "</div>" +
               pct + (res ? " " + res : "");
    }

    function renderMatrix(spec) {
        var m = spec.matrix;
        var head = "<tr><th></th>" + m.xs.map(function (x, i) {
            return "<th class='num'>" + DD.esc(x) + "<div class='denom'>" + DD.fmt(m.colsum[i]) + "</div></th>";
        }).join("") + "</tr>";
        var body = m.ys.map(function (y, r) {
            return "<tr><th style='text-align:left;white-space:nowrap'>" + DD.esc(y) +
                "<div class='denom'>" + DD.fmt(m.rowsum[r]) + "</div></th>" +
                m.cells[r].map(function (c, k) {
                    return "<td class='num' style='background:" + residColor(c.resid) +
                        ";cursor:" + (c.n ? "pointer" : "default") + "'" +
                        (c.n ? " data-x='" + DD.esc(m.xs[k]) + "' data-y='" + DD.esc(y) +
                               "' data-spec='" + DD.esc(spec.id) + "'" : "") +
                        " title='" + DD.esc(m.xs[k] + " × " + y) + "：n=" + c.n +
                        (c.resid !== null && !c.small ? "，Pearson 殘差 " + c.resid : "") + "'>" +
                        cellHtml(c) + "</td>";
                }).join("") + "</tr>";
        }).join("");
        return "<div class='table-wrap'><table class='data'><thead>" + head +
            "</thead><tbody>" + body + "</tbody></table></div>";
    }

    function renderBars(spec) {
        var mx = Math.max.apply(null, spec.rows.map(function (r) { return r.n; }));
        return "<div class='table-wrap'><table class='data'>" +
            "<thead><tr><th>k</th><th class='num'>region</th><th>分布</th>" +
            "<th class='num'>唯一解</th><th class='num'>唯一解率</th></tr></thead><tbody>" +
            spec.rows.map(function (r) {
                var w = (r.n / mx * 100).toFixed(2);
                var rate = r.rate === null
                    ? "<span class='denom'>n&lt;30 不報比率</span>"
                    : "<b>" + r.rate.toFixed(1) + "%</b>";
                return "<tr><td class='num'>" + r.k + "</td>" +
                    "<td class='num'>" + DD.fmt(r.n) + "</td>" +
                    "<td><span style='display:block;height:9px;background:var(--soft)'>" +
                    "<span style='display:block;height:100%;width:" + w + "%;background:var(--accent2)'></span>" +
                    "</span></td>" +
                    "<td class='num'>" + DD.fmt(r.unique) + "</td>" +
                    "<td class='num'>" + rate + "</td></tr>";
            }).join("") + "</tbody></table></div>";
    }

    function panel(spec) {
        var body = spec.kind === "bars" ? renderBars(spec) : renderMatrix(spec);
        return "<div class='section' style='box-shadow:none;margin-bottom:.9rem'>" +
            "<header><div><h3>" + DD.esc(spec.title) + "</h3>" +
            "<div class='denom'>" + DD.esc(spec.denom) + "</div></div>" +
            "<span class='src src-" + DD.esc(spec.src) + "'>◆ " + DD.esc(spec.srcLabel) + "</span></header>" +
            body +
            (spec.note ? "<div class='callout'>" + spec.note + "</div>" : "") +
            (spec.caveat ? "<div class='denom' style='margin-top:.4rem'>" + spec.caveat + "</div>" : "") +
            "</div>";
    }

    function renderCooccur() {
        var host = document.getElementById("view-cooccur");
        var inner = host.querySelector(".capability-off");
        var html = "<div class='callout warn'><b>怎麼讀。</b>格子裡第一行是計數，" +
            "第二行是欄百分比，<code>r</code> 是 Pearson 殘差（觀察−期望）/√期望 —— " +
            "|r| ≥ 2 才上色，因為更小的偏離跟隨機沒兩樣。" +
            "<b>點任一格</b>可把該交集套成篩選。</div>";
        html += A.cooccur.map(panel).join("");
        if (inner) inner.outerHTML = html; else host.insertAdjacentHTML("beforeend", html);

        host.querySelectorAll("td[data-x]").forEach(function (td) {
            td.onclick = function () {
                var spec = A.cooccur.find(function (s) { return s.id === td.dataset.spec; });
                if (!spec || !spec.filterMap) return;
                var fm = spec.filterMap;
                if (fm.x) DD.setFilter(fm.x, [td.dataset.x]);
                if (fm.y) DD.setFilter(fm.y, [td.dataset.y]);
                document.querySelector(".seg button[data-view='genome']").click();
            };
        });
    }

    function renderSelfcheck() {
        var host = document.getElementById("view-selfcheck");
        var s = A.selfcheck.summary;
        var mark = { PASS: "✅", FAIL: "❌", SKIP: "⊘" };
        var color = { PASS: "var(--accent)", FAIL: "var(--danger)", SKIP: "var(--muted)" };

        var html = "<div class='metrics' style='margin-bottom:.9rem'>" +
            "<div class='metric'><div class='label'>守恆等式通過</div>" +
            "<div class='value' style='color:var(--accent)'>" + s.pass + "</div>" +
            "<div class='note'>分母 " + s.total + " 項</div></div>" +
            "<div class='metric'><div class='label'>不成立</div>" +
            "<div class='value' style='color:" + (s.fail ? "var(--danger)" : "var(--muted)") + "'>" +
            s.fail + "</div><div class='note'>非 0 即為資料問題</div></div>" +
            "<div class='metric'><div class='label'>無法檢查</div>" +
            "<div class='value' style='color:var(--muted)'>" + s.skip + "</div>" +
            "<div class='note'>缺對應資料層</div></div></div>";

        html += "<div class='callout warn'><b>「無法檢查」是第三態。</b>" +
            "該項需要的資料層不在，既不是通過也不是失敗 —— " +
            "把它算成通過會讓人以為驗過了。</div>";

        html += "<div class='table-wrap'><table class='data'><thead><tr>" +
            "<th>#</th><th>守恆等式</th><th>狀態</th><th>實際</th></tr></thead><tbody>" +
            A.selfcheck.checks.map(function (c) {
                return "<tr><td class='mono'>" + DD.esc(c.id) + "</td>" +
                    "<td>" + DD.esc(c.title) +
                    (c.need ? "<div class='denom'>需要能力 <code>" + DD.esc(c.need) + "</code></div>" : "") +
                    "</td>" +
                    "<td style='color:" + color[c.status] + ";white-space:nowrap'>" +
                    mark[c.status] + " " + DD.esc(c.statusLabel) + "</td>" +
                    "<td class='denom'>" + DD.esc(c.detail) + "</td></tr>";
            }).join("") + "</tbody></table></div>";

        html += "<h3 style='margin-top:1.1rem'>覆蓋率與缺口成因</h3>" +
            "<div class='denom'>只給比率會讓人以為缺的是隨機的。這裡逐項說明缺的是什麼。</div>" +
            "<div class='table-wrap' style='margin-top:.4rem'><table class='data'><thead><tr>" +
            "<th>項目</th><th class='num'>分子</th><th class='num'>分母</th>" +
            "<th class='num'>比率</th><th>缺的是什麼</th></tr></thead><tbody>" +
            A.selfcheck.coverage.map(function (c) {
                var r = c.den ? (c.num / c.den * 100).toFixed(1) + "%" : "—";
                return "<tr><td>" + DD.esc(c.title) + "</td>" +
                    "<td class='num'>" + DD.fmt(c.num) + "</td>" +
                    "<td class='num'>" + DD.fmt(c.den) + "</td>" +
                    "<td class='num'>" + r + "</td>" +
                    "<td class='denom'>" + DD.esc(c.gap) + "</td></tr>";
            }).join("") + "</tbody></table></div>";

        // 選擇偏誤前哨置頂於此區塊之上
        var sent = A.sentinel;
        var sentHtml = sent.available
            ? panel(sent)
            : "<div class='capability-off'><b>" + DD.esc(sent.title) + "</b><br>" +
              DD.esc(sent.reason) + "</div>";
        html = "<h3>選擇偏誤前哨</h3><div class='denom' style='margin-bottom:.4rem'>" +
            "所有甲基結論的前提檢查 —— 先看這張再看下面任何東西。</div>" +
            sentHtml + "<hr style='border:0;border-top:1px solid var(--line);margin:1.2rem 0'>" + html;

        var inner = host.querySelector(".capability-off");
        if (inner) inner.outerHTML = html; else host.insertAdjacentHTML("beforeend", html);
    }

    // 讓共出現的格子能反過來套篩選
    DD.setFilter = function (dimId, keys) {
        var f = DD.state.filters[dimId];
        if (!f) return;
        f.mode = "include";
        f.on = new Set(keys);
        DD.refresh();
    };

    renderCooccur();
    renderSelfcheck();
})();
