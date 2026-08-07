/* 位點詳情面板：點一個 sSNV 之後看到的東西。
 *
 * 分片是惰性載入的 —— 點到哪條染色體才抓那片 L2。載入失敗必須看得見，
 * 留白畫面會讓人以為「這個位點就是沒資料」。
 */
(function () {
    "use strict";
    var DD = window.__DD;
    if (!DD) return;

    var loading = {};   // chrom -> true 表示正在載

    function shardReady(chrom) {
        return DD.L2 && DD.L2[chrom];
    }

    function ensureShard(chrom, cb) {
        if (shardReady(chrom)) { cb(true); return; }
        if (loading[chrom]) { setTimeout(function () { ensureShard(chrom, cb); }, 120); return; }
        loading[chrom] = true;
        DD.loadShard("data/L2." + chrom + ".js", function (ok) {
            loading[chrom] = false;
            cb(ok && shardReady(chrom));
        });
    }

    /* 分片是逐染色體的陣列，但 CSR 存的是全域 region 索引。
       用 region_id 反查而非位置索引，避免兩邊編號體系不一致。 */
    function findRegions(chrom, pos) {
        var rows = (DD.L2 && DD.L2[chrom]) || [];
        var out = [];
        for (var i = 0; i < rows.length; i++) {
            var ap = rows[i].ap || [];
            for (var j = 0; j < ap.length; j++) {
                if (ap[j] === pos) { out.push(rows[i]); break; }
            }
        }
        return out;
    }

    function afTable(row, focus) {
        var af = row.af || [];
        if (!af.length) return "<div class='denom'>此 region 沒有 af_coverage。</div>";
        var rows = af.map(function (a) {
            var tot = a[1] + a[2];
            var pct = tot ? (a[1] / tot * 100).toFixed(1) + "%" : "—";
            var hit = a[0] === focus;
            return "<tr" + (hit ? " style='background:#e4efe9'" : "") + ">" +
                "<td class='num'>" + DD.fmt(a[0]) + (hit ? " ◀" : "") + "</td>" +
                "<td class='num'>" + DD.fmt(a[1]) + "</td>" +
                "<td class='num'>" + DD.fmt(a[2]) + "</td>" +
                "<td class='num'>" + pct + "</td></tr>";
        }).join("");
        return "<div class='table-wrap'><table class='data'>" +
            "<thead><tr><th class='num'>位置</th><th class='num'>ALT</th>" +
            "<th class='num'>REF</th><th class='num'>ALT 比例</th></tr></thead>" +
            "<tbody>" + rows + "</tbody></table></div>" +
            "<div class='denom'>ALT 比例的分母是該位點的 ALT+REF read 數（exact-PS 投影後），" +
            "<b>不是 VAF、不是 CCF</b>。</div>";
    }

    /* Read state matrix：A=看到 ALT、R=看到 REF、X=該位點未被此 read 覆蓋。
       full = 跨完整個 block 的 read；partial = 只覆蓋一部分（就是 X 的來源）。 */
    function readStateMatrix(row) {
        var pat = row.pat;
        if (!pat) {
            return "<div class='capability-off'><b>不可用。</b>" +
                "缺 MLHP 層（<code>populations_by_hp</code> / <code>subread_groups_by_hp</code>），" +
                "無法得知哪些 read 觀察到哪個 pattern。</div>";
        }
        var entries = [];
        Object.keys(pat.full || {}).forEach(function (hp) {
            Object.keys(pat.full[hp]).forEach(function (p) {
                entries.push({ src: "完整覆蓋", hp: hp, pattern: p, n: pat.full[hp][p] });
            });
        });
        Object.keys(pat.part || {}).forEach(function (hp) {
            Object.keys(pat.part[hp]).forEach(function (p) {
                entries.push({ src: "部分覆蓋", hp: hp, pattern: p, n: pat.part[hp][p] });
            });
        });
        if (!entries.length) {
            return "<div class='capability-off'>此 region 沒有通過門檻的 pattern。</div>";
        }
        entries.sort(function (a, b) { return b.n - a.n; });
        var pos = (row.ap || []).slice().sort(function (a, b) { return a - b; });
        var head = pos.map(function (p, i) {
            return "<th title='" + DD.fmt(p) + "'>S" + (i + 1) + "</th>";
        }).join("");
        var body = entries.map(function (e) {
            var cells = String(e.pattern).split("").map(function (ch) {
                return "<td class='state-cell " + DD.esc(ch) + "'>" + DD.esc(ch) + "</td>";
            }).join("");
            return "<tr><td>" + DD.esc(e.src) + "</td><td>HP" + DD.esc(e.hp) + "</td>" +
                cells + "<td class='num'>" + DD.fmt(e.n) + "</td></tr>";
        }).join("");
        var tot = entries.reduce(function (a, e) { return a + e.n; }, 0);
        return "<div class='table-wrap'><table class='state-matrix'>" +
            "<thead><tr><th>來源</th><th>HP</th>" + head + "<th>reads</th></tr></thead>" +
            "<tbody>" + body + "</tbody></table></div>" +
            "<div class='denom' style='margin-top:.35rem'>" +
            "<b class='state-cell A' style='padding:0 .25rem'>A</b> 看到 ALT　" +
            "<b class='state-cell R' style='padding:0 .25rem'>R</b> 看到 REF　" +
            "<b class='state-cell X' style='padding:0 .25rem'>X</b> 此 read 未覆蓋該位點<br>" +
            "共 " + DD.fmt(tot) + " 條 read×pattern；完整覆蓋 " + DD.fmt(pat.nFull || 0) + " 條。" +
            "<b>X 是覆蓋不足，不是資料品質差</b> —— ONT read 沒跨完整個 block 就會有 X。</div>";
    }

    function methylGraphBlock(chrom, pos) {
        var rec = (DD.L4 && DD.L4[chrom]) ? DD.L4[chrom][String(pos)] : null;
        if (!rec) {
            return "<div class='capability-off'>此位點沒有 ISM 資料，無法畫群關係圖。</div>";
        }
        return DD.methylGraph ? DD.methylGraph(rec) :
            "<div class='capability-off'>甲基群圖模組未載入。</div>";
    }

    function statusRow(row) {
        var cells = [
            ["unit_status", row.us], ["family_status", row.fs],
            ["objective_status", row.os], ["read_af_status", row.rs]
        ].map(function (p) {
            return "<div><div class='denom'>" + DD.esc(p[0]) + "</div><code>" +
                DD.esc(p[1] || "—") + "</code></div>";
        }).join("");
        return "<div class='grid-2' style='grid-template-columns:repeat(auto-fit,minmax(150px,1fr));" +
            "gap:.5rem;margin:.5rem 0'>" + cells + "</div>";
    }

    function treeHeader(row) {
        var tie = row.tu
            ? "<span class='tag' style='color:var(--accent)'>唯一最佳樹</span>"
            : "<span class='tag' style='color:var(--amber)'>同分並列 " + DD.esc(row.tc || "?") + " 棵</span>";
        return "<div style='display:flex;gap:.4rem;align-items:center;flex-wrap:wrap;margin-bottom:.35rem'>" +
            "<b>代表樹</b>" + tie +
            "<span class='tag' style='color:var(--muted)'>候選樹總數 " + DD.esc(row.tt || "?") + "</span>" +
            "<span class='src src-topology'>◆ topology</span></div>" +
            (row.tu ? "" :
             "<div class='callout warn' style='margin:.3rem 0'>這個 region 有 " +
             DD.esc(row.tc || "?") + " 棵同分最佳樹。下方畫的是其中<b>一棵代表</b>，" +
             "不是唯一解 —— 突變順序在這裡<b>不可當結論</b>。</div>");
    }

    function neighbourhood(row, all) {
        // 四類鄰域連結，分開標明各自的證據強度
        var sameUnit = [], samePS = [], nearby = [];
        for (var i = 0; i < all.length; i++) {
            var r = all[i];
            if (r.id === row.id) continue;
            if (r.u === row.u) sameUnit.push(r);
            else if (r.ps === row.ps) samePS.push(r);
            else if (row.ap.length && r.ap.length &&
                     Math.abs(r.ap[0] - row.ap[0]) <= 100000) nearby.push(r);
        }
        function block(title, list, color, note) {
            if (!list.length) return "";
            return "<div style='border-left:3px solid " + color + ";padding-left:.6rem;margin:.45rem 0'>" +
                "<b style='font-size:.85rem'>" + title + "（" + list.length + "）</b>" +
                "<div class='denom'>" + note + "</div>" +
                list.slice(0, 8).map(function (r) {
                    return "<div class='mono' style='font-size:.76rem'>" + DD.esc(r.id) +
                        " · k=" + r.k + " · " + (r.tu ? "唯一解" : "多解") + "</div>";
                }).join("") +
                (list.length > 8 ? "<div class='denom'>另有 " + (list.length - 8) + " 筆</div>" : "") +
                "</div>";
        }
        var html = block("同一 unit 的其他 block", sameUnit, "var(--accent)",
                         "同一個 strict read-linkage 單元被切成多塊 —— 有 read 層證據") +
                   block("同一 phase set 的其他 HP", samePS, "var(--accent2)",
                         "同一個 PS 的另一條單倍型 —— 有 read 層證據") +
                   block("座標鄰近 ±100 kb 的其他 region", nearby, "var(--muted)",
                         "⚠ 僅空間鄰近，<b>沒有 read 連結證據</b>，不可推論同源");
        return html || "<div class='denom'>此 region 在 ±100 kb 內沒有其他 region。</div>";
    }

    DD.renderDetail = function (i, _regionIdx) {
        var host = document.getElementById("detail");
        var chrom = DD.CHROMS[DD.chromOf[i]], pos = DD.posOf[i];
        host.innerHTML = "<div class='empty'>載入 " + DD.esc(chrom) + " 分片…</div>";

        ensureShard(chrom, function (ok) {
            if (!ok) {
                host.innerHTML = "<div class='capability-off'><b>分片載入失敗。</b>" +
                    "無法顯示 <code>" + DD.esc(chrom) + ":" + DD.fmt(pos) + "</code> 的區域詳情。" +
                    "請確認 <code>data/L2." + DD.esc(chrom) + ".js</code> 與 index.html 在同一資料夾。</div>";
                return;
            }
            var regions = findRegions(chrom, pos);
            if (!regions.length) {
                host.innerHTML = "<div class='capability-off'>分片已載入，但在 " +
                    DD.esc(chrom) + " 找不到含 " + DD.fmt(pos) + " 的 region。" +
                    "這代表 L1 索引與 L2 分片不同步 —— 屬於建置錯誤，請重跑 generator。</div>";
                return;
            }
            var all = DD.L2[chrom];
            var idx = 0;

            function paint() {
                var row = regions[idx];
                var multi = regions.length > 1
                    ? "<div class='callout warn'>這個 sSNV 落在 <b>" + regions.length +
                      " 個 region</b>。以下顯示第 " + (idx + 1) + " 個；" +
                      "<b>不會替你挑</b>，請用下方按鈕切換。<div style='margin-top:.35rem'>" +
                      regions.map(function (r, j) {
                          return "<button type='button' class='tag' data-j='" + j + "' " +
                              "style='cursor:pointer;color:" + (j === idx ? "var(--accent)" : "var(--muted)") +
                              ";margin-right:.3rem'>" + DD.esc(r.id.split("|").slice(1).join("|")) + "</button>";
                      }).join("") + "</div></div>"
                    : "";

                host.innerHTML =
                    "<div style='display:flex;gap:.5rem;align-items:baseline;flex-wrap:wrap'>" +
                    "<h3 class='mono'>" + DD.esc(chrom) + ":" + DD.fmt(pos) + "</h3>" +
                    "<span class='tag' style='color:var(--muted)'>k=" + row.k + "</span>" +
                    "<span class='tag' style='color:var(--muted)'>PS=" + DD.esc(row.ps) +
                    " HP=" + DD.esc(row.hp) + "</span></div>" +
                    "<div class='denom mono' style='word-break:break-all'>" + DD.esc(row.id) + "</div>" +
                    multi + statusRow(row) +
                    "<details open><summary>演化分支（代表樹）</summary><div class='details-body'>" +
                    treeHeader(row) + DD.drawTree(row, { edgeMode: "repOnly" }) +
                    "</div></details>" +
                    "<details open><summary>Locus 排列 — 位點間距與各自的 ALT 比例</summary>" +
                    "<div class='details-body'>" +
                    (DD.locusStrip ? DD.locusStrip(row, pos) : "") +
                    "<details><summary>同一份資料的表格版</summary><div class='details-body'>" +
                    afTable(row, pos) + "</div></details></div></details>" +
                    "<details open><summary>Read state matrix — 哪些 read 看到哪個 R/A/X pattern</summary>" +
                    "<div class='details-body'>" + readStateMatrix(row) + "</div></details>" +
                    "<details><summary>鄰域連結</summary><div class='details-body'>" +
                    neighbourhood(row, all) + "</div></details>" +
                    "<details open><summary>甲基分群 — 沿哪個軸？</summary><div class='details-body'>" +
                    (DD.methylPanel ? DD.methylPanel(chrom, pos) :
                     "<div class='capability-off'>甲基模組未載入。</div>") +
                    "</div></details>" +
                    "<details open><summary>甲基群關係 — 哪個群內有幾群、哪兩群之間有差異</summary>" +
                    "<div class='details-body'>" + methylGraphBlock(chrom, pos) + "</div></details>" +
                    "<details open><summary>甲基矩陣圖 — read × CpG 與 read × read 距離</summary>" +
                    "<div class='details-body'><div id='methylFig'></div></div></details>";

                if (DD.methylFigure) {
                    var fh = host.querySelector("#methylFig");
                    if (fh) DD.methylFigure(chrom, pos, fh);
                }

                host.querySelectorAll("button[data-j]").forEach(function (b) {
                    b.onclick = function () { idx = +b.dataset.j; paint(); };
                });
            }
            paint();
        });
    };

    // 全基因組圖上「有 ISM 資料的位點加圈」需要這支；P4 接上前一律回 false
    DD.hasISM = DD.hasISM || function () { return false; };
})();
