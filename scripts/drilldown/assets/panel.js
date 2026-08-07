/* 甲基雙面板圖的顯示層。
 *
 * PNG 本身是 1 px = 1 個資料格（那是唯一能在 genome-scale 撐住的做法 ——
 * inline SVG 每格一個 <rect> 實測是 1,533 KB/位點，PNG 只要 6–19 KB）。
 * 文字放大會糊，所以軸標與圖例走疊在上面的一層 SVG。
 */
(function () {
    "use strict";
    var DD = window.__DD;
    if (!DD) return;

    // boot 是 dashboard.js IIFE 內的區域變數，跨檔取不到 —— 重新從 DOM 解析
    // （曾直接寫 boot.palette 而整個詳情面板拋 "boot is not defined"）
    var boot = JSON.parse(document.getElementById("bootData").textContent);

    var loadingL5 = {};

    function ensureL5(chrom, cb) {
        if (DD.L5 && DD.L5[chrom]) { cb(true); return; }
        if (loadingL5[chrom]) { setTimeout(function () { ensureL5(chrom, cb); }, 120); return; }
        loadingL5[chrom] = true;
        DD.loadShard("data/L5." + chrom + ".js", function (ok) {
            loadingL5[chrom] = false;
            cb(ok && DD.L5 && !!DD.L5[chrom]);
        });
    }

    var BAR_COLOR = {
        HP: "HP1 淺藍 / HP1-1 深藍 / HP2 淺紫 / HP2-1 深紫 / HP3 teal / 未定 灰",
        ALT: "ALT 紅 / REF 琥珀 / 未知 淺灰",
        "T/N": "tumor 橘 / normal 綠",
        Strand: "+ 深灰 / − 淺灰（刻意非語意色）",
        lineage: "依階層標籤 HP2-1-1 上色（tree-aware）",
        cluster: "甲基自身分群（◆ 儀表板自算）"
    };

    function overlay(info) {
        /* 軸標用 HTML 百分比定位，不放進 SVG —— preserveAspectRatio='none'
           會把文字一起橫向拉伸（曾把「甲基 110 CpG」拉成整條寬）。 */
        var nb = info.bars.length, SB = info.barW, gap = info.gap, w = info.w;
        var pct = function (x) { return (x / w * 100).toFixed(3) + "%"; };
        var out = ["<div class='panel-marks'>"];
        out.push("<span style='left:" + pct(0) + "'>側欄</span>");
        out.push("<span style='left:" + pct(nb * SB) + "'>甲基 " + info.cpgs + " CpG</span>");
        if (info.hasDist) {
            out.push("<span style='left:" + pct(nb * SB + info.cpgs + gap) + "'>距離 " +
                     info.reads + "×" + info.reads + "</span>");
        }
        out.push("</div>");
        // 分隔線仍用 SVG（線可以被拉伸，不影響可讀性）
        var lines = ["<svg class='panel-ov' viewBox='0 0 " + w + " " + info.h +
                     "' preserveAspectRatio='none'>"];
        [nb * SB, nb * SB + info.cpgs].forEach(function (x) {
            lines.push("<line x1='" + x + "' y1='0' x2='" + x + "' y2='" + info.h +
                       "' stroke='rgba(255,253,247,.35)' stroke-width='.6'/>");
        });
        lines.push("</svg>");
        return lines.join("") + out.join("");
    }

    var PAL = boot.palette;

    function scaleBar(sc) {
        var stops = sc.stops.map(function (s) {
            return DD.esc(s.color) + " " + (s.t * 100).toFixed(0) + "%";
        }).join(",");
        return "<div class='scale'>" +
            "<div class='scale-title'>" + DD.esc(sc.title) + "</div>" +
            "<div class='scale-bar' style='background:linear-gradient(to right," + stops + ")'></div>" +
            "<div class='scale-ends'><span>" + DD.esc(sc.lo) + "</span>" +
            "<span>" + DD.esc(sc.hi) + "</span></div>" +
            "<div class='denom'><i class='swatch' style='background:" + DD.esc(sc.na) +
            "'></i>" + DD.esc(sc.naLabel) + "</div></div>";
    }

    function trackKey(t) {
        return "<div class='track-key'><b>" + DD.esc(t.title) + "</b>" +
            t.items.map(function (i) {
                return "<span><i class='swatch' style='background:" + DD.esc(i.color) +
                    "'></i>" + DD.esc(i.label) + "</span>";
            }).join("") + "</div>";
    }

    function legend(info) {
        var shown = {};
        info.bars.forEach(function (b) { shown[b] = 1; });

        var keys = PAL
            ? PAL.tracks.filter(function (t) { return shown[t.id]; }).map(trackKey).join("")
            : "<div class='denom' style='color:var(--danger)'>色盤未匯出，圖例降級 —— " +
              "generator 端 composite.palette_export() 失敗，請看 build log。</div>";
        var scales = PAL
            ? "<div class='scales'>" + PAL.scales.filter(function (s) {
                  return s.id !== "dist" || info.hasDist;
              }).map(scaleBar).join("") + "</div>"
            : "";
        var dend = (PAL && info.hasDendro)
            ? "<div class='track-key'><b>UPGMA 樹狀圖</b>" +
              "<span><i class='swatch' style='background:" + DD.esc(PAL.dendro.grey) +
              "'></i>跨群分支</span>" +
              "<span class='denom'>" + DD.esc(PAL.dendro.note) + "</span></div>"
            : "";

        var lin = info.lineageTotal
            ? "lineage 軌 join 到 <b>" + DD.fmt(info.lineageHit) + " / " +
              DD.fmt(info.lineageTotal) + "</b> 條 read（其餘塗灰 —— 那些 read 不在任何 lineage block）"
            : "";
        var clu = (info.clustered !== undefined)
            ? "UPGMA 分群涵蓋 <b>" + DD.fmt(info.clustered) + " / " + DD.fmt(info.reads) +
              "</b> 條 read，切成 <b>" + DD.esc(info.clusterK) + "</b> 群" +
              "（ISM 的 optimal_k = " + DD.esc(info.optimalK) + "）；" +
              "<b>未涵蓋的 " + DD.fmt(info.reads - info.clustered) +
              " 條排在最下方、無樹枝、cluster 軌塗灰</b> —— 它們不是被丟掉，是分群階段沒收進去"
            : "";

        var scope = info.tumorOnly === true
            ? "<b style='color:var(--danger)'>目前顯示：只看 tumor read（" +
              DD.fmt(info.reads) + " / " + DD.fmt(info.readsAll || info.reads) +
              " 條）</b> —— 排除了 normal read，germline ASM 的對照基線不在圖上。<br>"
            : (info.readsAll || info.tumorOnly)
                ? "<b>目前顯示：全部 read（含 normal）</b> —— " +
                  "normal 是 germline ASM 的對照基線；若只關心腫瘤可切到「只看 tumor」。<br>"
                : "";
        return "<div class='panel-legend'>" + keys + dend + scales + "</div>" +
            "<div class='denom' style='line-height:1.7'>" + scope +
            "<b>read 依 UPGMA 葉序排列</b>（leaf_order.txt），不是依 HP —— " +
            "依 HP 排看不出甲基自己怎麼分。<br>" +
            (lin ? lin + "<br>" : "") + (clu ? clu + "<br>" : "") +
            "側欄由左至右：" + DD.esc(info.bars.join(" │ ")) + "　" +
            "<span class='src src-derived'>◆ cluster 軌與樹狀圖著色為儀表板自算</span></div>";
    }

    DD.methylFigure = function (chrom, pos, host) {
        host.innerHTML = "<div class='denom'>載入面板…</div>";
        ensureL5(chrom, function (ok) {
            var info = (ok && DD.L5[chrom]) ? DD.L5[chrom][String(pos)] : null;
            if (!info) {
                host.innerHTML = "<div class='capability-off'><b>此位點沒有預產面板。</b>" +
                    "面板是 ISM 算完後由 <code>--bake-panels</code> 產生的；" +
                    "本次 build 沒有涵蓋這個位點（或該位點 read/CpG 數不足以合成）。" +
                    "<br>重跑 generator 時加 <code>--bake-panels all</code> 可涵蓋全部。</div>";
                return;
            }
            /* 保持原始長寬比。先前 width:100% 把 180×296 的圖橫向拉滿容器，
               長寬比全毀 —— 甲基矩陣的「一格 = 一個 read×CpG」語意也跟著失真。
               改成預設 fit（等比縮到容器內）+ 可切換放大倍率 + 拖移。 */
            var variants = { all: info };
            if (info.tumorOnly) variants.T = info.tumorOnly;
            var curKey = "all";

            function paint() {
            var cur = variants[curKey];
            var ar = (cur.w / cur.h).toFixed(4);
            var info = cur;                       // 以下沿用既有變數名
            host.innerHTML =
                "<div class='panel-tools'>" +
                (variants.T
                 ? "<button type='button' data-v='all' aria-pressed='" + (curKey === "all") +
                   "'>全部 read</button>" +
                   "<button type='button' data-v='T' aria-pressed='" + (curKey === "T") +
                   "'>只看 tumor</button><span class='denom'>│</span>"
                 : "") +
                "<span class='denom'>" + info.w + " × " + info.h + " px（1 px = 1 格）</span>" +
                "<button type='button' data-z='fit' aria-pressed='true'>符合寬度</button>" +
                "<button type='button' data-z='1' aria-pressed='false'>1×</button>" +
                "<button type='button' data-z='2' aria-pressed='false'>2×</button>" +
                "<button type='button' data-z='4' aria-pressed='false'>4×</button>" +
                "<button type='button' data-z='8' aria-pressed='false'>8×</button>" +
                "<span class='denom'>放大後可拖移</span></div>" +
                "<div class='panel-scroll' data-mode='fit'>" +
                "<div class='panel-fig' style='aspect-ratio:" + ar + "'>" +
                "<img src='" + DD.esc(info.file) + "' alt='甲基與距離雙面板' " +
                "width='" + info.w + "' height='" + info.h + "'>" +
                overlay(info) + "</div></div>" + legend(info);

            host.querySelectorAll("button[data-v]").forEach(function (b) {
                b.onclick = function () { curKey = b.dataset.v; paint(); };
            });

            var scroll = host.querySelector(".panel-scroll");
            var fig = host.querySelector(".panel-fig");
            host.querySelectorAll("button[data-z]").forEach(function (b) {
                b.onclick = function () {
                    host.querySelectorAll("button[data-z]").forEach(function (x) {
                        x.setAttribute("aria-pressed", String(x === b));
                    });
                    var z = b.dataset.z;
                    if (z === "fit") {
                        scroll.dataset.mode = "fit";
                        fig.style.width = "";
                    } else {
                        scroll.dataset.mode = "zoom";
                        fig.style.width = (info.w * parseInt(z, 10)) + "px";
                    }
                };
            });
            // 拖移（放大時捲動容器會超出視窗）
            var drag = null;
            scroll.addEventListener("mousedown", function (e) {
                if (scroll.dataset.mode !== "zoom") return;
                drag = { x: e.clientX, y: e.clientY, l: scroll.scrollLeft, t: scroll.scrollTop };
                scroll.style.cursor = "grabbing";
                e.preventDefault();
            });
            window.addEventListener("mousemove", function (e) {
                if (!drag) return;
                scroll.scrollLeft = drag.l - (e.clientX - drag.x);
                scroll.scrollTop = drag.t - (e.clientY - drag.y);
            });
            window.addEventListener("mouseup", function () {
                drag = null; scroll.style.cursor = "";
            });
            }
            paint();
        });
    };
})();
