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

    function legend(info) {
        var BAR_SWATCH = {
            HP: "#60a5fa", ALT: "#dc2626", "T/N": "#f97316",
            Strand: "#334155", lineage: "#db2777", cluster: "#0d9488"
        };
        var swatches = info.bars.map(function (b) {
            return "<span class='bar-key' title='" + DD.esc(BAR_COLOR[b] || "") + "'>" +
                "<i class='swatch' style='background:" + (BAR_SWATCH[b] || "#667069") + "'></i>" +
                DD.esc(b) + "</span>";
        }).join("");
        var lin = info.lineageTotal
            ? "lineage 軌 join 到 <b>" + DD.fmt(info.lineageHit) + " / " +
              DD.fmt(info.lineageTotal) + "</b> 條 read（其餘塗灰 —— 那些 read 不在任何 lineage block）"
            : "";
        var clu = (info.clustered !== undefined)
            ? "分群涵蓋 <b>" + DD.fmt(info.clustered) + " / " + DD.fmt(info.reads) +
              "</b> 條 read（未涵蓋者排在最下方、cluster 軌塗灰）；" +
              "切成 <b>" + DD.esc(info.clusterK) + "</b> 群（ISM 的 optimal_k = " +
              DD.esc(info.optimalK) + "）"
            : "";
        return "<div class='legend'>" + swatches + "</div>" +
            "<div class='denom' style='line-height:1.7'>" +
            "<b>read 依甲基距離的葉序排列</b>，不是依 HP —— 依 HP 排看不出甲基自己怎麼分。<br>" +
            "左半 = read × CpG 的 β（藍 0 → 白 .5 → 紅 1，灰 = 無覆蓋）；" +
            (info.hasDist ? "右半 = read × read 距離（暗 = 近、亮 = 遠）。<br>" : "無距離矩陣。<br>") +
            (lin ? lin + "<br>" : "") + (clu ? clu + "<br>" : "") +
            "側欄由左至右：" + DD.esc(info.bars.join(" │ ")) + "　" +
            "<span class='src src-derived'>◆ cluster 軌為儀表板自算</span></div>";
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
            host.innerHTML =
                "<div class='panel-fig'>" +
                "<img src='" + DD.esc(info.file) + "' alt='甲基與距離雙面板' " +
                "width='" + info.w + "' height='" + info.h + "'>" +
                overlay(info) + "</div>" + legend(info);
        });
    };
})();
