/* 甲基面板：回答「甲基有沒有分群、沿哪個軸分群」。
 *
 * 三個必須守住的判讀紀律：
 *
 *   1. 「未檢定」是第三態。*Valid=false、群數 <2、或該 run 根本沒算這個軸 ——
 *      三者都不是「不顯著」。用斜線底紋與獨立文案標示，絕不塗成白底。
 *
 *   2. 甲基不參與拓撲推論。這裡的一切都是事後 characterize；
 *      以甲基推 lineage 再用 lineage 解釋甲基會構成循環。
 *
 *   3. 母體不是全部 sSNV。只有通過 ISM 自己的 coverage / CpG gate 那批有資料，
 *      而那個 gate 與 CpG 密度、覆蓋深度相關 —— 不是隨機缺失。
 */
(function () {
    "use strict";
    var DD = window.__DD;
    if (!DD) return;

    var boot = JSON.parse(document.getElementById("bootData").textContent);
    var ISM = boot.ism;
    var ALPHA = 0.05;

    var STATE = { SIG: "顯著", NS: "不顯著", NT: "未檢定" };

    /* 三態判定。回傳 {state, why} —— why 說明為什麼未檢定，
       因為「未檢定」的原因不只一種，混在一起會失去資訊。 */
    function axisState(rec, ax) {
        if (!rec) return { state: "NT", why: "此位點不在 ISM 產物中" };
        var p = rec[ax.id + "_p"];
        var valid = rec[ax.id + "_valid"];
        var n = rec[ax.id + "_n"];
        if (valid === false) return { state: "NT", why: "檢定本身判定為不可用（*Valid=false）" };
        if (n !== null && n !== undefined && n < 2) {
            return { state: "NT", why: "分組數 " + n + " < 2，無從比較" };
        }
        if (p === null || p === undefined) return { state: "NT", why: "此 run 沒有算這個軸" };
        return { state: p <= ALPHA ? "SIG" : "NS", why: "" };
    }

    var COLOR = { SIG: "var(--accent)", NS: "var(--muted)", NT: "var(--muted)" };

    function fmtP(p) {
        if (p === null || p === undefined) return "—";
        if (p < 1e-4) return p.toExponential(1);
        return p.toFixed(4);
    }

    function axisTable(rec) {
        if (!ISM) {
            return "<div class='capability-off'><b>甲基層不可用。</b>" +
                "缺 ISM run，無法回答「甲基沿哪個軸分群」。</div>";
        }
        var rows = ISM.axes.map(function (ax) {
            var st = axisState(rec, ax);
            var p = rec ? rec[ax.id + "_p"] : null;
            var eff = rec ? rec[ax.id + "_eff"] : null;
            var n = rec ? rec[ax.id + "_n"] : null;
            var bg = st.state === "NT"
                ? "background:repeating-linear-gradient(45deg,transparent,transparent 5px," +
                  "rgba(102,112,105,.07) 5px,rgba(102,112,105,.07) 10px)"
                : (st.state === "SIG" ? "background:#e4efe9" : "");
            return "<tr style='" + bg + "'>" +
                "<td>" + DD.esc(ax.title) + "</td>" +
                "<td class='num'>" + fmtP(p) + "</td>" +
                "<td class='num'>" + (eff === null || eff === undefined ? "—" : eff.toFixed(3)) + "</td>" +
                "<td class='num'>" + (n === null || n === undefined ? "—" : n) + "</td>" +
                "<td style='color:" + COLOR[st.state] + ";white-space:nowrap'>" +
                (st.state === "NT" ? "<i>" + STATE[st.state] + "</i>" : STATE[st.state]) + "</td>" +
                "<td class='denom'>" + DD.esc(st.why) + "</td></tr>";
        }).join("");

        var missing = ISM.missingAxes && ISM.missingAxes.length
            ? "<div class='callout stop'>這個 ISM run 沒有算：<b>" +
              ISM.missingAxes.map(DD.esc).join("、") +
              "</b>。那些軸顯示為「未檢定」，<b>不等於不顯著</b>。</div>"
            : "";

        var sig = ISM.axes.filter(function (ax) { return axisState(rec, ax).state === "SIG"; });
        var tested = ISM.axes.filter(function (ax) { return axisState(rec, ax).state !== "NT"; });
        var verdict;
        if (!tested.length) {
            verdict = "<div class='callout warn'><b>此位點無任何軸可檢定。</b>" +
                "不能因此說「甲基沒有分群」—— 那是沒測，不是測了沒有。</div>";
        } else if (!sig.length) {
            verdict = "<div class='callout'>已檢定 " + tested.length + " 個軸，<b>沒有一個顯著</b>" +
                "（判準 p ≤ " + ALPHA + "）。在這個位點，甲基的分群與這些標籤軸都對不上。</div>";
        } else {
            verdict = "<div class='callout'>甲基分群與 <b>" +
                sig.map(function (a) { return DD.esc(a.title); }).join("、") +
                "</b> 對得上（p ≤ " + ALPHA + "；已檢定 " + tested.length + " 軸中的 " + sig.length + " 個）。" +
                "<b>這是關聯不是因果</b> —— 甲基沒有參與拓撲推論，這裡是事後 characterize。</div>";
        }

        return verdict + missing +
            "<div class='table-wrap'><table class='data'><thead><tr>" +
            "<th>分組軸</th><th class='num'>p</th><th class='num'>效果量</th>" +
            "<th class='num'>群數</th><th>判定</th><th>備註</th></tr></thead>" +
            "<tbody>" + rows + "</tbody></table></div>" +
            "<div class='denom' style='margin-top:.4rem'>" +
            "斜線底紋 = 未檢定（第三態，<b>不等於不顯著</b>）。" +
            "窗 ±" + DD.fmt(ISM.windowBp || 0) + " bp —— " +
            "窗大小決定 CpG 數，直接影響分群強度，<b>不同窗的 run 不可互比</b>。" +
            "</div>";
    }

    function readsPanel(rec) {
        if (!rec) return "";
        var n = rec.NumReads, cpg = rec.NumCpGs;
        var warn = String(rec.ClusterDispersionWarn || "").toLowerCase() === "true"
            ? "<div class='callout warn'>離散度警告：各群的組內離散程度差異大。" +
              "PERMANOVA 對此敏感，可能產生偽陽性 —— 甲基自身分群那一列的 p 值要保守看。</div>"
            : "";
        return "<div class='metrics' style='margin:.5rem 0'>" +
            "<div class='metric'><div class='label'>read 數</div><div class='value'>" +
            DD.esc(n || "—") + "</div><div class='note'>此位點視窗內</div></div>" +
            "<div class='metric'><div class='label'>CpG 數</div><div class='value'>" +
            DD.esc(cpg || "—") + "</div><div class='note'>甲基矩陣的欄數</div></div>" +
            "</div>" + warn;
    }

    /* 對外：給 detail.js 呼叫 */
    DD.methylPanel = function (chrom, pos) {
        var rec = (DD.L4 && DD.L4[chrom]) ? DD.L4[chrom][String(pos)] : null;
        if (!ISM) {
            return "<div class='capability-off'><b>甲基層不可用。</b>缺 ISM run。</div>";
        }
        if (!rec) {
            return "<div class='capability-off'><b>此位點沒有 ISM 資料。</b>" +
                "它不在 ISM 的產物清單中 —— 可能被 ISM 自己的 coverage / CpG gate 濾掉。" +
                "<br><b>這不代表此處甲基沒有分群</b>，只代表沒測。</div>";
        }
        return readsPanel(rec) + axisTable(rec);
    };

    // 全基因組圖的「有 ISM 資料加圈」
    if (ISM && ISM.hitDir) {
        var hit = new Set(ISM.hitDir);
        DD.hasISM = function (i) { return hit.has(i); };
    }

    // 篩選維度：甲基沿哪個軸顯著（由 L4 決定，需分片載入後才有值，
    // 所以這裡只註冊 valueFn，維度本身由 Python 端在 boot 決定要不要出現）
    DD.valueFns = DD.valueFns || {};
})();
