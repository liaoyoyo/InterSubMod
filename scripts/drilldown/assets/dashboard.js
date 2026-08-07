/* 多層下鑽儀表板 — 互動核心。
 *
 * 三條設計原則，違反任一條就會誤導判讀：
 *
 *   1. 預設全顯示。每個篩選維度初始「全部勾選」，所以不篩掉任何東西。
 *      分母橫幅永遠在場，一旦有維度非全選就變色並列出動了哪些維度。
 *
 *   2. 圖層與篩選正交。圖層只改畫法、不改 filtered 集合；UI 用不同控件形狀
 *      （篩選=checkbox、圖層=pill toggle）避免混淆「藏起來」與「篩掉」。
 *
 *   3. 命中不取最近。一個像素常疊了多個 sSNV（chr1 在 1400px 寬下是 178 kb/px），
 *      取最近會靜默丟掉其他候選，所以命中 >1 時彈消歧清單讓使用者自己選。
 */
(function () {
    "use strict";

    var DD = window.__DD || (window.__DD = {});
    var boot = JSON.parse(document.getElementById("bootData").textContent);
    var L1 = boot.l1;
    var CAPS = boot.capabilities;
    var DIMS = boot.dims;

    /* ── L1 解碼：delta → 絕對座標 ─────────────────────────────── */
    var N = L1.n;
    var chromOf = L1.chrom;
    var kOf = L1.k;
    var flagOf = L1.flags;
    var AXIS = L1.axis || [];
    var posOf = new Int32Array(N);
    (function () {
        var prev = -1, acc = 0;
        for (var i = 0; i < N; i++) {
            if (chromOf[i] !== prev) { acc = L1.dpos[i]; prev = chromOf[i]; }
            else { acc += L1.dpos[i]; }
            posOf[i] = acc;
        }
    })();

    var CHROMS = L1.chroms, CHROMLEN = L1.chromLen;
    var MAXLEN = Math.max.apply(null, CHROMLEN);

    function regionsOf(i) {
        var a = L1.csrOff[i], b = L1.csrOff[i + 1], out = [];
        for (var j = a; j < b; j++) out.push(L1.csrVal[j]);
        return out;
    }
    function hasFlag(i, bit) { return (flagOf[i] >> bit) & 1; }

    /* ── 狀態 ───────────────────────────────────────────────────── */
    // 與 payload.py 的 AXIS_BITS 一致：1 HP / 2 ALT / 4 lineage / 8 有測但不顯著
    var AXIS_COLOR = { 1: "#285f8f", 2: "#a94336", 4: "#6b5592", 8: "#c9c8bd" };

    var state = {
        view: "genome",
        filters: {},          // dimId -> {mode:"include"|"exclude", on:Set}
        layers: {},           // layerId -> bool
        filtered: null,       // Int32Array of L1 index
        cur: -1,              // 目前選中的 sSNV index
        curRegion: -1,
        zoomChrom: -1         // -1 = 全基因組；否則單染色體模式
    };

    // 初始：每個維度全勾 → 不篩掉任何東西
    DIMS.forEach(function (d) {
        state.filters[d.id] = { mode: "include", on: new Set(d.keys.map(function (k) { return k.v; })) };
    });
    boot.layers.forEach(function (l) { state.layers[l.id] = l.on !== false; });

    /* ── 維度取值 ───────────────────────────────────────────────── */
    var VALUE_FN = {
        k: function (i) { return String(Math.min(kOf[i], 8)); },
        chrom: function (i) { return CHROMS[chromOf[i]]; },
        unique: function (i) { return hasFlag(i, 1) ? "all" : (hasFlag(i, 0) ? "some" : "none"); },
        multiRegion: function (i) { return hasFlag(i, 3) ? "yes" : "no"; },
        ranked: function (i) { return hasFlag(i, 2) ? "yes" : "no"; },
        axis: function (i) {
            var c = AXIS[i] || 0;
            if (c === 0) return "0";
            if (c === 8) return "8";
            // 多軸同時顯著要單獨成一類，否則會被算進每一個單軸而重複計數
            var bits = (c & 1 ? 1 : 0) + (c & 2 ? 1 : 0) + (c & 4 ? 1 : 0);
            return bits > 1 ? "multi" : String(c);
        }
    };
    // drop-in 註釋的取值是 Python 端算好的逐 sSNV 陣列
    var ANNOT = boot.annotValues || {};

    function valueOf(dimId, i) {
        var f = VALUE_FN[dimId];
        if (f) return f(i);
        var a = ANNOT[dimId];
        if (a) return a[i] || "no";
        var ext = DD.valueFns && DD.valueFns[dimId];
        return ext ? ext(i) : "";
    }

    function passes(i, skipDim) {
        for (var d = 0; d < DIMS.length; d++) {
            var dim = DIMS[d];
            if (dim.id === skipDim || dim.disabled) continue;
            var f = state.filters[dim.id];
            var has = f.on.has(valueOf(dim.id, i));
            if (f.mode === "include" ? !has : has) return false;
        }
        return true;
    }

    function recompute() {
        var out = [];
        for (var i = 0; i < N; i++) if (passes(i)) out.push(i);
        state.filtered = Int32Array.from(out);
    }

    function activeDims() {
        return DIMS.filter(function (d) {
            if (d.disabled) return false;
            var f = state.filters[d.id];
            return f.mode === "exclude" ? f.on.size > 0 : f.on.size !== d.keys.length;
        });
    }

    /* ── 分母橫幅 ───────────────────────────────────────────────── */
    function renderScope() {
        var el = document.getElementById("scope");
        var act = activeDims();
        var n = state.filtered.length;
        el.className = "scope" + (act.length ? " filtered" : "");
        var parts = ["<span>顯示 <span class='n'>" + fmt(n) + "</span> / " + fmt(N) +
                     " <span class='ub ub-ssnv'>sSNV</span></span>"];
        if (act.length) {
            parts.push("<span class='dims'>已套用 " + act.length + " 個維度：" +
                act.map(function (d) {
                    var f = state.filters[d.id];
                    return esc(d.title) + " " + (f.mode === "exclude" ? "排除" : "") +
                           f.on.size + "/" + d.keys.length;
                }).join("、") + "</span>");
            parts.push("<button type='button' id='resetFilters'>全部還原</button>");
        } else {
            parts.push("<span class='dims'>未套用任何篩選 —— 預設全顯示</span>");
        }
        el.innerHTML = parts.join("");
        var r = document.getElementById("resetFilters");
        if (r) r.onclick = function () {
            DIMS.forEach(function (d) {
                state.filters[d.id] = { mode: "include", on: new Set(d.keys.map(function (k) { return k.v; })) };
            });
            refresh();
        };
    }

    /* ── 篩選 UI ────────────────────────────────────────────────── */
    function renderDims() {
        var host = document.getElementById("dims");
        host.innerHTML = DIMS.map(function (d) {
            var f = state.filters[d.id];
            // 每個 key 的 n 是「其他維度已篩選後」的數量 —— 讓使用者知道勾了會剩多少
            var counts = {};
            for (var x = 0; x < N; x++) {
                if (!passes(x, d.id)) continue;
                var v = valueOf(d.id, x);
                counts[v] = (counts[v] || 0) + 1;
            }
            var keys = d.keys.map(function (k) {
                var n = counts[k.v] || 0;
                return "<label class='key" + (n ? "" : " zero") + "'>" +
                    "<input type='checkbox' data-dim='" + esc(d.id) + "' data-key='" + esc(k.v) + "'" +
                    (f.on.has(k.v) ? " checked" : "") + (d.disabled ? " disabled" : "") + ">" +
                    (k.color ? "<i class='swatch' style='background:" + esc(k.color) + "'></i>" : "") +
                    "<span class='lab' title='" + esc(k.label || k.v) + "'>" + esc(k.label || k.v) + "</span>" +
                    "<span class='n'>" + fmt(n) + "</span></label>";
            }).join("");
            var off = d.disabled
                ? "<div class='denom' style='color:var(--danger)'>不可用 — " + esc(d.disabledReason || "缺對應資料層") + "</div>"
                : "";
            return "<div class='dim" + (d.disabled ? " disabled" : "") + "'>" +
                "<header><h3>" + esc(d.title) + "</h3>" +
                "<span class='src src-" + esc(d.src || "topology") + "'>◆ " + esc(d.srcLabel || "topology") + "</span>" +
                "<span class='tools'>" +
                "<button type='button' data-act='all' data-dim='" + esc(d.id) + "'>全選</button>" +
                "<button type='button' data-act='none' data-dim='" + esc(d.id) + "'>清空</button>" +
                "<button type='button' data-act='invert' data-dim='" + esc(d.id) + "'>反選</button>" +
                "<button type='button' data-act='mode' data-dim='" + esc(d.id) + "' aria-pressed='" +
                (f.mode === "exclude") + "' title='包含 ⇄ 排除'>排除</button>" +
                "</span></header>" + off +
                "<div class='keys'>" + keys + "</div></div>";
        }).join("");

        host.querySelectorAll("input[type=checkbox]").forEach(function (cb) {
            cb.onchange = function () {
                var f = state.filters[cb.dataset.dim];
                if (cb.checked) f.on.add(cb.dataset.key); else f.on.delete(cb.dataset.key);
                refresh();
            };
        });
        host.querySelectorAll("button[data-act]").forEach(function (b) {
            b.onclick = function () {
                var d = DIMS.find(function (x) { return x.id === b.dataset.dim; });
                var f = state.filters[d.id];
                if (b.dataset.act === "all") f.on = new Set(d.keys.map(function (k) { return k.v; }));
                else if (b.dataset.act === "none") f.on = new Set();
                else if (b.dataset.act === "invert") {
                    var nx = new Set();
                    d.keys.forEach(function (k) { if (!f.on.has(k.v)) nx.add(k.v); });
                    f.on = nx;
                } else if (b.dataset.act === "mode") {
                    f.mode = f.mode === "include" ? "exclude" : "include";
                }
                refresh();
            };
        });
    }

    /* ── 圖層 ───────────────────────────────────────────────────── */
    function renderLayers() {
        var host = document.getElementById("layers");
        host.innerHTML = boot.layers.map(function (l) {
            return "<label class='layer' data-on='" + (state.layers[l.id] ? 1 : 0) + "' title='" +
                esc(l.hint || "") + "'><input type='checkbox' data-layer='" + esc(l.id) + "'" +
                (state.layers[l.id] ? " checked" : "") + (l.disabled ? " disabled" : "") + ">" +
                esc(l.title) + "</label>";
        }).join("");
        host.querySelectorAll("input[data-layer]").forEach(function (cb) {
            cb.onchange = function () {
                state.layers[cb.dataset.layer] = cb.checked;
                cb.closest(".layer").dataset.on = cb.checked ? 1 : 0;
                drawGenome();            // 圖層只重畫，不動 filtered
            };
        });
    }

    /* ── 全基因組畫布 ───────────────────────────────────────────── */
    var cv = document.getElementById("genomeCanvas");
    var ctx = cv.getContext("2d");
    var geom = null;      // {left,top,rowH,plotW,scaleLen}
    var rowIdx = null;    // 每列的排序索引，供二分命中

    function buildRowIndex() {
        // filtered 已依 (chrom,pos) 排序（L1 本身就是排序的），所以逐列切片即可
        rowIdx = [];
        for (var c = 0; c < CHROMS.length; c++) rowIdx.push({ xs: [], ids: [] });
        var f = state.filtered;
        for (var j = 0; j < f.length; j++) {
            var i = f[j], c2 = chromOf[i];
            if (state.zoomChrom >= 0 && c2 !== state.zoomChrom) continue;
            rowIdx[c2].xs.push(posOf[i]);
            rowIdx[c2].ids.push(i);
        }
        for (var c3 = 0; c3 < rowIdx.length; c3++) {
            rowIdx[c3].xs = Float64Array.from(rowIdx[c3].xs);
            rowIdx[c3].ids = Int32Array.from(rowIdx[c3].ids);
        }
    }

    function drawGenome() {
        var dpr = Math.min(window.devicePixelRatio || 1, 2);
        var W = cv.clientWidth, H = cv.clientHeight;
        cv.width = W * dpr; cv.height = H * dpr;
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, W, H);

        var single = state.zoomChrom >= 0;
        var rows = single ? [state.zoomChrom] : CHROMS.map(function (_, i) { return i; });
        var left = 62, right = 64, top = 22, bottom = 26;   // right 需容納右側計數標籤，否則最長染色體的數字會被裁掉
        var plotW = W - left - right;
        var rowH = (H - top - bottom) / rows.length;
        var scaleLen = single ? CHROMLEN[state.zoomChrom] : MAXLEN;
        geom = { left: left, top: top, rowH: rowH, plotW: plotW, scaleLen: scaleLen, rows: rows };

        ctx.font = "11px " + getComputedStyle(document.body).getPropertyValue("--mono");
        ctx.textBaseline = "middle";

        for (var r = 0; r < rows.length; r++) {
            var c = rows[r], y = top + r * rowH + rowH / 2;
            var trackW = plotW * (CHROMLEN[c] / scaleLen);

            ctx.fillStyle = "#667069"; ctx.textAlign = "right";
            ctx.fillText(CHROMS[c], left - 8, y);

            ctx.fillStyle = "#e7e5da";
            ctx.fillRect(left, y - 6, trackW, 12);

            var xs = rowIdx[c].xs;

            // 密度帶：全顯示在擁擠處靠單點無法讀，底下補一條 1 Mb 分箱灰階
            if (state.layers.density && xs.length) {
                var BIN = 1e6, nb = Math.ceil(CHROMLEN[c] / BIN);
                var bins = new Uint32Array(nb), mx = 0;
                for (var b0 = 0; b0 < xs.length; b0++) {
                    var bi = Math.min(nb - 1, Math.floor(xs[b0] / BIN));
                    bins[bi]++; if (bins[bi] > mx) mx = bins[bi];
                }
                var bw = Math.max(plotW * (BIN / scaleLen), 0.7);
                for (var b = 0; b < nb; b++) {
                    if (!bins[b]) continue;
                    ctx.fillStyle = "rgba(102,112,105," + (0.10 + 0.5 * Math.sqrt(bins[b] / mx)).toFixed(3) + ")";
                    ctx.fillRect(left + plotW * (b * BIN / scaleLen), y + 7, bw, 4);
                }
            }

            // sSNV 點層
            var ids = rowIdx[c].ids;
            for (var i2 = 0; i2 < xs.length; i2++) {
                var x = left + plotW * (xs[i2] / scaleLen);
                var id = ids[i2];
                var uniq = hasFlag(id, 1);
                ctx.fillStyle = uniq ? "rgba(23,107,88,.72)" : "rgba(182,110,32,.72)";
                ctx.fillRect(x - 0.9, y - 5, 1.8, 10);
            }

            /* 甲基軸軌道：每條染色體下方一條細軌，顏色 = 該位點沿哪個軸顯著。
               先前用「每個點加藍圈」的做法在 19,849 點下糊成一片，
               而且回答不了「依據哪個指標分群」—— 顏色編碼兩者都解決。 */
            if (state.layers.axisTrack && AXIS.length) {
                for (var t = 0; t < xs.length; t++) {
                    var code = AXIS[ids[t]] || 0;
                    if (!code) continue;
                    var col = AXIS_COLOR[code];
                    if (!col) {
                        var b = (code & 1 ? 1 : 0) + (code & 2 ? 1 : 0) + (code & 4 ? 1 : 0);
                        col = b > 1 ? "#17231e" : AXIS_COLOR[8];
                    }
                    ctx.fillStyle = col;
                    ctx.fillRect(left + plotW * (xs[t] / scaleLen) - 0.85, y + 11, 1.7, 6);
                }
            }

            if (state.cur >= 0 && chromOf[state.cur] === c) {
                var cx = left + plotW * (posOf[state.cur] / scaleLen);
                ctx.strokeStyle = "#17231e"; ctx.lineWidth = 1.6;
                ctx.beginPath(); ctx.arc(cx, y, 7, 0, 6.284); ctx.stroke();
            }

            ctx.fillStyle = "#98a09a"; ctx.textAlign = "left";
            ctx.fillText(fmt(xs.length), left + trackW + 5, y);
        }

        document.getElementById("genomeHint").textContent =
            (single ? CHROMS[state.zoomChrom] + " 單染色體模式（點左側染色體名可返回全景）"
                    : "全基因組 · 22 條共用 chr1 的 bp 尺度") +
            " · 上排點：綠=唯一解 / 琥珀=非唯一" +
            (AXIS.length ? " · 下排細軌：藍=HP 顯著 / 紅=ALT 顯著 / 紫=lineage 顯著 / " +
             "黑=多軸 / 淺灰=有測但不顯著（不畫=無 ISM 資料）" : "") +
            " · 點擊選取，同一像素有多個時會列出候選";
    }

    /* 二分搜尋 + 雙向擴張，收集全部命中而非只取最近 */
    function hitsAt(px, py) {
        if (!geom) return [];
        var r = Math.floor((py - geom.top) / geom.rowH);
        if (r < 0 || r >= geom.rows.length) return [];
        var c = geom.rows[r];
        var yc = geom.top + r * geom.rowH + geom.rowH / 2;
        if (Math.abs(py - yc) > Math.min(geom.rowH / 2, 14)) return [];

        var xs = rowIdx[c].xs, ids = rowIdx[c].ids;
        if (!xs.length) return [];
        var R = 7;
        var bpPerPx = geom.scaleLen / geom.plotW;
        var target = (px - geom.left) * bpPerPx;
        var tol = R * bpPerPx;

        var lo = 0, hi = xs.length;
        while (lo < hi) { var mid = (lo + hi) >> 1; if (xs[mid] < target - tol) lo = mid + 1; else hi = mid; }
        var out = [];
        for (var j = lo; j < xs.length && xs[j] <= target + tol; j++) out.push(ids[j]);
        return out;
    }

    function onCanvasClick(ev) {
        var rect = cv.getBoundingClientRect();
        var px = ev.clientX - rect.left, py = ev.clientY - rect.top;

        if (px < geom.left - 4) {                      // 點染色體名 → 縮放切換
            var r = Math.floor((py - geom.top) / geom.rowH);
            if (r >= 0 && r < geom.rows.length) {
                state.zoomChrom = state.zoomChrom >= 0 ? -1 : geom.rows[r];
                buildRowIndex(); drawGenome(); syncHash();
            }
            return;
        }
        var hits = hitsAt(px, py);
        closePop();
        if (!hits.length) return;
        if (hits.length === 1) { select(hits[0]); return; }
        showPop(px, py, hits);
    }

    function showPop(px, py, hits) {
        var wrap = document.getElementById("genomeWrap");
        var el = document.createElement("div");
        el.className = "pop"; el.id = "hitPop";
        el.style.left = Math.min(px + 10, wrap.clientWidth - 280) + "px";
        el.style.top = Math.min(py + 10, wrap.clientHeight - 200) + "px";
        el.innerHTML = "<h4>此處有 " + hits.length + " 個 sSNV — 選一個（不替你挑最近的）</h4>" +
            hits.map(function (i) {
                var rs = regionsOf(i);
                return "<button type='button' data-i='" + i + "'>" +
                    "<div class='mono'>" + CHROMS[chromOf[i]] + ":" + fmt(posOf[i]) + "</div>" +
                    "<div class='meta'>k=" + kOf[i] + " · " + rs.length + " region · " +
                    (hasFlag(i, 1) ? "唯一解" : hasFlag(i, 0) ? "部分唯一" : "非唯一") + "</div></button>";
            }).join("");
        wrap.appendChild(el);
        el.querySelectorAll("button").forEach(function (b) {
            b.onclick = function () { select(+b.dataset.i); closePop(); };
        });
    }
    function closePop() { var p = document.getElementById("hitPop"); if (p) p.remove(); }

    /* ── 虛擬清單：不截斷 ───────────────────────────────────────── */
    var ROWH = 30;
    function renderList() {
        var box = document.getElementById("vlist");
        var inner = document.getElementById("vlistInner");
        var f = state.filtered;
        inner.style.height = (f.length * ROWH) + "px";
        var top = box.scrollTop, h = box.clientHeight;
        var a = Math.max(0, Math.floor(top / ROWH) - 12);
        var b = Math.min(f.length, Math.ceil((top + h) / ROWH) + 12);
        var html = [];
        for (var j = a; j < b; j++) {
            var i = f[j], rs = regionsOf(i);
            html.push("<div class='vrow' role='option' data-i='" + i + "' aria-selected='" +
                (i === state.cur) + "' style='top:" + (j * ROWH) + "px;height:" + ROWH + "px'>" +
                "<span class='pos'>" + CHROMS[chromOf[i]] + ":" + fmt(posOf[i]) + "</span>" +
                "<span class='tags'>" +
                "<span class='tag' style='color:var(--muted)'>k" + kOf[i] + "</span>" +
                (rs.length > 1 ? "<span class='tag' style='color:var(--purple)'>" + rs.length + "區</span>" : "") +
                "<span class='tag' style='color:" + (hasFlag(i, 1) ? "var(--accent)" : "var(--amber)") + "'>" +
                (hasFlag(i, 1) ? "唯一" : hasFlag(i, 0) ? "部分" : "多解") + "</span>" +
                "</span></div>");
        }
        inner.innerHTML = html.join("");
        inner.querySelectorAll(".vrow").forEach(function (el) {
            el.onclick = function () { select(+el.dataset.i); };
        });
        document.getElementById("listCount").textContent = fmt(f.length) + " 筆";
    }

    /* ── 選取與麵包屑 ───────────────────────────────────────────── */
    function select(i) {
        state.cur = i;
        var rs = regionsOf(i);
        state.curRegion = rs.length ? rs[0] : -1;
        drawGenome(); renderList(); renderCrumb();
        if (DD.renderDetail) DD.renderDetail(i, state.curRegion);
        syncHash();
    }

    function renderCrumb() {
        var el = document.getElementById("crumb");
        var parts = ["<button type='button' data-go='genome'>" + esc(boot.sample) + "</button>"];
        if (state.zoomChrom >= 0) {
            parts.push("<span class='caret'>›</span><button type='button' data-go='chrom'>" +
                CHROMS[state.zoomChrom] + "</button>");
        }
        if (state.cur >= 0) {
            var i = state.cur;
            if (state.zoomChrom < 0) {
                parts.push("<span class='caret'>›</span><button type='button' data-go='chrom' data-c='" +
                    chromOf[i] + "'>" + CHROMS[chromOf[i]] + "</button>");
            }
            parts.push("<span class='caret'>›</span><span>" + fmt(posOf[i]) + "</span>");
        }
        el.innerHTML = parts.join("");
        el.querySelectorAll("button").forEach(function (b) {
            b.onclick = function () {
                if (b.dataset.go === "genome") { state.zoomChrom = -1; state.cur = -1; }
                else if (b.dataset.go === "chrom") {
                    state.zoomChrom = b.dataset.c !== undefined ? +b.dataset.c : state.zoomChrom;
                }
                buildRowIndex(); drawGenome(); renderList(); renderCrumb(); syncHash();
            };
        });
    }

    /* ── URL hash：可貼、可書籤、支援上一頁 ─────────────────────── */
    var hashLock = false;
    function syncHash() {
        if (hashLock) return;
        var p = [];
        if (state.zoomChrom >= 0) p.push("c=" + CHROMS[state.zoomChrom]);
        if (state.cur >= 0) p.push("p=" + CHROMS[chromOf[state.cur]] + ":" + posOf[state.cur]);
        var act = activeDims();
        if (act.length) {
            p.push("f=" + act.map(function (d) {
                var f = state.filters[d.id];
                return d.id + (f.mode === "exclude" ? "!" : "~") + Array.from(f.on).join("|");
            }).join(";"));
        }
        history.replaceState(null, "", p.length ? "#" + p.join("&") : location.pathname);
    }
    function readHash() {
        var h = location.hash.replace(/^#/, "");
        if (!h) return;
        hashLock = true;
        h.split("&").forEach(function (kv) {
            var eq = kv.indexOf("="); if (eq < 0) return;
            var k = kv.slice(0, eq), v = kv.slice(eq + 1);
            if (k === "c") { var ci = CHROMS.indexOf(v); if (ci >= 0) state.zoomChrom = ci; }
            else if (k === "f") {
                v.split(";").forEach(function (spec) {
                    var m = spec.match(/^([^~!]+)([~!])(.*)$/); if (!m) return;
                    var f = state.filters[m[1]]; if (!f) return;
                    f.mode = m[2] === "!" ? "exclude" : "include";
                    f.on = new Set(m[3] ? m[3].split("|") : []);
                });
            }
        });
        hashLock = false;
        var pm = h.match(/(?:^|&)p=(chr[^:]+):(\d+)/);
        if (pm) {
            var ci2 = CHROMS.indexOf(pm[1]), want = +pm[2];
            for (var i = 0; i < N; i++) if (chromOf[i] === ci2 && posOf[i] === want) { state.cur = i; break; }
        }
    }

    /* ── 工具 ───────────────────────────────────────────────────── */
    function fmt(n) { return (n | 0).toLocaleString("en-US"); }
    function esc(s) {
        return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
            return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
        });
    }

    function refresh() {
        recompute(); buildRowIndex();
        renderScope(); renderDims(); drawGenome(); renderList(); renderCrumb();
        syncHash();
    }

    /* ── 對外介面：讓 P2+ 的模組掛進來而不用改這支 ───────────────── */
    DD.state = state;
    DD.CHROMS = CHROMS;
    DD.posOf = posOf; DD.chromOf = chromOf; DD.kOf = kOf;
    DD.regionsOf = regionsOf; DD.hasFlag = hasFlag;
    DD.fmt = fmt; DD.esc = esc;
    DD.refresh = refresh; DD.select = select;
    DD.valueFns = DD.valueFns || {};

    /* ── 啟動 ───────────────────────────────────────────────────── */
    readHash();
    renderLayers();
    refresh();
    // 初始選取必須延到所有模組載完 —— renderDetail 定義在 detail.js，
    // 那支在本檔之後才載入。setTimeout(0) 也不夠：瀏覽器允許在 <script> 標籤
    // 之間讓出，計時器可能早於 detail.js 執行就觸發（深連結因此靜默失效過）。
    // DOMContentLoaded 保證所有 inline script 都跑完了。
    if (state.cur >= 0) {
        var bootSelect = function () { select(state.cur); };
        if (document.readyState === "loading") {
            document.addEventListener("DOMContentLoaded", bootSelect);
        } else {
            setTimeout(bootSelect, 0);
        }
    }

    cv.addEventListener("click", onCanvasClick);
    document.getElementById("vlist").addEventListener("scroll", renderList);
    window.addEventListener("resize", function () { drawGenome(); renderList(); });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") closePop(); });

    // 分片載入失敗必須看得見 —— 留白畫面會讓人以為「這個樣本就是沒資料」
    DD.loadShard = function (url, onDone) {
        var s = document.createElement("script");
        s.src = url;
        s.onload = function () { if (onDone) onDone(true); };
        s.onerror = function () {
            var box = document.getElementById("shardError");
            if (box) {
                box.style.display = "";
                box.innerHTML = "<b>分片載入失敗：</b><code>" + esc(url) + "</code><br>" +
                    "請確認 index.html 與 data/ 在同一個資料夾。若是從壓縮檔解開，" +
                    "檢查 data/ 目錄有沒有一起解出來。";
            }
            if (onDone) onDone(false);
        };
        document.head.appendChild(s);
    };
})();
