#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ISM case-heatmap DISPLAY STANDARD — single source of truth (SoT).

ALL workstation renderers import this so sidebars / palette / colormaps never drift.
Standard fixed sidebars (order, from reads.tsv): [cluster?] HP | ALT | T/N | Strand
  - cluster : only for clustering analyses (unsupervised group id); DISTINCT palette (teal/pink), never orange/green
  - HP      : hp           (HP1/HP1-1/HP2/HP2-1/HP3/unphase) germline-haplotype backbone, blue/purple family
  - ALT     : alt_support  (ALT carries somatic variant / REF) red / amber
  - T/N     : is_tumor      tumor orange / normal green  (omit when tumor-only = uniform)
  - Strand  : strand        + / -  (artifact check) slate
Colormaps: methylation beta = RdBu_r (0 blue / .5 white / 1 red / NaN grey);
           read×read distance = magma (dark=near / bright=far).
SNV marker: magenta (#ff00aa) line drawn BETWEEN flanking CpGs (variant is rarely a CpG).

Ref standard: memory feedback_ism_case_heatmap_standard_sidebars (2026-06-16) + this SoT (2026-06-18).
"""
import bisect
import colorsys
import numpy as np

# ----------------------------------------------------------------- PALETTE (SoT)
HP_COL = {'1': '#60a5fa', '1-1': '#1e3a8a',   # HP1 family light/dark BLUE
          '2': '#c084fc', '2-1': '#6b21a8',   # HP2 family light/dark PURPLE
          '3': '#14b8a6', '0': '#9ca3af', '': '#d1d5db'}  # HP3 teal / unphase grey
HP_NAME = {'1': 'HP1', '1-1': 'HP1-1(somatic)', '2': 'HP2', '2-1': 'HP2-1(somatic)',
           '3': 'HP3(unphase somatic)', '0': 'unphase'}
TN_COL = {'1': '#f97316', '0': '#22c55e'}                 # tumor orange / normal green
ALT_COL = {'ALT': '#dc2626', 'REF': '#fbbf24', '': '#e5e7eb', 'NA': '#e5e7eb', 'UNKNOWN': '#e5e7eb'}  # ALT red / REF amber / unknown grey
STRAND_COL = {'+': '#334155', '-': '#94a3b8', '': '#e5e7eb'}                     # slate (non-semantic hue)
# cluster ids: DISTINCT from every semantic bar (no orange/green/blue/purple/red/amber)
CLUSTER_COL = ['#db2777', '#0d9488', '#0891b2', '#65a30d', '#7c3aed', '#475569']
SNV_COL = '#facc15'        # SNV line/arrow: yellow, thin
SEP_COL = '#c0c0c0'        # cluster separator: silver
NA_GREY = '#d1d5db'
HP_ORDER = ['1', '1-1', '2', '2-1', '3', '0', '']

# ----------------------------------------------------------------- 演化樹 tree-aware cluster color (SoT)
# phylo coarse_label = 遞迴二元樹路徑 ('1'/'1-1'/'1-2'/'1-1-2'...). 取代扁平 CLUSTER_COL 當「演化樹群色」,
# 讓樹結構用顏色顯示且與 HP/ALT/T-N 脫離。配色三維(deterministic, 同 label 跨 region 固定):
#   HUE: 第一次分裂(segs[1]) 跳到「不同色相錨點」(明顯對比, root1 左=magenta/右=cyan); 更深=錨點鄰域微幅漂移(左−/右+逐層減半)→表親可分仍歸該支。
#   L  : 隨深度遞減(深=暗); 冗餘且 luminance 故色盲安全; 兄弟 -2 微亮。  S: 固定高飽和。  例外(other/outlier/edge/'-')=灰。
# 護欄: 錨點只用 cluster-reserved {magenta,rose,cyan},脫離 HP 藍(.58)/紫(.78)/teal(.47)+ALT 紅/REF 琥珀/T-N 橘綠。唯一近鄰 cyan↔HP3-teal(HP3罕見、不同欄)。
# ⚠ 對比天花板: 扣 HP+保留軸後真正明顯不同色相只剩 ~3 區(magenta/rose/cyan); 兩大支 magenta vs cyan=護欄內最大對比。詳 docs/references/20260623_ism_observation_panel_and_visual_convention_spec_01.md
CLUSTER_TREE_ROOT = {'1': 0.86, '2': 0.95, '3': 0.50}
CLUSTER_TREE_ANCHOR = {('1', '1'): 0.86, ('1', '2'): 0.50, ('2', '1'): 0.95, ('2', '2'): 0.78,
                       ('3', '1'): 0.50, ('3', '2'): 0.86}
CLUSTER_TREE_L = {1: 0.62, 2: 0.55, 3: 0.46, 4: 0.38, 5: 0.31}
CLUSTER_TREE_SAT = 0.68
CLUSTER_TREE_GREY = '#9ca3af'


def cluster_tree_color(label):
    """phylo coarse_label → hex (tree-aware, 取代扁平 CLUSTER_COL). other/outlier/edge/'-'/'' → 中性灰。"""
    if label in ('other', 'outlier', 'edge', '-', '', None):
        return CLUSTER_TREE_GREY
    segs = str(label).split('-')
    root = segs[0]
    if len(segs) == 1:
        hue = CLUSTER_TREE_ROOT.get(root, 0.86)
    else:
        hue = CLUSTER_TREE_ANCHOR.get((root, segs[1]), CLUSTER_TREE_ROOT.get(root, 0.86))
        w = 0.05
        for s in segs[2:]:
            w *= 0.5
            hue += w if s == '2' else -w
    depth = len(segs)
    L = CLUSTER_TREE_L.get(depth, 0.25)
    if depth > 1 and segs[-1] == '2':
        L += 0.04
    r, g, b = colorsys.hls_to_rgb(hue % 1.0, max(0.0, min(0.80, L)), CLUSTER_TREE_SAT)
    return '#%02x%02x%02x' % (round(r * 255), round(g * 255), round(b * 255))


def hp_rank(hp):
    return HP_ORDER.index(hp) if hp in HP_ORDER else len(HP_ORDER)


def tn_of(v):
    return '1' if str(v) in ('1', 'true', 'True') else '0'


# ----------------------------------------------------------------- colormaps
def rdbu_hex(b, quant=None):
    """methylation beta -> RdBu_r hex (0 blue, .5 white, 1 red); NaN -> grey."""
    if b != b:
        return NA_GREY
    if quant:
        b = round(b * quant) / quant
    b = max(0.0, min(1.0, b))
    if b < 0.5:
        t = b / 0.5; c0, c1 = (33, 102, 172), (247, 247, 247)
    else:
        t = (b - 0.5) / 0.5; c0, c1 = (247, 247, 247), (178, 24, 43)
    return '#%02x%02x%02x' % tuple(int(c0[i] + (c1[i] - c0[i]) * t) for i in range(3))


_MAGMA = [(0.0, (13, 8, 38)), (0.25, (80, 18, 82)), (0.5, (165, 44, 96)),
          (0.75, (229, 92, 75)), (1.0, (254, 194, 135))]


def magma_hex(t, quant=None):
    """distance ratio [0,1] -> magma hex (0 dark/near, 1 bright/far); NaN -> dark grey."""
    if t != t:
        return '#374151'
    if quant:
        t = round(t * quant) / quant
    t = max(0.0, min(1.0, t))
    for k in range(len(_MAGMA) - 1):
        t0, c0 = _MAGMA[k]; t1, c1 = _MAGMA[k + 1]
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0
            return '#%02x%02x%02x' % tuple(int(c0[i] + (c1[i] - c0[i]) * f) for i in range(3))
    return '#fec287'


def mpl_cmaps():
    """matplotlib colormaps for PNG renderers (meth RdBu_r + NaN grey ; dist magma)."""
    import matplotlib.pyplot as plt
    meth = plt.cm.RdBu_r.copy(); meth.set_bad(NA_GREY)
    return meth, 'magma'


def _set_cjk():
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False


def legend_groups(present_labels, kind='meth', n_cluster=0):
    """categorical groups [(category, [(name,hex)])] for the sidebars actually shown."""
    groups = []
    if 'cluster' in present_labels:
        groups.append(('cluster 無監督群', [(f'群{k+1}', CLUSTER_COL[k % len(CLUSTER_COL)]) for k in range(max(2, n_cluster))]))
    if 'HP' in present_labels:
        groups.append(('HP germline 單倍型', [(HP_NAME[g], HP_COL[g]) for g in ['1', '1-1', '2', '2-1', '3', '0']]))
    if 'ALT' in present_labels:
        groups.append(('Allele 帶變異?', [('ALT 帶變異', ALT_COL['ALT']), ('REF 參考', ALT_COL['REF']), ('UNKNOWN', ALT_COL['UNKNOWN'])]))
    if 'T/N' in present_labels:
        groups.append(('T/N 樣本', [('tumor 腫瘤', TN_COL['1']), ('normal 正常', TN_COL['0'])]))
    if 'Strand' in present_labels:
        groups.append(('Strand 股向', [('+ 正股', STRAND_COL['+']), ('− 負股', STRAND_COL['-'])]))
    return groups


def grouped_legend(ax, present_labels, kind='meth', n_cluster=0, metric='NHD'):
    """compact category-grouped legend (one tight row per sidebar category) + a real
    colormap GRADIENT bar (meth β 0/0.5/1 distinct + NaN ; dist 近→遠)."""
    from matplotlib.patches import Rectangle
    ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    groups = legend_groups(present_labels, kind, n_cluster)

    def width_of(name):   # CJK chars are ~2x wider than Latin -> weight them (generous spacing)
        w = sum(0.020 if ord(c) > 0x2E7F else 0.011 for c in name)
        return 0.021 + w + 0.018

    rows = list(groups)

    def _count_lines():   # account for wrapping so line-height is right (no bottom overflow)
        n = 0
        for cat, members in rows:
            x = 0.215; n += 1
            for name, _ in members:
                if x + width_of(name) > 0.90:
                    n += 1; x = 0.215
                x += width_of(name)
        return n + 1   # + colormap gradient line
    nline = _count_lines()
    lh = 1.0 / nline
    yline = 0
    for cat, members in rows:
        y = 1 - (yline + 0.5) * lh
        ax.text(0.004, y, cat + "：", fontsize=6.3, fontweight='bold', va='center', ha='left')
        x = 0.215
        for name, hexv in members:
            if x + width_of(name) > 0.90:
                yline += 1; y = 1 - (yline + 0.5) * lh; x = 0.215
            ax.add_patch(Rectangle((x, y - 0.26 * lh), 0.015, 0.52 * lh, transform=ax.transAxes,
                                   facecolor=hexv, edgecolor='#0008', lw=0.4, clip_on=False))
            ax.text(x + 0.019, y, name, fontsize=6.1, va='center', ha='left')
            x += width_of(name)
        yline += 1
    # gradient colorbar line
    y = 1 - (yline + 0.5) * lh
    if kind == 'meth':
        ax.text(0.004, y, "甲基 β：", fontsize=6.3, fontweight='bold', va='center', ha='left')
        gx0, gw = 0.215, 0.26; ncell = 26
        for i in range(ncell):
            ax.add_patch(Rectangle((gx0 + gw * i / ncell, y - 0.26 * lh), gw / ncell + 0.001, 0.52 * lh,
                                   transform=ax.transAxes, facecolor=rdbu_hex(i / (ncell - 1)), edgecolor='none', clip_on=False))
        for frac, lab in [(0.0, "0 未甲基"), (0.5, "0.5"), (1.0, "1 甲基化")]:
            ax.text(gx0 + gw * frac, y - 0.42 * lh, lab, fontsize=5.6, va='top',
                    ha=('left' if frac == 0 else 'right' if frac == 1 else 'center'), color='#334155')
        xg = gx0 + gw + 0.03
        ax.add_patch(Rectangle((xg, y - 0.26 * lh), 0.015, 0.52 * lh, transform=ax.transAxes,
                               facecolor=NA_GREY, edgecolor='#0008', lw=0.4, clip_on=False))
        ax.text(xg + 0.019, y, "NaN 未覆蓋", fontsize=6.1, va='center', ha='left')
        xg2 = xg + 0.20
        ax.plot([xg2, xg2], [y - 0.28 * lh, y + 0.28 * lh], color=SNV_COL, lw=1.0, transform=ax.transAxes, clip_on=False)
        ax.text(xg2 + 0.012, y, "SNV(在CpG之間)", fontsize=6.1, va='center', ha='left')
    else:
        ax.text(0.004, y, f"距離 {metric}：", fontsize=6.3, fontweight='bold', va='center', ha='left')
        gx0, gw = 0.215, 0.30; ncell = 26
        for i in range(ncell):
            ax.add_patch(Rectangle((gx0 + gw * i / ncell, y - 0.26 * lh), gw / ncell + 0.001, 0.52 * lh,
                                   transform=ax.transAxes, facecolor=magma_hex(i / (ncell - 1)), edgecolor='none', clip_on=False))
        ax.text(gx0, y - 0.42 * lh, "0 近/相似", fontsize=5.6, va='top', ha='left', color='#334155')
        ax.text(gx0 + gw, y - 0.42 * lh, "遠/相異", fontsize=5.6, va='top', ha='right', color='#334155')
        ax.text(gx0 + gw + 0.04, y, "對角暗塊=該群成立 · 右側樹=UPGMA(依群著色)", fontsize=6.0, va='center', ha='left', color='#334155')


def _link_colors(Z, leaf_cluster):
    """color each dendrogram link by cluster if all its leaves share one cluster, else grey."""
    n = len(Z) + 1
    leaves_under = {i: {i} for i in range(n)}
    color = {}
    cl_ids = sorted(set(leaf_cluster.values()))
    cmap = {c: CLUSTER_COL[k % len(CLUSTER_COL)] for k, c in enumerate(cl_ids)}
    for i, (a, b, _, _) in enumerate(Z):
        node = n + i
        leaves_under[node] = leaves_under[int(a)] | leaves_under[int(b)]
        cls = {leaf_cluster.get(le) for le in leaves_under[node]}
        color[node] = cmap[next(iter(cls))] if len(cls) == 1 else '#64748b'
    return color


def mpl_panel(kind, data, sidebars, outpath, info=None, cpgs=None, pos=None, refalt=None,
              dendro_Z=None, leaf_cluster=None, n_cluster=0, tn_split=None, dpi=120):
    """ONE panel PNG + bottom grouped legend, for HTML <img> display. GENERIC (no case title;
    locus/V live in HTML). kind='meth' read×CpG β(RdBu_r); kind='dist' SQUARE read×read (magma)
    + optional right UPGMA dendrogram coloured by cluster. info={n_reads,n_cpg,window}."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _set_cjk()
    meth_cmap, dist_cmap = mpl_cmaps()
    info = info or {}
    nsb = len(sidebars)
    cl_seq = next((h for l, h in sidebars if l == 'cluster'), sidebars[-1][1])
    nrow = data.shape[0]
    HFIG = 4.8
    HLEG = 0.6
    nr = info.get('n_reads', '?'); win = info.get('window', '?'); metric = info.get('metric', 'NHD')

    if kind == 'meth':
        ncol = data.shape[1]
        wr = [0.045] * nsb + [1.0]
        fig = plt.figure(figsize=(max(4.6, 0.045 * ncol + 0.25 * nsb + 1.4), HFIG))
        gs = fig.add_gridspec(2, nsb + 1, width_ratios=wr, height_ratios=[1, HLEG], wspace=0.025, hspace=0.05)
        for c, (lab, hexes) in enumerate(sidebars):
            _sb(fig.add_subplot(gs[0, c]), hexes, lab)
        a = fig.add_subplot(gs[0, nsb])
        a.imshow(data, aspect="auto", cmap=meth_cmap, vmin=0, vmax=1, interpolation="nearest")
        sx = snv_fractional_x(cpgs, pos)
        if sx is not None:
            a.axvline(sx, color=SNV_COL, lw=1.2, zorder=4)
            a.annotate("", xy=(sx, -0.5), xytext=(sx, -nrow * 0.05 - 1.2),
                       arrowprops=dict(arrowstyle="-|>", color=SNV_COL, lw=1.3), annotation_clip=False, zorder=5)
            lbl = f"SNV {refalt[0]}>{refalt[1]}" if refalt and refalt[0] and refalt[1] else "SNV"
            a.text(sx, -nrow * 0.05 - 1.5, lbl, ha="center", va="bottom", fontsize=7.5,
                   color="#a16207", fontweight="bold", clip_on=False)
        _cluster_sep(a, cl_seq, SEP_COL, 0.9)
        if tn_split:   # thick tumor|normal boundary (tumor×normal overlay)
            a.axhline(tn_split - 0.5, color="black", lw=2.2)
            a.text(-0.8, tn_split * 0.5, "tumor", rotation=90, ha="right", va="center", fontsize=6.5, color="#f97316")
            a.text(-0.8, tn_split + (data.shape[0] - tn_split) * 0.5, "normal", rotation=90, ha="right", va="center", fontsize=6.5, color="#22c55e")
        a.set_xticks([]); a.set_yticks([])
        fig.suptitle(f"{info.get('label', '甲基 read×CpG')} · {ncol} CpG · {nr} reads · 取樣窗 ±{win}bp", fontsize=7.5, color="#475569", y=0.985)
    else:
        have_tree = dendro_Z is not None and leaf_cluster is not None
        wr = [0.045] * nsb + ([1.0, 0.4] if have_tree else [1.0])
        sum_wr = sum(wr)
        main_h = HFIG * 1.0 / (1.0 + HLEG)          # main-row height (in)
        figw = main_h * sum_wr / 1.0                # makes the heatmap column (wr=1.0) a SQUARE box
        fig = plt.figure(figsize=(figw, HFIG))
        gs = fig.add_gridspec(2, len(wr), width_ratios=wr, height_ratios=[1, HLEG], wspace=0.02, hspace=0.05)
        for c, (lab, hexes) in enumerate(sidebars):
            _sb(fig.add_subplot(gs[0, c]), hexes, lab)
        a = fig.add_subplot(gs[0, nsb])
        vmax = max(0.5, float(np.nanmax(data)) if np.isfinite(np.nanmax(data)) else 0.5)
        a.imshow(data, aspect="auto", cmap=dist_cmap, vmin=0, vmax=vmax, interpolation="nearest")  # fills square box
        _cluster_sep(a, cl_seq, '#ffffff', 0.7, alpha=0.6)
        if tn_split:   # tumor|normal boundary on BOTH axes (square read×read)
            a.axhline(tn_split - 0.5, color="black", lw=1.8)
            a.axvline(tn_split - 0.5, color="black", lw=1.8)
        a.set_xticks([]); a.set_yticks([])
        if have_tree:
            from scipy.cluster.hierarchy import dendrogram
            dax = fig.add_subplot(gs[0, nsb + 1])
            lc = _link_colors(dendro_Z, leaf_cluster)
            dendrogram(dendro_Z, orientation='right', ax=dax, no_labels=True,
                       link_color_func=lambda nid: lc.get(nid, '#64748b'))
            dax.invert_yaxis()
            dax.margins(y=0)
            dax.set_xticks([]); dax.set_yticks([])
            for s in dax.spines.values():
                s.set_visible(False)
            dax.set_title("UPGMA", fontsize=6.5)
        fig.suptitle(f"{info.get('dist_label', '距離 read×read ' + metric)} · {nr} reads · 取樣窗 ±{win}bp · 暗近亮遠", fontsize=7.5, color="#475569", y=0.985)
    axl = fig.add_subplot(gs[1, :])
    grouped_legend(axl, [s[0] for s in sidebars], kind=kind, n_cluster=n_cluster, metric=metric)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight"); plt.close(fig)
    return outpath


def _cluster_sep(a, cl_seq, color, lw, alpha=1.0):
    prev = None
    for j in range(len(cl_seq)):
        if prev is not None and cl_seq[j] != prev:
            a.axhline(j - 0.5, color=color, lw=lw, alpha=alpha)
        prev = cl_seq[j]


def _sb(ax, hexes, lab=None, horiz=False):
    from matplotlib.colors import to_rgb
    rgb = np.array([to_rgb(h) for h in hexes])
    ax.imshow(rgb[np.newaxis] if horiz else rgb[:, None], aspect="auto", interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    if lab and not horiz:   # vertical title above the thin column -> never overlaps neighbours
        ax.text(0.5, 1.004, lab, transform=ax.transAxes, rotation=90, ha='center', va='bottom', fontsize=6)
    for s in ax.spines.values():
        s.set_visible(False)


def mpl_dual_panel(meth_mat, dist_mat, sidebars, cpgs, pos, title, outpath, n_cluster=0, dpi=110):
    """PNG: LEFT read×CpG methylation + RIGHT read×read distance + bottom GROUPED legend.
    sidebars = ordered [(label,[hex per read])] from sidebar_specs(). For HTML <img> display."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _set_cjk()
    meth_cmap, dist_cmap = mpl_cmaps()
    nsb = len(sidebars); nrow, ncol = meth_mat.shape
    wr = [0.05] * nsb + [1.0, 0.16] + [0.05] * nsb + [1.0]
    fig = plt.figure(figsize=(9.6, 5.0))
    gs = fig.add_gridspec(2, len(wr), width_ratios=wr, height_ratios=[1, 0.42], wspace=0.04, hspace=0.02)
    cl_seq = sidebars[0][1]
    # meth
    c = 0
    for lab, hexes in sidebars:
        _sb(fig.add_subplot(gs[0, c]), hexes, lab); c += 1
    axm = fig.add_subplot(gs[0, c]); c += 1
    axm.imshow(meth_mat, aspect="auto", cmap=meth_cmap, vmin=0, vmax=1, interpolation="nearest")
    sx = snv_fractional_x(cpgs, pos)
    if sx is not None:
        axm.axvline(sx, color=SNV_COL, lw=2.0, zorder=4)
        axm.plot([sx], [-0.6], marker="v", color=SNV_COL, markersize=8, clip_on=False, zorder=5)
    prev = None
    for j in range(nrow):
        if prev is not None and cl_seq[j] != prev:
            axm.axhline(j - 0.5, color="black", lw=0.8)
        prev = cl_seq[j]
    axm.set_xticks([]); axm.set_yticks([]); axm.set_xlabel(f"{ncol} CpG · SNV在CpG之間", fontsize=7)
    axm.set_title("甲基 read×CpG", fontsize=8.5)
    fig.add_subplot(gs[0, c]).axis("off"); c += 1   # gap
    # dist
    for lab, hexes in sidebars:
        _sb(fig.add_subplot(gs[0, c]), hexes, lab); c += 1
    axd = fig.add_subplot(gs[0, c])
    vmax = max(0.5, float(np.nanmax(dist_mat)) if np.isfinite(np.nanmax(dist_mat)) else 0.5)
    axd.imshow(dist_mat, aspect="auto", cmap=dist_cmap, vmin=0, vmax=vmax, interpolation="nearest")
    prev = None
    for j in range(nrow):
        if prev is not None and cl_seq[j] != prev:
            axd.axhline(j - 0.5, color="white", lw=0.6, alpha=0.6)
        prev = cl_seq[j]
    axd.set_xticks([]); axd.set_yticks([]); axd.set_xlabel("read×read 距離 · 暗=近/亮=遠", fontsize=7)
    axd.set_title("距離 read×read (對角塊=群成立)", fontsize=8.5)
    # grouped legend across bottom
    axl = fig.add_subplot(gs[1, :])
    grouped_legend(axl, [s[0] for s in sidebars], n_cluster)
    fig.suptitle(title, fontsize=9.5, y=1.0)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight"); plt.close(fig)
    return outpath


# ----------------------------------------------------------------- standard sidebars
def sidebar_specs(reads, ids, cluster_of=None, include_tn=True, include_strand=True):
    """Return ordered [(label, [hex per id])]. Order (left→right, heatmap on the right):
    HP | ALT | T/N? | Strand? | cluster?  — cluster INNERMOST (closest to heatmap, 2026-06-18).
    reads = {id: {hp,alt_support,is_tumor,strand}}."""
    specs = [('HP', [HP_COL.get(reads[i].get('hp', ''), NA_GREY) for i in ids]),
             ('ALT', [ALT_COL.get(reads[i].get('alt_support', ''), '#e5e7eb') for i in ids])]
    if include_tn:
        specs.append(('T/N', [TN_COL.get(tn_of(reads[i].get('is_tumor', '1')), '#e5e7eb') for i in ids]))
    if include_strand:
        specs.append(('Strand', [STRAND_COL.get(reads[i].get('strand', ''), '#e5e7eb') for i in ids]))
    if cluster_of is not None:
        cl_ids = sorted({cluster_of[i] for i in ids if i in cluster_of})
        cmap = {c: CLUSTER_COL[k % len(CLUSTER_COL)] for k, c in enumerate(cl_ids)}
        specs.append(('cluster', [cmap.get(cluster_of.get(i), NA_GREY) for i in ids]))
    return specs


def snv_fractional_x(cpgs, pos):
    """fractional column index of the SNV BETWEEN flanking CpGs (variant rarely a CpG)."""
    if not cpgs:
        return None
    i = bisect.bisect_left(cpgs, int(pos))
    if i == 0:
        return -0.5
    if i >= len(cpgs):
        return len(cpgs) - 0.5
    l, r = cpgs[i - 1], cpgs[i]
    return (i - 1) + ((int(pos) - l) / (r - l) if r > l else 0.5)


def legend_entries(include_cluster=False, n_cluster=0, include_tn=True):
    """(label, hex) pairs for a standard legend."""
    out = []
    if include_cluster:
        for k in range(min(n_cluster, len(CLUSTER_COL))):
            out.append((f'cluster{k+1}', CLUSTER_COL[k]))
    for g in ['1', '1-1', '2', '2-1', '3', '0']:
        out.append((HP_NAME[g], HP_COL[g]))
    out += [('ALT(變異)', ALT_COL['ALT']), ('REF', ALT_COL['REF'])]
    if include_tn:
        out += [('tumor', TN_COL['1']), ('normal', TN_COL['0'])]
    out += [('strand +', STRAND_COL['+']), ('strand −', STRAND_COL['-'])]
    return out


# ----------------------------------------------------------------- SVG dual-panel (inline, portable)
def svg_dual_panel(meth_mat, dist_mat, sidebars, cpgs, pos, cell_h=None):
    """Inline SVG: LEFT read×CpG methylation (RdBu) + RIGHT read×read distance (magma),
    sharing the read order; standard sidebars on the left of each panel + cluster top bar on dist.
    sidebars = ordered [(label,[hex per read])] from sidebar_specs(). Renders in ANY html viewer."""
    nrow, ncol = meth_mat.shape
    CH = cell_h or max(2, min(5, 300 // max(nrow, 1)))
    CW = 4; DW = CH; SBW = 9; GAP = 18; LBL = 24
    nsb = len(sidebars)
    sbtot = nsb * SBW + 2
    ax0 = sbtot                                   # meth grid x
    Aw = ax0 + ncol * CW
    bsb = Aw + GAP
    bx0 = bsb + sbtot
    W = bx0 + nrow * DW + 4
    TOP = SBW + 2
    gtop = TOP + 2
    H = gtop + nrow * CH + LBL
    vmax = max(0.5, float(np.nanmax(dist_mat)) if np.isfinite(np.nanmax(dist_mat)) else 0.5)
    P = [f'<svg viewBox="0 0 {W} {H}" width="100%" preserveAspectRatio="xMidYMid meet" '
         f'shape-rendering="crispEdges" font-family="system-ui" style="background:#11161f">']
    # left sidebars for BOTH panels (shared read order)
    for s, (lab, cols) in enumerate(sidebars):
        xa = s * SBW; xb = bsb + s * SBW
        for j in range(nrow):
            y = gtop + j * CH
            P.append(f'<rect x="{xa}" y="{y}" width="{SBW-1}" height="{CH}" fill="{cols[j]}"/>')
            P.append(f'<rect x="{xb}" y="{y}" width="{SBW-1}" height="{CH}" fill="{cols[j]}"/>')
    # dist top bar = cluster (or HP if no cluster) to show both axes share order
    topbar = sidebars[0][1]
    for j in range(nrow):
        P.append(f'<rect x="{bx0 + j*DW}" y="{gtop-SBW}" width="{DW}" height="{SBW-1}" fill="{topbar[j]}"/>')
    # meth grid (run-length, quantized)
    for j in range(nrow):
        y = gtop + j * CH; row = meth_mat[j]; c = 0
        while c < ncol:
            col = rdbu_hex(row[c], 12); run = 1
            while c + run < ncol and rdbu_hex(row[c + run], 12) == col:
                run += 1
            P.append(f'<rect x="{ax0 + c*CW}" y="{y}" width="{run*CW}" height="{CH}" fill="{col}"/>')
            c += run
    # dist grid
    for j in range(nrow):
        y = gtop + j * CH; row = dist_mat[j]; c = 0
        while c < nrow:
            col = magma_hex(row[c] / vmax if row[c] == row[c] else float('nan'), 16); run = 1
            while c + run < nrow and magma_hex(row[c+run] / vmax if row[c+run] == row[c+run] else float('nan'), 16) == col:
                run += 1
            P.append(f'<rect x="{bx0 + c*DW}" y="{y}" width="{run*DW}" height="{CH}" fill="{col}"/>')
            c += run
    # separators where first sidebar (cluster/HP) changes
    grp = sidebars[0][1]; prev = None
    for j in range(nrow):
        if prev is not None and grp[j] != prev:
            yy = gtop + j * CH
            P.append(f'<line x1="{ax0}" y1="{yy}" x2="{Aw}" y2="{yy}" stroke="#000" stroke-width="0.7"/>')
            P.append(f'<line x1="{bx0}" y1="{yy}" x2="{bx0+nrow*DW}" y2="{yy}" stroke="#fff" stroke-width="0.5" opacity="0.6"/>')
        prev = grp[j]
    # SNV on meth panel
    sx = snv_fractional_x(cpgs, pos)
    if sx is not None:
        x = ax0 + sx * CW
        P.append(f'<line x1="{x:.1f}" y1="{gtop-2}" x2="{x:.1f}" y2="{gtop+nrow*CH}" stroke="{SNV_COL}" stroke-width="1.4"/>')
        P.append(f'<text x="{x:.1f}" y="{gtop-4}" fill="{SNV_COL}" font-size="7" text-anchor="middle">SNV</text>')
    barlabels = '｜'.join(s[0] for s in sidebars)
    P.append(f'<text x="0" y="{H-12}" fill="#94a3b8" font-size="6.5">{barlabels}</text>')
    P.append(f'<text x="{ax0}" y="{H-4}" fill="#cbd5e1" font-size="7.5">甲基 read×CpG ({ncol}CpG·{nrow}reads)</text>')
    P.append(f'<text x="{bx0}" y="{H-4}" fill="#cbd5e1" font-size="7.5">距離 read×read (對角塊=群·暗近亮遠)</text>')
    P.append('</svg>')
    return ''.join(P)
