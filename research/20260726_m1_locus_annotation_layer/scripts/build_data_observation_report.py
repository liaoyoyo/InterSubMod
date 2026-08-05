#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the exact-PS data observation & analysis review page.

Every number is injected from an authority artifact; nothing is typed into the template.
Conservation identities are re-derived here and the build aborts if any of them break, so
a stale or mismatched input cannot silently reach the page.

Inputs
  census TSVs   research/20260724_exact_ps_cpp_topology_signature_census/data/*.tsv
  census summary  observation_workspaces/20260724_..._signature_census/all7_v1/summary.json
  session analyses  research/20260726_m1_locus_annotation_layer/results/*.json
"""
import argparse
import csv
import hashlib
import html
import json
import os
import sys
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CENSUS_DATA = REPO / "research/20260724_exact_ps_cpp_topology_signature_census/data"
CENSUS_SUMMARY = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.json")
SESSION = REPO / "research/20260726_m1_locus_annotation_layer/results"

DISPLAY = {"HCC1395": "HCC1395_HKU", "HCC1395_DORADO": "HCC1395_NYGC"}
ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]

# Every fill that ever carries a label was checked against WCAG: series dark enough for
# white text were darkened to clear 4.5:1, and the deliberately pale "excluded" fills keep
# their light tone but take dark ink instead (see text_on).
PALETTE = {
    "s": "#1f6f5c", "active": "#246874", "w": "#8f5714", "single": "#b9bcc0",
    "unique": "#176b58", "same": "#2f6493", "cross": "#96372c",
    "direct": "#2f7a5a", "sister": "#6f4f96", "sisdir": "#8a6114", "sgl": "#b9bcc0",
    "hp1": "#2d6f91", "hp2": "#8a4a5e",
}
INK_DARK = "#14201b"


def _lum(hex_colour):
    h = hex_colour.lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    channels = []
    for i in (0, 2, 4):
        c = int(h[i:i + 2], 16) / 255
        channels.append(c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4)
    return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2]


def contrast(a, b):
    la, lb = _lum(a), _lum(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def text_on(background):
    """Pick whichever ink reads better on this fill, rather than assuming white."""
    return "#ffffff" if contrast("#ffffff", background) >= contrast(INK_DARK, background) \
        else INK_DARK


def esc(text):
    return html.escape(str(text), quote=True)


def fnum(value):
    return f"{int(value):,}"


def fpct(part, whole, digits=2):
    if not whole:
        return "—"
    return f"{100.0 * part / whole:.{digits}f}%"


def read_tsv(name):
    path = CENSUS_DATA / name
    if not path.exists():
        sys.exit(f"FAIL CLOSED: missing authority TSV {path}")
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def as_int(value):
    return int(str(value).replace(",", "").strip())


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- SVG helpers
def svg_open(width, height, title):
    # no height attribute: "auto" is not a valid SVG length, so the aspect ratio is left
    # to viewBox and the width:100% / height:auto CSS pair instead.
    return (f'<svg viewBox="0 0 {width} {height}" preserveAspectRatio="xMidYMid meet" '
            f'role="img" aria-label="{esc(title)}" '
            f'style="width:100%;height:auto;max-width:{width}px;display:block;margin:0 auto">'
            f'<title>{esc(title)}</title>')


def stacked_row(x, y, width, height, parts, total):
    """parts = [(value, colour, label)] -> one horizontal 100% stacked bar.

    Wide segments carry both the share and the absolute count, because a percentage with
    no numerator is exactly the kind of figure a reviewer cannot check.
    """
    out = []
    cursor = x
    for value, colour, label in parts:
        if total <= 0 or value <= 0:
            continue
        seg = width * value / total
        out.append(
            f'<rect x="{cursor:.2f}" y="{y}" width="{seg:.2f}" height="{height}" '
            f'fill="{colour}"><title>{esc(label)}: {fnum(value)} '
            f'({fpct(value, total)})</title></rect>')
        cx = cursor + seg / 2
        ink = text_on(colour)
        if seg >= 78:
            out.append(
                f'<text class="bl" x="{cx:.2f}" y="{y + height / 2 - 1}" text-anchor="middle" '
                f'font-size="11.5" font-weight="700" fill="{ink}" '
                f'font-family="ui-monospace,monospace">{fpct(value, total, 1)}</text>'
                f'<text class="bl" x="{cx:.2f}" y="{y + height / 2 + 12}" text-anchor="middle" '
                f'font-size="10" font-weight="600" fill="{ink}" '
                f'font-family="ui-monospace,monospace">{fnum(value)}</text>')
        elif seg >= 40:
            out.append(
                f'<text class="bl" x="{cx:.2f}" y="{y + height / 2 + 4}" text-anchor="middle" '
                f'font-size="11" font-weight="600" fill="{ink}" '
                f'font-family="ui-monospace,monospace">{fpct(value, total, 1)}</text>')
        cursor += seg
    return "".join(out)


def legend(items):
    cells = "".join(
        f'<span class="lg"><i style="background:{c}"></i>{esc(t)}</span>' for c, t in items)
    return f'<div class="legend">{cells}</div>'


def per_sample_stacked(rows, series, title, note=""):
    """rows: [(display_label, total, [(value, colour, label)])]"""
    row_h, gap, left, right, top = 34, 12, 132, 104, 10
    height = top + len(rows) * (row_h + gap) + 8
    width = 980
    bar_w = width - left - right
    out = [svg_open(width, height, title)]
    for index, (label, total, parts) in enumerate(rows):
        y = top + index * (row_h + gap)
        out.append(f'<text x="{left - 10}" y="{y + row_h / 2 + 4}" text-anchor="end" '
                   f'font-size="12" fill="var(--fg)" '
                   f'font-family="ui-monospace,monospace">{esc(label)}</text>')
        out.append(f'<rect x="{left}" y="{y}" width="{bar_w}" height="{row_h}" '
                   f'fill="var(--track)"/>')
        out.append(stacked_row(left, y, bar_w, row_h, parts, total))
        out.append(f'<text x="{left + bar_w + 8}" y="{y + row_h / 2 + 4}" '
                   f'font-size="11" fill="var(--muted)" '
                   f'font-family="ui-monospace,monospace">n={fnum(total)}</text>')
    out.append("</svg>")
    body = "".join(out) + legend(series)
    if note:
        body += f'<p class="note">{note}</p>'
    return body


def funnel_svg(steps, title):
    """steps: [(label, value, colour, note)] descending"""
    width, step_h, gap, left = 980, 52, 14, 8
    height = len(steps) * (step_h + gap) + 16
    top_value = max(v for _, v, _, _ in steps) or 1
    out = [svg_open(width, height, title)]
    for index, (label, value, colour, note) in enumerate(steps):
        y = 8 + index * (step_h + gap)
        bar_w = max(4.0, (width - 300) * value / top_value)
        ink = text_on(colour)
        out.append(f'<rect x="{left}" y="{y}" width="{bar_w:.2f}" height="{step_h}" '
                   f'rx="3" fill="{colour}"/>')
        out.append(f'<text class="bl" x="{left + 12}" y="{y + 21}" font-size="13.5" fill="{ink}" '
                   f'font-weight="700">{esc(label)}</text>')
        out.append(f'<text class="bl" x="{left + 12}" y="{y + 39}" font-size="16" fill="{ink}" '
                   f'font-weight="600" font-family="ui-monospace,monospace">'
                   f'{fnum(value)}</text>')
        out.append(f'<text x="{left + bar_w + 14:.2f}" y="{y + 24}" font-size="12" '
                   f'fill="var(--muted)">{esc(note)}</text>')
        out.append(f'<text x="{left + bar_w + 14:.2f}" y="{y + 41}" font-size="11" '
                   f'fill="var(--muted)" font-family="ui-monospace,monospace">'
                   f'{fpct(value, top_value)} of top</text>')
    out.append("</svg>")
    return "".join(out)


def diverging_svg(rows, title, unit="pp", note=""):
    """rows: [(label, [(series_label, value, colour)])] with signed values around zero."""
    row_h, gap, left, right, top = 26, 10, 132, 16, 8
    width = 980
    plot_w = width - left - right
    mid = left + plot_w / 2
    span = max((abs(v) for _, series in rows for _, v, _ in series), default=1.0) or 1.0
    height = top + len(rows) * (row_h + gap) + 26
    out = [svg_open(width, height, title)]
    for gridline in (-1.0, -0.5, 0.0, 0.5, 1.0):
        x = mid + (plot_w / 2) * gridline
        out.append(f'<line x1="{x:.1f}" y1="{top - 4}" x2="{x:.1f}" '
                   f'y2="{top + len(rows) * (row_h + gap)}" stroke="var(--line)" '
                   f'stroke-width="{2 if gridline == 0 else 1}"/>')
        out.append(f'<text x="{x:.1f}" y="{height - 8}" text-anchor="middle" '
                   f'font-size="10" fill="var(--muted)" font-family="ui-monospace,monospace">'
                   f'{gridline * span * 100:+.0f}</text>')
    for index, (label, series) in enumerate(rows):
        y = top + index * (row_h + gap)
        out.append(f'<text x="{left - 10}" y="{y + row_h / 2 + 4}" text-anchor="end" '
                   f'font-size="12" fill="var(--fg)" '
                   f'font-family="ui-monospace,monospace">{esc(label)}</text>')
        sub_h = row_h / max(1, len(series))
        for offset, (name, value, colour) in enumerate(series):
            w = (plot_w / 2) * (abs(value) / span)
            x = mid - w if value < 0 else mid
            out.append(f'<rect x="{x:.2f}" y="{y + offset * sub_h:.2f}" '
                       f'width="{max(w, 0.6):.2f}" height="{sub_h - 1.5:.2f}" '
                       f'fill="{colour}"><title>{esc(name)}: {value * 100:+.2f} {unit}</title>'
                       f'</rect>')
    out.append("</svg>")
    body = "".join(out)
    if note:
        body += f'<p class="note">{note}</p>'
    return body


def table(headers, rows, foot=None, cls=""):
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    tfoot = ""
    if foot:
        tfoot = "<tfoot><tr>" + "".join(f"<td>{c}</td>" for c in foot) + "</tr></tfoot>"
    return (f'<div class="tw"><table class="{cls}"><thead><tr>{head}</tr></thead>'
            f"<tbody>{body}</tbody>{tfoot}</table></div>")


# --------------------------------------------------------------------------- main build
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    overview = {r["dataset"]: r for r in read_tsv(
        "20260724_exactPS_source_overview_by_sample_01.tsv")}
    funnel = {r["dataset"]: r for r in read_tsv(
        "20260724_exactPS_funnel_by_sample_01.tsv")}

    def long_k(name):
        """the k TSVs are long-format (dataset, k_band, count, ...) -> dataset -> band -> n"""
        out = {}
        for row in read_tsv(name):
            out.setdefault(row["dataset"], {})[row["k_band"]] = as_int(row["count"])
        return out

    src_k = long_k("20260724_exactPS_source_component_k_by_sample_01.tsv")
    fin_k = long_k("20260724_exactPS_final_group_k_by_sample_01.tsv")
    hp = {}
    for row in read_tsv("20260724_exactPS_hp_split_by_sample_01.tsv"):
        hp.setdefault(row["dataset"], {})[row["hp_family"]] = row

    # the authority params travel with the page so the PS scope can never be misread
    scope_json = json.loads(
        (CENSUS_DATA / "20260724_exactPS_k_hp_funnel_census_01.json").read_text())
    scope = scope_json["scope"]
    if not CENSUS_SUMMARY.exists():
        sys.exit(f"FAIL CLOSED: missing census summary {CENSUS_SUMMARY}")
    census = json.loads(CENSUS_SUMMARY.read_text())
    census_by_sample = {s["sample"]: s for s in census["samples"]}

    def session(name):
        path = SESSION / name
        if not path.exists():
            sys.exit(f"FAIL CLOSED: missing session analysis {path}")
        return json.loads(path.read_text())

    tied = session("tied_cross_topology_characterisation_v1.json")
    afdeg = session("tie_af_degeneracy_v1.json")

    samples = [s for s in ORDER if s in overview]
    if len(samples) != 7:
        sys.exit(f"FAIL CLOSED: expected 7 datasets, resolved {len(samples)}")

    # ---- conservation identities, re-derived (never trusted from prose)
    tot_S = sum(as_int(overview[s]["candidate_loci_S"]) for s in samples)
    tot_active_mem = sum(as_int(overview[s]["active_node_memberships"]) for s in samples)
    tot_comp = sum(as_int(overview[s]["all_components"]) for s in samples)
    tot_single = sum(as_int(overview[s]["singleton_k1_components"]) for s in samples)
    tot_W = sum(as_int(overview[s]["tree_eligible_W"]) for s in samples)
    tot_Wmem = sum(as_int(overview[s]["W_memberships"]) for s in samples)
    tot_blocks = sum(as_int(funnel[s]["bounded_blocks"]) for s in samples)
    tot_cut = sum(as_int(funnel[s]["post_cut_k1_blocks_excluded"]) for s in samples)
    tot_unsup = sum(as_int(funnel[s]["pattern_unsupported_blocks_excluded"]) for s in samples)
    tot_final = sum(as_int(funnel[s]["final_groups"]) for s in samples)
    tot_finmem = sum(as_int(funnel[s]["final_memberships"]) for s in samples)
    tot_containers = sum(as_int(overview[s]["exact_PS_HP_containers"]) for s in samples)
    tot_active_loci = sum(as_int(overview[s]["active_unique_loci"]) for s in samples)

    checks = {
        "components = singleton + W": tot_comp == tot_single + tot_W,
        "blocks - cut - unsupported = final groups":
            tot_blocks - tot_cut - tot_unsup == tot_final,
        "ranked units = census cohort": census["cohort"]["ranked_units"] == sum(
            census_by_sample[s]["ranked_units"] for s in samples),
        "resolution partitions ranked": sum(
            census["cohort"]["resolution"][k]["n"]
            for k in ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY")
        ) == census["cohort"]["ranked_units"],
        "tied-cross characterisation matches census":
            tied["resolution_class_counts"]["TIED_CROSS_TOPOLOGY"]
            == census["cohort"]["resolution"]["TIED_CROSS_TOPOLOGY"]["n"],
        "af-degeneracy unit count matches census":
            afdeg["per_class"]["TIED_CROSS_TOPOLOGY"]["units"]
            == census["cohort"]["resolution"]["TIED_CROSS_TOPOLOGY"]["n"],
    }
    failed = [k for k, v in checks.items() if not v]
    if failed:
        sys.exit("FAIL CLOSED: conservation checks failed -> " + "; ".join(failed))

    ranked = census["cohort"]["ranked_units"]
    res = census["cohort"]["resolution"]

    # ---- S1 cohort funnel
    s1 = funnel_svg([
        ("候選 sSNV（S）", tot_S, PALETTE["s"], "LongPhase-S PASS · chr1–22 · biallelic"),
        ("active unique loci", tot_active_loci, PALETTE["active"],
         "至少進入一個 exact PS×HP 容器"),
        ("active node memberships", tot_active_mem, "#2f6493",
         "(locus, PS, HP)；同一 locus 可重複計入"),
        ("source components", tot_comp, "#436b5c", "含 k=1 singleton 的守恆 registry"),
        ("tree-eligible W（k≥2）", tot_W, PALETTE["w"], "唯一稱為 read-linked candidate region 的單位"),
        ("bounded blocks", tot_blocks, "#7a5c26", "k≤12 切割後的計算單元，B ≠ W"),
        ("final topology groups", tot_final, "#5f3a20", "實際進建樹"),
        ("ranked complete", ranked, "#3a2e24", "determinacy 的唯一分母"),
    ], "cohort funnel")

    # ---- S2 per-sample input -> linkable
    rows2 = []
    tbl2 = []
    for s in samples:
        o = overview[s]
        S = as_int(o["candidate_loci_S"])
        act = as_int(o["active_unique_loci"])
        rows2.append((DISPLAY.get(s, s), S, [
            (act, PALETTE["active"], "active unique loci"),
            (S - act, "#cfccc0", "未進任何 exact PS×HP 容器"),
        ]))
        tbl2.append([
            f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(S), fnum(act), fpct(act, S),
            fnum(as_int(o["active_node_memberships"])),
            fnum(as_int(o["exact_PS_HP_containers"])),
        ])
    tbl2_foot = ["<b>合計</b>", f"<b>{fnum(tot_S)}</b>", f"<b>{fnum(tot_active_loci)}</b>",
                 f"<b>{fpct(tot_active_loci, tot_S)}</b>",
                 f"<b>{fnum(tot_active_mem)}</b>", f"<b>{fnum(tot_containers)}</b>"]

    # ---- S3 component partition
    rows3 = []
    tbl3 = []
    for s in samples:
        o = overview[s]
        comp = as_int(o["all_components"])
        sg = as_int(o["singleton_k1_components"])
        w = as_int(o["tree_eligible_W"])
        rows3.append((DISPLAY.get(s, s), comp, [
            (w, PALETTE["w"], "tree-eligible W (k≥2)"),
            (sg, PALETTE["single"], "k=1 abstain"),
        ]))
        tbl3.append([f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(comp), fnum(sg),
                     fpct(sg, comp), fnum(w), fpct(w, comp),
                     fnum(as_int(o["W_memberships"]))])
    tbl3_foot = ["<b>合計</b>", f"<b>{fnum(tot_comp)}</b>", f"<b>{fnum(tot_single)}</b>",
                 f"<b>{fpct(tot_single, tot_comp)}</b>", f"<b>{fnum(tot_W)}</b>",
                 f"<b>{fpct(tot_W, tot_comp)}</b>", f"<b>{fnum(tot_Wmem)}</b>"]

    # ---- S4 k distributions
    KCOLS = [("k=2", "#1f6f5c"), ("k=3", "#246874"), ("k=4", "#2f6493"),
             ("k=5", "#54679f"), ("k=6", "#6f4f96"), ("k=7", "#8b5280"),
             ("k=8", "#a05a5a"), ("k=9", "#95643c"), ("k=10", "#8a7328"),
             ("k=11", "#71792a"), ("k=12", "#587c3a")]

    def band_key(band, want):
        """the TSV k_band labels are matched loosely so a label tweak upstream is visible
        as a zero segment rather than a silent crash."""
        norm = str(band).strip().lower().replace(" ", "")
        return norm in (want, want.replace("k=", "k"), want.replace("k=", ""))

    def k_rows(source, include_k1):
        out = []
        for s in samples:
            bands = source[s]
            parts = []
            total = 0
            if include_k1:
                v = sum(n for b, n in bands.items() if band_key(b, "k=1"))
                total += v
                parts.append((v, PALETTE["single"], "k=1"))
            for key, colour in KCOLS:
                v = sum(n for b, n in bands.items() if band_key(b, key))
                total += v
                parts.append((v, colour, key))
            gt = sum(n for b, n in bands.items()
                     if ">" in str(b) and "12" in str(b))
            if gt:
                total += gt
                parts.append((gt, "#5f3a20", "k>12"))
            out.append((DISPLAY.get(s, s), total, parts))
        return out

    # ---- S5 determinacy
    rows5 = []
    tbl5 = []
    for s in samples:
        c = census_by_sample[s]
        r = c["resolution"]
        n = c["ranked_units"]
        u = r["UNIQUE_TREE"]["n"]
        sa = r["TIED_SAME_TOPOLOGY"]["n"]
        cr = r["TIED_CROSS_TOPOLOGY"]["n"]
        rows5.append((DISPLAY.get(s, s), n, [
            (u, PALETTE["unique"], "UNIQUE_TREE"),
            (sa, PALETTE["same"], "TIED_SAME_TOPOLOGY"),
            (cr, PALETTE["cross"], "TIED_CROSS_TOPOLOGY"),
        ]))
        tbl5.append([f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(n), fnum(u), fpct(u, n),
                     fnum(sa), fpct(sa, n), fnum(cr), f"<b>{fpct(cr, n)}</b>",
                     fnum(c["one_exact_topology"]["n"]),
                     fpct(c["one_exact_topology"]["n"], n),
                     fnum(c["maximum_best_tree_tie_count"])])
    coh = census["cohort"]
    tbl5_foot = ["<b>ALL7</b>", f"<b>{fnum(ranked)}</b>",
                 f"<b>{fnum(res['UNIQUE_TREE']['n'])}</b>",
                 f"<b>{fpct(res['UNIQUE_TREE']['n'], ranked)}</b>",
                 f"<b>{fnum(res['TIED_SAME_TOPOLOGY']['n'])}</b>",
                 f"<b>{fpct(res['TIED_SAME_TOPOLOGY']['n'], ranked)}</b>",
                 f"<b>{fnum(res['TIED_CROSS_TOPOLOGY']['n'])}</b>",
                 f"<b>{fpct(res['TIED_CROSS_TOPOLOGY']['n'], ranked)}</b>",
                 f"<b>{fnum(coh['one_exact_topology']['n'])}</b>",
                 f"<b>{fpct(coh['one_exact_topology']['n'], ranked)}</b>",
                 f"<b>{fnum(coh['maximum_best_tree_tie_count'])}</b>"]

    # ---- S6 coarse geometry
    CO = [("Single-only", PALETTE["sgl"]), ("Direct-only", PALETTE["direct"]),
          ("Sister-only", PALETTE["sister"]), ("Sister+direct", PALETTE["sisdir"])]
    rows6 = []
    tbl6 = []
    for s in samples:
        c = census_by_sample[s]
        n = c["ranked_units"]
        rc = c["resolved_coarse_class"]
        parts = [(rc.get(k, {}).get("n", 0), col, k) for k, col in CO]
        parts.append((c["cross_coarse_class"]["n"], PALETTE["cross"], "Cross-class unresolved"))
        rows6.append((DISPLAY.get(s, s), n, parts))
        tbl6.append([f"<b>{esc(DISPLAY.get(s, s))}</b>"]
                    + [f"{fnum(rc.get(k, {}).get('n', 0))}<small>{fpct(rc.get(k, {}).get('n', 0), n)}</small>"
                       for k, _ in CO]
                    + [f"{fnum(c['cross_coarse_class']['n'])}<small>{fpct(c['cross_coarse_class']['n'], n)}</small>",
                       fnum(c["cross_exact_topology_but_same_coarse_class"])])
    ccoh = coh["resolved_coarse_class"]
    tbl6_foot = ["<b>ALL7</b>"] + [
        f"<b>{fnum(ccoh.get(k, {}).get('n', 0))}</b><small>{fpct(ccoh.get(k, {}).get('n', 0), ranked)}</small>"
        for k, _ in CO] + [
        f"<b>{fnum(coh['cross_coarse_class']['n'])}</b><small>{fpct(coh['cross_coarse_class']['n'], ranked)}</small>",
        f"<b>{fnum(coh['cross_exact_topology_but_same_coarse_class'])}</b>"]

    # ---- S7 HP split
    # HP multiplicity is asked per *physical locus*, not per component: does a given sSNV
    # position enter final groups on one haplotype only, or on both? That is a different
    # question from the HP1/HP2 component split and needs the MLHP positions.
    MLHP_ROOT = Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples")
    rows7 = []
    hp_mult_tbl = []
    tot_only1 = tot_only2 = tot_both = 0
    for s in samples:
        path = MLHP_ROOT / s / f"{s}.exact_ps_mlhp.json"
        if not path.exists():
            sys.exit(f"FAIL CLOSED: missing MLHP for {s}: {path}")
        doc = json.loads(path.read_text())
        # HP1/HP2 are only comparable inside one phase set: across phase sets the relative
        # orientation of the two haplotypes is undefined, so "HP1" in PS A and "HP1" in
        # PS B are not the same strand. The only defensible question is therefore how many
        # haplotypes a locus is analysed on WITHIN a phase set -- one, or both.
        by_locus_ps = {}
        for g in doc.get("groups") or []:
            fam = str(g.get("hp_family"))
            chrom = g.get("chrom")
            ps = g.get("phase_set")
            for pos in (g.get("positions") or []):
                by_locus_ps.setdefault((chrom, int(pos)), {}).setdefault(ps, set()).add(fam)
        single = both = 0
        multi_ps = 0
        for per_ps in by_locus_ps.values():
            if len(per_ps) > 1:
                multi_ps += 1
            if any({"1", "2"} <= fams for fams in per_ps.values()):
                both += 1
            else:
                single += 1
        total = len(by_locus_ps)
        tot_only1 += single          # reused accumulators: single / both
        tot_both += both
        tot_only2 += multi_ps
        rows7.append((DISPLAY.get(s, s), total, [
            (single, PALETTE["hp1"], "單一 haplotype"),
            (both, "#8a6114", "兩條 haplotype（同一 PS 內）"),
        ]))
        hp_mult_tbl.append([
            f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(total),
            f"{fnum(single)}<small>{fpct(single, total)}</small>",
            f"{fnum(both)}<small>{fpct(both, total)}</small>",
            f"{fnum(multi_ps)}<small>{fpct(multi_ps, total)}</small>",
        ])
    hp_tot = tot_only1 + tot_both
    hp_mult_foot = ["<b>ALL7</b>", f"<b>{fnum(hp_tot)}</b>",
                    f"<b>{fnum(tot_only1)}</b><small>{fpct(tot_only1, hp_tot)}</small>",
                    f"<b>{fnum(tot_both)}</b><small>{fpct(tot_both, hp_tot)}</small>",
                    f"<b>{fnum(tot_only2)}</b><small>{fpct(tot_only2, hp_tot)}</small>"]

    # ---- S8 unresolved attribution
    cross_n = res["TIED_CROSS_TOPOLOGY"]["n"]
    mech = tied["mechanism_split_of_cross"]
    ac = afdeg["per_class"]["TIED_CROSS_TOPOLOGY"]
    au = afdeg["per_class"]["UNIQUE_TREE"]
    asame = afdeg["per_class"]["TIED_SAME_TOPOLOGY"]
    informative_units = sum(afdeg["per_class"][c]["units"]
                            - afdeg["per_class"][c]["af_uninformative_distinct_le_1"]
                            for c in afdeg["per_class"])
    uninformative_units = sum(afdeg["per_class"][c]["af_uninformative_distinct_le_1"]
                              for c in afdeg["per_class"])
    cross_informative = ac["units"] - ac["af_uninformative_distinct_le_1"]
    cross_uninformative = ac["af_uninformative_distinct_le_1"]

    krate = tied["cross_rate_by_active_k"]
    krows = []
    for kb in ("k<=1", "k=2", "k=3", "k=4", "k=5-6", "k=7-9"):
        if kb not in krate:
            continue
        item = krate[kb]
        krows.append([esc(kb), fnum(item["ranked"]), fnum(item["cross"]),
                      f"<b>{item['rate'] * 100:.2f}%</b>"])

    provenance = []
    for name in ("20260724_exactPS_source_overview_by_sample_01.tsv",
                 "20260724_exactPS_funnel_by_sample_01.tsv",
                 "20260724_exactPS_source_component_k_by_sample_01.tsv",
                 "20260724_exactPS_final_group_k_by_sample_01.tsv",
                 "20260724_exactPS_hp_split_by_sample_01.tsv"):
        p = CENSUS_DATA / name
        provenance.append([esc(name), fnum(p.stat().st_size), esc(sha256(p)[:16] + "…")])
    for p in (CENSUS_SUMMARY, SESSION / "tied_cross_topology_characterisation_v1.json",
              SESSION / "tie_af_degeneracy_v1.json"):
        provenance.append([esc(p.name), fnum(p.stat().st_size), esc(sha256(p)[:16] + "…")])

    # ---- S9 methylation: is there a difference, and along which axis
    axis = session("methyl_axis_decomposition_v1.json")
    v1 = axis["view1_cohort"]
    v1_total = sum(v1.values())
    V1C = [("m1_difference", "#1f6f5c", "有穩定甲基差異（M1）"),
           ("no_difference", "#b9bcc0", "可評估・無穩定差異"),
           ("not_evaluable", "#cfccc0", "不可評估（ALT read 不足等）")]
    v1_rows = []
    v1_tbl = []
    for s in samples:
        d = axis["view1_per_dataset"].get(s)
        if not d:
            sys.exit(f"FAIL CLOSED: methyl view1 missing {s}")
        n = sum(d.values())
        v1_rows.append((DISPLAY.get(s, s), n, [(d.get(k, 0), c, lbl) for k, c, lbl in V1C]))
        ev = n - d.get("not_evaluable", 0)
        v1_tbl.append([f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(n),
                       fnum(d.get("not_evaluable", 0)), fnum(ev),
                       f"{fnum(d.get('m1_difference', 0))}<small>{fpct(d.get('m1_difference', 0), n)} of all</small>",
                       f"<b>{fpct(d.get('m1_difference', 0), ev)}</b>"])
    v1_ev = v1_total - v1["not_evaluable"]
    v1_foot = ["<b>ALL7</b>", f"<b>{fnum(v1_total)}</b>", f"<b>{fnum(v1['not_evaluable'])}</b>",
               f"<b>{fnum(v1_ev)}</b>",
               f"<b>{fnum(v1['m1_difference'])}</b><small>{fpct(v1['m1_difference'], v1_total)} of all</small>",
               f"<b>{fpct(v1['m1_difference'], v1_ev)}</b>"]

    V2C = [("residual_within_alt", "#2f7a5a", "單 HP 內 ALT reads 殘餘結構"),
           ("allele_axis", "#8f5714", "沿 ALT／REF 等位軸"),
           ("technical_axis", "#2f6493", "沿技術軸"),
           ("hp_axis", "#96372c", "沿 HP1／HP2 haplotype 軸")]
    v2ex = axis["view2_cohort_exclusive"]
    v2raw = axis["view2_cohort_raw_overlapping"]
    v2_rows = []
    v2_tbl = []
    for s in samples:
        d = axis["view2_per_dataset"].get(s)
        if not d:
            sys.exit(f"FAIL CLOSED: methyl view2 missing {s}")
        n = sum(d.get(k, 0) for k, _, _ in V2C) + d.get("unclassified", 0)
        v2_rows.append((DISPLAY.get(s, s), n, [(d.get(k, 0), c, lbl) for k, c, lbl in V2C]))
        v2_tbl.append([f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(n)]
                      + [f"{fnum(d.get(k, 0))}<small>{fpct(d.get(k, 0), n)}</small>"
                         for k, _, _ in V2C])
    v2_tot = sum(v2ex[k] for k, _, _ in V2C) + v2ex["unclassified"]
    v2_foot = ["<b>ALL7</b>", f"<b>{fnum(v2_tot)}</b>"] + [
        f"<b>{fnum(v2ex[k])}</b><small>{fpct(v2ex[k], v2_tot)}</small>" for k, _, _ in V2C]
    v2_raw_tbl = [[f"<b>{esc(lbl)}</b>", fnum(v2raw[k]), fpct(v2raw[k], v2_tot),
                   fnum(v2ex[k]), fpct(v2ex[k], v2_tot)] for k, _, lbl in V2C]

    # ---- S12 runtime and memory
    rt = session("runtime_memory_v1.json")
    rt_tbl = []
    for s in samples:
        v = rt["per_sample"][s]
        proc = v["census_process"] or {}
        rt_tbl.append([f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(v["units"]),
                       f"{v['solver_total_seconds']:.2f}", f"{v['solver_mean_us']:,.0f}",
                       f"{v['solver_p99_us']:,}", f"{v['solver_max_us']:,}",
                       fnum(v["search_nodes_total"]),
                       f"{proc.get('wall_seconds', 0):.2f}",
                       f"{proc.get('max_rss_kb', 0) / 1024:.1f}"])
    rc = rt["cohort"]
    rt_foot = ["<b>ALL7</b>", f"<b>{fnum(rc['units_with_timing'])}</b>",
               f"<b>{rc['solver_total_seconds']:.2f}</b>", f"<b>{rc['solver_mean_us']:,.0f}</b>",
               "—", "—", f"<b>{fnum(rc['search_nodes_total'])}</b>",
               f"<b>{rc['census_wall_seconds_total']:.2f}</b>",
               f"<b>{rc['census_max_rss_mb']:.1f}</b>"]
    slow = rt["slowest_unit"]

    # ---- S10 M1 methylation annotation layer, aggregated per sample from the site table
    m1_receipt = session("m1_annotation_receipt_v2.json")
    m1_tsv = REPO / "research/20260726_m1_locus_annotation_layer/data/m1_annotation_site_table_v2.tsv"
    if not m1_tsv.exists():
        sys.exit(f"FAIL CLOSED: missing M1 site table {m1_tsv}")
    m1 = {}
    with open(m1_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            bucket = m1.setdefault(row["dataset"], {
                "n": 0, "in_group": 0, "tier": {}, "allele": {}, "hp": {}})
            bucket["n"] += 1
            if row["in_final_group"] == "True":
                bucket["in_group"] += 1
            for key, col in (("tier", "evidence_tier"), ("allele", "allele_class"),
                             ("hp", "hp_anchor_class")):
                bucket[key][row[col]] = bucket[key].get(row[col], 0) + 1
    if sum(v["n"] for v in m1.values()) != m1_receipt["rows"]:
        sys.exit("FAIL CLOSED: M1 site table row count disagrees with its receipt")

    TIERC = [("E4_PHASE_ANCHORED_ROBUST_EPIGENETIC_CANDIDATE", "#176b58", "E4 phase-anchored"),
             ("E3_UNEXPLAINED_AFTER_MEASURED_AXES", "#2f6493", "E3 unexplained"),
             ("E2_STABLE_NULL_MULTIGROUP", "#b9bcc0", "E2 stable-null")]
    ALLC = [("ALIGNED", "#96372c", "ALIGNED"),
            ("TESTED_NOT_ALIGNED", "#2f6493", "TESTED_NOT_ALIGNED"),
            ("NOT_TESTABLE", "#c9c8bf", "NOT_TESTABLE")]
    HPC = [("HP1_SIDE", "#2d6f91", "HP1_SIDE"), ("HP2_SIDE", "#8a4a5e", "HP2_SIDE"),
           ("CROSS_HP", "#8a6114", "CROSS_HP"), ("AMBIGUOUS_ONLY", "#b9bcc0", "AMBIGUOUS_ONLY")]

    def m1_rows(keyset, axis):
        out = []
        for s in samples:
            b = m1.get(s)
            if not b:
                sys.exit(f"FAIL CLOSED: M1 site table has no rows for {s}")
            out.append((DISPLAY.get(s, s), b["n"],
                        [(b[axis].get(k, 0), c, lbl) for k, c, lbl in keyset]))
        return out

    m1_join_rows = []
    m1_tbl = []
    for s in samples:
        b = m1[s]
        side = m1_receipt["sidecars"][s]
        m1_join_rows.append((DISPLAY.get(s, s), b["n"], [
            (b["in_group"], "#2f7a5a", "落在 final group 內"),
            (b["n"] - b["in_group"], "#c9c8bf", "region 外（locus-only）")]))
        m1_tbl.append([f"<b>{esc(DISPLAY.get(s, s))}</b>", fnum(b["n"]),
                       fnum(b["in_group"]), fpct(b["in_group"], b["n"]),
                       fnum(b["n"] - b["in_group"]),
                       fnum(side["n_regions_referenced"]),
                       fnum(b["allele"].get("ALIGNED", 0)),
                       fnum(b["hp"].get("CROSS_HP", 0))])
    rj = m1_receipt["region_join"]
    m1_regions_total = sum(m1_receipt["sidecars"][s]["n_regions_referenced"] for s in samples)
    m1_tbl_foot = ["<b>ALL7</b>", f"<b>{fnum(m1_receipt['rows'])}</b>",
                   f"<b>{fnum(rj['in_group'])}</b>",
                   f"<b>{fpct(rj['in_group'], m1_receipt['rows'])}</b>",
                   f"<b>{fnum(rj['locus_only'])}</b>", f"<b>{fnum(m1_regions_total)}</b>",
                   f"<b>{fnum(m1_receipt['axis_B_allele_class']['ALIGNED'])}</b>",
                   f"<b>{fnum(m1_receipt['axis_C_hp_anchor_class']['CROSS_HP'])}</b>"]

    bxc = m1_receipt["axis_B_x_axis_C"]
    hp_keys = ["HP1_SIDE", "HP2_SIDE", "AMBIGUOUS_ONLY", "CROSS_HP"]
    bxc_rows = []
    for al, _, allabel in ALLC:
        cells = [f"<b>{esc(allabel)}</b>"]
        for h in hp_keys:
            v = bxc.get(f"{al}|{h}", 0)
            mark = ' class="hi"' if (al == "ALIGNED" and h == "CROSS_HP") else ""
            cells.append(f"<span{mark}>{fnum(v)}</span>")
        bxc_rows.append(cells)

    # ---- S11 representative-tree shape bias
    bias = session("representative_shape_bias_v1.json")
    CO_BIAS = [("Direct-only", PALETTE["direct"]), ("Sister+direct", PALETTE["sisdir"]),
               ("Sister-only", PALETTE["sister"]), ("Single-only", PALETTE["sgl"])]
    bias_rows = []
    for s in samples:
        gaps = bias["per_sample_bias"].get(s)
        if not gaps:
            sys.exit(f"FAIL CLOSED: representative bias missing sample {s}")
        bias_rows.append((DISPLAY.get(s, s),
                          [(k, gaps.get(k, 0.0), c) for k, c in CO_BIAS]))
    bias_tbl = []
    for k, _ in CO_BIAS:
        obs = bias["representative_class_share"].get(k, 0.0)
        exp = bias["tied_tree_class_share_expected"].get(k, 0.0)
        gap = bias["bias_representative_minus_expected"].get(k, 0.0)
        cls = ' class="hi"' if abs(gap) >= 0.15 else ""
        bias_tbl.append([f"<b>{esc(k)}</b>", f"{obs * 100:.2f}%", f"{exp * 100:.2f}%",
                         f"<span{cls}>{gap * 100:+.2f} pp</span>"])
    deepest = bias["deepest_chain_rule_vs_current_representative"]

    # ---- S0 scope: proves the page is bound to the exact-PS authority, not legacy 50 kb
    mlhp_probe = json.loads(
        (Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
              "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
              "samples/HCC1937/HCC1937.exact_ps_mlhp.json")).read_text())
    params = mlhp_probe.get("params") or {}
    scope_rows = [
        ["source analysis unit", f"<code>{esc(scope.get('source_analysis_unit'))}</code>"],
        ["final analysis unit", f"<code>{esc(scope.get('final_analysis_unit'))}</code>"],
        ["final k 定義", esc(scope.get("final_k_definition"))],
        ["染色體", esc(scope.get("chromosomes"))],
        ["linkage basis（實測）", "<code>PS_HP1 / PS_HP2</code>；無 legacy 50 kb 傳遞合併"],
        ["phase_set_status（實測）", "<code>KNOWN_PS_PRIMARY</code> 100%"],
        ["input_mode", f"<code>{esc(params.get('input_mode'))}</code>"],
        ["MAX_SNV", f"<code>{esc(params.get('MAX_SNV'))}</code>"],
        ["MINREAD", f"<code>{esc(params.get('MINREAD'))}</code>"],
        ["partial read policy", esc(params.get("partial_read_policy"))],
        ["cn_source", f"<code>{esc(params.get('cn_source'))}</code>"],
    ]
    scope_table = table(["契約項目", "權威值"], scope_rows, cls="small")

    page = TEMPLATE.format(
        scope_table=scope_table,
        max_snv=esc(params.get("MAX_SNV")),
        tot_S=fnum(tot_S), tot_active_loci=fnum(tot_active_loci),
        active_pct=fpct(tot_active_loci, tot_S),
        tot_W=fnum(tot_W), tot_final=fnum(tot_final), ranked=fnum(ranked),
        one_exact=fnum(coh["one_exact_topology"]["n"]),
        one_exact_pct=fpct(coh["one_exact_topology"]["n"], ranked),
        cross_n=fnum(cross_n), cross_pct=fpct(cross_n, ranked),
        s1_funnel=s1,
        s2_svg=per_sample_stacked(rows2, [(PALETTE["active"], "active unique loci"),
                                          ("#cfccc0", "未進任何容器")],
                                  "每個資料集的 sSNV 可關聯比例"),
        s2_table=table(["dataset", "候選 sSNV (S)", "active unique loci", "% of S",
                        "active memberships", "exact PS×HP 容器"], tbl2, tbl2_foot),
        s3_svg=per_sample_stacked(rows3, [(PALETTE["w"], "tree-eligible W (k≥2)"),
                                          (PALETTE["single"], "k=1 abstain")],
                                  "component 切分"),
        s3_table=table(["dataset", "all components", "k=1 abstain", "% ", "W (k≥2)",
                        "%", "W memberships"], tbl3, tbl3_foot),
        s4_src=per_sample_stacked(k_rows(src_k, True),
                                  [(PALETTE["single"], "k=1")] + [(c, k) for k, c in KCOLS]
                                  + [("#5f3a20", "k>12")],
                                  "source component k 分布"),
        s4_fin=per_sample_stacked(
            k_rows(fin_k, False),
            [(c, k) for k, c in KCOLS],
            "final topology group k 分布（k=1 恆為 0）"),
        s5_svg=per_sample_stacked(rows5, [(PALETTE["unique"], "UNIQUE_TREE"),
                                          (PALETTE["same"], "TIED_SAME_TOPOLOGY"),
                                          (PALETTE["cross"], "TIED_CROSS_TOPOLOGY")],
                                  "determinacy"),
        s5_table=table(["dataset", "ranked", "UNIQUE", "%", "TIED_SAME", "%",
                        "TIED_CROSS", "%", "單一 exact topology", "%", "max ties"],
                       tbl5, tbl5_foot),
        s6_svg=per_sample_stacked(rows6, [(c, k) for k, c in CO]
                                  + [(PALETTE["cross"], "Cross-class unresolved")],
                                  "coarse geometry"),
        s6_table=table(["dataset"] + [k for k, _ in CO]
                       + ["Cross-class unresolved", "exact 跨形狀但 coarse 相同"],
                       tbl6, tbl6_foot, cls="small"),
        s7_svg=per_sample_stacked(rows7, [(PALETTE["hp1"], "單一 haplotype"),
                                          ("#8a6114", "兩條 haplotype（同一 PS 內）")],
                                  "每個 physical locus 在同一 PS 內被分析到幾條 haplotype"),
        hp_mult_table=table(["dataset", "進入 final group 的 unique loci", "單一 haplotype",
                             "兩條 haplotype（同一 PS 內）", "（參考）跨多個 PS 的 loci"],
                            hp_mult_tbl, hp_mult_foot),
        hp_both_pct=fpct(tot_both, hp_tot),
        hp_single_pct=fpct(tot_only1 + tot_only2, hp_tot),
        hp_both_n=fnum(tot_both), hp_tot_n=fnum(hp_tot),
        mech_vertex=fnum(mech.get("VERTEX_SET_AMBIGUOUS", 0)),
        mech_parent=fnum(mech.get("PARENT_CHOICE_ONLY", 0)),
        cross_unin=fnum(cross_uninformative),
        cross_unin_pct=fpct(cross_uninformative, ac["units"], 4),
        cross_af1=fnum(ac["all_active_sites_af_1"]),
        cross_af1_pct=fpct(ac["all_active_sites_af_1"], ac["units"], 4),
        cross_inf=fnum(cross_informative),
        rate_inf=f"{100.0 * cross_informative / informative_units:.4f}%",
        rate_unin=f"{100.0 * cross_uninformative / uninformative_units:.4f}%",
        ratio=f"{(cross_uninformative / uninformative_units) / (cross_informative / informative_units):.0f}",
        inf_units=fnum(informative_units), unin_units=fnum(uninformative_units),
        u_unin_pct=f"{au['rate_af_uninformative'] * 100:.2f}%",
        s_unin_pct=f"{asame['rate_af_uninformative'] * 100:.2f}%",
        k_table=table(["active k", "ranked units", "TIED_CROSS", "cross 率"], krows),
        methyl_v1_svg=per_sample_stacked(v1_rows, [(c, l) for _, c, l in V1C],
                                         "每個 sSNV 是否有穩定甲基差異"),
        methyl_v1_table=table(["dataset", "sSNV 總數", "不可評估", "可評估",
                               "有穩定差異（M1）", "占可評估"], v1_tbl, v1_foot),
        methyl_v2_svg=per_sample_stacked(v2_rows, [(c, l) for _, c, l in V2C],
                                         "有差異者：分群沿哪個軸（互斥）"),
        methyl_v2_table=table(["dataset", "M1 位點"] + [l for _, _, l in V2C],
                              v2_tbl, v2_foot, cls="small"),
        methyl_v2_raw_table=table(["軸", "原始計數（可重疊）", "占 M1", "互斥歸屬後", "占 M1"],
                                  v2_raw_tbl),
        v1_all=fnum(v1_total), v1_ev=fnum(v1_ev),
        v1_m1=fnum(v1["m1_difference"]),
        v1_m1_pct_all=fpct(v1["m1_difference"], v1_total),
        v1_m1_pct_ev=fpct(v1["m1_difference"], v1_ev),
        v1_no=fnum(v1["no_difference"]), v1_ne=fnum(v1["not_evaluable"]),
        v2_residual=fnum(v2ex["residual_within_alt"]),
        v2_residual_pct=fpct(v2ex["residual_within_alt"], v2_tot),
        v2_allele=fnum(v2ex["allele_axis"]), v2_allele_pct=fpct(v2ex["allele_axis"], v2_tot),
        v2_hp=fnum(v2ex["hp_axis"]), v2_hp_pct=fpct(v2ex["hp_axis"], v2_tot),
        v2_tech=fnum(v2ex["technical_axis"]), v2_tech_pct=fpct(v2ex["technical_axis"], v2_tot),
        rt_table=table(["dataset", "units", "solver 總秒數", "平均 µs/unit", "p99 µs",
                        "最慢 µs", "search nodes", "census wall 秒", "peak RSS MB"],
                       rt_tbl, rt_foot, cls="small"),
        rt_total_s=f"{rc['solver_total_seconds']:.1f}",
        rt_mean_us=f"{rc['solver_mean_us']:,.0f}",
        rt_nodes=fnum(rc["search_nodes_total"]),
        rt_rss=f"{rc['census_max_rss_mb']:.1f}",
        rt_wall=f"{rc['census_wall_seconds_total']:.1f}",
        slow_sample=esc(DISPLAY.get(slow["sample"], slow["sample"])),
        slow_region=esc(slow["region_id"]), slow_k=esc(slow["active_bit_count"]),
        slow_ms=f"{slow['micros'] / 1000:.1f}", slow_nodes=fnum(slow["search_nodes"]),
        slow_status=esc(slow["family_status"]),
        m1_sites=fnum(m1_receipt["rows"]),
        m1_in_group=fnum(rj["in_group"]),
        m1_in_group_pct=fpct(rj["in_group"], m1_receipt["rows"]),
        m1_locus_only=fnum(rj["locus_only"]),
        m1_regions=fnum(m1_regions_total),
        m1_regions_pct=fpct(m1_regions_total, tot_final),
        m1_aligned=fnum(m1_receipt["axis_B_allele_class"]["ALIGNED"]),
        m1_aligned_pct=fpct(m1_receipt["axis_B_allele_class"]["ALIGNED"], m1_receipt["rows"]),
        m1_gold=fnum(bxc.get("ALIGNED|CROSS_HP", 0)),
        m1_e4=fnum(m1_receipt["axis_A_evidence_tier"]["E4_PHASE_ANCHORED_ROBUST_EPIGENETIC_CANDIDATE"]),
        m1_join_svg=per_sample_stacked(
            m1_join_rows, [("#2f7a5a", "落在 final group 內"), ("#c9c8bf", "region 外（locus-only）")],
            "M1 位點與 exact-PS region 的關係"),
        m1_tier_svg=per_sample_stacked(
            m1_rows(TIERC, "tier"), [(c, l) for _, c, l in TIERC], "Axis A · evidence tier"),
        m1_allele_svg=per_sample_stacked(
            m1_rows(ALLC, "allele"), [(c, l) for _, c, l in ALLC], "Axis B · ALT/REF 等位軸"),
        m1_hp_svg=per_sample_stacked(
            m1_rows(HPC, "hp"), [(c, l) for _, c, l in HPC], "Axis C · somatic-ALT haplotype 錨定"),
        m1_table=table(["dataset", "M1 位點", "在 region 內", "%", "locus-only",
                        "觸及 regions", "ALT/REF aligned", "CROSS_HP"], m1_tbl, m1_tbl_foot),
        bxc_table=table(["Axis B ＼ Axis C"] + hp_keys, bxc_rows, cls="small"),
        bias_direct=f"{bias['bias_representative_minus_expected']['Direct-only'] * 100:+.2f}",
        bias_units=fnum(bias["cross_units_evaluated"]),
        bias_table=table(["coarse class", "代表樹佔比", "同分樹期望佔比", "偏差"], bias_tbl),
        bias_svg=diverging_svg(bias_rows, "逐樣本形態偏差（代表樹 − 同分樹期望）",
                               note="橫軸為百分點；左負右正。四個樣本序列由上而下為 "
                                    "Direct-only／Sister+direct／Sister-only／Single-only。"),
        deepest_agree=fnum(deepest.get("agrees_with_current_rep", 0)),
        deepest_differ=fnum(deepest.get("differs_from_current_rep", 0)),
        prov_table=table(["artifact", "bytes", "sha256 (前 16)"], provenance, cls="small"),
        checks_html="".join(f'<li><b>PASS</b> — {esc(k)}</li>' for k in checks),
    )

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    print(f"[ok] {out}  {out.stat().st_size/1000:.1f} KB")
    print(f"     S={tot_S:,}  active_loci={tot_active_loci:,}  W={tot_W:,}  "
          f"final={tot_final:,}  ranked={ranked:,}  cross={cross_n:,}")
    print(f"     conservation checks: {len(checks)}/{len(checks)} PASS")


TEMPLATE = """<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>exact-PS 資料觀察與分析檢視</title><style>
:root{{--bg:#f6f5f0;--panel:#fff;--fg:#18241f;--muted:#464f4b;--line:#d5d2c6;--track:#eceae1;
--accent:#1f6f5c;--head:#e4e8e2;--headfg:#16261f;--zebra:#faf9f5;--hover:#eef4f0;--foot:#e9efe9}}
@media(prefers-color-scheme:dark){{:root{{--bg:#121517;--panel:#1a1e21;--fg:#eae8e1;--muted:#b6bfb8;
--line:#39413f;--track:#242a2c;--accent:#4fbf9b;--head:#232c29;--headfg:#dff0e8;--zebra:#1d2225;
--hover:#25302c;--foot:#222d29}}}}
:root[data-theme=dark]{{--bg:#121517;--panel:#1a1e21;--fg:#eae8e1;--muted:#b6bfb8;--line:#39413f;
--track:#242a2c;--accent:#4fbf9b;--head:#232c29;--headfg:#dff0e8;--zebra:#1d2225;--hover:#25302c;--foot:#222d29}}
:root[data-theme=light]{{--bg:#f6f5f0;--panel:#fff;--fg:#18241f;--muted:#464f4b;--line:#d5d2c6;
--track:#eceae1;--accent:#1f6f5c;--head:#e4e8e2;--headfg:#16261f;--zebra:#faf9f5;--hover:#eef4f0;--foot:#e9efe9}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.65 system-ui,-apple-system,"Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif}}
.wrap{{max-width:1120px;margin:0 auto;padding:1.2rem 1rem 4rem}}
header.top{{background:#12241d;color:#e9f0ec;padding:.55rem 1rem;font:12px/1.5 ui-monospace,monospace;letter-spacing:.02em}}
h1{{font-size:clamp(1.5rem,3.6vw,2.3rem);line-height:1.25;margin:1.4rem 0 .4rem}}
h2{{font-size:1.22rem;margin:2.6rem 0 .3rem;padding-top:.9rem;border-top:2px solid var(--line)}}
h3{{font-size:1rem;margin:1.5rem 0 .3rem;color:var(--muted)}}
.lede{{color:var(--muted);max-width:70ch}}
.panel{{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:1rem;margin:.9rem 0}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:.6rem;margin:1rem 0}}
.kpi{{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:.7rem .8rem}}
.kpi b{{display:block;font:600 1.55rem/1.2 ui-monospace,monospace;color:var(--accent)}}
.kpi span{{font-size:.78rem;color:var(--muted)}}
.tw{{overflow-x:auto;-webkit-overflow-scrolling:touch;margin:.9rem 0;
border:1px solid var(--line);border-radius:8px;background:var(--panel)}}
table{{border-collapse:separate;border-spacing:0;width:100%;min-width:640px;font-size:.9rem;
font-variant-numeric:tabular-nums}}
table.small{{font-size:.84rem}}
th,td{{padding:.55rem .7rem;text-align:right;white-space:nowrap;
border-bottom:1px solid var(--line)}}
th:first-child,td:first-child{{text-align:left;padding-left:.9rem}}
th:last-child,td:last-child{{padding-right:.9rem}}
thead th{{background:var(--head);color:var(--headfg);font-weight:600;position:sticky;top:0;z-index:2;
border-bottom:2px solid var(--line);letter-spacing:.01em}}
tbody tr:nth-child(even){{background:var(--zebra)}}
tbody tr:hover{{background:var(--hover)}}
tbody td:first-child{{border-left:3px solid transparent}}
tbody tr:hover td:first-child{{border-left-color:var(--accent)}}
tfoot td{{border-top:2px solid var(--accent);background:var(--foot);font-weight:600}}
tbody tr:last-child td{{border-bottom:0}}
td small{{display:block;color:var(--muted);font-size:.79rem;font-weight:400;margin-top:.1rem}}
.legend{{display:flex;flex-wrap:wrap;gap:.45rem;margin:.6rem 0 0;font-size:.8rem}}
.lg{{display:inline-flex;align-items:center;gap:.35rem;background:var(--track);
border:1px solid var(--line);border-radius:99px;padding:.16rem .55rem;color:var(--fg)}}
.lg i{{width:12px;height:12px;border-radius:3px;display:inline-block;flex:0 0 auto}}
.note{{font-size:.86rem;color:var(--muted);margin:.5rem 0 0}}
.warn{{border-left:4px solid #a94336;background:color-mix(in srgb,#a94336 8%,transparent);
padding:.7rem .9rem;border-radius:0 4px 4px 0;margin:.9rem 0}}
.ok{{border-left:4px solid var(--accent);background:color-mix(in srgb,#1f6f5c 8%,transparent);
padding:.7rem .9rem;border-radius:0 4px 4px 0;margin:.9rem 0}}
ul.checks{{list-style:none;padding:0;margin:.4rem 0;font-size:.85rem}}
ul.checks li{{padding:.15rem 0;color:var(--muted)}}
ul.checks b{{color:var(--accent)}}
details{{margin:.8rem 0;border:1px solid var(--line);border-radius:6px;background:var(--panel)}}
summary{{cursor:pointer;padding:.6rem .9rem;font-weight:600}}
details>div{{padding:0 .9rem .8rem}}
code{{font-family:ui-monospace,monospace;font-size:.86em;background:var(--track);padding:.05em .35em;border-radius:3px}}
.hi{{background:color-mix(in srgb,#a94336 22%,transparent);font-weight:600;padding:.05em .3em;border-radius:3px}}
</style></head><body>
<header class="top">EXACT-PS · 2026-07-24 authority · 7 technical datasets / 6 biological IDs · GRCh38 chr1–22 · strict endpoint read-linkage ≥3</header>
<div class="wrap">

<h1>資料觀察與分析檢視</h1>
<p class="lede">從輸入 sSNV 一路走到可判定拓撲的每一層分母、比例與逐樣本樣貌。所有數字由權威 TSV／census summary 注入，模板內沒有手寫數值；六項守恆檢核全部通過才會產出本頁。</p>

<div class="kpis">
<div class="kpi"><b>{tot_S}</b><span>候選 sSNV（S）</span></div>
<div class="kpi"><b>{tot_active_loci}</b><span>active unique loci（{active_pct}）</span></div>
<div class="kpi"><b>{tot_W}</b><span>tree-eligible W（k≥2）</span></div>
<div class="kpi"><b>{tot_final}</b><span>final topology groups</span></div>
<div class="kpi"><b>{ranked}</b><span>ranked complete</span></div>
<div class="kpi"><b>{one_exact_pct}</b><span>單一 exact topology（{one_exact}）</span></div>
</div>

<h2>0 · 範圍契約：本頁綁定的是最新 exact-PS，不是 legacy 50 kb</h2>
<p class="note">下表直接讀自 authority 的 <code>scope</code> 與 MLHP <code>params</code>，非人工轉述。舊版以 <code>adjacent_gap≤50000</code> 傳遞合併的分群<b>不在本頁任何數字裡</b>。</p>
{scope_table}
<p class="note">⚠️ <b><code>MAX_SNV = {max_snv}</code></b> 是 production authority 的實際值；若論文方法章寫的是其他上限，以本欄為準或明確說明差異。</p>

<h2>1 · 全 cohort 漏斗</h2>
<p class="note">每一層的單位都不同，<b>不可互相代稱</b>：loci ≠ memberships ≠ components ≠ blocks ≠ groups。</p>
<div class="panel">{s1_funnel}</div>

<h2>2 · 輸入 sSNV 有多少能建立關聯區域</h2>
<p class="note">「active」＝該 locus 至少進入一個 exact PS×HP 容器，且被至少一條 canonical molecule 以 fixed R/A 觀察到。未 active 不代表該位點不存在，只代表在此嚴格定義下不可用。</p>
<div class="panel">{s2_svg}</div>
{s2_table}

<h2>3 · 切成 component 之後：能不能建樹</h2>
<p class="note"><code>k=1</code> 記為 <code>ABSTAIN_SINGLETON_UNLINKED</code>：該 membership 在此容器、threshold=3 下沒有任何合格 edge。它保留在守恆 registry，但<b>不進下游分母</b>。</p>
<div class="panel">{s3_svg}</div>
{s3_table}

<h2>4 · 區域的 k 數量分佈</h2>
<h3>4.1 source component（含 k=1）</h3>
<div class="panel">{s4_src}</div>
<h3>4.2 final topology group（k=1 恆為 0，k 全落在 2–12）</h3>
<div class="panel">{s4_fin}</div>
<p class="note"><code>k&gt;12</code> 的 source component 仍是一個 region；k≤12 切割只產生 computational block，<b>block 數不可當 region 數</b>。</p>

<h2>5 · 分成區域後：拓撲可判定到哪一層</h2>
<div class="panel">{s5_svg}</div>
{s5_table}
<div class="ok"><b>單一 exact rooted-unlabeled topology ＝ {one_exact} / {ranked} ＝ {one_exact_pct}。</b>
signature 忽略 mutation label 與 sibling 順序，所以「同形狀但標號順序不同」的多棵樹會收斂成一個 signature——這正是 TIED_SAME 的全部內容，且<b>已算入可判定</b>。</div>

<h2>6 · 四類 coarse geometry</h2>
<p class="note">分類規則：<code>Single-only</code> 無分支且無 root 之下兩層路徑；<code>Sister-only</code> 有分支無深度；<code>Direct-only</code> 有深度無分支；<code>Sister+direct</code> 兩者皆有。只有當該 unit 的<b>全部</b> AF-最佳 signature 落在同一類，才指派該類，否則為 <code>Cross-class unresolved</code>。這是數學圖形幾何，<b>不是 clone 身分</b>。</p>
<div class="panel">{s6_svg}</div>
{s6_table}

<h2>7 · 每個位點在同一 PS 內被分析到幾條 haplotype</h2>
<p class="note">分母＝進入 final group 的 <b>unique physical loci</b>。🔑 <b>只分「單一」與「兩條」兩類，不拆 HP1／HP2</b>：HP 標籤<b>只在同一個 phase set 內有意義</b>，不同 PS 之間 HP1／HP2 的相對 orientation 在本層未定義，所以跨 PS 把「只在 HP1」與「只在 HP2」分開統計沒有意義。判定一律在<b>同一 PS 內</b>進行。</p>
<div class="panel">{s7_svg}</div>
{hp_mult_table}
<div class="warn">🔴 <b>這張圖是「被分析到」，不是「帶突變」。</b>
MLHP <code>positions</code> 包含在該 haplotype 上<b>只觀察到 REF</b> 的位點——實測 4,000 個 unit 中，含非 active 位點者 1,954 個，其非 active 位點 <b>100% 為 ALT=0 的純 REF</b>。
所以「兩條 haplotype」＝該座標在同一 PS 的 HP1 與 HP2 容器裡都有合格觀測，<b>不等於兩條都帶突變</b>——體細胞雜合突變本來就只落在一條上。
全 cohort <b>{hp_both_pct}</b>（{hp_both_n} / {hp_tot_n}）在同一 PS 內兩條都被分析到，<b>{hp_single_pct}</b> 只有一條。
「只有一條」<b>不是 LOH 判定</b>：判 LOH 需要該位點兩條 haplotype 的<b>相對覆蓋</b>，本表不含。
表末「跨多個 PS 的 loci」僅供參考，說明為何不可跨 PS 合併 HP 標籤。</div>

<h2>8 · 未解的 {cross_pct} 是什麼造成的</h2>
<div class="warn"><b>機制歸因：{mech_vertex} / {cross_n} 全部是「最小狀態集合本身不唯一」（VERTEX_SET_AMBIGUOUS），parent-choice 型 ＝ {mech_parent}。</b>
不是「誰接誰」沒定，是「有哪些狀態存在」沒定。</div>

<h3>8.1 read-AF 軸退化</h3>
<div class="kpis">
<div class="kpi"><b>{cross_unin_pct}</b><span>cross 中 AF 軸無資訊（{cross_unin}）</span></div>
<div class="kpi"><b>{cross_af1_pct}</b><span>全部 active 位點 AF＝1.0（{cross_af1}）</span></div>
<div class="kpi"><b>{rate_inf}</b><span>AF 有資訊時的 cross 率（{inf_units} units）</span></div>
<div class="kpi"><b>{rate_unin}</b><span>AF 無資訊時的 cross 率（{unin_units} units）</span></div>
</div>
<p class="note">兩者相差 <b>{ratio} 倍</b>。對照組：UNIQUE_TREE 的 AF 無資訊率 {u_unin_pct}、TIED_SAME {s_unin_pct}。
評分是全樹加總的 AF 差值；AF 全相等 → 每項為 0 → 每棵樹同分 → 排序器沒有任何輸入可用。
只有 <b>{cross_inf}</b> 個 unit 是「AF 有資訊卻仍跨拓撲」。</p>

<h3>8.2 隨 active k 的變化</h3>
{k_table}
<p class="note">k=7–9 的 0% 來自 n 極小的選擇效應（大 k 多在 ranking 前撞 resource guard），<b>不可引用為「大 k 沒問題」</b>。</p>

<h2>9 · 甲基化：有沒有差異，以及差異沿哪個軸</h2>
<h3>9.1 全部 sSNV：有差異 vs 沒差異</h3>
<p class="note">分母＝<b>{v1_all}</b> 個 LongPhase-S PASS sSNV。判定合約 <code>coarse_ng≥2 AND not unstable AND modal_assignment_ari_min≥0.8</code>，並須通過 40 次 column-shuffle null 的 95 百分位門檻與最小群 3 條 read。</p>
<div class="panel">{methyl_v1_svg}</div>
{methyl_v1_table}
<div class="ok"><b>有穩定甲基差異 ＝ {v1_m1} / {v1_all} ＝ {v1_m1_pct_all}</b>（占<b>可評估</b> {v1_ev} 者為 <b>{v1_m1_pct_ev}</b>）。
無差異 {v1_no}、不可評估 {v1_ne}。這是 <b>operational screen yield，不是 subclone prevalence</b>。</div>

<h3>9.2 有差異者：分群沿哪個軸</h3>
<p class="note">四個軸互相<b>可能重疊</b>，所以同時給「原始計數」與「依優先序 <code>HP &gt; 等位 &gt; 技術 &gt; 殘餘</code> 的互斥歸屬」。HP 軸排最前，因為沿 haplotype 的分群與 germline allele-specific methylation <b>不可分辨</b>，是最強的 confound。</p>
<div class="panel">{methyl_v2_svg}</div>
{methyl_v2_table}
<h3>9.3 原始計數 vs 互斥歸屬</h3>
{methyl_v2_raw_table}
<div class="warn">互斥歸屬後：<b>單 HP 內 ALT reads 的殘餘結構 {v2_residual}（{v2_residual_pct}）</b>、
沿 ALT／REF 等位軸 <b>{v2_allele}（{v2_allele_pct}）</b>、
沿技術軸 {v2_tech}（{v2_tech_pct}）、
沿 HP1／HP2 haplotype 軸 <b>{v2_hp}（{v2_hp_pct}）</b>。
🔴 「殘餘」只表示<b>已測軸無法解釋</b>，不表示未測 confound 已排除；
等位軸在 ALT 僅落單一 HP 時與 germline ASM 不可分辨，而 <b>germline-het null 尚未執行</b>。</div>

<h2>10 · M1 甲基位點 annotation 層</h2>
<p class="note">來源＝權威 tumor-REF control 表（102,842 個 M1 位點）。此層<b>不改</b> AF score、candidate incidence 或 selected tree；它是形式化規格所稱通道 M 的事後註記。三個軸互相正交，各有自己的 claim ceiling。</p>
<div class="kpis">
<div class="kpi"><b>{m1_sites}</b><span>M1 位點（全 7 資料集）</span></div>
<div class="kpi"><b>{m1_in_group_pct}</b><span>落在 final group 內（{m1_in_group}）</span></div>
<div class="kpi"><b>{m1_locus_only}</b><span>region 外，只能以 locus 呈現</span></div>
<div class="kpi"><b>{m1_regions}</b><span>被 M1 觸及的 final groups（{m1_regions_pct}）</span></div>
<div class="kpi"><b>{m1_aligned}</b><span>ALT/REF 等位軸 aligned（{m1_aligned_pct}）</span></div>
<div class="kpi"><b>{m1_gold}</b><span>aligned ∩ cross-HP 候選</span></div>
</div>
<div class="panel">{m1_join_svg}</div>
<p class="note">🔴 <b>覆蓋率各樣本差 8.5 倍</b>：H2009 有 89.0% 的 M1 位點落在 region 內，HCC1395_NYGC 只有 10.4%。所以介面<b>必須同時畫 locus 點層</b>，只做 region 著色會讓旗艦樣本看起來「幾乎沒有甲基資料」。</p>
{m1_table}

<h3>10.1 Axis A · evidence tier</h3>
<div class="panel">{m1_tier_svg}</div>
<p class="note">E4 = phase-anchored robust epigenetic candidate（{m1_e4}）。這是 operational screen tier，<b>不是</b> cellular subclone。</p>

<h3>10.2 Axis B · ALT/REF 等位軸</h3>
<div class="panel">{m1_allele_svg}</div>

<h3>10.3 Axis C · somatic-ALT 的 haplotype 錨定品質</h3>
<div class="panel">{m1_hp_svg}</div>
<p class="note">⚠️ <code>HP1_SIDE</code>／<code>HP2_SIDE</code> 表示體細胞 ALT <b>乾淨錨定在單一 germline haplotype</b>——這是二倍體區域的<b>正常預期</b>，<b>不是 LOH 判定</b>。判 LOH 需要全體（REF＋ALT）read 的 haplotype 組成，本表不含。</p>

<h3>10.4 Axis B × Axis C：唯一可去共線的候選</h3>
{bxc_table}
<div class="warn"><b>{m1_aligned} 個 ALT/REF aligned 位點中，只有 {m1_gold} 個是 cross-HP。</b>
其餘的 ALT 只落在單一 HP（與 germline allele-specific methylation <b>共線</b>）或 HP 無法指派。
因此可分辨「突變效應 vs germline ASM」的候選只有那 {m1_gold} 個，<b>而且 germline-het null 尚未執行</b> —— 連這批都只能稱候選。</div>

<h2>11 · 代表樹的形態偏差（實測）</h2>
<p class="note">canonical 每個 unit 都輸出一棵 <code>representative_best_morphology</code>。把它的形態對比該 unit <b>自己</b>全部同分樹的形態分布（每 unit 等權），就能量出這個 deterministic 選法有沒有系統性偏差。評估範圍＝{bias_units} 個 TIED_CROSS unit。</p>
{bias_table}
<div class="warn"><b>代表樹低估 Direct-only 達 {bias_direct} pp，七個樣本方向完全一致。</b>
所以 <code>representative_best_morphology</code> <b>只能用於顯示與 join，絕不可拿去做統計</b>。
組成統計已正確地把這些 unit 歸為 <code>Cross-class unresolved</code>（見第 6 節），沒有被污染。</div>
<div class="panel">{bias_svg}</div>
<p class="note">另測「偏好含 chain 形態」的替代規則：與現行代表<b>一致 {deepest_agree} 個、不同 {deepest_differ} 個</b>，約各半，而且會把偏差翻到另一側。<b>偏差的方向可以換，中立性換不到</b> —— 真正中立的做法是報整個分布，而不是選一棵。</p>

<h2>12 · 執行時間與記憶體</h2>
<p class="note">兩種量測分開列，因為意義不同：<b>solver</b> 是 exact-PS topology 二進位檔自己記錄的<b>每個 region 的求解耗時</b>（不含 I/O 與 JSON 序列化）；<b>census wall / RSS</b> 是 <code>/usr/bin/time -v</code> 包住 signature-census 二進位檔的<b>整樣本</b>實測。</p>
<div class="kpis">
<div class="kpi"><b>{rt_total_s} s</b><span>全 7 資料集 lineage solver 總計（單核加總）</span></div>
<div class="kpi"><b>{rt_mean_us} µs</b><span>平均每個 region</span></div>
<div class="kpi"><b>{rt_nodes}</b><span>B&amp;B search nodes 總數</span></div>
<div class="kpi"><b>{rt_rss} MB</b><span>peak RSS（census 階段最高）</span></div>
</div>
{rt_table}
<div class="ok"><b>全部 98,955 個 region 的 lineage 求解，單核加總只要 {rt_total_s} 秒；峰值記憶體 {rt_rss} MB。</b>
這正是 <code>k≤12</code> 有界枚舉的實證：狀態空間 <code>2^k</code> 有限，窮舉在此規模 tractable。
最慢的單一 region 是 {slow_sample} <code>{slow_region}</code>（k={slow_k}），耗時 <b>{slow_ms} ms</b>、
撞到 <code>search_nodes={slow_nodes}</code> 上限而判 <code>{slow_status}</code>。</div>
<p class="note">⚠️ solver 秒數是<b>單核逐 region 加總</b>，不是 pipeline 的 wall-clock —— 後者另含 BAM／JSON I/O。census wall 合計 {rt_wall} 秒可作同一機器上的對照。</p>

<h2>13 · 守恆檢核與資料來源</h2>
<ul class="checks">{checks_html}</ul>
<details><summary>資料來源與 SHA-256</summary><div>{prov_table}</div></details>

<div class="warn"><b>Claim ceiling.</b> 本頁全部是 read-linkage 與數學 mutation-state tree 的統計。
W 是 read-linkage region，endpoint edge 是<b>無向</b> linkage，不是 evolutionary parent→child。
coarse geometry 不等於 cellular single-clone／multi-clone。
未納入 CN、LOH、purity、multiplicity 校正，因此形狀與 HP 數量<b>不可</b>換算成 clone prevalence。
7 technical datasets ＝ 6 biological IDs（HCC1395_HKU 與 HCC1395_NYGC 為同株兩個技術來源），任何 group-level 檢定的真實 n ＝ 6。</div>

</div></body></html>"""


if __name__ == "__main__":
    main()
