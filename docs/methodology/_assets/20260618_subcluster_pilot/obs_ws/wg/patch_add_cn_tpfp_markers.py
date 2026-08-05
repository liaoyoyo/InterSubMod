#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
patch_add_cn_tpfp_markers.py — add CN-count + TP/FP markers to the WG dashboard.

Reads enriched_cn_tpfp.json (produced by join_cn_tpfp_markers.py — real coordinate
joins, §13 no fabrication), merges its fields into the embedded D cases, and injects:
  - badges per card : TP/FP · k_CN·CN(SAVANA) · SEQC2 dir+CN
  - selectable filters (control bar) : TP/FP · k_CN · SEQC2
  - matches() filter clauses + listener wiring
  - modal full CN detail line
  - badge CSS colours

Usage: python3 patch_add_cn_tpfp_markers.py            # -> <dash>.patched2
       python3 patch_add_cn_tpfp_markers.py --inplace  # backup + overwrite
"""
import json, os, sys, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "20260618_HCC1395_observation_dashboard_FULL.standalone.html")
ENR = os.path.join(HERE, "enriched_cn_tpfp.json")

CSS = (".b-tp{background:#064e3b;color:#6ee7b7}.b-fp{background:#450a0a;color:#fca5a5}"
       ".b-cn{background:#0c2a4d;color:#7dd3fc}.b-seq{background:#2a1a40;color:#d8b4fe}")

FILTERS = (
 '<div class="cg"><label>TP/FP</label><select id="ftp"><option value="all">全部</option>'
 '<option value="TP">TP</option><option value="FP">FP</option></select></div>\n'
 ' <div class="cg"><label>k_CN(SAVANA)</label><select id="fkcn"><option value="all">全部</option>'
 '<option value="1">1 (LOH)</option><option value="2">2 (het)</option></select></div>\n'
 ' <div class="cg"><label>SEQC2</label><select id="fseq"><option value="all">全部</option>'
 '<option value="gain">gain</option><option value="loss">loss</option>'
 '<option value="loh">loh</option><option value="neutral">neutral</option></select></div>\n ')

BADGE_VARS = (
 " const tpfp=(c.tp_fp==='FP')?'<span class=\"badge b-fp\" title=\"kism false-positive somatic\">FP</span>'"
 ":'<span class=\"badge b-tp\" title=\"kism true-positive somatic（本 dashboard 全為 TP）\">TP</span>';\n"
 " const cnb=(c.k_cn!=null)?`<span class=\"badge b-cn\" title=\"SAVANA：k_CN=CN推定allele狀態數"
 "(minorCN≥0.5→2 het / <0.5→1 LOH)·總CN ${c.savana_cn}(minor ${c.savana_minor})\">k_CN ${c.k_cn}·CN ${c.savana_cn}</span>`:'';\n"
 " const seqcb=c.seqc2_dir?`<span class=\"badge b-seq\" title=\"SEQC2 truth CN/方向\">SEQC2 ${c.seqc2_dir}"
 "${(c.seqc2_cn!=null&&c.seqc2_cn!=='')?(' '+c.seqc2_cn):''}</span>`:'';\n")

MATCH_CLAUSES = (
 "\n const tpf=$('ftp')?$('ftp').value:'all'; if(tpf!=='all'&&(c.tp_fp||'TP')!==tpf)return false;"
 "\n const kc=$('fkcn')?$('fkcn').value:'all'; if(kc!=='all'&&String(c.k_cn)!==kc)return false;"
 "\n const sq=$('fseq')?$('fseq').value:'all'; if(sq!=='all'&&(c.seqc2_dir||'')!==sq)return false;")

MODAL_DETAIL = (
 '<span>TP/FP <b>${c.tp_fp||\'TP\'}</b></span>'
 '<span>k_CN <b>${c.k_cn??\'NA\'}</b></span>'
 '<span>SAVANA CN <b>${c.savana_cn??\'NA\'}</b>(minor ${c.savana_minor??\'NA\'})</span>'
 '<span>SEQC2 <b>${c.seqc2_dir||\'NA\'} ${c.seqc2_cn??\'\'}</b></span>')


def patch(s, enr):
    if 'id="ftp"' in s:
        raise SystemExit('ABORT: already patched (id="ftp" present).')

    # --- 0. merge enriched fields into embedded D cases ---
    tag = '<script id="d" type="application/json">'
    a = s.find(tag) + len(tag)
    b = s.find("</script>", a)
    D = json.loads(s[a:b])
    miss = 0
    for c in D["cases"]:
        e = enr.get(c["id"])
        if e:
            c.update(e)
        else:
            miss += 1
    new_d = json.dumps(D, ensure_ascii=False)
    s = s[:a] + new_d + s[b:]

    def repl(hay, old, new, label, n=1):
        cnt = hay.count(old)
        if cnt != n:
            raise SystemExit(f"ABORT: anchor [{label}] found {cnt} times (expected {n}).")
        return hay.replace(old, new)

    # --- 1. CSS ---
    s = repl(s, "</style>", CSS + "</style>", "css")
    # --- 2. control-bar filters after CN(fc) group ---
    fc_end = '<option value="Low">Low</option><option value="Normal">Normal</option></select></div>'
    s = repl(s, fc_end, fc_end + "\n " + FILTERS, "ctrl/filters")
    # --- 3. badge vars (insert before ` return `) + append badges to template ---
    cn_def_end = 'title="CN ${c.cn}">CN:${c.cn}</span>`:\'\';\n return '
    s = repl(s, cn_def_end, cn_def_end.replace("\n return ", "\n" + BADGE_VARS + " return "), "badge/vars")
    s = repl(s, "${susp}${loh}${cn}`;}", "${susp}${loh}${cn}${tpfp}${cnb}${seqcb}`;}", "badge/template")
    # --- 4. matches() clauses after cnf ---
    cnf = "const cnf=$('fc').value; if(cnf!=='all'&&c.cn!==cnf)return false;"
    s = repl(s, cnf, cnf + MATCH_CLAUSES, "matches/clauses")
    # --- 5. listener array ---
    s = repl(s, "'ts','fq'].forEach", "'ts','fq','ftp','fkcn','fseq'].forEach", "listener")
    # --- 6. modal detail ---
    modal_anchor = '<span>CN <b>${c.cn}</b></span><span>grade <b>${c.grade}</b></span></div>'
    s = repl(s, modal_anchor,
             '<span>CN <b>${c.cn}</b></span><span>grade <b>${c.grade}</b></span>' + MODAL_DETAIL + '</div>',
             "modal/detail")
    return s, miss


def main():
    enr = json.load(open(ENR, encoding="utf-8"))
    s = open(F, encoding="utf-8").read()
    out, miss = patch(s, enr)
    print(f"merged enriched into D; cases missing enrichment: {miss}")
    if "--inplace" in sys.argv:
        bak = F + ".premarkers_bak"
        if not os.path.exists(bak):
            shutil.copy2(F, bak)
        open(F, "w", encoding="utf-8").write(out)
        print(f"[ok] patched IN PLACE ({len(out)//1024} KB)  backup: {bak}")
    else:
        op = F + ".patched2"
        open(op, "w", encoding="utf-8").write(out)
        print(f"[ok] wrote {op} ({len(out)//1024} KB)")


if __name__ == "__main__":
    main()
