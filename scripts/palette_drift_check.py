#!/usr/bin/env python3
"""[palette_drift_check] 唯讀稽核: ISM renderer 是否違反「import ism_heatmap_std、禁本地硬編色盤」規約。
偵測 (a) .py 本地定義 HP/ALT 語意色 dict 或 cluster PALETTE (b) .html 圖例嵌死 fill=#hex 當 HP/群。
canonical 真值 = import ism_heatmap_std(checker 本身複用 SoT,不硬編)。分級 WRONG(值≠canonical,正在錯色)/DUP(值=canonical 但本地副本=漂移風險)/CLUSTER。
唯讀,不改檔。用法: python3 palette_drift_check.py [DIR ...]  (預設掃 subcluster_pilot scripts + asset HTML)"""
import ast
import glob
import json
import os
import re
import sys

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts")  # canonical SoT 位置(renderer 皆由此 import)
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
try:
    import ism_heatmap_std as H
    CANON = {"HP": dict(H.HP_COL), "ALT": dict(H.ALT_COL), "TN": dict(H.TN_COL), "CLUSTER": list(H.CLUSTER_COL)}
except Exception as e:
    print(f"FATAL: cannot import ism_heatmap_std for canonical truth: {e}")
    sys.exit(1)

# 語意軸的 canonical hex(任何本地 PALETTE 用到這些=撞語意軸)
SEMANTIC_HEX = {H.HP_COL["1"], H.HP_COL["2"], H.HP_COL["1-1"], H.HP_COL["2-1"],
                H.ALT_COL["ALT"], H.ALT_COL["REF"], H.TN_COL["1"], H.TN_COL["0"]}
# 已知常見錯色(歷史 bug): HP1 用 ALT 紅 / HP2 用某藍 / 群色用紅綠
KNOWN_BAD = {"#dc2626", "#2563eb", "#16a34a", "#f59e0b", "#7c3aed"}
HEX = re.compile(r"^#[0-9a-fA-F]{6}$")

PILOT = os.path.join(os.path.dirname(__file__), "..", "docs", "methodology", "_assets", "20260618_subcluster_pilot")


def lit(node):
    """ast 常數值(py3.8+ ast.Constant);非常數回 None。"""
    return node.value if isinstance(node, ast.Constant) else None


def scan_py(path):
    """回 {imports_H, dicts:[(kind,line,actual,issues)], palettes:[(name,line,hexes,collide)], parse_ok}。"""
    src = open(path, encoding="utf-8", errors="replace").read()
    imports_H = "ism_heatmap_std" in src
    out = {"imports_H": imports_H, "dicts": [], "palettes": [], "parse_ok": True}
    try:
        tree = ast.parse(src)
    except Exception:
        out["parse_ok"] = False
        return out
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        v = node.value
        ln = node.lineno
        # dict literal
        if isinstance(v, ast.Dict):
            keys = [lit(k) for k in v.keys]
            vals = [lit(x) for x in v.values]
            kv = {k: x for k, x in zip(keys, vals) if isinstance(k, str)}
            ks = set(kv)
            if "1-1" in ks or "2-1" in ks:  # HP somatic 子標籤是 haplotag 專屬
                # 必須是「色」dict(值為 hex),排除 label-remap dict(如 LABMAP {"1-1":"1-1"})
                if any(isinstance(kv.get(k), str) and HEX.match(kv[k]) for k in ("1", "2", "1-1", "2-1")):
                    issues = []
                    for k in ("1", "2", "1-1", "2-1"):
                        if k in kv and isinstance(kv[k], str) and HEX.match(kv[k]) and kv[k].lower() != CANON["HP"].get(k, "").lower():
                            issues.append(f"HP{k}={kv[k]}≠{CANON['HP'].get(k)}")
                    out["dicts"].append(("HP", ln, kv, issues))
            elif {"ALT", "REF"} <= ks:
                if any(isinstance(kv.get(k), str) and HEX.match(kv[k]) for k in ("ALT", "REF")):
                    issues = []
                    for k in ("ALT", "REF"):
                        if isinstance(kv.get(k), str) and HEX.match(kv[k]) and kv[k].lower() != CANON["ALT"].get(k, "").lower():
                            issues.append(f"{k}={kv[k]}≠{CANON['ALT'].get(k)}")
                    out["dicts"].append(("ALT", ln, kv, issues))
        # list literal that looks like a cluster palette
        elif isinstance(v, ast.List) and node.targets and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            hexes = [lit(e) for e in v.elts]
            hexes = [h for h in hexes if isinstance(h, str) and HEX.match(h)]
            if len(hexes) >= 3 and re.search(r"PAL|COLOR|CMAP|PALETTE", name, re.I):
                collide = sorted(set(h for h in hexes if h in SEMANTIC_HEX or h in KNOWN_BAD))
                out["palettes"].append((name, ln, hexes, collide))
    return out


def scan_html(path):
    """圖例 fill=#hex 緊跟 HP/群/REF/ALT 標籤 → 比 canonical。"""
    src = open(path, encoding="utf-8", errors="replace").read()
    hits = []
    for m in re.finditer(r'fill="(#[0-9a-fA-F]{6})"[^>]*/>\s*<text[^>]*>([^<]{0,12})', src):
        hexv, label = m.group(1), m.group(2).strip()
        ln = src[:m.start()].count("\n") + 1
        exp = None
        if label in ("HP1",): exp = H.HP_COL["1"]
        elif label in ("HP2",): exp = H.HP_COL["2"]
        elif label in ("HP1-1",): exp = H.HP_COL["1-1"]
        elif label in ("HP2-1",): exp = H.HP_COL["2-1"]
        elif label in ("REF",): exp = H.ALT_COL["REF"]
        elif label in ("ALT",): exp = H.ALT_COL["ALT"]
        if exp and hexv.lower() != exp.lower():
            hits.append((ln, label, hexv, exp))
    return hits


def main():
    dirs = sys.argv[1:] or [os.path.join(PILOT, "scripts")]
    py_files = []
    for d in dirs:
        py_files += sorted(glob.glob(os.path.join(d, "*.py")))
    html_files = sorted(glob.glob(os.path.join(PILOT, "*.html")) + glob.glob(os.path.join(PILOT, "obs_ws", "**", "*.html"), recursive=True))

    viol, dup = [], []
    cluster_collide, na, clean, parse_fail = [], [], [], []
    for f in py_files:
        r = scan_py(f)
        b = os.path.basename(f)
        if not r["parse_ok"]:
            parse_fail.append(b); continue
        wrong = [(k, ln, iss) for (k, ln, _, iss) in r["dicts"] if iss]
        dups = [(k, ln) for (k, ln, _, iss) in r["dicts"] if not iss]
        coll = [(n, ln, c) for (n, ln, _, c) in r["palettes"] if c]
        if wrong:
            viol.append((b, r["imports_H"], wrong))
        elif dups:
            dup.append((b, r["imports_H"], dups))
        elif coll:
            cluster_collide.append((b, r["imports_H"], coll))
        elif r["imports_H"]:
            clean.append(b)
        else:
            na.append(b)

    html_viol = []
    for f in html_files:
        hits = scan_html(f)
        if hits:
            html_viol.append((os.path.relpath(f, PILOT), hits))

    # ---- report ----
    print("=" * 78)
    print("PALETTE-DRIFT AUDIT (read-only)  canonical = scripts/ism_heatmap_std.py")
    print(f"scanned: {len(py_files)} .py + {len(html_files)} .html")
    print("=" * 78)
    print(f"\n🔴 VIOLATIONS — 本地語意色 dict 值≠canonical (正在錯色) : {len(viol)}")
    for b, imp, wrong in viol:
        tag = " [WORST: import H 卻仍硬編]" if imp else ""
        for k, ln, iss in wrong:
            print(f"   {b}:{ln}  {k}dict  {'; '.join(iss)}{tag}")
    print(f"\n🟠 DUPLICATES — 本地語意色 dict 值=canonical 但仍本地副本(漂移風險) : {len(dup)}")
    for b, imp, dups in dup:
        tag = " [import H]" if imp else ""
        print(f"   {b}  " + ", ".join(f"{k}@L{ln}" for k, ln in dups) + tag)
    print(f"\n🟡 CLUSTER-PALETTE 撞語意軸色 : {len(cluster_collide)}")
    for b, imp, coll in cluster_collide:
        for n, ln, c in coll:
            print(f"   {b}:{ln}  {n} 用到語意色 {c}")
    print(f"\n🔴 HTML 圖例嵌死錯色 : {len(html_viol)}")
    for rel, hits in html_viol:
        for ln, label, hexv, exp in hits:
            print(f"   {rel}:{ln}  {label}={hexv}≠{exp}")
    print(f"\n🟢 CLEAN(import H, 無本地語意 dict): {len(clean)}   ⚪ N/A(非標準側欄 renderer): {len(na)}   parse-fail: {len(parse_fail)}")
    total_bad = len(viol) + len(dup) + len(cluster_collide) + len(html_viol)
    print(f"\nSUMMARY: 需處理 = {len(viol)} WRONG + {len(dup)} DUP + {len(cluster_collide)} CLUSTER + {len(html_viol)} HTML = {total_bad}")

    rep = {"viol": [[b, imp, w] for b, imp, w in viol], "dup": [[b, imp, d] for b, imp, d in dup],
           "cluster_collide": [[b, imp, c] for b, imp, c in cluster_collide],
           "html_viol": [[r, h] for r, h in html_viol], "clean_n": len(clean), "na_n": len(na),
           "parse_fail": parse_fail, "scanned_py": len(py_files), "scanned_html": len(html_files)}
    outp = os.path.join(os.path.dirname(__file__), "..", "docs", "methodology", "_assets", "20260618_subcluster_pilot", "palette_drift_audit.json")
    json.dump(rep, open(outp, "w"), ensure_ascii=False, indent=1)
    print(f"\nJSON: {os.path.relpath(outp)}")


if __name__ == "__main__":
    main()
