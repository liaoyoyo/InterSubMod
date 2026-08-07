#!/usr/bin/env python3
"""註釋 drop-in 資料夾（擴充能力）。

使用者把 BED / TSV 丟進 <out>/annotations/（或 --annot-dir 指定的夾），
generator 自動偵測、自動註冊成篩選維度與畫布圖層 —— **不用改任何程式碼**。
這沿用 ISM C++ `--loh-bed` 的既成樣板，只是把「一個固定旗標」放大成「一個資料夾」。

支援兩種格式（副檔名判斷，內容再驗）：

  區間檔  .bed / .bed.gz            chrom, start, end[, name]
          → 維度值 = 落在區間內 / 不在。第 4 欄相異值 ≤24 時展開成子值。
          例：CNV 區段、LOH 區段、脆弱位點

  位點檔  .tsv / .csv / .txt(.gz)   需含 chrom + pos（欄名容忍多種寫法）
          → 維度值 = 該位點的某一欄。例：gene/drug 註釋、外部真值

檔名（去副檔名）就是維度名稱，所以 `cgc_genes.bed` 會出現一個叫 cgc_genes 的篩選維度。
"""
from __future__ import annotations

import bisect
import csv
import gzip
import os
from pathlib import Path

from capability import Capability, Registry, finalize, make, probe, probe_linkage

CHROM_ALIASES = ("chrom", "chr", "#chrom", "#chr", "chromosome", "Chr", "CHROM")
POS_ALIASES = ("pos", "position", "start", "Pos", "POS", "coordinate")
MAX_SUBVALUES = 24          # 第 4 欄相異值超過這個數就只留 overlap/no-overlap


def _open(p):
    return gzip.open(p, "rt", newline="") if str(p).endswith(".gz") else open(p, "r", newline="")


def _norm_chrom(c):
    c = str(c).strip()
    return c if c.startswith("chr") else ("chr" + c if c and c not in (".", "-") else "")


def _load_bed(path):
    """回傳 {chrom: (starts[], ends[], names[])}，starts 已排序。"""
    per = {}
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                continue
            c = _norm_chrom(f[0])
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3].strip() if len(f) > 3 else ""
            per.setdefault(c, []).append((s, e, name))
    out = {}
    for c, ivs in per.items():
        ivs.sort()
        out[c] = ([i[0] for i in ivs], [i[1] for i in ivs], [i[2] for i in ivs])
    return out


# 區段型 TSV（SAVANA CNA / 一般 CN segment）的欄名。偵測到就當區間檔處理，
# 使用者不必先跑 scripts/analysis/savana_to_smcnbed.py 轉成 BED。
SEG_CHROM = ("chromosome", "chrom", "chr", "#chrom", "seqnames")
SEG_START = ("start", "startPos", "begin")
SEG_END = ("end", "endPos", "stop")
SEG_CN = ("copyNumber", "copy_number", "cn", "CN", "total_cn")
SEG_MINOR = ("minorAlleleCopyNumber", "minor_cn", "minorCN", "minor_allele_cn")

# 與 savana_to_smcnbed.py 的預設一致，讓兩條路徑得到同樣的分類
GAIN_CN, LOSS_CN, LOH_MINOR = 3.0, 1.0, 0.5


def _load_segments(path):
    """SAVANA / CN segment TSV → {chrom: (starts, ends, labels)}。

    分類與 scripts/analysis/savana_to_smcnbed.py 相同：
      copyNumber >= 3          → gain
      copyNumber <= 1          → loss
      minorAlleleCopyNumber < 0.5 → loh（可與 gain/loss 併存，取較特異者）
      其餘                      → neutral
    """
    with _open(path) as fh:
        head = fh.readline().rstrip("\r\n")
    delim = "\t" if head.count("\t") >= head.count(",") else ","
    cols = [c.strip() for c in head.split(delim)]
    ck = next((c for c in cols if c in SEG_CHROM), None)
    sk = next((c for c in cols if c in SEG_START), None)
    ek = next((c for c in cols if c in SEG_END), None)
    cn = next((c for c in cols if c in SEG_CN), None)
    mn = next((c for c in cols if c in SEG_MINOR), None)
    if not (ck and sk and ek and cn):
        return None, cols

    per = {}
    with _open(path) as fh:
        for r in csv.DictReader(fh, delimiter=delim):
            c = _norm_chrom(r.get(ck))
            try:
                st, en = int(float(r[sk])), int(float(r[ek]))
                v = float(r[cn])
            except (TypeError, ValueError, KeyError):
                continue
            try:
                m = float(r[mn]) if mn and r.get(mn) not in (None, "", "NA", "nan") else None
            except (TypeError, ValueError):
                m = None
            if v >= GAIN_CN:
                lab = "gain"
            elif v <= LOSS_CN:
                lab = "loss"
            else:
                lab = "neutral"
            if m is not None and m < LOH_MINOR:
                lab = "loh" if lab == "neutral" else lab + "+loh"
            if c:
                per.setdefault(c, []).append((st, en, lab))
    out = {}
    for c, ivs in per.items():
        ivs.sort()
        out[c] = ([i[0] for i in ivs], [i[1] for i in ivs], [i[2] for i in ivs])
    return out, cols


def _load_points(path):
    """回傳 ({(chrom,pos): {col: val}}, 欄名清單)。"""
    with _open(path) as fh:
        sample = fh.read(8192)
    delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    rows, cols = {}, []
    with _open(path) as fh:
        rd = csv.DictReader(fh, delimiter=delim)
        cols = [c.strip() for c in (rd.fieldnames or [])]
        ck = next((c for c in cols if c in CHROM_ALIASES), None)
        pk = next((c for c in cols if c in POS_ALIASES), None)
        if not ck or not pk:
            return None, cols
        for r in rd:
            c = _norm_chrom(r.get(ck))
            try:
                p = int(float(r.get(pk)))
            except (TypeError, ValueError):
                continue
            if c:
                rows[(c, p)] = {k: (v or "").strip() for k, v in r.items()
                                if k not in (ck, pk) and v is not None}
    return rows, cols


def _overlaps(idx, chrom, pos):
    """回傳命中的 name（或 "" 表示命中但無名），沒命中回 None。"""
    t = idx.get(chrom)
    if not t:
        return None
    starts, ends, names = t
    i = bisect.bisect_right(starts, pos) - 1
    # BED 半開區間 [start, end)；往前多看幾個以處理巢狀區間
    while i >= 0 and i >= bisect.bisect_right(starts, pos) - 12:
        if starts[i] <= pos < ends[i]:
            return names[i]
        i -= 1
    return None


def load(reg: Registry, annot_dir, l1) -> Capability:
    cap = make(reg, "annotations", "註釋 drop-in 資料夾",
               enables=["filter:annot", "layer:annotBand"])
    d = Path(annot_dir) if annot_dir else None
    if not d or not d.is_dir():
        cap.state = "ABSENT"
        probe(cap, "S0", False,
              f"註釋資料夾不存在：{d}　"
              f"（建好後把 .bed / .tsv 丟進去即可，不用改程式）")
        cap.payload = {"dims": [], "dir": str(d) if d else ""}
        return cap
    probe(cap, "S0", True, str(d))

    files = sorted(p for p in d.iterdir()
                   if p.is_file() and not p.name.startswith(".")
                   and p.name.lower().endswith((".bed", ".bed.gz", ".tsv", ".tsv.gz",
                                                ".csv", ".csv.gz", ".txt", ".txt.gz")))
    if not files:
        cap.state = "ABSENT"
        probe(cap, "S1", False, f"{d} 裡沒有可辨識的檔案（.bed/.tsv/.csv/.txt，可 .gz）")
        cap.payload = {"dims": [], "dir": str(d)}
        return cap

    # 還原 L1 的絕對座標
    chroms = l1["chroms"]
    abspos, acc, prev = [], 0, -1
    for i in range(l1["n"]):
        if l1["chrom"][i] != prev:
            acc = l1["dpos"][i]
            prev = l1["chrom"][i]
        else:
            acc += l1["dpos"][i]
        abspos.append(acc)

    dims, notes = [], []
    for f in files:
        name = f.name
        for suf in (".gz",):
            if name.endswith(suf):
                name = name[: -len(suf)]
        stem = os.path.splitext(name)[0]
        dim_id = "annot_" + "".join(ch if ch.isalnum() else "_" for ch in stem)[:40]

        try:
            if name.lower().endswith(".bed"):
                idx = _load_bed(f)
                if not idx:
                    notes.append(f"{f.name}：解析不到區間，略過")
                    continue
                vals, names_seen = [], set()
                for i in range(l1["n"]):
                    nm = _overlaps(idx, chroms[l1["chrom"][i]], abspos[i])
                    if nm is None:
                        vals.append("no")
                    else:
                        vals.append(nm or "yes")
                        if nm:
                            names_seen.add(nm)
                hit = sum(1 for v in vals if v != "no")
                if len(names_seen) > MAX_SUBVALUES:
                    vals = ["no" if v == "no" else "yes" for v in vals]
                    names_seen = set()
                keys = [{"v": "no", "label": "不在區間內", "color": "#e7e5da"}]
                if names_seen:
                    pal = ["#176b58", "#285f8f", "#b66e20", "#a94336", "#6b5592",
                           "#2f8f76", "#5f8f6b", "#7d4a12"]
                    keys += [{"v": n, "label": n, "color": pal[i % len(pal)]}
                             for i, n in enumerate(sorted(names_seen))]
                else:
                    keys.append({"v": "yes", "label": "落在區間內", "color": "#176b58"})
                dims.append({"id": dim_id, "title": stem, "src": "derived",
                             "srcLabel": "drop-in", "keys": keys, "values": vals,
                             "kind": "bed", "file": f.name, "hit": hit})
                notes.append(f"{f.name}：BED，{sum(len(v[0]) for v in idx.values()):,} 個區間，"
                             f"命中 {hit:,}/{l1['n']:,} 個 sSNV")
            else:
                # 先試區段格式（SAVANA CNA / CN segment），再退回位點表
                seg, segcols = _load_segments(f)
                if seg:
                    vals, labs = [], set()
                    for i in range(l1["n"]):
                        nm = _overlaps(seg, chroms[l1["chrom"][i]], abspos[i])
                        if nm is None:
                            vals.append("no")
                        else:
                            vals.append(nm)
                            labs.add(nm)
                    hit = sum(1 for v in vals if v != "no")
                    palette = {"gain": "#a94336", "loss": "#285f8f", "loh": "#b66e20",
                               "gain+loh": "#7d2d24", "loss+loh": "#1c4160",
                               "neutral": "#c9c8bd"}
                    keys = [{"v": "no", "label": "無區段覆蓋", "color": "#e7e5da"}]
                    keys += [{"v": n, "label": n, "color": palette.get(n, "#6b5592")}
                             for n in sorted(labs)]
                    dims.append({"id": dim_id, "title": stem, "src": "derived",
                                 "srcLabel": "drop-in · CN segment", "keys": keys,
                                 "values": vals, "kind": "segment", "file": f.name,
                                 "hit": hit})
                    notes.append(
                        f"{f.name}：CN segment（SAVANA 式），"
                        f"{sum(len(v[0]) for v in seg.values()):,} 個區段，"
                        f"命中 {hit:,}/{l1['n']:,} 個 sSNV；"
                        f"分類 {'、'.join(sorted(labs))}"
                        f"（門檻 gain≥{GAIN_CN} / loss≤{LOSS_CN} / loh minor<{LOH_MINOR}，"
                        f"與 savana_to_smcnbed.py 一致）")
                    continue

                rows, cols = _load_points(f)
                if rows is None:
                    notes.append(f"{f.name}：找不到 chrom/pos 欄（欄名 {', '.join(cols[:6])}…），略過")
                    continue
                # 挑一個相異值適中的欄當維度值
                best, best_vals = None, None
                for c in (cols or []):
                    vs = {r.get(c, "") for r in rows.values()}
                    vs.discard("")
                    if 1 < len(vs) <= MAX_SUBVALUES:
                        best, best_vals = c, vs
                        break
                vals = []
                for i in range(l1["n"]):
                    r = rows.get((chroms[l1["chrom"][i]], abspos[i]))
                    if r is None:
                        vals.append("no")
                    else:
                        vals.append((r.get(best, "") or "yes") if best else "yes")
                hit = sum(1 for v in vals if v != "no")
                keys = [{"v": "no", "label": "無此註釋", "color": "#e7e5da"}]
                pal = ["#176b58", "#285f8f", "#b66e20", "#a94336", "#6b5592"]
                for i, v in enumerate(sorted(best_vals) if best_vals else ["yes"]):
                    keys.append({"v": v, "label": v if best else "有此註釋",
                                 "color": pal[i % len(pal)]})
                dims.append({"id": dim_id, "title": stem, "src": "derived",
                             "srcLabel": "drop-in", "keys": keys, "values": vals,
                             "kind": "points", "file": f.name, "hit": hit,
                             "column": best})
                notes.append(f"{f.name}：位點表 {len(rows):,} 列"
                             + (f"，取欄 {best}" if best else "")
                             + f"，命中 {hit:,}/{l1['n']:,} 個 sSNV")
        except Exception as exc:                       # noqa: BLE001
            notes.append(f"{f.name}：讀取失敗 {type(exc).__name__}: {exc}")

    probe(cap, "S1", True, f"{len(files)} 個檔案 → {len(dims)} 個維度")
    for n in notes:
        probe(cap, "S2", True, n)
    if dims:
        probe_linkage(cap, sum(d["hit"] for d in dims), l1["n"] * len(dims),
                      "各註釋命中 sSNV 數合計 / (sSNV × 註釋數)")
    if not dims:
        cap.state = "PARTIAL"
        cap.reason = "資料夾裡有檔案但都無法解析成維度（見上方逐檔說明）"
    cap.payload = {"dims": dims, "dir": str(d)}
    return finalize(cap)
