#!/usr/bin/env python3
"""多層下鑽觀察儀表板 — CLI。

從全基因組 sSNV 一路下鑽到單一位點，每層都能開關圖層、勾選篩選、看維度共出現，
並在最細的一層把「該區域的演化分支」與「甲基的分群狀況」擺在一起。

硬核心只有 topology.jsonl；其餘資料層缺了就把對應面板降級並在頁面上寫明原因，
不會整份拒繪 —— 這與 build_exact_ps_layered_workstation.py 的哲學刻意相反。

用法:
    build_drilldown_dashboard.py --sample HCC1395 --probe-only
    build_drilldown_dashboard.py --sample HCC1395 --out <DIR> [--figs-mode copy]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "drilldown"))
sys.path.insert(0, str(HERE / "drilldown" / "sources"))
sys.path.insert(0, str(HERE / "drilldown" / "emit"))
sys.path.insert(0, str(HERE / "drilldown" / "panels"))

import capability as cap_mod                                   # noqa: E402
import topology as src_topology                                # noqa: E402
import ism as src_ism                                          # noqa: E402
import mlhp as src_mlhp                                        # noqa: E402
import annotation as src_annot                                 # noqa: E402
import payload as emit_payload                                 # noqa: E402
import shell as emit_shell                                     # noqa: E402
import selfcheck as mod_selfcheck                              # noqa: E402
import cooccurrence as mod_cooccur                             # noqa: E402
import bake as mod_bake                                        # noqa: E402
from capability import EXIT_REFUSE, Refuse, Registry           # noqa: E402

# 站點預設路徑。全部可用 CLI 覆寫；找不到就是對應能力 ABSENT，不是致命錯誤。
DEF_TOPOLOGY_ROOT = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
                     "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples")


DEF_ISM_ROOT = "/bip7_disk/liaoyoyo2001/ism_lineage_v1"


def _loci_with_ism(l1, ism_cap):
    """有 ISM 產物目錄的 sSNV，依 (chrom, pos) 排序。"""
    hit = ism_cap.payload["hit_dir"]
    chroms = l1["chroms"]
    out, acc, prev = [], 0, -1
    for i in range(l1["n"]):
        if l1["chrom"][i] != prev:
            acc = l1["dpos"][i]
            prev = l1["chrom"][i]
        else:
            acc += l1["dpos"][i]
        if i in hit:
            out.append((chroms[l1["chrom"][i]], acc))
    return out


def build_registry(args) -> Registry:
    reg = Registry()
    topo = src_topology.load(reg, args.topology, args.topology_receipt)
    if topo.usable:
        l1 = topo.payload["l1"]
        rids = {r["id"] for r in topo.payload["regions"]}
        src_mlhp.load(reg, args.mlhp, rids)
        src_ism.load(reg, args.ism_root, l1, l1["chroms"])
        src_annot.load(reg, args.annot_dir, l1)
    return reg


def print_matrix(reg: Registry) -> None:
    print(f"\n{'能力':<26}{'層級':<10}{'狀態':<10}{'停在':<6}{'連結率':>18}")
    print("-" * 74)
    for row in reg.matrix():
        link = row["linkage"]
        link_s = (f"{link['numerator']:,}/{link['denominator']:,} "
                  f"({link['rate'] * 100:.1f}%)") if link else "—"
        print(f"{row['title'][:24]:<26}{row['tier']:<10}{row['state_label']:<10}"
              f"{row['failed_stage'] or '—':<6}{link_s:>18}")
        for p in row["probes"]:
            print(f"    {'✓' if p['ok'] else '✗'} {p['stage']}  {p['detail']}")
        if row["reason"] and not all(p["ok"] for p in row["probes"]):
            print(f"    → {row['reason']}")
    print("-" * 74)
    print(reg.summary_line())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--topology", help="topology.jsonl（硬核心；預設由 --sample 推導）")
    ap.add_argument("--topology-receipt")
    ap.add_argument("--mlhp", help="MLHP json（read pattern；預設由 --sample 推導）")
    ap.add_argument("--out", help="輸出目錄；--probe-only 時可省略")
    ap.add_argument("--ism-root", default=DEF_ISM_ROOT,
                    help="ISM run 根目錄（擴充能力；缺則甲基面板降級）")
    ap.add_argument("--annot-dir",
                    help="註釋 drop-in 資料夾：把 .bed / .tsv 丟進去就會自動變成篩選維度，"
                         "不用改程式。預設 <out>/annotations/")
    ap.add_argument("--bake-panels", default="0",
                    help="產甲基雙面板 PNG：0（不產）/ N（前 N 個）/ all（全部）")
    ap.add_argument("--lineage-assign", help="read_lineage_assignments.tsv.gz（側欄 lineage 軌）")
    ap.add_argument("--lineage-paths", help="unit_lineage_paths.tsv.gz")
    ap.add_argument("--probe-only", action="store_true", help="只探測並印能力矩陣，不產頁")
    ap.add_argument("--figs-mode", choices=["copy", "link", "ref"], default="copy",
                    help="copy（預設，輸出夾自足可搬走）/ link（symlink，搬移會斷）/ ref（相對路徑）")
    args = ap.parse_args()

    if not args.topology:
        args.topology = f"{DEF_TOPOLOGY_ROOT}/{args.sample}/{args.sample}.topology.jsonl"
    if not args.topology_receipt:
        args.topology_receipt = args.topology.replace(".jsonl", ".receipt.json")
    if not args.mlhp:
        args.mlhp = f"{DEF_TOPOLOGY_ROOT}/{args.sample}/{args.sample}.exact_ps_mlhp.json"
    if not args.probe_only and not args.out:
        ap.error("需要 --out（或加 --probe-only）")
    if not args.annot_dir and args.out:
        args.annot_dir = str(Path(args.out) / "annotations")

    try:
        reg = build_registry(args)
        reg.require_core()
    except Refuse as exc:
        print(f"REFUSE: {exc}", file=sys.stderr)
        return EXIT_REFUSE

    print_matrix(reg)

    if args.probe_only:
        return 0

    import datetime
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    topo = reg.get("topology")

    ismc = reg.get("ism_dirs")
    emit_payload.fill_axis_codes(topo.payload["l1"], ismc)
    annotc = reg.get("annotations")
    annot_dims = (annotc.payload or {}).get("dims", []) if annotc else []
    dims = emit_payload.build_dims(topo, emit_payload.axis_dim(ismc) + [
        {k: v for k, v in d.items() if k != "values"} for d in annot_dims])
    # drop-in 維度的取值是逐 sSNV 的陣列，直接給前端查表
    boot_annot = {d["id"]: d["values"] for d in annot_dims}
    boot = emit_payload.build_boot(args.sample, reg, dims)
    boot["annotValues"] = boot_annot
    boot["ism"] = ({
        "axes": ismc.payload["axes"],
        "missingAxes": ismc.payload["missing_axes"],
        "windowBp": ismc.payload["window_bp"],
        "hitDir": sorted(ismc.payload["hit_dir"]),
    } if (ismc and ismc.usable and ismc.payload) else None)
    chroms = topo.payload["l1"]["chroms"]
    shards = emit_payload.write_shards(out, topo, chroms,
                                       reg.get("ism_dirs"), reg.get("mlhp"))

    # 甲基雙面板：ISM 算完後產圖，HTML 以相對路徑連結
    panel_info = None
    if args.bake_panels not in ("0", "", None) and ismc and ismc.usable:
        loci = _loci_with_ism(topo.payload["l1"], ismc)
        limit = 0 if args.bake_panels == "all" else int(args.bake_panels)
        lm = mod_bake.load_lineage_map(args.lineage_assign, args.lineage_paths)
        print(f"\n產甲基面板：{len(loci):,} 個候選位點"
              f"{'（限前 ' + str(limit) + ' 個）' if limit else '（全部）'}"
              f"　lineage map {len(lm):,} 條 read")
        panel_info = mod_bake.bake(ismc.payload["root"], out, loci, lm, cell_h=2, limit=limit)
        panel_info["shards"] = mod_bake.write_manifest(out, panel_info, chroms)
        print(f"　產出 {panel_info['made']:,} 張 / 略過 {panel_info['skipped']:,}"
              f"　共 {panel_info['bytes'] / 2**20:.1f} MB")
        boot["panels"] = {"made": panel_info["made"], "cellH": panel_info["cellH"]}

    meta = {
        "regions": len(topo.payload["regions"]),
        "cap_summary": reg.summary_line(),
        "built_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    }
    analysis = {
        "cooccur": mod_cooccur.build(reg),
        "sentinel": mod_cooccur.build_selection_sentinel(reg),
        "selfcheck": mod_selfcheck.run(reg),
    }
    size = emit_shell.write(out / "index.html", boot, reg.matrix(), meta, analysis)

    inputs = [r.as_dict() for r in reg.all_file_refs()]
    (out / "SELFCHECK.md").write_text(
        mod_selfcheck.to_markdown(analysis["selfcheck"], args.sample, inputs), encoding="utf-8")

    receipt = {
        "sample": args.sample,
        "panels": ({k: v for k, v in panel_info.items() if k != "panels"}
                   if panel_info else None),
        "built_at": meta["built_at"],
        "capabilities": reg.matrix(),
        "selfcheck": analysis["selfcheck"]["summary"],
        "inputs": inputs,
        "figs_mode": args.figs_mode,
        "shards": shards,
        "index_bytes": size,
    }
    (out / "receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")

    sc = analysis["selfcheck"]["summary"]
    print(f"\n自檢  {sc['pass']} 通過 / {sc['fail']} 不成立 / {sc['skip']} 無法檢查（共 {sc['total']}）")

    shard_total = sum(v["bytes"] for v in shards.values())
    shard_max = max(shards.items(), key=lambda kv: kv[1]["bytes"]) if shards else ("-", {"bytes": 0})
    print(f"\nindex.html   {size / 1024:>9,.0f} KB")
    print(f"分片 {len(shards):>2} 片   {shard_total / 1024:>9,.0f} KB"
          f"（最大 {shard_max[0]} {shard_max[1]['bytes'] / 1024:,.0f} KB）")
    print(f"receipt      {out / 'receipt.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
