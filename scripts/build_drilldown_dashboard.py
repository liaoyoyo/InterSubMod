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
import platform
import shlex
import subprocess
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
import lca_ab as src_lca_ab                                    # noqa: E402
import lineage_paths as src_lineage_paths                        # noqa: E402
import strict_edges as src_edges                               # noqa: E402
import payload as emit_payload                                 # noqa: E402
import shell as emit_shell                                     # noqa: E402
import selfcheck as mod_selfcheck                              # noqa: E402
import cooccurrence as mod_cooccur                             # noqa: E402
import bake as mod_bake                                        # noqa: E402
from capability import EXIT_REFUSE, Refuse, Registry           # noqa: E402

# 站點路徑一律外部提供，不寫死在原始碼裡 —— 寫死等於別人 clone 下來跑不動。
# 三種來源，優先序由高到低：
#   1. CLI 旗標
#   2. 環境變數 DD_TOPOLOGY_ROOT / DD_ISM_ROOT / DD_LINEAGE_ROOT
#   3. 專案根的 drilldown.paths.json（gitignored，各站點自己維護）
#
# 這是 LongLineage 的 scripts/ci/check_repo_hygiene.sh 強制的慣例
# （該 CI 會 exit 2 阻擋寫死的絕對路徑）；InterSubMod 沒有那道 gate，
# 所以本檔自己遵守。
import os as _os

_CFG_NAME = "drilldown.paths.json"


def _load_site_config(start: Path) -> dict:
    """往上找 drilldown.paths.json。找不到回空 dict，不是錯誤。"""
    d = start.resolve()
    for _ in range(5):
        f = d / _CFG_NAME
        if f.is_file():
            try:
                return json.loads(f.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                print(f"  [warn] {f} 無法解析，忽略：{exc}", file=sys.stderr)
                return {}
        if d.parent == d:
            break
        d = d.parent
    return {}


_SITE = _load_site_config(HERE)


def site_path(key: str, env: str, default=None, sample: str = ""):
    """取站點路徑：CLI > 環境變數 > 設定檔 > default(None)。

    多樣本站點可用 ``{"HCC1395": "...", "COLO829": "..."}`` mapping，
    或在字串中使用 ``{sample}`` placeholder。舊的純字串設定仍相容。
    """
    raw = _os.environ.get(env) or _SITE.get(key) or default
    if isinstance(raw, dict):
        raw = raw.get(sample) if sample else None
    if isinstance(raw, str) and sample:
        raw = raw.replace("{sample}", sample)
    return raw


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


def _missing_hint(key: str, env: str, flag: str) -> str:
    return (f"找不到 {key}。三種提供方式擇一：\n"
            f"  1. CLI:  {flag} <PATH>\n"
            f"  2. 環境變數:  export {env}=<PATH>\n"
            f"  3. 設定檔:  在專案根建 {_CFG_NAME}，內容如\n"
            f'       {{"topology_root": "...", "ism_root": "...", "lineage_root": "..."}}')


def build_registry(args) -> Registry:
    reg = Registry()
    topo = src_topology.load(reg, args.topology, args.topology_receipt,
                             expected_sample=args.sample)
    if topo.usable:
        l1 = topo.payload["l1"]
        rids = {r["id"] for r in topo.payload["regions"]}
        src_mlhp.load(reg, args.mlhp, rids, expected_sample=args.sample)
        src_ism.load(reg, args.ism_root, l1, l1["chroms"],
                     expected_sample=args.sample)
        src_annot.load(reg, args.annot_dir, l1)
        src_lca_ab.load(reg, args.lca_pre, args.lca_post,
                        expected_sample=args.sample)
        src_lineage_paths.load(reg, args.paths_dir, rids,
                               expected_sample=args.sample)
        src_edges.load(reg, args.chrom_root, topo.payload["regions"], reg.get("mlhp"),
                       expected_sample=args.sample)
    return reg


def _generator_provenance() -> dict:
    """產生器身分證據。dirty 只檢查本產生器相關原始碼，避免被 repo 內
    無關的研究產物干擾。"""
    repo = HERE.parent

    def _run(argv):
        try:
            return subprocess.run(argv, cwd=repo, capture_output=True, text=True,
                                  timeout=15, check=False)
        except (OSError, subprocess.TimeoutExpired):
            return None

    rev = _run(["git", "rev-parse", "HEAD"])
    source_paths = ["scripts/build_drilldown_dashboard.py", "scripts/drilldown"]
    dirty_work = _run(["git", "diff", "--quiet", "--", *source_paths])
    dirty_index = _run(["git", "diff", "--cached", "--quiet", "--", *source_paths])
    dirty = None
    if dirty_work is not None and dirty_index is not None:
        dirty = dirty_work.returncode != 0 or dirty_index.returncode != 0
    return {
        "git_commit": rev.stdout.strip() if rev and rev.returncode == 0 else "",
        "dirty": dirty,
        "script_path": str(Path(__file__).resolve()),
        "script_sha256": cap_mod.sha256_file(Path(__file__).resolve()),
        "python": platform.python_version(),
    }


def _output_inventory(out: Path) -> dict:
    """輸出類型 count/bytes；不做全檔 hash，避免對數萬張 PNG/JS 的巨大 IO。

    receipt.json 本身排除，否則「receipt 記錄自己的 bytes」會形成無限遞迴。
    """
    kinds = {}
    total_count = total_bytes = 0
    for path in out.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(out).as_posix()
        if rel == "receipt.json":
            continue
        if rel == "index.html":
            kind = "html"
        elif rel == "SELFCHECK.md":
            kind = "selfcheck"
        elif rel.startswith("panels/") and rel.endswith(".T.png"):
            kind = "panel_tumor_png"
        elif rel.startswith("panels/") and rel.endswith(".png"):
            kind = "panel_base_png"
        elif rel.startswith("igv/") and rel.endswith(".js"):
            kind = "igv_js"
        elif rel.startswith("data/") and rel.endswith(".js"):
            kind = "data_js"
        elif rel.startswith("annotations/"):
            kind = "annotation"
        else:
            kind = "other"
        try:
            size = path.stat().st_size
        except OSError:
            continue
        ent = kinds.setdefault(kind, {"count": 0, "bytes": 0})
        ent["count"] += 1
        ent["bytes"] += size
        total_count += 1
        total_bytes += size
    return {"excludes": ["receipt.json"], "total_count": total_count,
            "total_bytes": total_bytes, "by_type": kinds}


def _validation_record(summary: dict, generator: dict, reg: Registry) -> dict:
    """Derive validation gates; never infer scientific validity from data presence."""
    malformed = sorted(cid for cid in reg.order if reg.caps[cid].state == cap_mod.MALFORMED)
    unhashed = sorted(
        ref.path for ref in reg.all_file_refs()
        if ref.exists and ref.size_bytes > 0 and not ref.sha256
    )
    gates = {
        "selfcheck_complete": summary.get("fail", 0) == 0 and summary.get("skip", 0) == 0,
        "sample_identity_and_schema": not malformed,
        "generator_clean": generator.get("dirty") is False,
        "input_file_hashes_complete": not unhashed,
        # This dashboard does not consume SEQC2 truth calls / HC BED / som.py output.
        "truth_set_validation": False,
    }
    blockers = []
    if malformed:
        blockers.append("MALFORMED_CAPABILITIES:" + ",".join(malformed))
    if summary.get("fail", 0):
        blockers.append(f"SELF_CHECK_FAIL:{summary['fail']}")
    if summary.get("skip", 0):
        blockers.append(f"SELF_CHECK_SKIP:{summary['skip']}")
    if generator.get("dirty") is not False:
        blockers.append("GENERATOR_NOT_CLEAN")
    if unhashed:
        blockers.append("UNHASHED_INPUT_FILES:" + ",".join(unhashed))
    blockers.append("NO_TRUTH_SET_BENCHMARK")

    if malformed:
        status = "INPUT_IDENTITY_OR_SCHEMA_FAILED"
    elif summary.get("fail", 0):
        status = "SELF_CHECK_FAILED"
    elif generator.get("dirty") is not False:
        status = "SOURCE_STATE_BLOCKED"
    elif unhashed:
        status = "PROVENANCE_INCOMPLETE"
    elif summary.get("skip", 0):
        status = "SELF_CHECK_INCOMPLETE"
    else:
        status = "SELF_CHECK_PASS_INTERNAL"
    return {
        "status": status,
        "gates": gates,
        "malformed_capabilities": malformed,
        "unhashed_input_files": unhashed,
        "blockers": blockers,
        "publishable_internal_qa": all((
            gates["selfcheck_complete"], gates["sample_identity_and_schema"],
            gates["generator_clean"], gates["input_file_hashes_complete"],
        )),
        "scientifically_validated": False,
    }


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
    ap.add_argument("--ism-root",
                    help="ISM run 根目錄（擴充能力；缺則甲基面板降級）")
    ap.add_argument("--chrom-root",
                    help="partition 的 chromosomes/ 目錄（strict endpoint edges；預設由 --sample 推導）")
    ap.add_argument("--lca-pre",
                    help="pre-LCA receipt 目錄（量化 LCA 增益）")
    ap.add_argument("--lca-post",
                    help="post-LCA receipt 目錄")
    ap.add_argument("--annot-dir",
                    help="註釋 drop-in 資料夾：把 .bed / .tsv 丟進去就會自動變成篩選維度，"
                         "不用改程式。預設 <out>/annotations/")
    ap.add_argument("--bake-panels", default="0",
                    help="產甲基雙面板 PNG：0（不產）/ N（前 N 個）/ all（全部）")
    ap.add_argument("--igv-mode", choices=["build", "skip"], default="build",
                    help="build（預設，產逐 region IGV JS）/ skip（僅供輕量 QA；receipt/UI 明示未產）")
    ap.add_argument("--lineage-assign", help="read_lineage_assignments.tsv.gz（側欄 lineage 軌）")
    ap.add_argument("--lineage-paths", help="unit_lineage_paths.tsv.gz")
    ap.add_argument("--paths-dir", help="放 *.unit_lineage_paths.tsv.gz 的目錄（自檢 C8 用）")
    ap.add_argument("--probe-only", action="store_true", help="只探測並印能力矩陣，不產頁")
    ap.add_argument("--figs-mode", choices=["copy"], default="copy",
                    help="copy（唯一已實作模式；輸出夾自足可搬走）")
    args = ap.parse_args()

    lin_root = site_path("lineage_root", "DD_LINEAGE_ROOT", sample=args.sample)
    if not args.ism_root:
        args.ism_root = site_path("ism_root", "DD_ISM_ROOT", sample=args.sample)
    if not args.lca_pre and lin_root:
        args.lca_pre = f"{lin_root}/bam_pre_lca_receipts"
    if not args.lca_post and lin_root:
        args.lca_post = f"{lin_root}/bam"
    if not args.lineage_assign and lin_root:
        args.lineage_assign = (f"{lin_root}/assign/{args.sample}"
                               ".all.read_lineage_assignments.tsv.gz")
    if not args.lineage_paths and lin_root:
        args.lineage_paths = f"{lin_root}/paths/{args.sample}.unit_lineage_paths.tsv.gz"
    if not getattr(args, "paths_dir", None) and lin_root:
        args.paths_dir = f"{lin_root}/paths"

    if not args.topology:
        root = site_path("topology_root", "DD_TOPOLOGY_ROOT", sample=args.sample)
        if not root:
            ap.error(_missing_hint("topology_root", "DD_TOPOLOGY_ROOT", "--topology"))
        args.topology = f"{root}/{args.sample}/{args.sample}.topology.jsonl"
    if not args.topology_receipt:
        args.topology_receipt = args.topology.replace(".jsonl", ".receipt.json")
    if not args.mlhp:
        root = site_path("topology_root", "DD_TOPOLOGY_ROOT", sample=args.sample)
        args.mlhp = f"{root}/{args.sample}/{args.sample}.exact_ps_mlhp.json" if root else ""
    if not args.chrom_root:
        root = site_path("topology_root", "DD_TOPOLOGY_ROOT", sample=args.sample)
        args.chrom_root = f"{root}/{args.sample}/chromosomes" if root else ""
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
    boot["assetPolicy"] = {"igv": args.igv_mode, "methylPanels": str(args.bake_panels)}
    boot["observationScope"] = emit_payload.observation_scope(topo.payload["l1"], ismc)
    sec = reg.get("strict_edges")
    boot["edges"] = ({
        "thresholds": sec.payload["thresholds"],
        "nEdge": sec.payload["n_edge"], "nPass": sec.payload["n_pass"],
        "brokenActive": sec.payload.get("broken_active", 0),
        "checked": sec.payload["checked"],
    } if (sec and sec.usable and sec.payload) else None)
    try:
        import composite as _comp
        boot["palette"] = _comp.palette_export()
    except Exception as exc:                    # noqa: BLE001
        print(f"  [warn] 色盤匯出失敗，圖例將降級：{type(exc).__name__}: {exc}")
        boot["palette"] = None
    boot["ism"] = ({
        "axes": ismc.payload["axes"],
        "missingAxes": ismc.payload["missing_axes"],
        "windowBp": ismc.payload["window_bp"],
        "hitDir": sorted(ismc.payload["hit_dir"]),
        "alpha": boot["observationScope"]["alpha"],
        "multiplicity": boot["observationScope"]["multiplicity"],
        "axisDefinitions": boot["observationScope"]["axis_definitions"],
    } if (ismc and ismc.usable and ismc.payload) else None)
    chroms = topo.payload["l1"]["chroms"]

    # IGV 式 read 對齊圖（平行產；逐 region 讀多個 ISM 窗，單執行緒太慢）
    igv_map = {}
    if ismc and ismc.usable and args.igv_mode == "build":
        import igv as mod_igv
        igv_map = mod_igv.build_many(ismc.payload["root"], topo.payload["regions"])

    def igv_fn(chrom, row):
        return igv_map.get(row["id"])

    shards = emit_payload.write_shards(out, topo, chroms,
                                       reg.get("ism_dirs"), reg.get("mlhp"), igv_fn)

    # IGV 逐 region 落檔（點擊才載）
    if igv_map:
        import json as _json
        igv_dir = out / "igv"
        n_b = 0
        for rid, v in igv_map.items():
            safe = "".join(ch if (ch.isalnum() or ch in "._-") else "_" for ch in rid)
            f = igv_dir / (safe[:120] + ".js")
            f.parent.mkdir(parents=True, exist_ok=True)
            body = ("window.__DD=window.__DD||{};window.__DD.IGV=window.__DD.IGV||{};"
                    f"window.__DD.IGV[{_json.dumps(rid)}]="
                    f"{emit_payload.json_for_html(v)};")
            f.write_text(body, encoding="utf-8")
            n_b += len(body.encode())
        print(f"IGV 圖    {len(igv_map):,} 個檔　{n_b / 2**20:.0f} MB（點擊才載入）")

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

    input_refs = reg.all_file_refs()
    # panel baking 另外直接讀 lineage assignment/path，它們不一定會經過
    # capability loader，因此只在實際 bake 時補進 receipt。
    if panel_info:
        known = {r.path for r in input_refs}
        for raw in (args.lineage_assign, args.lineage_paths):
            if raw and str(raw) not in known:
                input_refs.append(cap_mod.file_ref(raw))
                known.add(str(raw))
    inputs = [r.as_dict() for r in input_refs]
    (out / "SELFCHECK.md").write_text(
        mod_selfcheck.to_markdown(analysis["selfcheck"], args.sample, inputs), encoding="utf-8")

    output_inventory = _output_inventory(out)
    generator = _generator_provenance()
    validation = _validation_record(analysis["selfcheck"]["summary"], generator, reg)
    receipt = {
        "schema_name": "intersubmod.drilldown_dashboard_receipt",
        "schema_version": "2.0.0",
        "sample": args.sample,
        "panels": ({k: v for k, v in panel_info.items() if k != "panels"}
                   if panel_info else None),
        "built_at": meta["built_at"],
        "generator": generator,
        "command": shlex.join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]),
        "claim_ceiling": "OBSERVATION_ONLY_NO_TRUTH_SET",
        "validation_status": validation["status"],
        "validation": validation,
        "capabilities": reg.matrix(),
        "selfcheck": analysis["selfcheck"]["summary"],
        "inputs": inputs,
        "figs_mode": args.figs_mode,
        "asset_policy": boot["assetPolicy"],
        "shards": shards,
        "index_bytes": size,
        "output_inventory": output_inventory,
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
