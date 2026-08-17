#!/usr/bin/env python3
"""把一個樣本的最終交付品收攏成單一目錄，並產生「讓人與 AI 都不會誤讀」的入口檔。

為什麼需要這支腳本
------------------
同一個樣本在磁碟上會同時存在多輪產物（v1 / v2 / v3 / 暫存 / 半途中斷的重跑）。
沒有標記時，下一個人 —— 或下一個 AI —— 沒有辦法分辨哪一份是結論、哪一份是過程
殘留。憑目錄名猜（「v3 應該比 v2 新吧」）在重跑中斷、命名不一致時就會猜錯，
而且錯得無聲無息。

這支腳本把「哪一份是交付品」寫成機械可驗證的事實：
  1. 交付目錄內有 DELIVERY.md，開頭就宣告信任度與 scope
  2. MANIFEST.json 逐檔記 sha256 + 角色 + 產生它的 commit
  3. verify.sh 可一鍵重驗，不必相信任何敘述
  4. 非交付路徑留 tombstone，明寫「這不是結果」

沿用既有三軸，不發明第四套（docs/data_specs/20260414_output資料信任度與生命週期_01.md）：
  信任度   CURRENT / PRE-FIX / SUPERSEDED / DEPRECATED
  生命週期 in_progress / validated / finalized
  刪除安全度 SAFE_DELETE / LOW_RISK / KEEP

用法
----
  build_final_delivery.py --sample HCC1395 --root /bip7_disk/.../HCC1395_final \\
      [--canonical-link <DIR>] [--tombstone <DIR> ...] [--dry-run]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

CHUNK = 16 * 1024 * 1024

# 每個角色一句話說明「這是什麼、可以拿來做什麼」。
# 交付目錄裡不該有任何檔案是「不知道為什麼在這」——沒登記的檔案會被列成 UNREGISTERED。
ROLES = {
    "bam": ("lineage_tagged_bam",
            "交付主體：每條 read 帶譜系標籤的座標排序 BAM，可直接進 IGV"),
    "bam_index": ("bam_index", "BAM 索引；缺它時 region query 會靜默回空"),
    "receipt": ("tag_bam_receipt", "逐染色體的產生紀錄（輸入路徑與各項計數）"),
    "assign": ("read_lineage_assignments", "重建配方：read × region 的 pattern 指派"),
    "paths": ("unit_lineage_paths", "重建配方：每個樹頂點的階層路徑"),
    "pre_lca": ("pre_lca_receipt", "LCA 增益 A/B 的對照組 receipt"),
    "vcf": ("somatic_vcf_input",
            "分析輸入：somatic sSNV 位點清單（全基因組 + 逐染色體分片）"),
    "ism": ("ism_methylation_run", "甲基分析輸出（逐位點目錄 + significance_summary）"),
    "dashboard": ("observation_dashboard", "互動觀察 HTML，自足可整包搬走"),
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def git_commit(repo: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=30)
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"],
                               capture_output=True, text=True, timeout=60)
        sha = out.stdout.strip() or "unknown"
        return sha + ("-dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def classify(rel: str) -> tuple[str, str]:
    """路徑 → (role_id, 說明)。未知者回 UNREGISTERED，讓它在 manifest 裡顯眼。"""
    p = rel.replace(os.sep, "/")
    if p.startswith("bam/"):
        if p.endswith(".bai") or p.endswith(".crai"):
            return ROLES["bam_index"]
        if p.endswith(".receipt.json"):
            return ROLES["receipt"]
        if p.endswith(".bam") or p.endswith(".cram"):
            return ROLES["bam"]
    if p.startswith("lineage/assign/"):
        return ROLES["assign"]
    if p.startswith("lineage/paths/"):
        return ROLES["paths"]
    if p.startswith("lineage/bam_pre_lca_receipts/"):
        return ROLES["pre_lca"]
    if p.startswith("inputs/vcf/"):
        return ROLES["vcf"]
    if p.startswith("ism/"):
        return ROLES["ism"]
    if p.startswith("dashboard/"):
        return ROLES["dashboard"]
    return ("UNREGISTERED", "未登記 —— 不知道這個檔為什麼在交付目錄裡")


def collect(root: Path, hash_limit_bytes: int) -> list[dict]:
    """走訪交付目錄。大檔仍算 sha256（正確性優先），但會回報進度。"""
    out = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in sorted(filenames):
            f = Path(dirpath) / fn
            if f.name in ("MANIFEST.json", "DELIVERY.md", "verify.sh", "REBUILD.md"):
                continue
            rel = str(f.relative_to(root))
            role, desc = classify(rel)
            size = f.stat().st_size
            rec = {"path": rel, "role": role, "description": desc, "bytes": size}
            # ISM/dashboard 是數萬個小檔，逐檔 sha256 會跑很久且對交接沒有增益，
            # 改記目錄層級的檔數與總量；BAM 這種「單一大檔就是交付主體」的一定要算。
            if size <= hash_limit_bytes or role in ("lineage_tagged_bam", "bam_index"):
                print(f"    sha256 {rel} ({size / 2**30:.1f} GiB)…", flush=True)
                rec["sha256"] = sha256_of(f)
            out.append(rec)
    return out


def summarize_dirs(entries: list[dict]) -> dict:
    agg: dict[str, dict] = {}
    for e in entries:
        top = e["path"].split("/")[0]
        a = agg.setdefault(top, {"files": 0, "bytes": 0})
        a["files"] += 1
        a["bytes"] += e["bytes"]
    return agg


DELIVERY_TMPL = """# {sample} · 最終交付品

> **這個目錄就是結論。** 同一個樣本在磁碟上曾有多輪產物（v1／v2／暫存／中斷的重跑），
> 那些都已刪除或標記為過程殘留。**要看結果、要引用數字、要進 IGV，都只用這裡的檔案。**

| | |
|---|---|
| 信任度 | **{trust}** — 依 `docs/data_specs/20260414_output資料信任度與生命週期_01.md §1` |
| 生命週期 | **{lifecycle}** |
| 產生時間 | {as_of} |
| 產生 commit | `{commit}` |
| 驗證方式 | `bash verify.sh`（重算 sha256 + read 數守恆 + symlink 有效性） |

## 目錄用途

| 路徑 | 是什麼 |
|---|---|
| `bam/` | **交付主體**。單一座標排序 BAM，每條 read 帶譜系標籤，可直接進 IGV |
| `inputs/` | **分析輸入**。VCF 位點清單放這裡而不是塞在某個 run 底下 —— 2026-08-13 曾因為 VCF 被放在 `ism_lineage_v1/vcf/`，清理舊 run 時連輸入一起刪掉 |
| `lineage/` | **重建配方**（約 1.5 GB）。有它就能從原始 BAM 重新產出 `bam/`，不必備份 267 GB |
| `ism/` | 甲基分析輸出（逐位點目錄 + `significance_summary.csv`） |
| `dashboard/` | 互動觀察 HTML。**整個資料夾一起搬**，`index.html` 需要同層的 `data/`／`panels/`／`igv/` |

## IGV 怎麼看

```
Color alignments by → tag → lg     譜系（HP2-1 與 HP2-1+ 同色，共 13 種）
Group alignments by → tag → ls     可信度（U 唯一／M 同分並列／P 部分覆蓋／A 拓撲未定）
```

只想看**已確認**的：用 `ls=U OR ln=1`（`ln` = 相容頂點數，等於 1 表示位置已定死）。
⚠ `HP` 是字串 `HP:Z`（九態 `1`/`2`/`1-1`/`2-1`/`3`），IGV 內建 haplotype 模式只認 `HP:i`，
所以要走 **Group by tag** 而非內建模式。

## 🔴 三件必須知道的事

1. **scope ≠ BAM 內容**。BAM 含 24 個 contig 共 {reads} 條（來源 BAM 的 100%），
   但**拓撲與甲基分析只做 chr1–22**。chrX／chrY 的 read 有 `HP`／`PS`、**無 `lv`** ——
   那是分析範圍，不是資料遺失。

2. **這裡不是備份**。`bam/` 只有一份實體。真正的保險是 `lineage/` 的重建配方
   ＋ canonical 未被動過的原始 BAM；照 `REBUILD.md` 重跑約 5 小時可得同一份。

3. **canonical 底下的是 symlink**。`.bam` 與 `.bai` 都必須連；只連 BAM 會讓 IGV
   **靜默失敗**（region query 回空且不報錯）。

## 給 AI 的說明

讀到這個檔案就停止猜測目錄版本號。`MANIFEST.json` 的 `artifacts[]` 是這個樣本
唯一的權威清單；任何不在其中的路徑都不是結論。若 `MANIFEST.json` 與本檔敘述不一致，
**以 `MANIFEST.json` 為準**（它是機械產生的，本檔是人寫的）。
"""

TOMBSTONE = """# ⚠ 這不是結果

本目錄是 **{sample}** 的**過程殘留**，信任度 **SUPERSEDED**
（`docs/data_specs/20260414_output資料信任度與生命週期_01.md §1`）。

**不要引用這裡的任何數字，也不要拿這裡的 BAM 進 IGV。**

最終交付品在：

    {final}

那裡的 `DELIVERY.md` 說明用途，`MANIFEST.json` 是權威清單，`verify.sh` 可一鍵重驗。

（標記時間 {as_of}）
"""


def upstream_ref(path: str | None, role: str, want_hash: bool) -> dict:
    """描述一個不在交付目錄內、但重建時必需的上游輸入。

    預設只記 size + mtime，不算 sha256 —— 來源 BAM 有 259 GiB，逐次雜湊要 ~10 分鐘。
    這是刻意的取捨：size+mtime 足以偵測「檔案被換掉」這個最常見的失效，
    要 byte 級保證時再加 --hash-upstream。REBUILD.md 會明寫用了哪一種。
    """
    if not path:
        return {"role": role, "path": None, "status": "NOT_RECORDED",
                "note": "產生 manifest 時未提供此路徑 —— 重建前必須自行確認"}
    f = Path(path)
    if not f.exists():
        return {"role": role, "path": str(f), "status": "MISSING",
                "note": "🔴 上游輸入已不存在 —— 這份交付品目前無法重建"}
    st = f.stat()
    rec = {"role": role, "path": str(f), "status": "PRESENT", "bytes": st.st_size,
           "mtime": datetime.fromtimestamp(st.st_mtime, timezone.utc)
                            .astimezone().isoformat(timespec="seconds")}
    if want_hash:
        print(f"    sha256 上游 {f.name} ({st.st_size / 2**30:.1f} GiB)…", flush=True)
        rec["sha256"] = sha256_of(f)
    return rec


REBUILD_TMPL = """# 重建 {sample} 的譜系標籤 BAM

> **為什麼記配方而不備份 267 GiB**：這份 BAM 是**衍生產物**，由下列輸入完全決定。
> 備份會漂移、會忘記是哪一版產的；配方不會。照本文重跑約 5 小時，可得到內容一致的同一份
> （`MANIFEST.json` 裡的 sha256 可逐檔比對）。

產生時間 `{as_of}` ／ InterSubMod commit `{commit}` ／ LongLineage commit `{ll_commit}`

## 輸入清單

{inputs_table}

{upstream_warn}
交付目錄內的 `lineage/assign/` 與 `lineage/paths/` 也是重建輸入，已隨交付品保存
（合計約 1.5 GB），它們的 sha256 在 `MANIFEST.json` 裡。

## 重建步驟

```bash
LL=/big7_disk/liaoyoyo2001/LongLineage
OUT=<新的輸出目錄>          # 不要直接覆蓋現有交付品
FINAL={root}

mkdir -p "$OUT"/{{bam,logs}}

# ① 逐染色體標記（24 條；每條約 1-6 分鐘，會自帶 .bai）
for c in chr{{1..22}} chrX chrY; do
  "$LL/build_lineage_migrate/bin/longlineage-tag-bam" \
    --in-bam    {src_bam} \
    --sidecar   {sidecar} \
    --assignments  "$FINAL/lineage/assign/{sample}.all.read_lineage_assignments.tsv.gz" \
    --lineage-paths "$FINAL/lineage/paths/{sample}.unit_lineage_paths.tsv.gz" \
    --topology  {topology} \
    --region "$c" --threads 16 \
    --out-bam "$OUT/bam/{sample}.$c.lineage_tagged.bam" \
    --receipt "$OUT/bam/{sample}.$c.tag_bam.receipt.json" \
    > "$OUT/logs/tag_bam.$c.log" 2>&1 || {{ echo "$c FAIL"; exit 1; }}
done

# ② 合併成單一 BAM（約 4 小時，I/O 為瓶頸）
parts=(); for c in chr{{1..22}} chrX chrY; do parts+=("$OUT/bam/{sample}.$c.lineage_tagged.bam"); done
samtools merge -@ 16 -f -o "$OUT/bam/{sample}.lineage_tagged.bam" "${{parts[@]}}"
samtools index -@ 16 "$OUT/bam/{sample}.lineage_tagged.bam"

# ③ 驗證守恆後才刪分檔（順序不可顛倒）
GOT=$(samtools idxstats "$OUT/bam/{sample}.lineage_tagged.bam" | awk '{{s+=$3+$4}} END{{print s}}')
[[ "$GOT" == "{reads}" ]] || {{ echo "read 數不符：$GOT ≠ {reads}，分檔全部保留"; exit 1; }}
for f in "${{parts[@]}}"; do rm -f "$f" "$f.bai"; done
```

## 驗證重建結果

```bash
bash InterSubMod/scripts/handoff/verify_final_delivery.sh "$OUT"
```

⚠ **重跑不保證 byte-identical**。`samtools merge` 的 BGZF 區塊切分會受執行緒數影響，
所以請比對**內容層**的不變量（read 數、各 contig 計數、標籤分布），
而不是整檔 sha256。`MANIFEST.json` 的 `reads_total` 與 receipt 各項計數是可靠的比對基準。

## 這份 BAM 的 scope

含 24 個 contig（chr1–22 ＋ chrX ＋ chrY）共 {reads} 條 read＝來源 BAM 的 100%。
但**拓撲與甲基分析只做 chr1–22**（`topology.jsonl` 不含性染色體），
所以 chrX／chrY 的 read 有 `HP`／`PS` 但**無 `lv`**。這是分析範圍，不是資料遺失。
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--root", required=True, help="最終交付目錄")
    ap.add_argument("--repo", default="/big7_disk/liaoyoyo2001/InterSubMod")
    ap.add_argument("--trust", default="CURRENT",
                    choices=["CURRENT", "PRE-FIX", "SUPERSEDED", "DEPRECATED"])
    ap.add_argument("--lifecycle", default="finalized",
                    choices=["in_progress", "validated", "finalized"])
    ap.add_argument("--hash-limit-mb", type=float, default=64.0,
                    help="超過此大小的非 BAM 檔跳過 sha256（BAM 一律算）")
    ap.add_argument("--tombstone", action="append", default=[],
                    help="要標記為過程殘留的目錄（可重複）")
    ap.add_argument("--src-bam", help="重建配方：LongPhase-S 輸出的來源 BAM")
    ap.add_argument("--sidecar", help="重建配方：read_tags sidecar")
    ap.add_argument("--topology", help="重建配方：topology.jsonl")
    ap.add_argument("--ll-repo", default="/big7_disk/liaoyoyo2001/LongLineage")
    ap.add_argument("--hash-upstream", action="store_true",
                    help="連上游大檔也算 sha256（來源 BAM 259 GiB 約需 10 分鐘）")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    root = Path(a.root)
    if not root.is_dir():
        print(f"交付目錄不存在：{root}", file=sys.stderr)
        return 2

    as_of = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    commit = git_commit(Path(a.repo))

    print(f"掃描 {root} …")
    entries = collect(root, int(a.hash_limit_mb * 2**20))
    dirs = summarize_dirs(entries)
    unreg = [e for e in entries if e["role"] == "UNREGISTERED"]

    bam = next((e for e in entries if e["role"] == "lineage_tagged_bam"), None)
    reads = "—"
    if bam:
        try:
            o = subprocess.run(["samtools", "idxstats", str(root / bam["path"])],
                               capture_output=True, text=True, timeout=600).stdout
            reads = f"{sum(int(l.split(chr(9))[2]) + int(l.split(chr(9))[3]) for l in o.strip().split(chr(10)) if l):,}"
        except (OSError, subprocess.SubprocessError, ValueError, IndexError):
            pass

    # 上游輸入：不在交付目錄內，但少了就無法重建。狀態要顯眼 ——
    # MISSING 代表這份交付品已經失去可重建性，那是必須立刻知道的事。
    upstream = [
        upstream_ref(a.src_bam, "source_tagged_bam", a.hash_upstream),
        upstream_ref(a.sidecar, "read_tags_sidecar", a.hash_upstream),
        upstream_ref(a.topology, "topology_jsonl", a.hash_upstream),
    ]
    ll_commit = git_commit(Path(a.ll_repo))

    manifest = {
        "schema_name": "intersubmod.final_delivery_manifest",
        "schema_version": "1.0.0",
        "sample": a.sample,
        "as_of": as_of,
        "generated_by_commit": commit,
        "trust_tier": a.trust,
        "lifecycle": a.lifecycle,
        "trust_tier_spec": "docs/data_specs/20260414_output資料信任度與生命週期_01.md §1",
        "root": str(root),
        "reads_total": reads,
        "directory_summary": {k: {"files": v["files"], "bytes": v["bytes"],
                                  "gib": round(v["bytes"] / 2**30, 2)}
                              for k, v in sorted(dirs.items())},
        "unregistered_count": len(unreg),
        "rebuildable": all(u["status"] == "PRESENT" for u in upstream),
        "upstream_inputs": upstream,
        "longlineage_commit": ll_commit,
        "artifacts": entries,
    }

    if a.dry_run:
        print(json.dumps({k: v for k, v in manifest.items() if k != "artifacts"},
                         ensure_ascii=False, indent=2))
        print(f"\n（dry-run）artifacts {len(entries)} 筆、未登記 {len(unreg)} 筆")
        return 0

    (root / "MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (root / "DELIVERY.md").write_text(
        DELIVERY_TMPL.format(sample=a.sample, trust=a.trust, lifecycle=a.lifecycle,
                             as_of=as_of, commit=commit, reads=reads), encoding="utf-8")

    rows = []
    for u in upstream:
        mark = {"PRESENT": "✅", "MISSING": "🔴", "NOT_RECORDED": "⚠"}[u["status"]]
        size = f"{u['bytes'] / 2**30:.2f} GiB" if u.get("bytes") else "—"
        h = u.get("sha256", "")
        rows.append(f"| {mark} `{u['role']}` | `{u.get('path') or '（未記錄）'}` | {size} | "
                    f"{'`' + h[:16] + '…`' if h else 'size+mtime'} |")
    inputs_table = ("| 狀態／角色 | 路徑 | 大小 | 完整性依據 |\n|---|---|---|---|\n"
                    + "\n".join(rows))
    warn = ""
    if not manifest["rebuildable"]:
        warn = ("\n> 🔴 **這份交付品目前無法重建** —— 上表有輸入標記為 MISSING／未記錄。\n"
                "> 在補回那些輸入之前，`bam/` 的 267 GiB 是不可再生的孤本。\n")
    elif not a.hash_upstream:
        warn = ("\n> ⚠ 上游輸入以 **size + mtime** 記錄，未算 sha256（來源 BAM 259 GiB，"
                "雜湊約需 10 分鐘）。\n> 要 byte 級保證時重跑本工具並加 `--hash-upstream`。\n")

    (root / "REBUILD.md").write_text(REBUILD_TMPL.format(
        sample=a.sample, as_of=as_of, commit=commit, ll_commit=ll_commit,
        inputs_table=inputs_table, upstream_warn=warn, root=str(root),
        src_bam=a.src_bam or "<來源 BAM>", sidecar=a.sidecar or "<sidecar>",
        topology=a.topology or "<topology.jsonl>", reads=reads.replace(",", "")),
        encoding="utf-8")

    for t in a.tombstone:
        d = Path(t)
        if d.is_dir():
            (d / "README_NOT_THE_RESULT.md").write_text(
                TOMBSTONE.format(sample=a.sample, final=root, as_of=as_of),
                encoding="utf-8")
            print(f"  tombstone → {d}")

    print(f"\n✓ MANIFEST.json  {len(entries)} 筆 artifact")
    for k, v in sorted(dirs.items()):
        print(f"    {k:<12}{v['files']:>7,} 檔  {v['bytes'] / 2**30:>8.2f} GiB")
    if unreg:
        print(f"\n⚠ 有 {len(unreg)} 個未登記檔案（交付目錄不該有來歷不明的東西）：")
        for e in unreg[:10]:
            print(f"    {e['path']}")
    print(f"\n✓ DELIVERY.md   信任度 {a.trust} · 生命週期 {a.lifecycle} · commit {commit[:12]}")
    st = "✅ 可重建" if manifest["rebuildable"] else "🔴 無法重建（上游輸入缺失）"
    print(f"✓ REBUILD.md    {st}")
    for u in upstream:
        print(f"    {u['status']:<14}{u['role']:<22}{u.get('path') or '（未提供）'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
