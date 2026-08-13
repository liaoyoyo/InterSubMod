#!/usr/bin/env python3
"""把 research/ 底下每個專案目錄機械地分成「最終產物」與「中間產物」。

為什麼需要這支腳本
------------------
research/ 已累積到 98 個專案目錄，其中絕大多數沒有任何狀態標記。憑感覺判斷
「這個還要不要留」既不可重現，也無法交接 —— 下一個人（或下一個 AI）會得到
不同答案。這支腳本把判斷寫成可執行的規則，讓同樣的輸入永遠得到同樣的分類。

它疊在既有三軸之上，不發明第四套分類法：
  信任度      CURRENT / PRE-FIX / SUPERSEDED / DEPRECATED
              → docs/data_specs/20260414_output資料信任度與生命週期_01.md §1
  生命週期    in_progress / validated / finalized
              → docs/CLAUDE.md
  刪除安全度  SAFE_DELETE / LOW_RISK / KEEP
              → 同上 §2

本腳本輸出的「處置 bucket」是上述三軸的**函數**，而非替代品。

決策樹（第一個 YES 就定案）
---------------------------
  Q0  量得出體積嗎？        → 量不出來（PATHOLOGICAL）走 MANUAL，需專案級處置
  Q1  有人引用它嗎？        → 【最終】依體積走 A（進 git）或 B（留磁碟＋索引）
  Q2  30 分鐘內重建得出來？  → 【中間】走 C（觀察期後刪）
  Q3  命中中間物 pattern？   → 【中間-立即】走 D（直接刪）
  四問皆否                  → MANUAL，進人工裁決佇列，不要猜

兩階段設計，以及為什麼必須這樣
------------------------------
第一版用 os.walk 走整個 research/，跑超過 10 分鐘還沒完 —— 因為
`paired_priority_bug_audit/` 這類目錄裡有數百萬個檔案（每個 region 一個子目錄，
各含 metadata / reads.tsv / methylation.csv / distance_matrix / *.png）。
所以拆成：

    --scan    限時量測每個目錄，寫入 inventory 快取（慢，偶爾跑一次）
    （預設）   讀快取 + git 資料產生 ledger（快，可常跑）

限時是刻意的設計而非權宜：**一個連體積都量不出來的目錄，本質上就不可能是
「進 git」的候選**，所以量測逾時本身就是有意義的分類訊號，記成 PATHOLOGICAL。

刻意設計的保守偏誤
------------------
- Q1 只要有**一個**外部引用就判「最終」。誤判成最終的代價是多留一份檔案；
  誤判成中間的代價是刪掉還有人在引用的東西。兩者不對稱，所以偏向留。
- 自我引用不算數（專案自己的 README 提到自己不構成「有人在用」）。
- DEPRECATED 永遠不進 D。刪掉它會讓「這個結論是錯的」失去物證，
  一律留 ARCHIVED.md redirect。

用法
----
    python3 scripts/handoff/classify_artifacts.py --scan          # 重建 inventory 快取
    python3 scripts/handoff/classify_artifacts.py -o ledger.tsv   # 產生 ledger
    python3 scripts/handoff/classify_artifacts.py --check ledger.tsv   # 漂移偵測

--check 會重新分類並與既有 TSV 比對，不一致就 exit 1。用它確認 ledger 沒有跟
現況脫節 —— 這正是本 repo 幾份 INDEX 壞掉的方式（寫完就再也沒人對過）。
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RESEARCH = REPO / "research"
DEFAULT_CACHE = REPO / "docs/handoff/20260813_實驗室交接_01/research_inventory.tsv"

SKIP_DIRS = {"_template"}

# 量測單一目錄的時間上限（秒）。超過即記為 PATHOLOGICAL。
DU_TIMEOUT = 10
FIND_TIMEOUT = 6

# 檔名層級的中間物特徵。命中即為「產生過程的副產品」，不是任何人要看的東西。
INTERMEDIATE_PATTERNS = [
    re.compile(r"\.tmp-"),
    re.compile(r"_debug\."),
    re.compile(r"verification-failure"),
    re.compile(r"\.bak$"),
    re.compile(r"\.pyc$"),
    re.compile(r"_v\d+\.py$"),            # 同一支腳本的編號版本
    re.compile(r"report_extract_input"),
    re.compile(r"\.stdout$|\.stderr$"),
]

# 引用偵測掃描的檔案型別。刻意排除 docs/methodology/_assets/ —— 那裡有 66 MB
# 的單一 HTML，掃它只會拖慢而不會產生有意義的引用關係。
GREP_PATHSPEC = [
    "*.md", "*.json", "*.tsv", "*.py", "*.sh", "*.html", "*.yml",
    ":!docs/methodology/_assets",
]

# 體積門檻：超過就不適合進 git。GitHub 對單檔 100 MB 硬拒、50 MB 警告；
# 這裡取遠低於它的值，因為本 repo 的 pack 已經 200 MB 以上。
GIT_SIZE_LIMIT = 10 * 1024 * 1024

BUCKET_MEANING = {
    "A": "keep-in-git（最終產物，體積可進版本控制）",
    "B": "keep-on-disk-indexed（最終產物，體積過大，留磁碟並在 manifest 登錄路徑）",
    "C": "archive-then-delete（中間產物，可重建，觀察期後刪）",
    "D": "delete-now（中間物 pattern，無保留價值）",
    "MANUAL": "人工裁決（規則不適用，不要猜）",
}


class CommandFailed(RuntimeError):
    """指令以非預期的方式失敗。

    刻意設計成例外而非回傳空字串 —— 本腳本的第一版就是因為 git grep 的 regex
    不合法（POSIX ERE 不支援非捕獲群組），錯誤被吞成空字串，讓「零外部引用」
    這個明顯錯誤的結論看起來像是真實發現。
    """


def run(cmd: list[str], ok_codes: tuple[int, ...] = (0,)) -> str:
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=1800)
    if r.returncode not in ok_codes:
        raise CommandFailed(
            f"{' '.join(cmd[:5])}… exit={r.returncode}\n{r.stderr.strip()[:400]}"
        )
    return r.stdout


def project_dirs() -> list[Path]:
    return [d for d in sorted(RESEARCH.iterdir())
            if d.is_dir() and not d.name.startswith(".") and d.name not in SKIP_DIRS]


# ---------------------------------------------------------------- scan 階段

def scan(cache_path: Path) -> None:
    """限時量測每個專案目錄，寫入 inventory 快取。

    輸出欄位：project / bytes / files / builder_scripts / intermediate_hits
    量測逾時的目錄，bytes 記為 PATHOLOGICAL，其餘欄位記 -1。
    """
    rows = []
    dirs = project_dirs()
    for i, d in enumerate(dirs, 1):
        rel = f"research/{d.name}"
        try:
            out = subprocess.run(["du", "-sb", rel], cwd=REPO, capture_output=True,
                                 text=True, timeout=DU_TIMEOUT)
            nbytes = int(out.stdout.split("\t")[0]) if out.returncode == 0 else None
        except (subprocess.TimeoutExpired, ValueError, IndexError):
            nbytes = None

        if nbytes is None:
            rows.append((d.name, "PATHOLOGICAL", -1, -1, -1))
            print(f"[{i}/{len(dirs)}] {d.name}: 量測逾時 → PATHOLOGICAL", file=sys.stderr)
            continue

        files = builders = hits = 0
        try:
            listing = subprocess.run(
                ["find", rel, "-type", "f"], cwd=REPO,
                capture_output=True, text=True, timeout=FIND_TIMEOUT)
            for line in listing.stdout.splitlines():
                files += 1
                base = os.path.basename(line)
                if base.endswith((".py", ".sh", ".mjs")):
                    builders += 1
                if any(p.search(base) for p in INTERMEDIATE_PATTERNS):
                    hits += 1
        except subprocess.TimeoutExpired:
            files = builders = hits = -1
        rows.append((d.name, str(nbytes), files, builders, hits))

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8") as fh:
        fh.write("project\tbytes\tfiles\tbuilder_scripts\tintermediate_hits\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    n_bad = sum(1 for r in rows if r[1] == "PATHOLOGICAL")
    print(f"寫入 {cache_path}（{len(rows)} 個專案，其中 {n_bad} 個 PATHOLOGICAL）")


def load_cache(cache_path: Path) -> dict[str, dict]:
    if not cache_path.exists():
        raise SystemExit(
            f"找不到 inventory 快取 {cache_path}\n"
            f"先跑：python3 {Path(__file__).relative_to(REPO)} --scan"
        )
    out = {}
    for ln in cache_path.read_text(encoding="utf-8").splitlines()[1:]:
        parts = ln.split("\t")
        if len(parts) != 5:
            continue
        name, b, f, bs, ih = parts
        out[name] = {
            "pathological": b == "PATHOLOGICAL",
            "bytes": 0 if b == "PATHOLOGICAL" else int(b),
            "files": int(f),
            "builder_scripts": int(bs),
            "intermediate_hits": int(ih),
        }
    return out


# ------------------------------------------------------------- 引用偵測

def collect_reference_counts(names: list[str]) -> dict[str, dict[str, set[str]]]:
    """統計每個專案目錄名被哪些 tracked 檔案提到，並區分引用的強度。

    回傳 {專案名: {"strong": {檔案…}, "weak": {檔案…}}}

    區分強弱是必要的，不是講究。research/ 裡有 `subclonal_reconstruction`、
    `feature_layered_observation`、`external_tools`、`data_registry` 這類
    **通用名稱**目錄 —— 它們同時也是這個 repo 的日常詞彙與 skill 名稱。
    只比對裸名的話，任何一篇談到「subclonal reconstruction」的文件都會被算成
    「有人在引用這個目錄」，於是所有東西都變成「最終產物」，ledger 就退化成
    「全部保留」而失去意義（第一版正是如此：C=1、D=0）。

      strong  該行含 `research/<名稱>` —— 明確指向這個目錄，不會誤判
      weak    只出現裸名 —— 可能只是在講這個概念

    把 98 個目錄名一次餵給 `git grep -F -f`（固定字串比對，不碰 regex），
    再在 Python 端分類。用一次 git grep 而非 98 次迴圈，是因為 tracked 檔案
    夠多，後者會慢到不堪用。實測單次約 10 秒。

    exit code 1 是「沒有任何命中」，屬正常；其餘非零一律拋出。
    """
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False,
                                     encoding="utf-8") as fh:
        fh.write("\n".join(names) + "\n")
        pattern_file = fh.name
    try:
        out = run(["git", "grep", "-I", "-F", "-f", pattern_file, "-n", "--"]
                  + GREP_PATHSPEC, ok_codes=(0, 1))
    finally:
        os.unlink(pattern_file)

    refs: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"strong": set(), "weak": set()})
    # 長名字優先比對，避免短名（如 loh_subclone_af）被含它的長名那行誤算
    ordered = sorted(names, key=len, reverse=True)
    for line in out.splitlines():
        parts = line.split(":", 2)
        if len(parts) < 3:
            continue
        src, _lineno, text = parts
        for name in ordered:
            if name not in text:
                continue
            kind = "strong" if f"research/{name}" in text else "weak"
            refs[name][kind].add(src)
    return refs


def tracked_stats() -> dict[str, tuple[int, int]]:
    """一次取得每個專案目錄下有幾個檔進了 git、合計多少 bytes。"""
    out = run(["git", "ls-files", "--", "research/"])
    agg: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for p in out.splitlines():
        parts = p.split("/", 2)
        if len(parts) < 3:
            continue
        proj = parts[1]
        agg[proj][0] += 1
        try:
            agg[proj][1] += (REPO / p).stat().st_size
        except OSError:
            pass
    return {k: (v[0], v[1]) for k, v in agg.items()}


def find_status_marker(d: Path) -> str:
    """從 README / 主要 .md 的 metadata 區塊抓狀態標記。

    支援兩種格式：HTML comment（AGENTS.md §16 明文規則）與 YAML frontmatter
    （少數 2026-08-06 產出的檔案誤用了它）。兩種都讀，輸出統一。
    """
    marker_re = re.compile(
        r"^\s*(?:狀態|status)\s*[:：]\s*([A-Za-z_\-]+)", re.IGNORECASE | re.MULTILINE)
    trust_re = re.compile(r"\b(CURRENT|PRE-FIX|SUPERSEDED|DEPRECATED)\b")
    candidates = []
    for pat in ("README*.md", "*.md", "docs/*.md"):
        candidates.extend(sorted(d.glob(pat))[:4])
    for f in candidates[:8]:
        try:
            head = f.read_text(encoding="utf-8", errors="replace")[:4000]
        except OSError:
            continue
        if (m := marker_re.search(head)):
            return m.group(1).lower()
        if (t := trust_re.search(head)):
            return t.group(1)
    return ""


def classify(row: dict) -> tuple[str, str]:
    """套決策樹，回傳 (bucket, 依據)。第一個 YES 就定案。"""
    # 保護規則優先於決策樹：DEPRECATED 不可靜默刪
    if row["status_marker"].upper() == "DEPRECATED":
        return "B", "DEPRECATED 保留物證（不可靜默刪，須留 ARCHIVED.md redirect）"

    # Q0 連體積都量不出來 → 專案級處置，不是逐檔判斷能解決的
    if row["pathological"]:
        return "MANUAL", (f"Q0 體積量測逾時（>{DU_TIMEOUT}s）＝ 檔案數量級失控，"
                          "需專案級處置而非逐檔分類")

    mb = row["bytes"] / 1e6
    # Q1 有人以路徑形式引用它嗎（strong ＝ 明確指向這個目錄）
    if row["refs_strong"] > 0:
        if row["bytes"] <= GIT_SIZE_LIMIT:
            return "A", f"Q1 有 {row['refs_strong']} 處路徑引用，{mb:.1f} MB 可進 git"
        return "B", f"Q1 有 {row['refs_strong']} 處路徑引用，但 {mb:.1f} MB 過大，留磁碟"

    # Q2 重建得出來嗎（有 tracked builder ＝ 有重建路徑）
    if row["builder_scripts"] > 0 and row["tracked_files"] > 0:
        return "C", (f"Q2 無路徑引用，但有 {row['builder_scripts']} 支 builder "
                     f"且 {row['tracked_files']} 檔已 tracked，可重建")

    # Q3 命中中間物 pattern。有裸名提及者不判 D —— 名字被人講過就還不到
    # 「直接刪」的把握，降一級到 C 讓它走觀察期。
    if row["intermediate_hits"] > 0:
        if row["refs_weak"] > 0:
            return "C", (f"Q3 命中 {row['intermediate_hits']} 個中間物 pattern，"
                         f"但有 {row['refs_weak']} 處裸名提及 → 降級走觀察期而非直接刪")
        return "D", f"Q3 命中 {row['intermediate_hits']} 個中間物 pattern 且無人提及"

    if row["refs_weak"] > 0:
        return "MANUAL", (f"僅 {row['refs_weak']} 處裸名提及、無路徑引用 —— "
                          "目錄名可能是通用詞，須人工確認是真引用還是同名巧合")
    return "MANUAL", "四問皆不適用，需人工裁決"


COLUMNS = [
    "project", "bucket", "rationale", "refs_strong", "refs_weak", "ref_sources",
    "files", "bytes", "tracked_files", "tracked_bytes",
    "builder_scripts", "intermediate_hits", "status_marker",
]


def build_rows(cache_path: Path) -> list[dict]:
    cache = load_cache(cache_path)
    dirs = project_dirs()
    names = [d.name for d in dirs]
    refs = collect_reference_counts(names)
    tracked = tracked_stats()

    rows = []
    for d in dirs:
        c = cache.get(d.name)
        if c is None:
            print(f"警告：{d.name} 不在 inventory 快取中，請重跑 --scan", file=sys.stderr)
            continue
        rel = f"research/{d.name}"
        r = refs.get(d.name, {"strong": set(), "weak": set()})
        # 自我引用不算數：專案自己的 README 提到自己，不構成「有人在用」
        strong = {s for s in r["strong"] if not s.startswith(rel + "/")}
        weak = {s for s in r["weak"] if not s.startswith(rel + "/")} - strong
        tf, tb = tracked.get(d.name, (0, 0))
        row = {
            "project": d.name,
            "pathological": c["pathological"],
            "files": c["files"],
            "bytes": c["bytes"],
            "builder_scripts": max(c["builder_scripts"], 0),
            "intermediate_hits": max(c["intermediate_hits"], 0),
            "tracked_files": tf,
            "tracked_bytes": tb,
            "refs_strong": len(strong),
            "refs_weak": len(weak),
            "status_marker": find_status_marker(d),
        }
        row["bucket"], row["rationale"] = classify(row)
        row["ref_sources"] = ";".join(sorted(strong)[:3]) or ";".join(sorted(weak)[:2])
        rows.append(row)
    return rows


def to_tsv(rows: list[dict]) -> str:
    lines = ["\t".join(COLUMNS)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in COLUMNS))
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", action="store_true", help="重建 inventory 快取（慢）")
    ap.add_argument("--cache", default=str(DEFAULT_CACHE), help="inventory 快取路徑")
    ap.add_argument("-o", "--output", help="寫入 ledger TSV 的路徑")
    ap.add_argument("--check", metavar="TSV", help="重新分類並比對，不一致 exit 1")
    args = ap.parse_args()

    if not RESEARCH.is_dir():
        print(f"找不到 {RESEARCH}", file=sys.stderr)
        return 2

    cache_path = Path(args.cache)
    if args.scan:
        scan(cache_path)
        return 0

    rows = build_rows(cache_path)
    tsv = to_tsv(rows)

    if args.check:
        old = Path(args.check).read_text(encoding="utf-8")
        if old == tsv:
            print(f"OK  {len(rows)} 個專案，分類與 {args.check} 一致")
            return 0
        old_map = {ln.split("\t")[0]: ln for ln in old.splitlines()[1:]}
        new_map = {ln.split("\t")[0]: ln for ln in tsv.splitlines()[1:]}
        for k in sorted(set(old_map) | set(new_map)):
            if old_map.get(k) != new_map.get(k):
                print(f"DRIFT  {k}")
                print(f"  舊: {old_map.get(k, '(不存在)')[:150]}")
                print(f"  新: {new_map.get(k, '(已移除)')[:150]}")
        return 1

    if args.output:
        Path(args.output).write_text(tsv, encoding="utf-8")
        print(f"寫入 {args.output}（{len(rows)} 個專案）")

    counts: dict[str, int] = defaultdict(int)
    sizes: dict[str, int] = defaultdict(int)
    tracked_n: dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r["bucket"]] += 1
        sizes[r["bucket"]] += r["bytes"]
        tracked_n[r["bucket"]] += r["tracked_files"]
    print(f"\n研究專案總數：{len(rows)}")
    print(f"{'bucket':<8}{'數量':>5}{'磁碟體積':>12}{'tracked檔':>10}   意義")
    for b in ["A", "B", "C", "D", "MANUAL"]:
        print(f"{b:<8}{counts[b]:>5}{sizes[b]/1e9:>11.2f}G{tracked_n[b]:>10}   "
              f"{BUCKET_MEANING[b]}")
    print(f"{'合計':<8}{len(rows):>5}{sum(sizes.values())/1e9:>11.2f}G"
          f"{sum(tracked_n.values()):>10}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
