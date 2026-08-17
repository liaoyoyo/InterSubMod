#!/usr/bin/env python3
"""檢查 README 記載的 LongLineage 版本，是否還等於它公開 main 的實際版本。

為什麼需要這支腳本
------------------
InterSubMod 的 README 會宣告搭檔 repository LongLineage 的能力狀態，而那個狀態是
**版本限定**的 —— 例如「public main 是否提供 longlineage-tag-bam」。這種宣告一旦
寫死就會腐爛：2026-08-17 實測時，README 記載的是 `5daf50f`（2026-07-19），
而 LongLineage 的公開 main 早已前進到 `583e03e`，中間隔了三個 PR、將近一個月。
更糟的是那段文字宣稱「public main 不支援 tagged BAM」，但當時 main 其實已經支援。

讀者無從察覺這種漂移 —— 兩邊看起來都很肯定，只是其中一邊已經不成立。
所以把「有沒有漂移」變成可執行的檢查，而不是靠人記得回來對。

判斷方式
--------
從 README 抽出被記載的 LongLineage commit，跟實際的公開 main 比對。
取得實際版本有兩條路徑，依序嘗試：

  1. 本機的 LongLineage repo（快，且不需網路）
  2. `git ls-remote`（無本機 repo 時的後備）

兩條都失敗時回報 SKIP 而非 FAIL —— 拿不到真值時不該假裝知道答案。

用法
----
    python3 scripts/handoff/check_companion_version.py
    python3 scripts/handoff/check_companion_version.py --quiet   # 只回 exit code

exit code：0 = 一致或無法取得真值　1 = 偵測到漂移
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LOCAL_COMPANION = Path("/big7_disk/liaoyoyo2001/LongLineage")
COMPANION_URL = "https://github.com/liaoyoyo/LongLineage.git"

README_FILES = ["README.md", "README.zh-TW.md"]

# README 用一個機器可讀的標記宣告搭檔版本：
#     <!-- companion-version: LongLineage public main = 583e03e -->
#
# 刻意不從散文抓。第一版用 regex 去比對「public main `<sha>`」這種句式，結果雙語
# README 一改寫法（改成「LongLineage **公開** `main` `583e03e`」）就整個對不上，
# 腳本回報 SKIP —— 而 SKIP 看起來像「沒問題」，等於偵測能力靜默消失。
# 用標記則散文可自由改寫，契約仍然成立。
RECORDED_RE = re.compile(
    r"<!--\s*companion-version:\s*LongLineage\s+public\s+main\s*=\s*([0-9a-f]{7,40})\s*-->"
)


def run(cmd: list[str], cwd: Path | None = None, timeout: int = 30) -> str | None:
    """跑指令回傳 stdout；任何失敗都回 None，由呼叫端決定語意。"""
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return r.stdout.strip() if r.returncode == 0 else None
    except (subprocess.TimeoutExpired, OSError):
        return None


def actual_companion_main() -> tuple[str | None, str]:
    """取得 LongLineage 公開 main 的實際 commit，回傳 (sha, 來源說明)。"""
    if LOCAL_COMPANION.is_dir():
        run(["git", "fetch", "--quiet", "origin"], cwd=LOCAL_COMPANION, timeout=90)
        sha = run(["git", "rev-parse", "origin/main"], cwd=LOCAL_COMPANION)
        if sha:
            return sha, f"本機 repo {LOCAL_COMPANION}"

    out = run(["git", "ls-remote", COMPANION_URL, "refs/heads/main"], timeout=90)
    if out:
        return out.split()[0], f"git ls-remote {COMPANION_URL}"

    return None, "取不到（本機 repo 不存在且 ls-remote 失敗）"


def recorded_versions() -> tuple[list[tuple[str, int, str]], list[str]]:
    """從每份 README 抽出被記載的搭檔版本。

    回傳 (命中清單, 缺標記的檔案清單)。

    缺標記要逐檔回報，不能只看「整體有沒有找到」—— 第一版就是那樣寫的，
    結果把 README.md 的標記刪掉、而 README.zh-TW.md 還留著時，腳本回報通過。
    雙語 README 各自對外，任一份漏掉都是真實的偵測缺口。
    """
    found: list[tuple[str, int, str]] = []
    missing: list[str] = []
    for name in README_FILES:
        path = REPO / name
        if not path.is_file():
            missing.append(f"{name}（檔案不存在）")
            continue
        hits_before = len(found)
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for m in RECORDED_RE.finditer(line):
                found.append((name, lineno, m.group(1)))
        if len(found) == hits_before:
            missing.append(name)
    return found, missing


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true", help="只回 exit code，不印明細")
    args = ap.parse_args()

    def say(*a):
        if not args.quiet:
            print(*a)

    recorded, missing = recorded_versions()
    if missing:
        # 這不是 SKIP 而是 FAIL：標記消失代表偵測機制本身壞了，
        # 而「沒有偵測到漂移」與「無法偵測漂移」必須是不同的結果。
        say("FAIL  下列檔案缺少 companion-version 標記：")
        for m in missing:
            say(f"        {m}")
        say("      預期格式：<!-- companion-version: LongLineage public main = <sha> -->")
        return 1

    actual, source = actual_companion_main()
    if actual is None:
        say(f"SKIP  無法取得 LongLineage 公開 main 的實際版本（{source}）。")
        say("      拿不到真值時不判定漂移。")
        return 0

    say(f"實際 LongLineage 公開 main：{actual[:12]}")
    say(f"  來源：{source}")
    say("")

    drifted = []
    for name, lineno, sha in recorded:
        ok = actual.startswith(sha)
        say(f"  {'OK  ' if ok else 'DRIFT'}  {name}:{lineno}  記載 {sha}")
        if not ok:
            drifted.append((name, lineno, sha))

    if not drifted:
        say(f"\n一致：{len(recorded)} 處記載全部等於實際的公開 main。")
        return 0

    say(f"\n偵測到 {len(drifted)} 處版本漂移。")
    say(f"把上列位置的 commit 改成 {actual[:7]}，並**同時檢查該處的能力陳述是否仍成立** ——")
    say("漂移的危險不在 hash 本身，而在附帶的「支援／不支援」敘述可能已經反轉。")
    say("2026-08-17 就發生過：README 說 public main 不提供 tagged BAM，但當時它已經提供了。")
    return 1


if __name__ == "__main__":
    sys.exit(main())
