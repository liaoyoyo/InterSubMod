#!/usr/bin/env python3
"""能力探測：決定這次 build 有哪些資料層可用，以及缺的要怎麼誠實交代。

這支的設計哲學與 build_exact_ps_layered_workstation.py 相反。那支是「缺任何一項就
整份拒繪」；這支只有 topology 是硬核心，其餘缺了就把面板降級並在頁面上寫明缺什麼、
探測停在哪一段。使用者要的是「有什麼顯示什麼」，不是「湊不齊就什麼都不給」。

探測分四段，任一段失敗就停在那段並記錄下來：
    S0 存在  → 路徑在不在
    S1 結構  → schema / header 對不對
    S2 體量  → 列數與 receipt 對不對得上
    S3 連結  → 與硬核心的 join 率（一律帶分母）
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

EXIT_REFUSE = 3

CORE = "CORE"
EXTENSION = "EXTENSION"

PRESENT = "PRESENT"
PARTIAL = "PARTIAL"
ABSENT = "ABSENT"
MALFORMED = "MALFORMED"

STATE_LABEL = {
    PRESENT: "可用",
    PARTIAL: "部分可用",
    ABSENT: "不存在",
    MALFORMED: "格式不符",
}


class Refuse(RuntimeError):
    """硬核心缺席。呼叫端印訊息後以 EXIT_REFUSE 離開，不產出半成品頁面。"""


@dataclass
class Probe:
    stage: str
    ok: bool
    detail: str


@dataclass
class Linkage:
    """與硬核心的連結率。分子分母都留著，因為只給比率會讓人看不出樣本量。"""
    numerator: int
    denominator: int
    method: str

    @property
    def rate(self) -> float:
        return self.numerator / self.denominator if self.denominator else 0.0

    def as_dict(self) -> dict:
        return {"numerator": self.numerator, "denominator": self.denominator,
                "rate": round(self.rate, 6), "method": self.method}


@dataclass
class Capability:
    id: str
    tier: str
    title: str
    state: str = ABSENT
    reason: str = ""
    paths: list = field(default_factory=list)
    probes: list = field(default_factory=list)
    linkage: Optional[Linkage] = None
    counts: dict = field(default_factory=dict)
    enables: list = field(default_factory=list)
    degrade_to: Optional[str] = None
    payload: Any = None          # 探測順便讀出來的資料，避免第二次 IO

    @property
    def usable(self) -> bool:
        return self.state in (PRESENT, PARTIAL)

    @property
    def failed_stage(self) -> Optional[str]:
        for p in self.probes:
            if not p.ok:
                return p.stage
        return None

    def as_dict(self) -> dict:
        return {
            "id": self.id, "tier": self.tier, "title": self.title,
            "state": self.state, "state_label": STATE_LABEL[self.state],
            "reason": self.reason,
            "failed_stage": self.failed_stage,
            "paths": [p.as_dict() if isinstance(p, FileRef) else str(p) for p in self.paths],
            "probes": [{"stage": p.stage, "ok": p.ok, "detail": p.detail} for p in self.probes],
            "linkage": self.linkage.as_dict() if self.linkage else None,
            "counts": self.counts,
            "enables": self.enables,
            "degrade_to": self.degrade_to,
        }


@dataclass
class FileRef:
    """一個輸入檔的身分證。頁面上要能一眼看出這份 HTML 是用哪些檔案產的。"""
    path: str
    exists: bool
    size_bytes: int = 0
    sha256: str = ""
    mtime: str = ""

    def as_dict(self) -> dict:
        return {"path": self.path, "exists": self.exists, "size_bytes": self.size_bytes,
                "sha256": self.sha256, "mtime": self.mtime}


def _iso(ts: float) -> str:
    import datetime
    return datetime.datetime.utcfromtimestamp(ts).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path, limit_bytes: int = 0) -> str:
    """整檔 SHA-256。limit_bytes > 0 時只雜湊前 N bytes（給 GB 級 BAM 用，
    並在 receipt 裡標明是部分雜湊，不假裝是全檔。）"""
    h = hashlib.sha256()
    read = 0
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            if limit_bytes:
                remain = limit_bytes - read
                if remain <= 0:
                    break
                chunk = chunk[:remain]
            h.update(chunk)
            read += len(chunk)
    return h.hexdigest()


def file_ref(path, hash_limit: int = 0) -> FileRef:
    p = Path(path)
    if not p.is_file():
        return FileRef(path=str(p), exists=False)
    stat = p.stat()
    return FileRef(path=str(p), exists=True, size_bytes=stat.st_size,
                   sha256=sha256_file(p, hash_limit), mtime=_iso(stat.st_mtime))


def open_maybe_gzip(path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", newline="")
    return open(p, "r", newline="")


def read_tsv_header(path) -> list:
    with open_maybe_gzip(path) as fh:
        line = fh.readline()
    return [c.strip().lstrip("﻿") for c in line.rstrip("\r\n").split("\t")]


def read_csv_header(path) -> list:
    with open_maybe_gzip(path) as fh:
        line = fh.readline()
    return [c.strip().lstrip("﻿") for c in line.rstrip("\r\n").split(",")]


def iter_tsv(path):
    with open_maybe_gzip(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield {k: (v.rstrip("\r") if isinstance(v, str) else v) for k, v in row.items()}


def first_jsonl(path) -> dict:
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                return json.loads(line)
    return {}


class Registry:
    """所有能力的集散地。除了存狀態，也負責把「這個面板為什麼不能用」講清楚。"""

    def __init__(self):
        self.caps: dict = {}
        self.order: list = []

    def add(self, cap: Capability) -> Capability:
        self.caps[cap.id] = cap
        self.order.append(cap.id)
        return cap

    def get(self, cap_id: str) -> Optional[Capability]:
        return self.caps.get(cap_id)

    def usable(self, cap_id: str) -> bool:
        c = self.caps.get(cap_id)
        return bool(c and c.usable)

    def require_core(self):
        for cid in self.order:
            c = self.caps[cid]
            if c.tier == CORE and not c.usable:
                raise Refuse(
                    f"硬核心能力 '{c.id}'（{c.title}）不可用 — 狀態 {c.state}"
                    f"，探測停在 {c.failed_stage or 'S0'}。\n  原因：{c.reason}"
                )

    def matrix(self) -> list:
        return [self.caps[cid].as_dict() for cid in self.order]

    def all_file_refs(self) -> list:
        out, seen = [], set()
        for cid in self.order:
            for ref in self.caps[cid].paths:
                if isinstance(ref, FileRef) and ref.path not in seen:
                    seen.add(ref.path)
                    out.append(ref)
        return out

    def summary_line(self) -> str:
        n = len(self.order)
        ok = sum(1 for cid in self.order if self.caps[cid].state == PRESENT)
        part = sum(1 for cid in self.order if self.caps[cid].state == PARTIAL)
        return f"能力 {ok}/{n} 完整、{part} 部分可用、{n - ok - part} 不可用"


def probe(cap: Capability, stage: str, ok: bool, detail: str) -> bool:
    """記一段探測結果。回傳 ok 讓呼叫端可以 `if not probe(...): return`。"""
    cap.probes.append(Probe(stage=stage, ok=ok, detail=detail))
    if not ok and not cap.reason:
        cap.reason = detail
    return ok


def probe_exists(cap: Capability, path, hash_limit: int = 0) -> bool:
    """S0：路徑存在嗎。順便建立 FileRef 供 receipt 使用。"""
    ref = file_ref(path, hash_limit)
    cap.paths.append(ref)
    if not ref.exists:
        cap.state = ABSENT
        probe(cap, "S0", False, f"檔案不存在：{ref.path}")
        return False
    probe(cap, "S0", True, f"{ref.size_bytes:,} bytes")
    return True


def probe_header(cap: Capability, actual: list, required: set, label: str) -> bool:
    """S1：header 必要欄位在不在。多出來的欄允許（上游加欄不該讓這裡壞掉），
    缺欄則降為 PARTIAL 並把缺的列出來。"""
    missing = sorted(required - set(actual))
    if missing:
        cap.state = MALFORMED
        probe(cap, "S1", False, f"{label} 缺必要欄位：{', '.join(missing)}")
        return False
    extra = len(set(actual) - required)
    probe(cap, "S1", True, f"{label} {len(actual)} 欄（必要 {len(required)} 全在，另有 {extra} 欄）")
    return True


def probe_schema(cap: Capability, rec: dict, expect_name: str) -> bool:
    got = rec.get("schema_name", "")
    if expect_name and got != expect_name:
        cap.state = MALFORMED
        probe(cap, "S1", False, f"schema_name 不符：期望 {expect_name}，實得 {got or '(無)'}")
        return False
    probe(cap, "S1", True, f"schema {got} v{rec.get('schema_version', '?')}")
    return True


def probe_count(cap: Capability, actual: int, expected: Optional[int], label: str) -> bool:
    """S2：列數與 receipt 對得上嗎。對不上不當錯誤，降 PARTIAL 並把兩個數字都留著。"""
    cap.counts[label] = actual
    if expected is None:
        probe(cap, "S2", True, f"{label} = {actual:,}（無 receipt 可比對）")
        return True
    cap.counts[label + "_expected"] = expected
    if actual != expected:
        cap.state = PARTIAL
        probe(cap, "S2", False, f"{label} 實得 {actual:,}，receipt 宣稱 {expected:,}（差 {actual - expected:+,}）")
        return False
    probe(cap, "S2", True, f"{label} = {actual:,}，與 receipt 一致")
    return True


def probe_linkage(cap: Capability, numerator: int, denominator: int, method: str) -> bool:
    """S3：與硬核心的 join 率。這裡不設通過門檻 —— 低就低，誠實顯示，
    由使用者決定那個面板可不可信。"""
    cap.linkage = Linkage(numerator=numerator, denominator=denominator, method=method)
    probe(cap, "S3", True, f"{numerator:,} / {denominator:,} = {cap.linkage.rate * 100:.2f}%（{method}）")
    return True


def finalize(cap: Capability, note: str = "") -> Capability:
    """探測跑完後定案。沒被降級過就是 PRESENT。"""
    if cap.state == ABSENT and any(p.ok for p in cap.probes):
        cap.state = PRESENT
    if cap.state not in (ABSENT, MALFORMED, PARTIAL):
        cap.state = PRESENT
    if note and not cap.reason:
        cap.reason = note
    return cap


def make(reg: Registry, cap_id: str, title: str, tier: str = EXTENSION,
         enables=None, degrade_to=None) -> Capability:
    return reg.add(Capability(id=cap_id, tier=tier, title=title,
                              enables=list(enables or []), degrade_to=degrade_to))


def guarded(fn: Callable[[], Capability]) -> Capability:
    """探測本身出錯（權限、編碼、壞檔）不該讓整個 build 掛掉 —— 那也是一種
    「這層不可用」，記成 MALFORMED 就好。"""
    cap = None
    try:
        return fn()
    except Refuse:
        raise
    except Exception as exc:  # noqa: BLE001 - 探測期任何例外都轉成能力不可用
        if cap is None:
            cap = Capability(id="unknown", tier=EXTENSION, title="unknown")
        cap.state = MALFORMED
        probe(cap, cap.probes[-1].stage if cap.probes else "S0", False,
              f"探測時發生例外：{type(exc).__name__}: {exc}")
        return cap
