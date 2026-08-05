#!/usr/bin/env python3
"""Bind verified computation evidence into the standalone teaching HTML."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


PLACEHOLDER = "__INTERSUBMOD_VERIFIED_EVIDENCE_JSON__"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    template = args.template.read_text(encoding="utf-8")
    if template.count(PLACEHOLDER) != 1:
        raise ValueError("template must contain the evidence placeholder exactly once")
    evidence = json.loads(args.evidence.read_text(encoding="utf-8"))
    if evidence.get("all_pass") is not True:
        raise ValueError("refusing to bind evidence whose all_pass is not true")
    compact = json.dumps(
        evidence,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).replace("</script", r"<\/script")
    rendered = template.replace(PLACEHOLDER, compact)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "evidence_sha256": sha256_file(args.evidence),
                "html_sha256": sha256_file(args.output),
                "output": str(args.output.resolve()),
                "size_bytes": args.output.stat().st_size,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
