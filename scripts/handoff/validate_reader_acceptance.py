#!/usr/bin/env python3
"""Fail-closed schema, provenance and semantic validation for a fresh reader receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

try:
    from jsonschema import Draft202012Validator, FormatChecker
    from jsonschema.exceptions import SchemaError
except ImportError:  # pragma: no cover - exercised by the CLI dependency error path.
    Draft202012Validator = None
    FormatChecker = None
    SchemaError = Exception


PACKAGE_RELATIVE = Path("docs/handoff/20260813_完整研究資料與軟體交接_01")
SCHEMA_RELATIVE = Path("schemas/reader-acceptance.schema.json")
QUESTIONS = [
    "Q_PROJECT",
    "Q_CONCLUSION",
    "Q_FINALITY",
    "Q_SOFTWARE_ROLES",
    "Q_MACHINES",
    "Q_VERIFY_CONTINUE",
]
PROMOTIONS = [
    "NO_CELLULAR_PROMOTION",
    "NO_ANCESTRY_PROMOTION",
    "NO_882579_ACCURACY_OR_PREVALENCE",
    "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE",
    "NO_FEATURE_AS_MAIN",
    "NO_LOCAL_AS_LIVE_PUBLISHED",
    "NO_BIP7_AS_BIP8",
    "NO_PYTHON_ROLE_CONFLATION",
]
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
FLAGS = re.IGNORECASE | re.DOTALL

# The receipt is intentionally not in this set.  A reader can evaluate a reachable
# source commit and add the receipt afterwards without making the receipt self-hash.
REQUIRED_MANIFEST_PATHS = {
    "00_INDEX.md",
    "20260813_研究結論時間與Finality_01.md",
    "20260813_軟體輸入輸出與研究流程_01.md",
    "20260813_bip7_bip8操作與驗證_01.md",
    "ai_context/CONTEXT.md",
    "ai_context/READER_ACCEPTANCE_PROMPT.md",
    "schemas/reader-acceptance.schema.json",
    "registries/artifact_registry.json",
    "registries/authority_superseded_crosswalk.json",
    "registries/claim_registry.json",
    "registries/machine_path_registry.json",
    "evidence/authority_manifest.json",
    "evidence/denominator_registry.tsv",
    "evidence/authority_replay_receipt.json",
    "evidence/longlineage_capability_matrix.md",
}

QUESTION_EVIDENCE = {
    "Q_PROJECT": {"00_INDEX.md", "20260813_軟體輸入輸出與研究流程_01.md"},
    "Q_CONCLUSION": {
        "20260813_研究結論時間與Finality_01.md",
        "evidence/denominator_registry.tsv",
    },
    "Q_FINALITY": {
        "20260813_研究結論時間與Finality_01.md",
        "registries/artifact_registry.json",
        "registries/authority_superseded_crosswalk.json",
    },
    "Q_SOFTWARE_ROLES": {
        "20260813_軟體輸入輸出與研究流程_01.md",
        "evidence/longlineage_capability_matrix.md",
    },
    "Q_MACHINES": {
        "20260813_bip7_bip8操作與驗證_01.md",
        "registries/machine_path_registry.json",
    },
    "Q_VERIFY_CONTINUE": {
        "00_INDEX.md",
        "evidence/authority_replay_receipt.json",
        "registries/claim_registry.json",
    },
}

PROMOTION_EVIDENCE = {
    "NO_CELLULAR_PROMOTION": {"20260813_研究結論時間與Finality_01.md"},
    "NO_ANCESTRY_PROMOTION": {"20260813_研究結論時間與Finality_01.md"},
    "NO_882579_ACCURACY_OR_PREVALENCE": {
        "20260813_研究結論時間與Finality_01.md",
        "evidence/denominator_registry.tsv",
    },
    "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE": {
        "00_INDEX.md",
        "evidence/authority_manifest.json",
    },
    "NO_FEATURE_AS_MAIN": {"evidence/longlineage_capability_matrix.md"},
    "NO_LOCAL_AS_LIVE_PUBLISHED": {"00_INDEX.md", "registries/claim_registry.json"},
    "NO_BIP7_AS_BIP8": {
        "20260813_bip7_bip8操作與驗證_01.md",
        "registries/machine_path_registry.json",
    },
    "NO_PYTHON_ROLE_CONFLATION": {"20260813_軟體輸入輸出與研究流程_01.md"},
}

# Each tuple is (human-readable concept, alternatives).  At least one alternative
# must occur for every concept.  Patterns intentionally cover concise English and
# Traditional-Chinese answers while retaining the scientific boundary.
QUESTION_RUBRICS: dict[str, list[tuple[str, tuple[str, ...]]]] = {
    "Q_PROJECT": [
        ("research handoff snapshot", (r"research handoff (?:snapshot|package)", r"研究交接(?:快照|包)")),
        ("ONT long reads", (r"(?:ont|nanopore).{0,30}long[- ]?reads?", r"(?:ont|奈米孔).{0,30}長讀")),
        ("sSNV read linkage", (r"(?:ssnv|somatic variant).{0,50}read linkage", r"體細胞突變.{0,50}read.{0,20}(?:linkage|連鎖)")),
        ("candidate reconstruction", (r"candidate (?:mutation[- ]state )?reconstruction", r"候選.{0,30}(?:mutation[- ]state|突變狀態).{0,20}重建")),
        ("methylation association", (r"methylation association", r"甲基化.{0,20}(?:association|關聯)")),
        ("not production/cellular lineage caller", (r"not.{0,50}(?:production release|cellular[- ]lineage caller)", r"不是.{0,50}(?:production|細胞譜系).{0,20}(?:release|caller|軟體)")),
    ],
    "Q_CONCLUSION": [
        ("zero confirmed cellular subclones", (r"confirmed cellular subclones?\s*(?:=|:|remain(?:s)?|is|are)?\s*0\b", r"(?:確認的?|confirmed).{0,20}細胞亞群.{0,15}(?:=|:|為|是)?\s*0\b")),
        ("zero confirmed linear ancestry", (r"confirmed linear ancestry(?: relationships?)?\s*(?:=|:|remain(?:s)?|is|are)?\s*0\b", r"(?:確認的?|confirmed).{0,20}線性祖先(?:關係)?.{0,15}(?:=|:|為|是)?\s*0\b")),
        ("methylation association only", (r"methylation.{0,40}association only", r"甲基化.{0,30}僅(?:為|是).{0,20}(?:association|關聯)")),
        ("frozen CN/LOH not integrated", (r"frozen.{0,70}(?:cn/loh|cn.{0,10}loh).{0,70}not integrated", r"(?:cn/loh|cn.{0,10}loh).{0,70}(?:凍結|frozen).{0,40}未整合")),
        ("88.2579 model-conditional graph shape", (r"88\.2579%?.{0,100}model[- ]conditional.{0,80}graph[- ]shape", r"88\.2579%?.{0,100}(?:模型條件|model).{0,80}(?:graph[- ]shape|圖形)")),
        ("88.2579 is not accuracy", (r"88\.2579%?.{0,140}(?:not|isn['’]?t).{0,35}accuracy", r"88\.2579%?.{0,140}不是.{0,35}(?:accuracy|正確率)")),
        ("88.2579 is not prevalence", (r"88\.2579%?.{0,180}(?:not|isn['’]?t).{0,35}prevalence", r"88\.2579%?.{0,180}不是.{0,35}(?:prevalence|盛行率|佔比)")),
    ],
    "Q_FINALITY": [
        ("FINAL_FOR_SCOPE", (r"\bfinal_for_scope\b",)),
        ("evidence_status", (r"\bevidence_status\b",)),
        ("finality field", (r"\bfinality\b",)),
        ("supersedes relation", (r"\bsupersedes\b",)),
        ("scope-bounded finality", (r"final.{0,60}(?:only|bounded).{0,50}scope", r"final.{0,50}僅.{0,30}scope")),
        ("not production/whole-research final", (r"not.{0,60}(?:production[- ]ready|whole research final)", r"不(?:代表|是).{0,60}(?:production|整體研究).{0,30}final")),
    ],
    "Q_SOFTWARE_ROLES": [
        ("LongPhase S/TO products", (r"longphase[- ]?(?:s|s/to|to).{0,100}(?:hp/ps|phased|recalibrated|tagged bam)",)),
        ("exact-PS/LongLineage candidate reconstruction", (r"(?:exact[- ]ps|longlineage).{0,100}(?:candidate famil|read assignment|candidate reconstruction)",)),
        ("InterSubMod read-level outputs", (r"intersubmod.{0,120}per[- ]region methylation.{0,100}(?:distance|read clustering|statistics)",)),
        ("two parallel provenance chains", (r"(?:two|dual|second|parallel).{0,50}provenance (?:chain|chains|branch|branches)", r"兩條.{0,30}provenance.{0,30}(?:chain|分支|鏈)")),
        ("Python science producer", (r"commit[- ]pinned python.{0,80}(?:research solver|analysis script).{0,60}science producer",)),
        ("presenter does not recompute science", (r"(?:validator|publication builder|html presenter).{0,100}(?:validated|validate|present).{0,100}(?:does not|must not|not).{0,50}recompute science", r"(?:validator|publication builder|html).{0,100}(?:驗證|呈現).{0,100}不.{0,30}重算 science")),
    ],
    "Q_MACHINES": [
        ("bip7 remains bounded/blocked", (r"bip7.{0,120}(?:blocked|fresh[- ]clone.{0,30}(?:pending|required)|receipt.{0,30}(?:bounded|preflight))", r"bip7.{0,120}(?:尚待|未完成|阻塞).{0,50}fresh[- ]clone")),
        ("bip8 is blocked without host receipt", (r"bip8.{0,120}(?:blocked|no.{0,50}host[- ]local receipt|has no.{0,50}receipt)", r"bip8.{0,120}(?:無|尚無|沒有).{0,50}(?:host|hostname|主機).{0,30}receipt")),
        ("host-specific receipts are not substitutable", (r"(?:host[- ]specific|host[- ]local).{0,80}receipt.{0,100}(?:cannot|must not|not).{0,60}(?:substitut|prove the other)", r"主機.{0,30}receipt.{0,80}不(?:能|得).{0,40}(?:取代|代表)")),
        ("doctor/build/smoke sequence", (r"doctor.{0,100}(?:build|test).{0,100}(?:synthetic|smoke)",)),
    ],
    "Q_VERIFY_CONTINUE": [
        ("19/19 authority byte replay", (r"19\s*/\s*19.{0,80}(?:byte|sha[- ]?256|hash).{0,50}(?:match|replay)", r"19\s*/\s*19.{0,80}(?:byte|sha[- ]?256|hash).{0,50}(?:相符|重播)")),
        ("not a science rerun", (r"19\s*/\s*19.{0,160}(?:not|does not).{0,60}(?:science rerun|science recomputation)", r"19\s*/\s*19.{0,160}不(?:是|代表).{0,50}(?:science|科學).{0,20}重跑")),
        ("registry validation", (r"registr(?:y|ies).{0,80}(?:validat|schema|unique)", r"registry.{0,80}(?:驗證|schema|唯一)")),
        ("blocking gates", (r"(?:release|publication|machine).{0,80}gates?.{0,60}(?:pass|blocked|fail[- ]closed)", r"gate.{0,80}(?:pass|blocked|阻塞)")),
        ("new cycle pins provenance", (r"new research cycle.{0,120}(?:pre[- ]decision audit).{0,160}(?:commit|input hash).{0,100}(?:scope|denominator)", r"新.{0,20}(?:research )?cycle.{0,120}(?:pre[- ]decision audit).{0,160}(?:commit|input hash).{0,100}(?:scope|denominator)")),
    ],
}

QUESTION_FORBIDDEN = {
    "Q_PROJECT": (r"(?:is|as) (?:a )?production[- ]ready", r"已.{0,20}production[- ]ready"),
    "Q_CONCLUSION": (
        r"confirmed cellular subclones?\s*(?:=|:|is|are)\s*[1-9]",
        r"confirmed linear ancestry(?: relationships?)?\s*(?:=|:|is|are)\s*[1-9]",
        r"88\.2579%?\s*(?:is|=|means|represents)\s*(?:the )?(?:biological )?(?:accuracy|prevalence)",
    ),
    "Q_FINALITY": (r"(?:filename|directory name|mtime).{0,50}(?:proves?|means?).{0,30}final", r"final_for_scope.{0,50}(?:means?|proves?).{0,30}production"),
    "Q_SOFTWARE_ROLES": (r"longlineage.{0,100}(?:topology|jsonl).{0,80}(?:is|becomes|directly).{0,30}intersubmod.{0,30}vcf", r"python(?:/html| and html)?.{0,30}only (?:presents?|renders?)"),
    "Q_MACHINES": (r"(?:bip7 and bip8|both machines).{0,40}(?:pass|complete|verified)",),
    "Q_VERIFY_CONTINUE": (r"19\s*/\s*19.{0,80}(?:proves?|means?).{0,50}(?:science rerun|scientifically correct|biological truth)",),
}

PROMOTION_RUBRICS: dict[str, list[tuple[str, tuple[str, ...]]]] = {
    "NO_CELLULAR_PROMOTION": [
        ("zero confirmed cellular subclones", QUESTION_RUBRICS["Q_CONCLUSION"][0][1]),
        ("candidate/read groups are not cellular truth", (r"(?:candidate|read|molecule).{0,80}(?:not|isn['’]?t).{0,40}(?:cellular truth|cellular subclone)", r"候選.{0,60}不是.{0,30}細胞.{0,20}真值")),
    ],
    "NO_ANCESTRY_PROMOTION": [
        ("zero confirmed ancestry", QUESTION_RUBRICS["Q_CONCLUSION"][1][1]),
        ("graph/read structure is not ancestry", (r"(?:graph shape|read dendrogram|local block).{0,80}(?:not|isn['’]?t).{0,40}(?:ancestry|phylogeny)", r"(?:圖形|read 樹|局部 block).{0,80}不是.{0,30}(?:祖先|譜系)")),
    ],
    "NO_882579_ACCURACY_OR_PREVALENCE": [
        ("88.2579 model-conditional graph shape", QUESTION_RUBRICS["Q_CONCLUSION"][4][1]),
        ("not accuracy", QUESTION_RUBRICS["Q_CONCLUSION"][5][1]),
        ("not prevalence", QUESTION_RUBRICS["Q_CONCLUSION"][6][1]),
    ],
    "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE": [
        ("seven technical datasets", (r"7 technical datasets?", r"7 個 technical datasets?")),
        ("six biological IDs", (r"6 biological ids?", r"6 個 biological ids?")),
        ("not seven biological replicates", (r"(?:not|isn['’]?t).{0,50}7 (?:independent )?biological (?:samples|replicates)", r"不是.{0,50}7 個(?:獨立)?生物(?:樣本|複本)")),
    ],
    "NO_FEATURE_AS_MAIN": [
        ("b9aaa12 preview", (r"b9aaa12.{0,80}(?:research preview|candidate)",)),
        ("non-production", (r"non[- ]production", r"非 production")),
        ("blocked phases", (r"p3/p4/p5/p7/p8.{0,30}blocked",)),
        ("feature is not main/release", (r"feature.{0,60}(?:not|isn['’]?t).{0,40}(?:main|release|production)", r"feature.{0,60}不(?:是|代表).{0,40}(?:main|release|production)")),
    ],
    "NO_LOCAL_AS_LIVE_PUBLISHED": [
        ("local source is not live publication", (r"local source.{0,80}(?:not|isn['’]?t).{0,50}live publication", r"local source.{0,80}不(?:是|等於).{0,40}(?:live|公開).{0,20}發布")),
        ("main Wiki Pages blocked", (r"main.{0,40}wiki.{0,40}pages.{0,80}blocked",)),
        ("publish and refetch required", (r"publish.{0,40}re[- ]?fetch", r"發布.{0,40}重新抓取")),
    ],
    "NO_BIP7_AS_BIP8": [
        ("bip7 cannot prove bip8", (r"bip7.{0,80}(?:cannot|does not|must not).{0,50}(?:prove|stand in for|substitute).{0,30}bip8", r"bip7.{0,80}不(?:能|得).{0,40}(?:代表|取代).{0,30}bip8")),
        ("bip8 blocked", (r"bip8.{0,80}blocked",)),
        ("host-local receipt needed", (r"host[- ]local receipt", r"主機.{0,30}receipt")),
    ],
    "NO_PYTHON_ROLE_CONFLATION": [
        ("Python solver can produce science", (r"commit[- ]pinned python.{0,80}(?:solver|analysis).{0,60}(?:can be|is).{0,40}science producer",)),
        ("validator/presenter role", (r"(?:validator|publication builder|html presenter).{0,100}(?:validates?|presents?|renders?)",)),
        ("presenter does not recompute science", (r"(?:does not|must not|cannot).{0,50}recompute science", r"不(?:能|得|會).{0,30}重算 science")),
    ],
}

PROMOTION_FORBIDDEN = {
    "NO_CELLULAR_PROMOTION": (r"confirmed cellular subclones?\s*(?:=|:|is|are)\s*[1-9]", r"(?:candidate|read group).{0,50}(?:proves?|confirms?).{0,30}cellular subclone"),
    "NO_ANCESTRY_PROMOTION": (r"confirmed linear ancestry(?: relationships?)?\s*(?:=|:|is|are)\s*[1-9]", r"(?:graph shape|read dendrogram).{0,50}(?:proves?|confirms?).{0,30}ancestry"),
    "NO_882579_ACCURACY_OR_PREVALENCE": QUESTION_FORBIDDEN["Q_CONCLUSION"][2:],
    "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE": (r"7 technical datasets?.{0,50}(?:are|=)\s*(?:exactly\s*)?7 (?:independent )?biological",),
    "NO_FEATURE_AS_MAIN": (r"b9aaa12.{0,50}(?:production[- ]ready|public main|released)",),
    "NO_LOCAL_AS_LIVE_PUBLISHED": (r"local source.{0,50}(?:is\s+(?!not\b)(?:already\s+)?|means?\s+)(?:live|published)",),
    "NO_BIP7_AS_BIP8": (r"bip7.{0,50}(?:proves?|means?).{0,30}bip8.{0,30}(?:pass|complete)",),
    "NO_PYTHON_ROLE_CONFLATION": (r"python(?:/html| and html)?.{0,30}only (?:presents?|renders?)", r"html.{0,50}recomputes? science"),
}

GENERIC_ANSWERS = {
    "a sufficiently complete bounded answer.",
    "the answer preserves the required boundary.",
    "see the documentation for details.",
    "the project is described in the package.",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def safe_path(package: Path, relative: object) -> Path:
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise ValueError(f"invalid package-relative path: {relative!r}")
    candidate = (package / relative).resolve()
    try:
        candidate.relative_to(package.resolve())
    except ValueError as error:
        raise ValueError(f"path escapes package: {relative}") from error
    if not candidate.is_file() or candidate.is_symlink():
        raise ValueError(f"path is not a regular package file: {relative}")
    return candidate


def run_git(repo: Path, *args: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True)


def discover_repo(package: Path) -> Path:
    probe = run_git(package, "rev-parse", "--show-toplevel")
    if probe.returncode != 0:
        raise ValueError("package is not inside a Git repository")
    return Path(probe.stdout.decode("utf-8").strip()).resolve()


def validate_schema(receipt: Any, schema_path: Path) -> list[str]:
    if Draft202012Validator is None or FormatChecker is None:
        return ["jsonschema dependency is unavailable; install requirements.txt"]
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        Draft202012Validator.check_schema(schema)
    except (OSError, UnicodeError, json.JSONDecodeError, SchemaError) as error:
        return [f"reader acceptance schema is invalid: {error}"]
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors: list[str] = []
    for issue in sorted(validator.iter_errors(receipt), key=lambda item: tuple(str(part) for part in item.absolute_path)):
        location = ".".join(str(part) for part in issue.absolute_path) or "$"
        errors.append(f"schema violation at {location}: {issue.message}")
    return errors


def missing_concepts(text: object, rubric: list[tuple[str, tuple[str, ...]]]) -> list[str]:
    if not isinstance(text, str):
        return [label for label, _ in rubric]
    return [label for label, alternatives in rubric if not any(re.search(pattern, text, FLAGS) for pattern in alternatives)]


def forbidden_matches(text: object, patterns: tuple[str, ...]) -> list[str]:
    if not isinstance(text, str):
        return []
    return [pattern for pattern in patterns if re.search(pattern, text, FLAGS)]


def is_generic(text: object) -> bool:
    if not isinstance(text, str):
        return True
    normalized = re.sub(r"\s+", " ", text.strip().lower())
    if normalized in GENERIC_ANSWERS:
        return True
    tokens = re.findall(r"[a-z0-9_./%'-]+|[\u4e00-\u9fff]", normalized)
    return len(tokens) < 10


def validate(
    receipt_path: Path,
    package: Path,
    *,
    require_head: bool = False,
) -> tuple[list[str], dict[str, Any]]:
    """Validate a receipt against current package bytes and a reachable source commit.

    ``require_head`` is retained for caller compatibility.  It no longer requires
    equality with HEAD: the tested commit may be any reachable ancestor, provided
    every manifest byte is identical in the commit and current package.
    """

    del require_head
    errors: list[str] = []
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return [f"receipt parse failed: {error}"], {"verdict": "FAIL", "errors": 1}

    schema_path = package / SCHEMA_RELATIVE
    errors.extend(validate_schema(receipt, schema_path))
    if not isinstance(receipt, dict):
        return errors or ["receipt must be an object"], {"verdict": "FAIL", "errors": len(errors) or 1}

    try:
        repo = discover_repo(package)
    except ValueError as error:
        errors.append(str(error))
        repo = None

    commit = receipt.get("tested_git_commit")
    commit_reachable = False
    commit_relation = "INVALID"
    if not isinstance(commit, str) or not COMMIT_RE.fullmatch(commit):
        errors.append("tested_git_commit must be a full 40-character lowercase SHA")
    elif repo is not None:
        exists = run_git(repo, "cat-file", "-e", f"{commit}^{{commit}}")
        if exists.returncode != 0:
            errors.append(f"tested_git_commit is not a reachable Git object: {commit}")
        else:
            ancestor = run_git(repo, "merge-base", "--is-ancestor", commit, "HEAD")
            if ancestor.returncode != 0:
                errors.append(f"tested_git_commit is not an ancestor of HEAD: {commit}")
            else:
                commit_reachable = True
                head = run_git(repo, "rev-parse", "HEAD").stdout.decode("utf-8").strip()
                commit_relation = "HEAD" if commit == head else "REACHABLE_ANCESTOR"

    manifest = receipt.get("package_source_manifest")
    manifest_by_path: dict[str, str] = {}
    manifest_current_matches = 0
    manifest_commit_matches = 0
    if not isinstance(manifest, list):
        errors.append("package_source_manifest must be a list")
        manifest = []
    for index, row in enumerate(manifest):
        if not isinstance(row, dict) or set(row) != {"path", "sha256"}:
            errors.append(f"invalid package_source_manifest row {index}")
            continue
        relative = row.get("path")
        expected = row.get("sha256")
        if not isinstance(relative, str):
            errors.append(f"manifest row {index} path is invalid")
            continue
        if relative in manifest_by_path:
            errors.append(f"duplicate source path: {relative}")
            continue
        if not isinstance(expected, str) or not SHA256_RE.fullmatch(expected):
            errors.append(f"manifest sha256 is invalid: {relative}")
            continue
        manifest_by_path[relative] = expected
        try:
            current_path = safe_path(package, relative)
        except ValueError as error:
            errors.append(str(error))
            continue
        current_hash = sha256(current_path)
        if current_hash != expected:
            errors.append(f"current package hash mismatch: {relative}")
        else:
            manifest_current_matches += 1
        if repo is not None and commit_reachable:
            try:
                repo_relative = current_path.relative_to(repo).as_posix()
            except ValueError:
                errors.append(f"manifest path is outside Git repository: {relative}")
                continue
            committed = run_git(repo, "show", f"{commit}:{repo_relative}")
            if committed.returncode != 0:
                errors.append(f"manifest path absent at tested_git_commit: {relative}")
            elif sha256_bytes(committed.stdout) != expected:
                errors.append(f"tested commit hash mismatch: {relative}")
            else:
                manifest_commit_matches += 1

    missing_required = sorted(REQUIRED_MANIFEST_PATHS - set(manifest_by_path))
    if missing_required:
        errors.append("package_source_manifest lacks required core paths: " + ", ".join(missing_required))

    questions = receipt.get("questions")
    if not isinstance(questions, list) or [row.get("question_id") for row in questions if isinstance(row, dict)] != QUESTIONS:
        errors.append("questions must contain the six ordered IDs exactly once")
        questions = []
    normalized_answers: list[str] = []
    semantic_concepts_checked = 0
    for row in questions:
        question_id = row.get("question_id")
        if row.get("status") != "PASS":
            errors.append(f"question not PASS: {question_id}")
        answer = row.get("answer")
        if is_generic(answer):
            errors.append(f"generic or insufficient answer: {question_id}")
        if isinstance(answer, str):
            normalized_answers.append(re.sub(r"\s+", " ", answer.strip().lower()))
        rubric = QUESTION_RUBRICS.get(question_id, [])
        semantic_concepts_checked += len(rubric)
        missing = missing_concepts(answer, rubric)
        if missing:
            errors.append(f"{question_id} missing required concepts: {', '.join(missing)}")
        forbidden = forbidden_matches(answer, QUESTION_FORBIDDEN.get(question_id, ()))
        if forbidden:
            errors.append(f"{question_id} contains forbidden claim(s): {len(forbidden)}")
        paths = row.get("evidence_paths")
        if not isinstance(paths, list):
            errors.append(f"question evidence_paths is not a list: {question_id}")
            continue
        unknown = sorted(set(paths) - set(manifest_by_path))
        if unknown:
            errors.append(f"{question_id} evidence path not in package_source_manifest: {', '.join(unknown)}")
        required_evidence = QUESTION_EVIDENCE.get(question_id, set())
        missing_evidence = sorted(required_evidence - set(paths))
        if missing_evidence:
            errors.append(f"{question_id} lacks required evidence paths: {', '.join(missing_evidence)}")
    if normalized_answers and len(normalized_answers) != len(set(normalized_answers)):
        errors.append("question answers must be independently written; duplicate generic answers found")

    checks = receipt.get("prohibited_promotions")
    if not isinstance(checks, list) or [row.get("check_id") for row in checks if isinstance(row, dict)] != PROMOTIONS:
        errors.append("prohibited_promotions must contain the eight ordered IDs exactly once")
        checks = []
    promotion_concepts_checked = 0
    for row in checks:
        check_id = row.get("check_id")
        if row.get("status") != "PASS":
            errors.append(f"prohibited promotion check not PASS: {check_id}")
        explanation = row.get("explanation")
        if is_generic(explanation):
            errors.append(f"generic or insufficient promotion explanation: {check_id}")
        rubric = PROMOTION_RUBRICS.get(check_id, [])
        promotion_concepts_checked += len(rubric)
        missing = missing_concepts(explanation, rubric)
        if missing:
            errors.append(f"{check_id} missing required explanation concepts: {', '.join(missing)}")
        forbidden = forbidden_matches(explanation, PROMOTION_FORBIDDEN.get(check_id, ()))
        if forbidden:
            errors.append(f"{check_id} contains forbidden promotion claim(s): {len(forbidden)}")
        paths = row.get("evidence_paths")
        if not isinstance(paths, list):
            errors.append(f"promotion evidence_paths is not a list: {check_id}")
            continue
        unknown = sorted(set(paths) - set(manifest_by_path))
        if unknown:
            errors.append(f"{check_id} evidence path not in package_source_manifest: {', '.join(unknown)}")
        required_evidence = PROMOTION_EVIDENCE.get(check_id, set())
        missing_evidence = sorted(required_evidence - set(paths))
        if missing_evidence:
            errors.append(f"{check_id} lacks required evidence paths: {', '.join(missing_evidence)}")

    if receipt.get("verdict") != "PASS":
        errors.append("receipt verdict is not PASS")
    summary = {
        "schema_version": receipt.get("schema_version"),
        "tested_git_commit": commit,
        "tested_commit_relation": commit_relation,
        "questions": len(questions),
        "question_semantic_concepts_checked": semantic_concepts_checked,
        "prohibited_promotions": len(checks),
        "promotion_semantic_concepts_checked": promotion_concepts_checked,
        "source_files": len(manifest_by_path),
        "required_source_files": len(REQUIRED_MANIFEST_PATHS),
        "current_hash_matches": manifest_current_matches,
        "tested_commit_hash_matches": manifest_commit_matches,
        "errors": len(errors),
        "verdict": "PASS" if not errors else "FAIL",
    }
    return errors, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--receipt", type=Path)
    parser.add_argument(
        "--allow-non-head",
        action="store_true",
        help="Deprecated compatibility flag; reachable ancestor commits are always accepted when all package hashes match.",
    )
    args = parser.parse_args()
    package = args.repo_root.resolve() / PACKAGE_RELATIVE
    receipt = args.receipt or package / "evidence/reader_acceptance_receipt.json"
    errors, summary = validate(receipt, package, require_head=False)
    print(json.dumps({**summary, "error_messages": errors}, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
