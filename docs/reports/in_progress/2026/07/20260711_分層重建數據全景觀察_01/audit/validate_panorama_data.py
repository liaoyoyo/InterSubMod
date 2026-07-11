#!/usr/bin/env python3
"""Validate panorama-report data conservation, provenance, and source disclosure."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote, urlparse


REPORT_STEM = "20260711_分層重建全景數據觀察_01"
SCRIPT_DIR = Path(__file__).resolve().parent
REPORT_DIR = SCRIPT_DIR.parent
DEFAULT_DATA = REPORT_DIR / f"{REPORT_STEM}.data.json"
DEFAULT_HTML = REPORT_DIR / f"{REPORT_STEM}.standalone.html"
DEFAULT_VALIDATION = REPORT_DIR / f"{REPORT_STEM}.validation.json"
DEFAULT_OUTPUT = SCRIPT_DIR / "after" / "data_quality.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fail-closed data and source-disclosure QA for the layered panorama HTML."
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--html", type=Path, default=DEFAULT_HTML)
    parser.add_argument("--validation", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    paths = {name: path.expanduser().resolve() for name, path in {
        "data": args.data,
        "html": args.html,
        "validation": args.validation,
        "output": args.output,
    }.items()}
    missing_inputs = [name for name in ("data", "html", "validation") if not paths[name].is_file()]
    if missing_inputs:
        raise SystemExit("missing input(s): " + ", ".join(missing_inputs))

    data = load_json(paths["data"])
    validation = load_json(paths["validation"])
    html_text = paths["html"].read_text(encoding="utf-8")
    checks = []

    def check(name: str, passed: bool, expected, actual) -> None:
        checks.append({"name": name, "passed": bool(passed), "expected": expected, "actual": actual})

    samples = data.get("samples", [])
    sample_names = [sample.get("sample") for sample in samples]
    biological_ids = [sample.get("biological_id") for sample in samples]
    check("dataset_rows", len(samples) == 7 and len(set(sample_names)) == 7, 7, {
        "rows": len(samples), "unique_names": len(set(sample_names))
    })
    check("biological_samples", len(set(biological_ids)) == 6, 6, len(set(biological_ids)))
    duplicate_biology = {key: value for key, value in Counter(biological_ids).items() if value > 1}
    check("technical_replicate_identity", duplicate_biology == {"HCC1395": 2}, {"HCC1395": 2}, duplicate_biology)

    for sample in samples:
        name = sample["sample"]
        funnel = sample["site_funnel"]
        regions = sample["regions"]
        units = sample["units"]
        topology = sample["region_candidate_topology"]
        check(
            f"{name}:site_funnel_conservation",
            funnel["universe"] == sum(funnel[key] for key in (
                "out_of_scope", "singleton", "cap_excluded", "read_unsupported", "retained"
            )),
            funnel["universe"],
            sum(funnel[key] for key in (
                "out_of_scope", "singleton", "cap_excluded", "read_unsupported", "retained"
            )),
        )
        check(
            f"{name}:region_bridge",
            regions["W_all_pre"] == regions["W_k1"] + regions["W_k_gt1"]
            and regions["W_k_gt1"] == regions["W_pre"]
            and regions["W_pre"] >= regions["W_ret"] >= regions["W_tree"] >= regions["W_primary"],
            "W_all_pre=W_k1+W_k_gt1; W_pre>=W_ret>=W_tree>=W_primary",
            {key: regions[key] for key in (
                "W_all_pre", "W_k1", "W_k_gt1", "W_pre", "W_ret", "W_tree", "W_primary"
            )},
        )
        check(
            f"{name}:hp_partition",
            regions["W_tree"] == regions["single_primary_HP"] + regions["multi_HP"] + regions["no_primary"],
            regions["W_tree"],
            regions["single_primary_HP"] + regions["multi_HP"] + regions["no_primary"],
        )
        check(
            f"{name}:l1_partition",
            units["primary"] == sum(units[key] for key in ("determined", "ambiguous", "solver_capped", "recurrence")),
            units["primary"],
            sum(units[key] for key in ("determined", "ambiguous", "solver_capped", "recurrence")),
        )
        check(
            f"{name}:topology_partition",
            regions["W_primary"] == sum(topology["topology_classes"].values()),
            regions["W_primary"],
            sum(topology["topology_classes"].values()),
        )
        check(
            f"{name}:hidden_partition",
            regions["W_primary"] == sum(topology["hidden_classes"].values()),
            regions["W_primary"],
            sum(topology["hidden_classes"].values()),
        )
        check(
            f"{name}:hp_h3_partition",
            regions["W_tree"] == sum(topology["hp_h3"].values()),
            regions["W_tree"],
            sum(topology["hp_h3"].values()),
        )
        check(
            f"{name}:k_partition",
            regions["W_ret"] == sum(sample["k_distribution"]["counts"].values()),
            regions["W_ret"],
            sum(sample["k_distribution"]["counts"].values()),
        )

    aggregate = {
        "site_universe": sum(sample["site_funnel"]["universe"] for sample in samples),
        "retained_sites": sum(sample["site_funnel"]["retained"] for sample in samples),
        "W_tree": sum(sample["regions"]["W_tree"] for sample in samples),
        "W_primary": sum(sample["regions"]["W_primary"] for sample in samples),
        "primary_units": sum(sample["units"]["primary"] for sample in samples),
        "complete_regions": sum(sample["region_candidate_topology"]["complete_regions"] for sample in samples),
        "incomplete_regions": sum(sample["region_candidate_topology"]["incomplete_regions"] for sample in samples),
    }
    check(
        "aggregate_complete_incomplete",
        aggregate["W_primary"] == aggregate["complete_regions"] + aggregate["incomplete_regions"],
        aggregate["W_primary"],
        aggregate["complete_regions"] + aggregate["incomplete_regions"],
    )

    production = data["clean_production"]
    latest_pass = sum(1 for item in production["latest_by_sample"].values() if item.get("status") == "PASS")
    check("production_role", "raw_all" in production["root"], "normalized raw-all root", production["root"])
    check("production_pass_count", production["n_pass"] == latest_pass, latest_pass, production["n_pass"])
    check(
        "production_completion_gate",
        production["state"] != "complete" or (
            production["n_pass"] == production["expected_datasets"]
            and production["aggregate_complete"]
            and production["sample_status_complete"]
        ),
        "complete only after 7/7 + aggregate + sample-status gates",
        {key: production[key] for key in (
            "state", "n_pass", "expected_datasets", "aggregate_complete", "sample_status_complete"
        )},
    )

    check("artifact_validation", validation.get("status") == "pass" and not validation.get("errors"), "pass", {
        "status": validation.get("status"), "errors": validation.get("errors")
    })
    check("html_hash", validation.get("output_html_sha256") == sha256(paths["html"]),
          validation.get("output_html_sha256"), sha256(paths["html"]))
    check("data_hash", validation.get("output_data_sha256") == sha256(paths["data"]),
          validation.get("output_data_sha256"), sha256(paths["data"]))

    check("inline_source_tooltips", 'class="source-tooltip inline-source"' not in html_text, 0,
          html_text.count('class="source-tooltip inline-source"'))
    check("sourced_value_annotations", html_text.count('class="sourced-value"') > 0, ">0",
          html_text.count('class="sourced-value"'))
    raw_links = re.findall(r'<a\b[^>]*href="([^"]+\.json)"[^>]*>(.*?)</a>', html_text, re.S | re.I)
    unique_hrefs = sorted(set(href for href, _ in raw_links))
    visible_json_link_labels = [
        re.sub(r"<[^>]+>", " ", label).strip()
        for href, label in raw_links
        if ".json" in re.sub(r"<[^>]+>", " ", label).lower()
    ]
    missing_links = []
    for href in unique_hrefs:
        if href.startswith("file:"):
            parsed = urlparse(href)
            target = Path(unquote(parsed.path))
        else:
            target = (paths["html"].parent / href).resolve()
        if not target.is_file():
            missing_links.append({"href": href, "resolved": str(target)})
    check("hidden_json_link_labels", not visible_json_link_labels, [], visible_json_link_labels)
    check("companion_json_links", len(unique_hrefs) >= 5 and not missing_links,
          {"minimum_unique": 5, "missing": 0}, {"unique": unique_hrefs, "missing": missing_links})

    failed = [item for item in checks if not item["passed"]]
    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "pass" if not failed else "fail",
        "inputs": {name: str(path) for name, path in paths.items() if name != "output"},
        "aggregate": aggregate,
        "production": {
            "root": production["root"],
            "state": production["state"],
            "pass": production["n_pass"],
            "expected": production["expected_datasets"],
            "aggregate_complete": production["aggregate_complete"],
        },
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "failed_checks": failed,
        "checks": checks,
    }
    paths["output"].parent.mkdir(parents=True, exist_ok=True)
    paths["output"].write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "checks_total", "checks_passed", "checks_failed", "aggregate", "production"
    )}, ensure_ascii=False, indent=2))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
