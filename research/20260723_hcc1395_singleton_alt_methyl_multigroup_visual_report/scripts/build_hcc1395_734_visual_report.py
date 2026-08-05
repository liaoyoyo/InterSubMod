#!/usr/bin/env python3
"""Build the source-backed HCC1395 734-site methyl multigroup HTML artifact.

The script derives a focused report from the already validated singleton report,
recomputes every HCC1395 denominator from the current source-authority-v7 audit,
and adds reproducible average-linkage (UPGMA) dendrograms for two M2-PASS examples.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import math
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORK = REPO / "research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report"
RESULTS = WORK / "results"
FIGURES = WORK / "figures"
BASE_ARTIFACT = (
    REPO
    / "research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json"
)
AUTHORITY = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "positional_singleton_methyl_multigroup_audit_v2_source_authority_v7"
)
SITE_AUDIT = AUTHORITY / "positional_singleton_site_audit.tsv.gz"
AUDIT_SUMMARY = AUTHORITY / "positional_singleton_audit_summary.json"
M2_PASS_CASES = AUTHORITY / "positional_singleton_m2_pass_cases.tsv"
ARTIFACT_OUT = WORK / "artifact.json"
SITE_TABLE_OUT = RESULTS / "20260723_HCC1395_M1穩定甲基多群734位點_01.tsv.gz"
BUILD_SUMMARY_OUT = RESULTS / "20260723_HCC1395_734位點視覺報告建置摘要_01.json"
RECOUNT_AUDIT = RESULTS / "20260723_HCC1395_734位點重算摘要_01.json"
STAGING_DB = RESULTS / "artifact_staging.sqlite"

EXPECTED = {
    "singleton_sites": 8279,
    "m1_evaluable": 8074,
    "m1_flagged": 734,
    "m2_pass": 2,
    "m2_not_evaluable": 732,
    "m2_not_run": 7545,
    "k_distribution": {2: 651, 3: 73, 4: 9, 5: 1},
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("true")


def parse_cluster_sizes(value: object) -> list[int]:
    parsed = json.loads(str(value))
    sizes = sorted((int(v) for v in parsed.values()), reverse=True)
    if not sizes or any(v < 1 for v in sizes):
        raise ValueError(f"Invalid cluster_sizes: {value}")
    return sizes


def html_figure(path: Path, alt: str, caption: str, source_id: str) -> dict:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    body = (
        '<figure style="margin:0;font-family:system-ui,sans-serif;color:#171411">'
        f'<img alt="{alt}" src="data:image/png;base64,{encoded}" '
        'style="display:block;width:100%;height:auto;border:0">'
        f'<figcaption style="margin-top:10px;font-size:13px;line-height:1.55;color:#565d66">{caption}</figcaption>'
        "</figure>"
    )
    return {"type": "html", "body": body, "sourceId": source_id}


def build_dendrogram(
    distance_cells: list[dict],
    locus: str,
    output_path: Path,
) -> dict:
    row_order = list(dict.fromkeys(row["row_read"] for row in distance_cells))
    column_order = list(dict.fromkeys(row["column_read"] for row in distance_cells))
    if row_order != column_order:
        raise AssertionError(f"{locus}: row and column label order differs")
    index = {label: idx for idx, label in enumerate(row_order)}
    matrix = np.full((len(row_order), len(row_order)), np.nan, dtype=float)
    groups: dict[str, str] = {}
    for cell in distance_cells:
        i = index[cell["row_read"]]
        j = index[cell["column_read"]]
        matrix[i, j] = float(cell["distance"])
        groups[cell["row_read"]] = cell["row_group"]
        groups[cell["column_read"]] = cell["column_group"]
    if not np.isfinite(matrix).all():
        raise AssertionError(f"{locus}: distance matrix has non-finite cells")
    if not np.allclose(matrix, matrix.T, atol=1e-6):
        raise AssertionError(f"{locus}: distance matrix is not symmetric")
    if not np.allclose(np.diag(matrix), 0.0, atol=1e-8):
        raise AssertionError(f"{locus}: distance matrix diagonal is not zero")

    condensed = squareform(matrix, checks=True)
    z_matrix = linkage(condensed, method="average", optimal_ordering=True)
    fig, ax = plt.subplots(figsize=(12.5, 10.5), constrained_layout=True)
    result = dendrogram(
        z_matrix,
        labels=row_order,
        orientation="right",
        leaf_font_size=8.5,
        color_threshold=0,
        above_threshold_color="#3F454D",
        link_color_func=lambda _: "#3F454D",
        ax=ax,
    )
    colors = {"Group A": "#2563EB", "Group B": "#D97706"}
    for label in ax.get_yticklabels():
        label.set_color(colors.get(groups[label.get_text()], "#3F454D"))
        label.set_fontweight("normal")
    ax.set_xlabel("Average-linkage merge height (Bernoulli methylation distance)", fontsize=10)
    ax.set_ylabel("Deterministic representative focal-ALT reads", fontsize=10)
    ax.set_title(
        f"{locus} representative focal-ALT read UPGMA",
        loc="left",
        fontsize=15,
        weight="bold",
        pad=34,
    )
    ax.text(
        0,
        1.015,
        "Neutral branches; blue/orange leaf labels are methyl Group A/B aliases",
        transform=ax.transAxes,
        fontsize=10,
        color="#565D66",
    )
    ax.grid(axis="x", color="#D7DCE2", linewidth=0.7, alpha=0.75)
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.savefig(output_path, dpi=180, facecolor="white")
    plt.close(fig)
    return {
        "locus": locus,
        "representative_reads": len(row_order),
        "leaf_order": result["ivl"],
        "max_merge_height": round(float(z_matrix[:, 2].max()), 6),
        "matrix_symmetric": True,
        "matrix_diagonal_zero": True,
        "figure": str(output_path.relative_to(REPO)),
        "sha256": sha256_file(output_path),
    }


def chart_source(source_id: str, description: str, metric_definitions: list[str]) -> dict:
    return {
        "id": source_id,
        "label": "HCC1395 source-authority-v7 focused recomputation",
        "path": "positional_singleton_site_audit.tsv.gz",
        "query": {
            "description": description,
            "engine": "Python 3.11 / pandas 2.3.3",
            "language": "Python",
            "sql": (
                "import pandas as pd\n"
                "audit = pd.read_csv('positional_singleton_site_audit.tsv.gz', sep='\\t')\n"
                "hcc = audit[(audit.dataset == 'HCC1395') & "
                "(audit.ssnv_branch == 'positional_singleton')].copy()\n"
                "flagged = hcc[(hcc.m1_status == 'FLAGGED') & "
                "hcc.stable_null_multigroup.astype(str).str.lower().eq('true')]\n"
            ),
            "tables_used": ["positional_singleton_site_audit.tsv.gz"],
            "filters": [
                "dataset=HCC1395",
                "ssnv_branch=positional_singleton",
                "M1 stable rows require m1_status=FLAGGED and stable_null_multigroup=true",
            ],
            "metric_definitions": metric_definitions,
        },
        "sha256": sha256_file(SITE_AUDIT),
    }


def create_staging_sources(datasets: dict[str, list[dict]]) -> dict[str, dict]:
    """Persist reviewed rows, execute the exact SELECTs, and return widget SQL sources."""
    table_specs = {
        "summary_metrics_sql": (
            "artifact_hcc1395_summary_metrics",
            "hcc1395_summary_metrics",
            [
                "singleton_sites",
                "m1_evaluable",
                "m1_evaluable_share",
                "m1_stable",
                "m1_all_share",
                "m1_evaluable_share_stable",
                "m2_pass",
            ],
            "HCC1395 headline denominator metrics",
            "singleton_audit",
        ),
        "funnel_sql": (
            "artifact_hcc1395_734_funnel",
            "hcc1395_734_funnel",
            ["stage", "count", "population", "denominator_note", "share_of_singletons"],
            "HCC1395 singleton-to-M1/M2 ordered gate counts",
            "singleton_audit",
        ),
        "k_distribution_sql": (
            "artifact_hcc1395_k_distribution",
            "hcc1395_k_distribution",
            ["k", "k_label", "count", "share", "cumulative_share", "denominator", "rank"],
            "HCC1395 stable methyl-group K distribution",
            "singleton_audit",
        ),
        "support_distribution_sql": (
            "artifact_hcc1395_support_distribution",
            "hcc1395_support_distribution",
            [
                "support_bin",
                "support_bin_order",
                "bin_lower",
                "bin_upper",
                "k",
                "k_label",
                "site_count",
                "bin_total",
                "share_of_734",
            ],
            "Binned grouped focal-ALT read support by methyl-group K",
            "singleton_audit",
        ),
        "m1_sites_sql": (
            "artifact_hcc1395_m1_sites",
            "hcc1395_m1_sites",
            [
                "locus",
                "truth",
                "k",
                "k_label",
                "alt_reads_before_peel",
                "core_reads",
                "grouped_reads",
                "other_after_peel_reads",
                "peeled_reads",
                "cluster_profile",
                "min_cluster_n",
                "min_group_fraction",
                "ari_min",
                "unstable",
                "null_calibrated_multigroup",
                "m2_status",
            ],
            "All 734 HCC1395 M1 stable methyl-multigroup audit rows",
            "singleton_audit",
        ),
        "stability_summary_sql": (
            "artifact_hcc1395_stability_summary",
            "hcc1395_stability_summary",
            [
                "sites_checked",
                "accepted_splits",
                "splits_above_null95",
                "split_pass_share",
                "median_between_within",
                "min_between_within",
                "median_observed_over_null95",
                "min_observed_over_null95",
                "assignment_ari_min",
                "modal_k_seed_fraction_min",
            ],
            "Assignment and 40-column-shuffle null stability summary",
            "recount_audit",
        ),
    }

    connection = sqlite3.connect(STAGING_DB)
    try:
        for _, (table_name, dataset_name, columns, _, _) in table_specs.items():
            frame = pd.DataFrame(datasets[dataset_name])[columns].copy()
            frame.insert(0, "_row_order", range(len(frame)))
            frame.to_sql(table_name, connection, if_exists="replace", index=False)
        connection.commit()
        for _, (table_name, dataset_name, columns, _, _) in table_specs.items():
            quoted = ", ".join(f'"{column}"' for column in columns)
            sql = f'SELECT {quoted} FROM "{table_name}" ORDER BY "_row_order" ASC'
            returned = connection.execute(sql).fetchall()
            if len(returned) != len(datasets[dataset_name]):
                raise AssertionError(f"{table_name}: SQLite row-count mismatch")
    finally:
        connection.close()

    db_sha = sha256_file(STAGING_DB)
    sources: dict[str, dict] = {}
    for source_id, (table_name, dataset_name, columns, description, upstream_source_id) in table_specs.items():
        quoted = ", ".join(f'"{column}"' for column in columns)
        sql = f'SELECT {quoted} FROM "{table_name}" ORDER BY "_row_order" ASC'
        sources[source_id] = {
            "id": source_id,
            "label": f"{description} via executed SQLite staging query",
            "path": "artifact_staging.sqlite",
            "query": {
                "description": (
                    "Actual query executed against the persisted artifact staging database; row count was exact-compared "
                    "with the reviewed snapshot dataset."
                ),
                "engine": f"SQLite {sqlite3.sqlite_version}",
                "language": "SQL",
                "sql": sql,
                "tables_used": [table_name, dataset_name],
                "filters": ["dataset=HCC1395", "ssnv_branch=positional_singleton"],
                "metric_definitions": [description],
            },
            "sha256": db_sha,
            "upstreamSourceId": upstream_source_id,
            "upstreamSha256": sha256_file(SITE_AUDIT if upstream_source_id == "singleton_audit" else RECOUNT_AUDIT),
        }
    return sources


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    with BASE_ARTIFACT.open() as handle:
        artifact = json.load(handle)

    audit = pd.read_csv(SITE_AUDIT, sep="\t", compression="gzip", low_memory=False)
    hcc = audit[(audit["dataset"] == "HCC1395") & (audit["ssnv_branch"] == "positional_singleton")].copy()
    evaluable = hcc[bool_series(hcc["m1_evaluable"])].copy()
    flagged = hcc[(hcc["m1_status"] == "FLAGGED") & bool_series(hcc["stable_null_multigroup"])].copy()

    observed = {
        "singleton_sites": int(len(hcc)),
        "m1_evaluable": int(len(evaluable)),
        "m1_flagged": int(len(flagged)),
        "m2_pass": int((hcc["m2_status"] == "PASS").sum()),
        "m2_not_evaluable": int((hcc["m2_status"] == "NOT_EVALUABLE").sum()),
        "m2_not_run": int((hcc["m2_status"] == "NOT_RUN").sum()),
    }
    for key, expected in EXPECTED.items():
        if key != "k_distribution" and observed[key] != expected:
            raise AssertionError(f"{key}: observed {observed[key]}, expected {expected}")

    flagged["k"] = pd.to_numeric(flagged["coarse_ng"], errors="raise").astype(int)
    flagged["ari_min"] = pd.to_numeric(flagged["modal_assignment_ari_min"], errors="raise")
    flagged["core_reads"] = pd.to_numeric(flagged["n_alt_after_peel"], errors="raise").astype(int)
    flagged["cluster_size_values"] = flagged["cluster_sizes"].map(parse_cluster_sizes)
    flagged["cluster_total"] = flagged["cluster_size_values"].map(sum)
    flagged["min_cluster_n"] = flagged["cluster_size_values"].map(min)
    flagged["min_group_fraction"] = flagged["min_cluster_n"] / flagged["cluster_total"]
    flagged["locus"] = (
        flagged["chrom"].astype(str)
        + ":"
        + flagged["pos"].astype(int).map(lambda value: f"{value:,}")
        + " "
        + flagged["ref"].astype(str)
        + ">"
        + flagged["alt"].astype(str)
    )
    flagged["cluster_profile"] = flagged["cluster_size_values"].map(lambda values: " / ".join(map(str, values)))
    flagged["k_label"] = flagged["k"].map(lambda value: f"K={value}")

    if bool_series(flagged["unstable"]).any():
        raise AssertionError("A flagged HCC1395 M1 row has unstable=true")
    if float(flagged["ari_min"].min()) < 0.8:
        raise AssertionError("A flagged HCC1395 M1 row has modal_assignment_ari_min < 0.8")
    if int(flagged["k"].min()) < 2:
        raise AssertionError("A flagged HCC1395 M1 row has K < 2")
    if int(flagged["min_cluster_n"].min()) < 3:
        raise AssertionError("A flagged HCC1395 M1 row has group size < 3")
    flagged["other_after_peel_reads"] = flagged["core_reads"] - flagged["cluster_total"]
    if (flagged["other_after_peel_reads"] < 0).any():
        raise AssertionError("A row has grouped cluster_total larger than n_alt_after_peel")

    k_counts = Counter(flagged["k"].tolist())
    if dict(sorted(k_counts.items())) != EXPECTED["k_distribution"]:
        raise AssertionError(
            f"K distribution: observed {dict(sorted(k_counts.items()))}, expected {EXPECTED['k_distribution']}"
        )

    site_rows = []
    for _, row in flagged.sort_values(["chrom", "pos", "ref", "alt"]).iterrows():
        site_rows.append(
            {
                "locus": row["locus"],
                "truth": str(row["truth_label"]),
                "k": int(row["k"]),
                "k_label": row["k_label"],
                "alt_reads_before_peel": int(row["n_alt_matrix"]),
                "core_reads": int(row["core_reads"]),
                "grouped_reads": int(row["cluster_total"]),
                "other_after_peel_reads": int(row["other_after_peel_reads"]),
                "peeled_reads": int(row["n_alt_peeled"]),
                "cluster_profile": row["cluster_profile"],
                "min_cluster_n": int(row["min_cluster_n"]),
                "min_group_fraction": round(float(row["min_group_fraction"]), 6),
                "ari_min": round(float(row["ari_min"]), 6),
                "unstable": False,
                "null_calibrated_multigroup": True,
                "m2_status": str(row["m2_status"]),
            }
        )

    pd.DataFrame(site_rows).to_csv(
        SITE_TABLE_OUT,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    singleton_rate = observed["m1_flagged"] / observed["singleton_sites"]
    evaluable_rate = observed["m1_flagged"] / observed["m1_evaluable"]
    k_rows = []
    cumulative = 0
    for rank, (k, count) in enumerate(sorted(k_counts.items()), start=1):
        cumulative += count
        k_rows.append(
            {
                "k": k,
                "k_label": f"K={k}",
                "count": count,
                "share": count / observed["m1_flagged"],
                "cumulative_share": cumulative / observed["m1_flagged"],
                "denominator": observed["m1_flagged"],
                "rank": rank,
            }
        )

    stability_rows = [
        {
            "locus": row["locus"],
            "k": int(row["k"]),
            "k_label": row["k_label"],
            "core_reads": int(row["core_reads"]),
            "grouped_reads": int(row["cluster_total"]),
            "ari_min": round(float(row["ari_min"]), 6),
            "min_cluster_n": int(row["min_cluster_n"]),
            "min_group_fraction": round(float(row["min_group_fraction"]), 6),
            "truth": str(row["truth_label"]),
            "m2_status": str(row["m2_status"]),
        }
        for _, row in flagged.sort_values(["k", "core_reads", "locus"]).iterrows()
    ]

    support_bin_edges = [5, 9, 19, 29, 39, 49, 59, 79, 99, math.inf]
    support_bin_labels = ["6–9", "10–19", "20–29", "30–39", "40–49", "50–59", "60–79", "80–99", "100+"]

    def support_bin(value: int) -> tuple[int, str, int, int | None]:
        lower = 6
        for order, (upper, label) in enumerate(zip(support_bin_edges[1:], support_bin_labels), start=1):
            if value <= upper:
                return order, label, lower, None if math.isinf(upper) else int(upper)
            lower = int(upper) + 1
        raise AssertionError(f"Could not bin grouped read support: {value}")

    support_counter: Counter[tuple[int, str, int, int | None, int]] = Counter()
    bin_totals: Counter[int] = Counter()
    for row in stability_rows:
        order, label, lower, upper = support_bin(int(row["grouped_reads"]))
        support_counter[(order, label, lower, upper, int(row["k"]))] += 1
        bin_totals[order] += 1
    support_distribution = []
    for (order, label, lower, upper, k), count in sorted(support_counter.items()):
        support_distribution.append(
            {
                "support_bin": label,
                "support_bin_order": order,
                "bin_lower": lower,
                "bin_upper": upper,
                "k": k,
                "k_label": f"K={k}",
                "site_count": count,
                "bin_total": bin_totals[order],
                "share_of_734": count / observed["m1_flagged"],
            }
        )
    if sum(row["site_count"] for row in support_distribution) != observed["m1_flagged"]:
        raise AssertionError("Grouped-read support bins do not sum to 734")

    metrics = [
        {
            "singleton_sites": observed["singleton_sites"],
            "m1_evaluable": observed["m1_evaluable"],
            "m1_evaluable_share": observed["m1_evaluable"] / observed["singleton_sites"],
            "m1_stable": observed["m1_flagged"],
            "m1_all_share": singleton_rate,
            "m1_evaluable_share_stable": evaluable_rate,
            "m2_pass": observed["m2_pass"],
        }
    ]

    with RECOUNT_AUDIT.open() as handle:
        recount = json.load(handle)
    null_support = recount["m1_stability_and_null_support"]
    rates_with_ci = recount["rates_with_wilson_95ci"]
    stability_summary = [
        {
            "sites_checked": int(null_support["sites_checked"]),
            "accepted_splits": int(null_support["coarse_passed_split_count"]),
            "splits_above_null95": int(null_support["coarse_passed_split_count"]),
            "split_pass_share": 1.0,
            "median_between_within": float(null_support["split_between_within_ratio"]["median"]),
            "min_between_within": float(null_support["split_between_within_ratio"]["min"]),
            "median_observed_over_null95": float(null_support["split_observed_over_null95_margin"]["median"]),
            "min_observed_over_null95": float(null_support["split_observed_over_null95_margin"]["min"]),
            "assignment_ari_min": float(null_support["modal_assignment_ari_min_distribution"]["min"]),
            "modal_k_seed_fraction_min": float(null_support["modal_k_fraction_quantiles"]["min"]),
        }
    ]

    focused_funnel = [
        {
            "stage": "Positional singleton",
            "count": observed["singleton_sites"],
            "population": "HCC1395",
            "denominator_note": "component_size=1",
            "share_of_singletons": 1.0,
        },
        {
            "stage": "M1 evaluable",
            "count": observed["m1_evaluable"],
            "population": "HCC1395",
            "denominator_note": "M1 input requirements met",
            "share_of_singletons": observed["m1_evaluable"] / observed["singleton_sites"],
        },
        {
            "stage": "M1 stable multigroup",
            "count": observed["m1_flagged"],
            "population": "HCC1395",
            "denominator_note": "M1 FLAGGED + stable_null_multigroup",
            "share_of_singletons": singleton_rate,
        },
        {
            "stage": "M2 determinate / PASS",
            "count": observed["m2_pass"],
            "population": "HCC1395",
            "denominator_note": "2 determinate, both PASS",
            "share_of_singletons": observed["m2_pass"] / observed["singleton_sites"],
        },
    ]

    focused_datasets = {
        "hcc1395_summary_metrics": metrics,
        "hcc1395_734_funnel": focused_funnel,
        "hcc1395_k_distribution": k_rows,
        "hcc1395_support_distribution": support_distribution,
        "hcc1395_stability_summary": stability_summary,
        "hcc1395_m1_sites": site_rows,
    }
    staging_sources = create_staging_sources(focused_datasets)

    base_datasets = artifact["snapshot"]["datasets"]
    dendrogram_specs = [
        (
            "chr14:86,272,476 A>T",
            "chr14_86272476_A_T_distance_cells",
            FIGURES / "01_chr14_86272476_A_T_representative_UPGMA.png",
            "chr14_raw",
            "chr14_upgma",
        ),
        (
            "chr22:47,466,517 A>G",
            "chr22_47466517_A_G_distance_cells",
            FIGURES / "02_chr22_47466517_A_G_representative_UPGMA.png",
            "chr22_raw",
            "chr22_upgma",
        ),
    ]
    dendrogram_results: list[dict] = []
    dendrogram_blocks: dict[str, dict] = {}
    for locus, dataset, output_path, source_id, block_id in dendrogram_specs:
        dendrogram_results.append(build_dendrogram(base_datasets[dataset], locus, output_path))
        block = html_figure(
            output_path,
            (
                f"{locus} representative focal-ALT read average-linkage UPGMA dendrogram; "
                "neutral branches and colored Group A/B leaf labels"
            ),
            (
                "以每群最多 15 個 deterministic medoid-nearest representative reads 建立；枝高是 Bernoulli "
                "甲基距離的 average-linkage 合併高度。這是 read methylation-distance dendrogram，"
                "不是突變祖先—後代演化樹，也不把 Group A/B 當作 HP。"
            ),
            source_id,
        )
        block["id"] = block_id
        dendrogram_blocks[block_id] = block

    now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    title = "HCC1395：734 個穩定 focal-ALT 甲基多群位點的統計與視覺證據"
    manifest = artifact["manifest"]
    manifest["title"] = title
    manifest["description"] = (
        "HCC1395 全 8,279 個 positional-singleton sSNV 的 M1 穩定甲基多群重算、734 位點逐筆表、"
        "群數與穩定度分布，以及兩個 M2-PASS 代表位點的 UPGMA/read-distance/read×CpG 視覺。"
    )
    manifest["generatedAt"] = now

    singleton_source = chart_source(
        "singleton_audit",
        "Recompute HCC1395 positional-singleton denominators, M1 stable flags, K distribution, and stability fields.",
        [
            "M1 stable multigroup = m1_status=FLAGGED AND stable_null_multigroup=true.",
            "Operational M1 rate uses 734/8,279; evaluable-conditional rate uses 734/8,074.",
            "K is coarse_ng, the number of stable focal-ALT methyl read groups, not clone count.",
            "modal_assignment_ari_min >=0.8 and unstable=false are required stability gates.",
        ],
    )
    for source_list in (manifest["sources"], artifact["sources"]):
        source_list[:] = [singleton_source if source["id"] == "singleton_audit" else source for source in source_list]
        if not any(source["id"] == "recount_audit" for source in source_list):
            source_list.append(
                {
                    "id": "recount_audit",
                    "label": "HCC1395 734-site assignment/null stability recount",
                    "path": "20260723_HCC1395_734位點重算摘要_01.json",
                    "query": {
                        "description": (
                            "Exact key-set reconciliation of 734 authority-v7 sites against 102,842 stable-assignment "
                            "objects and census of all 828 accepted recursive splits."
                        ),
                        "engine": "Python 3.11 streaming JSONL + pandas",
                        "language": "Python",
                        "sql": (
                            "import json\n"
                            "with open('20260723_HCC1395_734位點重算摘要_01.json') as fh:\n"
                            "    recount = json.load(fh)\n"
                            "stability = recount['m1_stability_and_null_support']\n"
                        ),
                        "tables_used": [
                            "positional_singleton_site_audit.tsv.gz",
                            "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
                        ],
                        "filters": ["dataset=HCC1395", "734 authority-v7 M1 stable site keys"],
                        "metric_definitions": [
                            "Accepted split requires observed between/within ratio above the interpolated 95th percentile of 40 column-shuffle nulls.",
                            "Empirical rank p is descriptive and is not relabeled as a classical p<0.05 test.",
                        ],
                    },
                    "sha256": sha256_file(RECOUNT_AUDIT),
                }
            )
        for source_id, source in staging_sources.items():
            source_list[:] = [existing for existing in source_list if existing["id"] != source_id]
            source_list.append(source)

    manifest["cards"] = [
        {
            "id": "singleton_sites_card",
            "description": "HCC1395 positional-singleton sSNV sites in the frozen whole-autosome audit.",
            "dataset": "hcc1395_summary_metrics",
            "sourceId": "summary_metrics_sql",
            "metrics": [{"label": "單一 sSNV 區域", "field": "singleton_sites", "format": "number"}],
        },
        {
            "id": "m1_evaluable_card",
            "description": "Sites with sufficient inputs to run the M1 methyl multigroup screen.",
            "dataset": "hcc1395_summary_metrics",
            "sourceId": "summary_metrics_sql",
            "metrics": [
                {"label": "M1 可評估", "field": "m1_evaluable", "format": "number"},
                {"label": "佔 8,279", "field": "m1_evaluable_share", "format": "percent"},
            ],
        },
        {
            "id": "m1_stable_card",
            "description": "Focal-ALT reads retain a null-calibrated, seed-stable methyl multigroup partition.",
            "dataset": "hcc1395_summary_metrics",
            "sourceId": "summary_metrics_sql",
            "metrics": [
                {"label": "M1 穩定多群", "field": "m1_stable", "format": "number"},
                {"label": "佔全部", "field": "m1_all_share", "format": "percent"},
                {"label": "佔可評估", "field": "m1_evaluable_share_stable", "format": "percent"},
            ],
        },
        {
            "id": "m2_pass_card",
            "description": "M2 determinate examples that remain compatible after the measured-axis checks.",
            "dataset": "hcc1395_summary_metrics",
            "sourceId": "summary_metrics_sql",
            "metrics": [{"label": "M2 PASS 代表位點", "field": "m2_pass", "format": "number"}],
        },
        {
            "id": "accepted_splits_card",
            "description": "All recursive splits retained across the 734 M1 sites and tested against column-shuffle null95.",
            "dataset": "hcc1395_stability_summary",
            "sourceId": "stability_summary_sql",
            "metrics": [
                {"label": "通過 null95 的 splits", "field": "splits_above_null95", "format": "number"},
                {"label": "通過率", "field": "split_pass_share", "format": "percent"},
            ],
        },
        {
            "id": "split_ratio_card",
            "description": "Between-group over within-group distance ratio across all accepted recursive splits.",
            "dataset": "hcc1395_stability_summary",
            "sourceId": "stability_summary_sql",
            "metrics": [
                {"label": "Between/within 中位數", "field": "median_between_within", "format": "number"},
                {"label": "最小值", "field": "min_between_within", "format": "number"},
            ],
        },
        {
            "id": "null_margin_card",
            "description": "Observed split ratio divided by its site-specific interpolated column-shuffle null95 threshold.",
            "dataset": "hcc1395_stability_summary",
            "sourceId": "stability_summary_sql",
            "metrics": [
                {"label": "Observed/null95 中位數", "field": "median_observed_over_null95", "format": "number"},
                {"label": "最小值", "field": "min_observed_over_null95", "format": "number"},
            ],
        },
        {
            "id": "seed_stability_card",
            "description": "Worst retained assignment agreement and modal-K recurrence across ten seeds.",
            "dataset": "hcc1395_stability_summary",
            "sourceId": "stability_summary_sql",
            "metrics": [
                {"label": "Assignment ARI 最小值", "field": "assignment_ari_min", "format": "number"},
                {"label": "Modal K 最少 seeds", "field": "modal_k_seed_fraction_min", "format": "percent"},
            ],
        },
    ]

    old_charts = {chart["id"]: chart for chart in manifest["charts"]}
    charts = [
        {
            "id": "hcc1395_734_funnel_chart",
            "title": "HCC1395 positional-singleton M1/M2 funnel",
            "subtitle": "8,279 → 8,074 → 734；M2 只有 2 個 determinate 且皆 PASS，732 個 M1 flags 為 NOT_EVALUABLE。",
            "type": "funnel",
            "intent": "funnel",
            "dataset": "hcc1395_734_funnel",
            "sourceId": "funnel_sql",
            "question": "How does the frozen HCC1395 singleton population narrow through M1 and M2 gates?",
            "rationale": "Ordered counts preserve every denominator and prevent M2 NOT_EVALUABLE from being read as a biological negative.",
            "showDescription": True,
            "encodings": {
                "x": {"field": "stage", "type": "ordinal", "label": "Gate stage"},
                "y": {"field": "count", "type": "quantitative", "label": "Sites"},
            },
            "palette": {"kind": "sequential"},
            "valueFormat": "number",
            "comparisonContext": {
                "grain": "HCC1395 positional-singleton site",
                "unit": "sites",
                "denominator": "stage-specific frozen populations",
            },
        },
        {
            "id": "hcc1395_k_distribution_chart",
            "title": "734 個 M1 穩定多群位點的 K 分布",
            "subtitle": "K=2 佔 651/734=88.69%；K 是 methyl read groups 數，不是 clone 數。",
            "type": "bar",
            "intent": "comparison",
            "dataset": "hcc1395_k_distribution",
            "sourceId": "k_distribution_sql",
            "question": "How many stable focal-ALT methyl groups are observed per HCC1395 M1 site?",
            "rationale": "A zero-based bar chart makes the dominant K=2 pattern and rare higher-K tail directly visible.",
            "showDescription": True,
            "encodings": {
                "x": {"field": "k_label", "type": "ordinal", "label": "Methyl groups (K)"},
                "y": {"field": "count", "type": "quantitative", "label": "Sites"},
            },
            "palette": {"kind": "sequential"},
            "valueFormat": "number",
            "comparisonContext": {
                "grain": "M1 stable HCC1395 site",
                "unit": "sites",
                "denominator": "734 M1 stable multigroup sites",
            },
        },
        {
            "id": "hcc1395_support_distribution_chart",
            "title": "734 位點的 grouped focal-ALT read support 分布",
            "subtitle": "各 support bin 依 K 堆疊；全部位點 assignment ARI=1.0，且每群至少 3 reads。",
            "type": "stackedBar",
            "intent": "distribution",
            "dataset": "hcc1395_support_distribution",
            "sourceId": "support_distribution_sql",
            "question": "How much grouped focal-ALT read support backs the retained sites, and how is K composed within each support bin?",
            "rationale": "Binned stacked counts reveal the support distribution without overplotting 734 individual points; the full site table preserves exact rows.",
            "showDescription": True,
            "encodings": {
                "x": {"field": "support_bin", "type": "ordinal", "label": "Grouped focal-ALT reads"},
                "y": {"field": "site_count", "type": "quantitative", "label": "Sites"},
                "color": {"field": "k_label", "type": "nominal", "label": "K"},
            },
            "palette": {"kind": "categorical"},
            "comparisonContext": {
                "grain": "M1 stable HCC1395 site",
                "unit": "sites",
                "denominator": "734 M1 stable multigroup sites",
            },
        },
        old_charts["chr14_86272476_A_T_distance_heatmap"],
        old_charts["chr22_47466517_A_G_distance_heatmap"],
    ]
    manifest["charts"] = charts

    old_tables = {table["id"]: table for table in manifest["tables"]}
    manifest["tables"] = [
        {
            "id": "hcc1395_734_sites_table",
            "title": "HCC1395 734 個 M1 穩定甲基多群位點",
            "subtitle": "逐位點可搜尋/排序；cluster profile 只呈現群大小，不暴露內部分群標籤。",
            "dataset": "hcc1395_m1_sites",
            "sourceId": "m1_sites_sql",
            "defaultSort": {"field": "core_reads", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "locus", "label": "Locus", "type": "text"},
                {"field": "truth", "label": "Truth", "type": "text"},
                {"field": "k", "label": "K", "type": "number"},
                {"field": "core_reads", "label": "Core ALT reads", "type": "number"},
                {"field": "grouped_reads", "label": "Grouped reads", "type": "number"},
                {"field": "other_after_peel_reads", "label": "Other after peel", "type": "number"},
                {"field": "peeled_reads", "label": "Peeled", "type": "number"},
                {"field": "cluster_profile", "label": "Cluster sizes", "type": "text"},
                {"field": "min_cluster_n", "label": "Min group n", "type": "number"},
                {"field": "min_group_fraction", "label": "Min group share", "type": "percent"},
                {"field": "ari_min", "label": "Min ARI", "type": "number"},
                {"field": "unstable", "label": "Unstable", "type": "text"},
                {"field": "m2_status", "label": "M2", "type": "text"},
            ],
        },
        old_tables["case_statistics_table"],
        old_tables["source_inventory_table"],
    ]

    old_blocks = {block["id"]: block for block in manifest["blocks"]}
    manifest["blocks"] = [
        {"id": "title", "type": "markdown", "body": f"# {title}"},
        {
            "id": "executive_summary",
            "type": "markdown",
            "sourceId": "singleton_audit",
            "body": (
                "## 重點結論\n\n"
                "- **734/8,279=8.87%**：以全部 HCC1395 positional-singleton sSNV 為分母。\n"
                "  Wilson 95% CI 為 **8.27–9.50%**。\n"
                "- **734/8,074=9.09%**：以 M1 可評估位點為分母；Wilson 95% CI 為 **8.48–9.74%**。\n"
                "- 734 個均為 `M1 FLAGGED + stable_null_multigroup=true`；assignment ARI 最小值實際為 **1.0**。"
                "合計 **828/828 accepted splits** 的 observed ratio 均高於各自 column-shuffle null95。\n"
                "- 群數以 **K=2：651/734=88.69%** 為主；K=3 為 73、K=4 為 9、K=5 為 1。\n"
                "- 這證明的是 **共同 focal ALT reads 內可重現的甲基分子多群訊號**；不能直接換算成 734 個 subclones。"
            ),
        },
        {
            "id": "headline_metrics",
            "type": "metric-strip",
            "cardIds": ["singleton_sites_card", "m1_evaluable_card", "m1_stable_card", "m2_pass_card"],
        },
        {
            "id": "evidence_boundary",
            "type": "markdown",
            "sourceId": "claim_contract",
            "body": (
                "## 734 的證據強度與界線\n\n"
                "**可以說：** 734 個位點通過 M1 的統計 null-calibration 與跨 seed assignment 穩定性門檻，"
                "是 focal-ALT read 內的 methyl-defined multigroup candidates。\n\n"
                "**不能說：** 734 個位點已證實各自代表不同細胞 clone/subclone、clone 數，或父子突變關係。"
                "UPGMA 的枝長是甲基距離，不是時間或突變祖先方向。"
            ),
        },
        {"id": "funnel_heading", "type": "markdown", "body": "## 從 8,279 到 734：分母與 gate"},
        {"id": "hcc1395_734_funnel", "type": "chart", "chartId": "hcc1395_734_funnel_chart"},
        {
            "id": "funnel_explanation",
            "type": "markdown",
            "sourceId": "singleton_audit",
            "body": (
                "### 圖怎麼讀\n\n"
                "8,279 是所有 positional-singleton 位點；其中 8,074 可執行 M1，734 通過穩定多群 gate。"
                "M2 對 734 個逐一檢查額外測量軸，但只有 2 個 determinate 且 PASS；另外 732 個是 "
                "NOT_EVALUABLE，**不是 FAIL，也不是沒有甲基多群**。"
            ),
        },
        {"id": "distribution_heading", "type": "markdown", "body": "## 734 個位點的群數與穩定度"},
        {
            "id": "stability_metrics",
            "type": "metric-strip",
            "cardIds": ["accepted_splits_card", "split_ratio_card", "null_margin_card", "seed_stability_card"],
        },
        {"id": "k_distribution", "type": "chart", "chartId": "hcc1395_k_distribution_chart"},
        {"id": "support_distribution", "type": "chart", "chartId": "hcc1395_support_distribution_chart"},
        {
            "id": "stability_explanation",
            "type": "markdown",
            "sourceId": "recount_audit",
            "body": (
                "### 穩定性不是只看一張樹\n\n"
                "每個點是一個位點。M1 先要求觀察到的群間/群內分離超過門檻，再以 column-shuffle null "
                "校準，最後在 10 個 seeds 檢查 modal K 與 assignment agreement。堆疊圖顯示保留下來的位點"
                "在 grouped read support bins 中的 K 組成；ARI 的最小值另由上方卡片顯示為 1.0。"
                "734 不是只憑視覺挑出的案例。\n\n"
                "**統計口徑提醒：** operational gate 是 observed ratio > interpolated null 95th percentile。"
                "因只有 40 shuffles，另算的 empirical rank p 最高可為 3/41=0.0732；所以不能把 828 個 splits "
                "改寫成『每個 classical p<0.05』。"
            ),
        },
        {"id": "site_table_heading", "type": "markdown", "body": "## 734 個位點逐筆查核"},
        {"id": "hcc1395_734_sites", "type": "table", "tableId": "hcc1395_734_sites_table", "layout": "full"},
        {
            "id": "upgma_teaching",
            "type": "markdown",
            "sourceId": "claim_contract",
            "body": (
                "## 代表位點：UPGMA、read-distance 與 read×CpG 怎麼一起看\n\n"
                "1. **UPGMA**：把 representative focal-ALT reads 依 Bernoulli 甲基距離作 average linkage；"
                "看群是否在距離空間分開。\n"
                "2. **read-distance heatmap**：檢查組內較近、組間較遠的 block pattern。\n"
                "3. **read×CpG heatmap**：確認差異來自多個 CpG 的圖樣，而非單一缺失值；灰色保留未呼叫 CpG。\n\n"
                "三圖相容才支持 read-level epigenetic partition；仍不等於 cellular lineage tree。"
            ),
        },
        old_blocks["case_statistics"],
        old_blocks["chr14_86272476_A_T_heading"],
        dendrogram_blocks["chr14_upgma"],
        old_blocks["chr14_86272476_A_T_distance"],
        old_blocks["chr14_86272476_A_T_methyl"],
        old_blocks["chr14_86272476_A_T_explanation"],
        old_blocks["chr22_47466517_A_G_heading"],
        dendrogram_blocks["chr22_upgma"],
        old_blocks["chr22_47466517_A_G_distance"],
        old_blocks["chr22_47466517_A_G_methyl"],
        old_blocks["chr22_47466517_A_G_explanation"],
        {
            "id": "method",
            "type": "markdown",
            "sourceId": "recount_audit",
            "body": (
                "## 方法與可重現門檻\n\n"
                "1. Scope 固定為 HCC1395、chr1–22、8,279 個 positional-singleton sSNV。\n"
                "2. read 間距離採 Bernoulli methylation distance；average linkage/UPGMA 只作甲基距離分群。\n"
                "3. M1 recursive split 要求每群至少 3 reads、between/within≥1.3，觀察值高於 40 次 "
                "column-shuffle null 的第 95 百分位（有效 null≥32）。\n"
                "4. 重複 10 seeds；modal K≥70%、modal assignment minimum ARI≥0.80、unstable=false。\n"
                "5. 本報告再從 source-authority-v7 逐列重算 8,279／8,074／734 與 K 分布，並輸出全部 734 列。"
            ),
        },
        {
            "id": "limitations",
            "type": "markdown",
            "sourceId": "claim_contract",
            "body": (
                "## 限制與正確結論\n\n"
                "- M1 的 734 是 **統計穩定 methyl multigroup candidates**，不是 biological prevalence。\n"
                "- 只有 2 個位點在目前 M2 測量軸下 determinate/PASS；其他 732 個資訊不足，不能判為反例。\n"
                "- 要升級為 clone/subclone 證據，仍需第二個獨立 sSNV marker、tumor-REF/matched-normal control、"
                "exact-locus CN/purity/multiplicity/CCF，以及跨技術或 single-cell/colony/spatial 正交驗證。\n"
                "- 因此目前最精確的結論是：**HCC1395 約 9% 可評估 positional-singleton 位點，能在共享 ALT 的 reads "
                "內偵測到可重現的甲基分子子結構。**"
            ),
        },
        old_blocks["source_inventory"],
    ]

    artifact["snapshot"]["generatedAt"] = now
    artifact["snapshot"]["status"] = "ready"
    artifact["snapshot"]["datasets"].update(focused_datasets)

    inventory = artifact["snapshot"]["datasets"]["public_source_inventory"]
    current_inventory = {
        "positional_singleton_audit_summary.json": AUDIT_SUMMARY,
        "positional_singleton_site_audit.tsv.gz": SITE_AUDIT,
        "positional_singleton_m2_pass_cases.tsv": M2_PASS_CASES,
    }
    for row in inventory:
        if row["file"] in current_inventory:
            path = current_inventory[row["file"]]
            row["sha256"] = sha256_file(path)
            row["size_bytes"] = path.stat().st_size
            row["role"] = f"source-authority-v7 {row['role']}"
    if not any(row["file"] == RECOUNT_AUDIT.name for row in inventory):
        inventory.append(
            {
                "file": RECOUNT_AUDIT.name,
                "role": "independent_hcc1395_assignment_null_recount",
                "sha256": sha256_file(RECOUNT_AUDIT),
                "size_bytes": RECOUNT_AUDIT.stat().st_size,
                "used_for": "734-key reconciliation, 828 accepted-split null/stability census, and Wilson intervals",
            }
        )

    artifact["package_info"] = {
        "deliveryMode": "portable-html",
        "mode": "report",
        "readOnly": True,
        "taskType": "comprehensive-validation",
    }

    with ARTIFACT_OUT.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    summary = {
        "generated_at": now,
        "task_type": "B comprehensive validation",
        "research_goals": ["G3", "G4"],
        "inputs": {
            "site_audit": str(SITE_AUDIT),
            "site_audit_sha256": sha256_file(SITE_AUDIT),
            "base_validated_artifact": str(BASE_ARTIFACT),
            "base_validated_artifact_sha256": sha256_file(BASE_ARTIFACT),
            "independent_recount": str(RECOUNT_AUDIT),
            "independent_recount_sha256": sha256_file(RECOUNT_AUDIT),
            "sqlite_staging": str(STAGING_DB),
            "sqlite_staging_sha256": sha256_file(STAGING_DB),
        },
        "observed": observed,
        "rates": {
            "m1_stable_over_all_singletons": singleton_rate,
            "m1_stable_over_evaluable": evaluable_rate,
        },
        "k_distribution": dict(sorted(k_counts.items())),
        "stability": {
            "minimum_ari": float(flagged["ari_min"].min()),
            "maximum_ari": float(flagged["ari_min"].max()),
            "minimum_group_n": int(flagged["min_cluster_n"].min()),
            "minimum_core_reads": int(flagged["core_reads"].min()),
            "median_core_reads": float(flagged["core_reads"].median()),
            "maximum_core_reads": int(flagged["core_reads"].max()),
            "minimum_grouped_reads": int(flagged["cluster_total"].min()),
            "sites_with_other_after_peel_reads": int((flagged["other_after_peel_reads"] > 0).sum()),
            "unstable_true": int(bool_series(flagged["unstable"]).sum()),
            "null_calibrated_true": int(bool_series(flagged["stable_null_multigroup"]).sum()),
        },
        "dendrograms": dendrogram_results,
        "outputs": {
            "artifact": str(ARTIFACT_OUT),
            "artifact_sha256": sha256_file(ARTIFACT_OUT),
            "site_table": str(SITE_TABLE_OUT),
            "site_table_sha256": sha256_file(SITE_TABLE_OUT),
        },
        "claim_ceiling": (
            "Stable focal-ALT read-level methyl multigroup candidates; not confirmed cellular clones, "
            "subclones, or mutation ancestry trees."
        ),
    }
    with BUILD_SUMMARY_OUT.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
