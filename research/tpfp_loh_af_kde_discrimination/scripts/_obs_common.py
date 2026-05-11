"""Common helpers for obs01-obs08 scripts.

Shared dataset loader + style constants so each obs script can stay focused on
its own plotting logic.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination")
DATA_DIR = PROJECT_DIR / "data"
FIG_NEW_DIR = PROJECT_DIR / "figures" / "new"
MASTER_TSV = DATA_DIR / "master.tsv.gz"

PALETTE = {
    "dark": "#1E2A44",
    "bg": "#F7F3EC",
    "accent": "#A85540",
    "green": "#009E73",
    "red": "#D55E00",
    "blue": "#5B8DB7",
    "grey": "#8A8A8A",
    "tp": "#2E7D5B",
    "fp": "#D55E00",
}

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "HCC1937",
    "HCC1954", "H2009", "H1437", "COLO829",
]

LOH_ORDER = ["None", "LOH_Noise", "LOH_Weak", "LOH_Subclone", "LOH_Strong"]
AF_ORDER = ["Extreme", "Intermediate", "Near-half"]
CN_ORDER = ["T0", "T1", "T2", "T3", "T4"]
NG_ORDER = [1, 2, 3, 4]
NR_ORDER = ["low", "mid", "high"]

DIM_ORDERS = {
    "LOH_Subtype": LOH_ORDER,
    "AF_class": AF_ORDER,
    "cn_tier_F": CN_ORDER,
    "HPFineNGroups": NG_ORDER,
    "nr_band": NR_ORDER,
}

DIM_SHORT_LABEL = {
    "LOH_Subtype": "LOH",
    "AF_class": "AF",
    "cn_tier_F": "CN",
    "HPFineNGroups": "NG",
    "nr_band": "NR",
}

# Sample-mode combos present in master TSV (ordered for per-sample grid)
# Grouped by mode: 6 TO sample-modes first, then 7 paired_full. 13 total → 4×4 grid (3 empty slots)
SAMPLE_MODE_ORDER: list[tuple[str, str]] = [
    # TO mode (6 samples)
    ("HCC1395", "to_pileup"),
    ("HCC1395_DORADO", "to_pileup"),
    ("HCC1937", "to_pileup"),
    ("HCC1954", "to_pileup"),
    ("H2009", "to_pileup"),
    ("H1437", "to_pileup"),
    # paired_full mode (7 samples)
    ("COLO829", "paired_full"),
    ("HCC1395", "paired_full"),
    ("HCC1395_DORADO", "paired_full"),
    ("HCC1937", "paired_full"),
    ("HCC1954", "paired_full"),
    ("H2009", "paired_full"),
    ("H1437", "paired_full"),
]

# Per sample-mode metadata for power-tier annotation (FP counts from master TSV)
# 2026-04-22 update: 5 new TO samples from archive_rerun_20260422 added
SAMPLE_FP_META: dict[tuple[str, str], dict] = {
    # TO mode
    ("HCC1395", "to_pileup"):        {"fp": 11606, "tier": "[HIGH]", "flag": None},
    ("HCC1395_DORADO", "to_pileup"): {"fp": 11572, "tier": "[HIGH]", "flag": None},
    ("HCC1937", "to_pileup"):        {"fp": 12032, "tier": "[HIGH]", "flag": None},
    ("HCC1954", "to_pileup"):        {"fp": 50218, "tier": "[HIGH]", "flag": None},
    ("H2009", "to_pileup"):          {"fp": 11989, "tier": "[HIGH]", "flag": None},
    ("H1437", "to_pileup"):          {"fp": 13442, "tier": "[HIGH]", "flag": None},
    # paired_full mode
    ("COLO829", "paired_full"):        {"fp": 2273, "tier": "[HIGH]", "flag": None},
    ("HCC1395", "paired_full"):        {"fp": 627, "tier": "[MID]", "flag": "FP<1000"},
    ("HCC1395_DORADO", "paired_full"): {"fp": 240, "tier": "[MID]", "flag": "FP<1000"},
    ("HCC1937", "paired_full"):        {"fp": 195, "tier": "[MID]", "flag": "FP<1000"},
    ("H2009", "paired_full"):          {"fp": 86,  "tier": "[LOW]", "flag": "FP<100"},
    ("HCC1954", "paired_full"):        {"fp": 29,  "tier": "[LOW]", "flag": "FP<100"},
    ("H1437", "paired_full"):          {"fp": 8,   "tier": "[LOW]", "flag": "FP<20 — not interpretable"},
}

# S1-S5 scheme definitions (mirrored from step6_tpfp_detailed.py)
def build_scheme_masks(df: pd.DataFrame) -> dict:
    cn_tier = df["CovM_used"].copy()
    # SEQC2-grounded boundaries
    edges = [0.65, 0.99, 1.33, 1.82]
    tiers = pd.cut(cn_tier, bins=[-np.inf] + edges + [np.inf],
                   labels=[f"T{i}" for i in range(len(edges) + 1)], include_lowest=True)
    loh = df["LOH_Subtype"].fillna("None")
    afc = df["AF_class"]
    near_half = afc == "Near-half"
    extreme = afc == "Extreme"
    intermediate = afc == "Intermediate"
    t2 = tiers == "T2"  # diploid proxy (CN=2)
    low_cn = tiers.isin(["T0", "T1", "T2"])
    high_cn = tiers.isin(["T3", "T4"])

    # v2 scheme masks — aligned with step6_tpfp_detailed.py SCHEME_DEFS
    S1 = (loh == "LOH_Strong") & extreme                                 # LOH_Strong + Extreme
    S2 = (loh.isin(["LOH_Subclone", "LOH_Weak"])) & intermediate         # Subclonal LOH + Intermediate AF
    S3 = (loh == "None") & near_half & tiers.astype(str).isin(["T1", "T2"])  # Diploid Het (CN=1-2)
    S4 = (loh == "None") & extreme                                        # NonLOH Extreme (no-discriminator)
    S5 = (S1 | S2 | S3) & (~S4)                                           # Combo union minus S4 override
    return {"S1": S1, "S2": S2, "S3": S3, "S4": S4, "S5": S5}


def load_master() -> pd.DataFrame:
    df = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip", low_memory=False)
    if "AF_class" not in df.columns:
        raise ValueError("AF_class missing from master TSV")
    return df


def ensure_fig_dir() -> Path:
    FIG_NEW_DIR.mkdir(parents=True, exist_ok=True)
    return FIG_NEW_DIR


def apply_style() -> None:
    # Use Noto Sans CJK JP (available on this box) as primary so 中文 glyphs render.
    # Fall back to DejaVu Sans for any Latin chars it may lack.
    plt.rcParams.update({
        "figure.facecolor": PALETTE["bg"],
        "axes.facecolor": PALETTE["bg"],
        "axes.edgecolor": PALETTE["dark"],
        "axes.labelcolor": PALETTE["dark"],
        "xtick.color": PALETTE["dark"],
        "ytick.color": PALETTE["dark"],
        "text.color": PALETTE["dark"],
        "font.family": "sans-serif",
        "font.sans-serif": ["Noto Sans CJK JP", "DejaVu Sans", "Arial"],
        "axes.unicode_minus": False,
    })


# ---------------------------------------------------------------------------
# L1-L3 + consistency helpers (added 2026-04-22 for obs09-obs12)
# ---------------------------------------------------------------------------

CN_EDGES = [0.65, 0.99, 1.33, 1.82]


def build_5d_key(df: pd.DataFrame) -> pd.DataFrame:
    """Return DataFrame with 5 categorical keys added as columns (inplace-free).

    Adds / normalises:
      - LOH_Subtype (NaN → "None")
      - AF_class (existing)
      - cn_tier_F (derived from CovM_used with SEQC2 edges)
      - HPFineNGroups (int, drops 0 rows via NaN)
      - nr_band (low<60, mid 60-120, high>=120)
    """
    out = df.copy()
    out["LOH_Subtype"] = out["LOH_Subtype"].fillna("None").astype(str)
    out["cn_tier_F"] = pd.cut(
        out["CovM_used"],
        bins=[-np.inf] + CN_EDGES + [np.inf],
        labels=[f"T{i}" for i in range(len(CN_EDGES) + 1)],
        include_lowest=True,
    ).astype(str)
    out["nr_band"] = out["NumReads"].apply(
        lambda x: "low" if x < 60 else ("mid" if x < 120 else "high")
    )
    out["HPFineNGroups"] = pd.to_numeric(out["HPFineNGroups"], errors="coerce").astype("Int64")
    return out


def compute_tp_rate_by_key(
    df: pd.DataFrame,
    keys: list[str],
    min_n: int = 20,
) -> pd.DataFrame:
    """Groupby `keys`, return tidy DataFrame with n / n_tp / tp_rate / wilson_lo / wilson_hi.

    Rows with n < min_n have tp_rate/wilson masked NaN (count kept).
    """
    agg = (
        df.groupby(keys, observed=True, dropna=False)
        .agg(n=("tp_label", "size"), n_tp=("tp_label", "sum"))
        .reset_index()
    )
    agg["tp_rate"] = agg["n_tp"] / agg["n"]
    # Wilson 95% CI
    z = 1.96
    p = agg["tp_rate"]
    n = agg["n"].astype(float)
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    halfwidth = (z / denom) * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    agg["wilson_lo"] = centre - halfwidth
    agg["wilson_hi"] = centre + halfwidth
    mask_low = agg["n"] < min_n
    for col in ("tp_rate", "wilson_lo", "wilson_hi"):
        agg.loc[mask_low, col] = np.nan
    return agg


def per_sample_grid_axes(
    n_rows: int = 4,
    n_cols: int = 4,
    figsize: tuple[float, float] = (20.0, 16.0),
    n_used: int | None = None,
) -> tuple[plt.Figure, list[plt.Axes]]:
    """Per-sample grid (row-major). Default 4x4=16 slots for 13 SAMPLE_MODE_ORDER combos.

    If ``n_used`` is provided, axes beyond that index are hidden (axis('off')).
    Returns the full flat axes list (length = n_rows * n_cols) so callers can
    zip with SAMPLE_MODE_ORDER and have trailing unused axes pre-hidden.
    """
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes_flat = list(axes.flatten()) if n_rows * n_cols > 1 else [axes]
    if n_used is not None:
        for i in range(n_used, len(axes_flat)):
            axes_flat[i].axis("off")
    return fig, axes_flat


def sample_label(sample: str, mode: str) -> str:
    """Short label for panel titles."""
    mode_short = "TO" if mode == "to_pileup" else "pf"
    return f"{sample}_{mode_short}"


def annotate_power_tier(
    ax: plt.Axes,
    sample: str,
    mode: str,
    baseline_tp: float | None = None,
    n_total: int | None = None,
) -> None:
    """Add power-tier + optional baseline + n badge to upper-right of a panel.

    2026-04-22 update: accepts baseline_tp and n_total so each panel self-documents
    its sample-level context (per the V2 plan §T11 requirement).
    Draws red frame when FP<100.
    """
    meta = SAMPLE_FP_META.get((sample, mode))
    if meta is None:
        return
    lines = [f"{meta['tier']} FP={meta['fp']:,}"]
    if baseline_tp is not None:
        lines.append(f"base={baseline_tp:.2f}")
    if n_total is not None:
        lines.append(f"n={n_total:,}")
    if meta["flag"]:
        lines.append(meta["flag"])
    ax.text(
        0.98, 0.98, "\n".join(lines), transform=ax.transAxes,
        ha="right", va="top", fontsize=7, color=PALETTE["dark"],
        bbox=dict(facecolor="white", edgecolor=PALETTE["grey"], alpha=0.85, pad=1.8),
    )
    if meta["fp"] < 100:
        for spine in ax.spines.values():
            spine.set_edgecolor(PALETTE["red"])
            spine.set_linewidth(1.5)


# ---------------------------------------------------------------------------
# V2 (2026-04-22) helpers — research-question suptitle + takeaway caption box
# + TO/paired grid split
# ---------------------------------------------------------------------------

SAMPLE_MODE_TO: list[tuple[str, str]] = [sm for sm in SAMPLE_MODE_ORDER if sm[1] == "to_pileup"]
SAMPLE_MODE_PAIRED: list[tuple[str, str]] = [sm for sm in SAMPLE_MODE_ORDER if sm[1] == "paired_full"]


def research_suptitle(
    fig: plt.Figure,
    question: str,
    context: str | None = None,
    y: float = 0.995,
    fontsize: int = 16,
) -> None:
    """Set a research-question style figure suptitle.

    The top line is the research question (larger, bold).
    The optional ``context`` line is a smaller plain-text sub-line (no LaTeX markup,
    so CJK renders correctly).
    """
    # Main question: larger, bold
    fig.suptitle(
        question,
        fontsize=fontsize,
        fontweight="bold",
        color=PALETTE["dark"],
        y=y,
    )
    if context:
        # Secondary context line just below suptitle; smaller + grey
        fig.text(
            0.5, y - 0.018,
            context,
            ha="center", va="top",
            fontsize=max(9, fontsize - 4),
            color=PALETTE["grey"],
            style="italic",
        )


def takeaway_caption(
    fig: plt.Figure,
    text: str,
    color: str = "green",
    y: float = 0.01,
    fontsize: int = 10,
) -> None:
    """Render a coloured takeaway caption box at the bottom of the figure.

    ``color``: "green" (支持/一致), "red" (駁斥/異常), "blue" (中性觀察), "grey" (待驗證).
    """
    bg = {
        "green": "#E6F3EC",
        "red": "#FBE5DC",
        "blue": "#E1ECF6",
        "grey": "#EAEAEA",
    }.get(color, "#EAEAEA")
    edge = {
        "green": PALETTE["green"],
        "red": PALETTE["red"],
        "blue": PALETTE["blue"],
        "grey": PALETTE["grey"],
    }.get(color, PALETTE["grey"])
    fig.text(
        0.5, y,
        text,
        ha="center", va="bottom",
        fontsize=fontsize,
        color=PALETTE["dark"],
        wrap=True,
        bbox=dict(
            facecolor=bg,
            edgecolor=edge,
            linewidth=1.2,
            pad=6,
            boxstyle="round,pad=0.4",
        ),
    )


def mode_subset_grid(
    mode: str,
    figsize: tuple[float, float] | None = None,
) -> tuple[plt.Figure, list[plt.Axes], list[tuple[str, str]]]:
    """Per-mode subset grid for TO-only (6 samples, 2x3) or paired-only (7 samples, 3x3).

    Returns (fig, flat_axes, sample_mode_list). Unused axes are pre-hidden.
    """
    if mode == "to_pileup":
        sm_list = SAMPLE_MODE_TO   # 6 samples
        n_rows, n_cols = 2, 3
        # Wider + taller for readability after 2026-04-22 V2 review
        default_fig = (22.0, 13.0)
    elif mode == "paired_full":
        sm_list = SAMPLE_MODE_PAIRED   # 7 samples
        n_rows, n_cols = 3, 3
        default_fig = (22.0, 17.0)
    else:
        raise ValueError(f"Unknown mode: {mode}")
    fs = figsize or default_fig
    fig, axes = plt.subplots(n_rows, n_cols, figsize=fs)
    axes_flat = list(axes.flatten())
    for i in range(len(sm_list), len(axes_flat)):
        axes_flat[i].axis("off")
    return fig, axes_flat, sm_list


def mode_display_name(mode: str) -> str:
    return "TO mode (6 samples)" if mode == "to_pileup" else "paired_full mode (7 samples)"
