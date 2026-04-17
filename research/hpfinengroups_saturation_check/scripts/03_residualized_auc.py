"""
B.1-1 質疑驗證：HPFineNGroups residualized AUC 多變數控制

核心問題：原聲稱 residualized AUC = 0.617（控制 NR+AF）；若加入 Coverage_Multiple、LOH、AF，
         HPFineNGroups 的獨立預測力是否仍 ≥ 0.60？

方法：
  M1 (Raw):       AUC(TP ~ NGroups)
  M2 (L1 NR):     殘差化 TP ~ NR 後 AUC(residual ~ NGroups)
  M2b (L1 AF):    殘差化 TP ~ AF 後 AUC
  M2c (L1 Cov):   殘差化 TP ~ Coverage 後 AUC
  M3 (L2 NR+AF):  殘差化 TP ~ NR+AF 後 AUC  (原聲稱 0.617)
  M4 (L3 All):    殘差化 TP ~ NR+AF+Coverage_Multiple+LOH 後 AUC

兩種殘差化策略：
  A) residual approach: y_residual = y - p_hat(X_control); AUC(y_residual ~ NGroups)
  B) Δ AUC approach: AUC(X_control+NGroups) - AUC(X_control alone)
  兩者都呈現

輸出：
  data/residualized_auc_to_nonloh.tsv
  data/residualized_auc_to_all.tsv
  data/residualized_auc_paired.tsv
  data/residualized_auc_per_sample.tsv
  figures/03_residualized_auc_bar.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
DATA_TSV = ROOT / "output" / "synthesis" / "observation_workspaces" / \
           "20260327_loh_round1_cross_sample_audit" / "all_region_rows.tsv.gz"
OUT = Path(__file__).resolve().parents[1]
DATA_OUT = OUT / "data"
FIG_OUT = OUT / "figures"


def load() -> pd.DataFrame:
    cols = ["sample", "mode", "truth_label", "NumReads", "HPFineNGroups",
            "to_loh_bed_hit", "caller_af", "Coverage_Multiple"]
    df = pd.read_csv(DATA_TSV, sep="\t", usecols=cols,
                     dtype={"sample": "category", "mode": "category",
                            "truth_label": "category", "to_loh_bed_hit": "category"})
    df["y"] = (df["truth_label"] == "TP").astype(int)
    df["is_loh"] = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1"]).astype(int)
    df = df.dropna(subset=["NumReads", "HPFineNGroups", "caller_af", "Coverage_Multiple"])
    df["NumReads"] = df["NumReads"].astype(float)
    df["HPFineNGroups"] = df["HPFineNGroups"].astype(float)
    df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
    df["Coverage_Multiple"] = pd.to_numeric(df["Coverage_Multiple"], errors="coerce")
    df = df.dropna(subset=["caller_af", "Coverage_Multiple"])
    return df


def compute_auc(y: np.ndarray, score: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, score))


def delta_auc(y: np.ndarray, X_control: np.ndarray, X_feat: np.ndarray) -> dict:
    """AUC(control+feat) - AUC(control alone)."""
    if len(np.unique(y)) < 2:
        return {"auc_control": float("nan"), "auc_full": float("nan"), "delta": float("nan")}
    model_c = LogisticRegression(max_iter=500, solver="lbfgs")
    model_c.fit(X_control, y)
    auc_c = compute_auc(y, model_c.predict_proba(X_control)[:, 1])
    X_full = np.concatenate([X_control, X_feat], axis=1)
    model_f = LogisticRegression(max_iter=500, solver="lbfgs")
    model_f.fit(X_full, y)
    auc_f = compute_auc(y, model_f.predict_proba(X_full)[:, 1])
    return {"auc_control": auc_c, "auc_full": auc_f, "delta": auc_f - auc_c}


def residual_auc(y: np.ndarray, X_control: np.ndarray, feat: np.ndarray) -> dict:
    """AUC of residual after removing control. Positive residual approach."""
    if len(np.unique(y)) < 2:
        return {"auc_raw": float("nan"), "auc_residual": float("nan")}
    auc_raw = compute_auc(y, feat)
    model = LogisticRegression(max_iter=500, solver="lbfgs")
    model.fit(X_control, y)
    p_hat = model.predict_proba(X_control)[:, 1]
    resid = y - p_hat  # signed residual
    # use feat to predict resid sign (correlation-based surrogate) → equivalent to Δ AUC but with feat alone
    # here we compute AUC of feat against binary "y>p_hat" which tests if feat adds above baseline
    y_above = (y > p_hat).astype(int)
    auc_resid = compute_auc(y_above, feat)
    return {"auc_raw": auc_raw, "auc_residual": auc_resid}


def run_panel(df: pd.DataFrame, label: str) -> dict:
    """Run M1-M4 panel on a dataframe."""
    if df.empty or len(np.unique(df["y"])) < 2:
        return {"label": label, "n": len(df)}
    y = df["y"].to_numpy()
    feat = df[["HPFineNGroups"]].to_numpy()
    # Control sets
    nr = df[["NumReads"]].to_numpy()
    af = df[["caller_af"]].to_numpy()
    cov = df[["Coverage_Multiple"]].to_numpy()
    loh = df[["is_loh"]].to_numpy()
    nr_af = df[["NumReads", "caller_af"]].to_numpy()
    all_ctrl = df[["NumReads", "caller_af", "Coverage_Multiple", "is_loh"]].to_numpy()
    # Compute
    M1 = compute_auc(y, df["HPFineNGroups"])
    M2_nr = delta_auc(y, nr, feat)
    M2_af = delta_auc(y, af, feat)
    M2_cov = delta_auc(y, cov, feat)
    M3_nr_af = delta_auc(y, nr_af, feat)
    M4_all = delta_auc(y, all_ctrl, feat)
    res_all = residual_auc(y, all_ctrl, df["HPFineNGroups"].to_numpy())
    return {
        "label": label, "n": len(df),
        "tp_rate": float(y.mean()),
        "M1_raw_auc": M1,
        "M2_NR_auc_control": M2_nr["auc_control"], "M2_NR_auc_full": M2_nr["auc_full"], "M2_NR_delta": M2_nr["delta"],
        "M2b_AF_auc_control": M2_af["auc_control"], "M2b_AF_auc_full": M2_af["auc_full"], "M2b_AF_delta": M2_af["delta"],
        "M2c_Cov_auc_control": M2_cov["auc_control"], "M2c_Cov_auc_full": M2_cov["auc_full"], "M2c_Cov_delta": M2_cov["delta"],
        "M3_NR_AF_auc_control": M3_nr_af["auc_control"], "M3_NR_AF_auc_full": M3_nr_af["auc_full"], "M3_NR_AF_delta": M3_nr_af["delta"],
        "M4_ALL_auc_control": M4_all["auc_control"], "M4_ALL_auc_full": M4_all["auc_full"], "M4_ALL_delta": M4_all["delta"],
        "M4_residual_above_baseline_auc": res_all["auc_residual"],
    }


def main() -> int:
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    print("[Load]", flush=True)
    df = load()
    print(f"  rows: {len(df):,}")

    panels = []
    # Scenario panels
    panels.append(run_panel(df[(df["mode"] == "to") & (df["is_loh"] == 0)], "TO NonLOH all"))
    panels.append(run_panel(df[(df["mode"] == "to") & (df["is_loh"] == 0) & (df["NumReads"] >= 80)], "TO NonLOH NR>=80"))
    panels.append(run_panel(df[(df["mode"] == "to") & (df["is_loh"] == 1)], "TO LOH all"))
    panels.append(run_panel(df[df["mode"] == "to"], "TO all"))
    panels.append(run_panel(df[(df["mode"] == "paired") & (df["is_loh"] == 0)], "Paired NonLOH all"))
    panels.append(run_panel(df[df["mode"] == "paired"], "Paired all"))
    panel_df = pd.DataFrame(panels)
    panel_df.to_csv(DATA_OUT / "residualized_auc_panels.tsv", sep="\t", index=False)
    print("\n[Scenario panels]")
    print(panel_df.to_string(index=False))

    # Per-sample (TO NonLOH)
    print("\n[Per-sample (TO NonLOH)]")
    ps_rows = []
    for sample in sorted(df["sample"].unique()):
        sub = df[(df["sample"] == sample) & (df["mode"] == "to") & (df["is_loh"] == 0)]
        res = run_panel(sub, f"TO NonLOH {sample}")
        res["sample"] = sample
        ps_rows.append(res)
    ps_df = pd.DataFrame(ps_rows)
    ps_df.to_csv(DATA_OUT / "residualized_auc_per_sample.tsv", sep="\t", index=False)
    print(ps_df[["sample", "n", "tp_rate", "M1_raw_auc", "M3_NR_AF_auc_full",
                 "M3_NR_AF_delta", "M4_ALL_auc_full", "M4_ALL_delta"]].to_string(index=False))

    # Plot: horizontal bar of M1, M3, M4 across scenarios
    fig, ax = plt.subplots(figsize=(12, 6))
    y_labels = panel_df["label"].tolist()
    yp = np.arange(len(y_labels))
    width = 0.25
    ax.barh(yp - width, panel_df["M1_raw_auc"], height=width, label="M1 Raw AUC", color="steelblue")
    ax.barh(yp, panel_df["M3_NR_AF_auc_full"], height=width, label="M3 +NR+AF full model AUC", color="orange")
    ax.barh(yp + width, panel_df["M4_ALL_auc_full"], height=width, label="M4 +NR+AF+Cov+LOH full model AUC", color="green")
    ax.axvline(0.5, color="gray", lw=0.5, ls="--")
    ax.axvline(0.60, color="orange", lw=0.8, ls="--", label="0.60 threshold")
    ax.axvline(0.617, color="red", lw=0.8, ls="--", label="0.617 (original claim)")
    ax.set_yticks(yp)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("AUC")
    ax.set_title("B.1-1: Residualized AUC — HPFineNGroups 多變數控制")
    ax.legend(loc="lower right")
    ax.set_xlim(0.4, 1.0)
    for i, row in panel_df.iterrows():
        ax.text(row["M4_ALL_auc_full"] + 0.005, i + width,
                f"Δ(M4)={row['M4_ALL_delta']:+.4f}", fontsize=7, va="center")
    plt.tight_layout()
    plt.savefig(FIG_OUT / "03_residualized_auc_bar.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    # Final verdict
    to_nonloh = panel_df[panel_df["label"] == "TO NonLOH all"].iloc[0]
    print("\n=== FINAL VERDICT (TO NonLOH, B.1-1) ===")
    print(f"M1 raw AUC = {to_nonloh['M1_raw_auc']:.4f}")
    print(f"M3 NR+AF full AUC = {to_nonloh['M3_NR_AF_auc_full']:.4f} (original claim ≈0.617)")
    print(f"M3 Δ vs control = {to_nonloh['M3_NR_AF_delta']:+.4f}")
    print(f"M4 NR+AF+Cov+LOH full AUC = {to_nonloh['M4_ALL_auc_full']:.4f}")
    print(f"M4 Δ vs control = {to_nonloh['M4_ALL_delta']:+.4f}")
    if to_nonloh['M4_ALL_delta'] < 0.005:
        print("VERDICT: INCREMENT NEGLIGIBLE — HPFineNGroups 訊號完全被 NR+AF+Cov+LOH 吸收")
    elif to_nonloh['M4_ALL_delta'] < 0.02:
        print("VERDICT: WEAK INCREMENT — 結論降級考量")
    else:
        print("VERDICT: ROBUST INCREMENT — 結論 16 ⭐3 保留")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
