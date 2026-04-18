#!/usr/bin/env python3
"""Run only M5, M6, M7 and synthesize verdict (M1-M4 already complete)."""
import sys, time
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")

from build_beyond_auc_validation import (
    run_m5_bootstrap, run_m6_nonlinear, run_m7_distribution, synthesize_verdict,
    TSV_DIR, FIG_DIR, ALL_FEATURES, ISM_FEATURES, CALLER_FEATURES,
)
from observation_common import (
    load_master_dataset, setup_plot_style, encode_truth_binary, ensure_dir,
)
import pandas as pd

ensure_dir(TSV_DIR)
ensure_dir(FIG_DIR)
setup_plot_style()

df = load_master_dataset()
df["truth_binary"] = encode_truth_binary(df)
print(f"Paired: {(df['mode']=='paired').sum():,}  TO: {(df['mode']=='to').sum():,}")

t0 = time.time()
m5 = run_m5_bootstrap(df)
m6_cv, m6_imp = run_m6_nonlinear(df)
m7 = run_m7_distribution(df)

# Load existing M1-M4 results for verdict
m1 = pd.read_csv(TSV_DIR / "m1_pr_auc_results.tsv", sep="\t")
m2 = pd.read_csv(TSV_DIR / "m2_residualization_results.tsv", sep="\t")
m3_promising = pd.read_csv(TSV_DIR / "m3_promising_subgroups.tsv", sep="\t")

verdict_df, conclusion = synthesize_verdict(m1, m2, m3_promising, m5, m6_cv, m6_imp, m7)
print(f"\nDone in {(time.time()-t0)/60:.1f} min")
