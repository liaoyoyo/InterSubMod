# figures/ 復原狀態（2026-05-11）

此目錄因 2026-05-11 git history rewrite 操作意外影響，原 56 PNG 全部消失。

## 9 個子目錄 + 對應重產 script（大多可重產）

| Sub-dir | Recovery script | Status |
|---------|----------------|--------|
| `20260423_B2_hcc1954_outlier/` | `InterSubMod/scripts/analysis/20260423_B2_hcc1954_outlier.py` | ✅ 可重產 |
| `20260423_phase3_synthesis/` | `InterSubMod/scripts/analysis/20260423_phase3_synthesis.py` | ✅ 可重產 |
| `20260425_obs18_af_split/` | （無對應 script） | ⚠ 需手動重產 |
| `20260426_hp_ratio_3d_auc/` | （無對應 script） | ⚠ 需手動重產 |
| `20260426_longphase_old_vs_new/` | （無對應 script） | ⚠ 需手動重產 |
| `20260426_methylation_3d_addon/` | （無對應 script） | ⚠ 需手動重產 |
| `20260430_v3f_ablation/` | `InterSubMod/scripts/analysis/v3f_ablation_*.py` | ✅ 可重產 |
| `20260501_latest_binary_f1/` | （無對應 script） | ⚠ 需手動重產 |
| `20260501_v5_force_path2only/` | `InterSubMod/scripts/analysis/v5_force_path2only_*.py` | ✅ 可重產 |

## 對應數據

數據還在 `output/` symlink → `/big7_disk/liaoyoyo2001/big7_disk_output/`（由 .gitignore 保護不入 git）。

## 注意

新 PNG 重產後**不會自動 commit**（pre-commit hook 阻擋 binary）。
這是業界規範：figures 由 script 重產，不存進 git。
