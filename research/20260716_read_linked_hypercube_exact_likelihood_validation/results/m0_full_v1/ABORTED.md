<!--
建立時間: 2026-07-16 06:39 +08:00
目標: 記錄低效M0執行中止原因，避免把空目錄誤認為完成結果
處理範圍: current-v5全7 datasets；未產生分析表
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_m0_likelihood_census.py
-->

# M0 full v1 — ABORTED

- **Input**：current-v5 7-dataset `layered_region_view_*.json`。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_m0_likelihood_census.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v1`
- **Actual**：約10分鐘後由本cycle主動送出SIGINT；trace停在逐候選EM iteration，未產生TSV／receipt，也未修改canonical。
- **原因**：EM在simplex boundary附近收斂過慢，不適合約20萬以上vertex-set fits。
- **處理**：v2改用具有analytic gradient的SLSQP simplex optimizer，並以EM作小型oracle交叉驗證；v1目錄保留作audit，不刪除。
