<!--
建立時間: 2026-07-13
目標: 驗證固定 seed 的兩次完整 compositional analysis 是否 byte-identical
處理範圍: 全部分析輸出檔
關聯檔案: reproducibility_hashes.tsv
-->

# Reproducibility receipt

> **PASS：18/18 files byte-identical。**

- First run：`/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/rerun_a`
- Second run：`/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/rerun_b`
- 比對方式：相對檔名集合＋每檔 SHA-256；包含 stochastic summaries、JSON、report與 provenance。
- Fixed seed與replicate數由各 run的 `provenance.json` 記錄。
