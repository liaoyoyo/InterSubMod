<!--
建立時間：2026-07-24
目標：標示合成測試與 oracle artifact 的權威版本
處理範圍：synthetic only；不可作 HCC1395 或跨樣本 evidence
關聯檔案：test_exact_ps_topology_af.py、synthetic_oracle_input.json
-->

# Synthetic oracle artifacts

- 權威輸入：`synthetic_oracle_input.json`
- 權威定稿輸出：`synthetic_oracle_output_v2.jsonl`
- 權威定稿 receipt：`synthetic_oracle_output_v2.receipt.json`

無版本尾碼的兩個 output 是 label contract 定稿前的編輯中快照；依 repo
no-delete 規範保留，但不得作 parity oracle。所有內容皆為 synthetic，
只驗證 schema、exact solver 與 AF factorization，不是生物結果。
