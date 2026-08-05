<!--
建立時間: 2026-07-18 00:00 +08:00
目標: 即時記錄 Hypercube exact-preserving prepared-base 與 structural cache 修補的設計決定、偏離、折衷與未決事項
處理範圍: 新 release source/tests/diagnostics；不修改 frozen v2 output
關聯檔案:
  - InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/pre-decision-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py
-->

# Hypercube exact-preserving solver/cache remediation — Implementation notes

> **Task type**：E Hotfix / performance-remediation probe  
> **服務目標**：G3 read-level方法可擴展性、G4 reproducibility、G5外部可驗證  
> **狀態**：implementation_complete／scale_probe_pending  
> **決策上限**：本輪只證明exact-preserving與局部工程收益；不宣稱正式全量已解鎖。

## 設計決定

1. `hypercube_exact.py`把固定 predecessor/group constraints準備一次，後續僅依序附加
   objective equality與no-good rows。
2. 保留SciPy `milp()`每次呼叫；不冒充persistent HiGHS model。
3. ranking cache只存在單一`run()`生命週期，dict以完整frozen tuple判等。
4. 只存`complete=true`且通過結構不變量重驗的enumeration；store/hit各deep-copy。
5. likelihood、BQ、error-grid與bootstrap每個unit重新計算，不進cache。
6. cache diagnostics只進receipt runtime區，不進unit/candidate/responsibility scientific
   semantic hashes。

## 偏離

- 無。若後續必須改candidate universe、objective、timeout語意或ranking分母，立即停止並另立spec。

## 折衷

- Python端可避免重複universe、group reduction、base sparse rows建構，但SciPy 1.13.1
  `milp()`仍會每輪建立HiGHS model與presolve；真正persistent backend留待`highspy`另案。
- cache key包含source SHA、solver time limit與max sets；可能降低命中率，但先換取可審核性。
- 不設定wall-time assertion；以call count、答案等價與counter守恆作穩定驗證。

## 未決

1. 33個長尾units的新source實際cache hit率與wall-time。
2. chr6與H2009 chr2雙pilot是否達wall≤4h、incomplete≤1%、exact-limit coverage≥90%。
3. 是否值得另案導入`highspy.Highs` persistent model、checkpoint/resume與per-unit deadline。

## 完成結果

- Source SHA：solver=`1def0de1952d...a32f77`；ranker=`b28d494563be...84f0a9`。
- Prepared-base probe：24 optima、ordered digest不變、build 25→1、solve仍25。
- Cache probe：3 lookups＝1 miss＋2 hits；solver invocations=1；incomplete stored=0。
- Targeted solver/ranker：61/61 PASS。
- Core solver/exhaustive/ranker/orchestration：93/93 PASS。
- Full research test discovery：341/341 PASS，136.479秒，exit 0。
- Static：5支Python檔`py_compile` exit 0；trailing-whitespace matches=0。
- Release：Runbook v2 SHA保持`8a347a...`；v3 SHA=`d278da...`且exact-11 3/3 PASS。
- Result：`PASS_FOR_PROBE`；真實scale仍須v3雙pilot。

## 修改前基準

- Git HEAD：`0ee2fa1`
- `hypercube_exact.py` SHA-256：
  `9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`
- `build_m2_patterns_and_rank.py` SHA-256：
  `c82210f25870d1405cc070c18096fb7d1c2b908fb4d8a7858aece7ac4b151d28`
- Targeted tests：51/51 PASS，2.638秒，exit 0。
- `AAAA,k=4,all_loci,max_sets=32`：24 candidates、objective=3、complete=true；
  `_build_problem` 25 calls、wall 0.266119秒。

## Step → Verify

1. Prepared-base refactor  
   → Verify：24-optima仍完整且ordered-ID digest相同；base build=1、MILP solve=25。
2. Complete-only process cache  
   → Verify：complete hit只solve一次；incomplete重算；mutation isolation與key separation通過。
3. Receipt diagnostics  
   → Verify：lookups=hits+misses、solver_invocations=misses、entries=stores-evictions。
4. Regression與oracle  
   → Verify：targeted、exhaustive、full ranking tests全PASS；`py_compile`與`git diff --check`通過。
5. 後續pilot gate  
   → Verify：新release雙pilot receipts PASS後，才可把scale verdict由PROBE升級。
