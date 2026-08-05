<!--
建立時間: 2026-07-18 00:00 +08:00
目標: 以新的 exact-preserving prepared-base/cache source 身分，先執行雙 pilot gates，再決定是否進入 M2 全量 extraction/ranking
處理範圍: Task Type B；chr1-chr22 × 7 technical datasets / 6 biological samples；154 extraction + 154 ranking tasks；雙 direct pilots
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式凍結全量執行與驗證協定_01.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md
  - InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/pre-decision-audit.md
前版文件:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md
-->

# M2 exact-preserving修補正式執行 Runbook v3

> **TL;DR**：v3只更新ranker與Hypercube solver身分，加入prepared-base與complete-only process cache；流程仍固定為 `zero-conflict bootstrap → PRE → exact-11 freeze → 雙 pilot 6 gates → extraction 12 gates → ranking 12 gates → POST → independent verifier → final numeric summary → presentation freeze → FINAL HTML → frozen browser QA`。目前scale verdict是`PROBE`；雙pilot未通過前禁止進full 154。每個正式gate仍是不可覆寫JSON＋SHA sidecar，任何不符即停止。（影響：高；exact-equivalence：PASS；formal actual-data：WAIT）

> **Task Type**：B — Comprehensive validation。Subset 只是 pilot，不能代替 chr1–22 × 7 datasets 的最後數字。服務 G3／G4／G5。

## 0. 結果可以與不可以代表什麼

| 項目 | 正式定義 |
|---|---|
| technical datasets | 7：COLO829、H1437、H2009、HCC1395、HCC1395_DORADO、HCC1937、HCC1954 |
| biological samples | 6；HCC1395 與 HCC1395_DORADO 是同一生物樣本的兩份 technical datasets |
| chromosomes | chr1–chr22；chrX 因 ploidy、PAR、X-inactivation、sex/CN/LOH 模型不同而不混入此 autosomal primary run |
| primary reconstruction unit | HP family × exact known PS × read-linked component × bridge threshold |
| exact solver scope | observed-ALT effective k≤12；k>12 為 local-only／abstain，不能冒充 global optimum |
| h* | minimum-extra-state count；可能包含 partial representative 與 connector，不是 hidden clone count |
| V | distinct globally minimum-extra vertex sets 數 |
| T | 所有 V rows 的 parent-choice counts 加總；因此 T≥V≥1 |
| claim ceiling | solver-certified regional mutation-state candidate sets，以及宣告的 molecule/error contract 下可被 reads 區分的候選；不是唯一真實 clone tree |

Resource attestation 證明的是 frozen producer 所持久化的 launch-time process/disk snapshot；在「無 hostile same-UID actor」威脅模型下提供 checksum provenance，不是 trusted execution 或 cryptographic process ancestry。

## 1. 固定路徑、attempt-exclusive logs 與輔助函式

每次 fresh 或 resume 都建立新的 attempt 目錄；因此 stdout、stderr、GNU time 與 verify logs 不覆寫舊紀錄。

```bash
set -euo pipefail
umask 022
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

REPO=/big7_disk/liaoyoyo2001/InterSubMod
PYTHON=/bip7_disk/liaoyoyo2001/miniconda3/bin/python3
OUT_BASE=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation
RUN_ROOT="$OUT_BASE/m2_frozen_release_v3"
ATTEMPT_PARENT="$OUT_BASE/m2_frozen_release_v3_attempts"
mkdir -p "$ATTEMPT_PARENT"
test -d "$ATTEMPT_PARENT" && test ! -L "$ATTEMPT_PARENT"
ATTEMPT_ID="$(date -u +%Y%m%dT%H%M%S.%NZ)-$$-$RANDOM"
ATTEMPT="$ATTEMPT_PARENT/$ATTEMPT_ID"
mkdir "$ATTEMPT"

CONTROL="$RUN_ROOT/control"
CONTRACT="$RUN_ROOT/release_contract"
PILOTS="$RUN_ROOT/pilots"
PILOT_DORADO="$PILOTS/HCC1395_DORADO_chr6"
PILOT_H2009="$PILOTS/H2009_chr2"
EXTRACTION="$RUN_ROOT/full_extraction"
RANKING="$RUN_ROOT/full_ranking"
PRESENTATION="$RUN_ROOT/presentation_contract"
HTML="$RUN_ROOT/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.html"
HTML_QA="$RUN_ROOT/html_qa_final"

ORIGIN_MANIFEST=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json
PRE="$CONTROL/pre_input_identity_receipt.json"
POST="$CONTROL/post_input_identity_receipt.json"
FINAL_VERIFY="$CONTROL/full_m2_independent_verification.json"
NUMERIC="$CONTROL/final_numeric_summary.json"

LIVE_SCRIPTS="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts"
LIVE_INPUT_VERIFY="$LIVE_SCRIPTS/verify_m2_input_identity.py"
LIVE_FREEZER="$LIVE_SCRIPTS/freeze_m2_release_contract.py"
LIVE_FULL_EXTRACT="$LIVE_SCRIPTS/run_full_m2_extraction.py"
LIVE_NUMERIC="$LIVE_SCRIPTS/build_final_numeric_summary.py"
LIVE_PRESENTATION_FREEZER="$LIVE_SCRIPTS/freeze_presentation_snapshot.py"

cd "$REPO"
test -x "$PYTHON"
test "$(sha256sum "$ORIGIN_MANIFEST" | awk '{print $1}')" = 16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45

authenticate_sidecar() {
  local target="$1"
  test -f "$target" && test ! -L "$target"
  test -f "$target.sha256" && test ! -L "$target.sha256"
  test "$(stat -c %a "$target")" = 444
  test "$(stat -c %h "$target")" = 1
  test "$(stat -c %a "$target.sha256")" = 444
  test "$(stat -c %h "$target.sha256")" = 1
  (cd "$(dirname "$target")" && sha256sum -c "$(basename "$target").sha256")
}

seal_direct_receipt() {
  local target="$1" sidecar="$1.sha256" digest covered extra
  test -f "$target" && test ! -L "$target"
  test -f "$sidecar" && test ! -L "$sidecar"
  test "$(stat -c %h "$target")" -eq 1
  test "$(stat -c %h "$sidecar")" -eq 1
  test "$(wc -l < "$sidecar")" -eq 1
  LC_ALL=C grep -Eq '^[0-9a-f]{64}  [^[:space:]]+$' "$sidecar"
  IFS=' ' read -r digest covered extra < "$sidecar"
  [[ "$digest" =~ ^[0-9a-f]{64}$ ]]
  test "$covered" = "$(basename "$target")"
  test -z "${extra:-}"
  (cd "$(dirname "$target")" && \
    sha256sum --strict --status -c "$(basename "$sidecar")")
  chmod 0444 -- "$target" "$sidecar"
  authenticate_sidecar "$target"
}

expect_rc() {
  local observed="$1" expected="$2" label="$3"
  if (( observed != expected )); then
    printf 'NO_GO %s observed_rc=%s expected_rc=%s\n' "$label" "$observed" "$expected" >&2
    return 91
  fi
}

seal_diagnostic_log() {
  local target="$1"
  chmod 0444 "$target"
  (cd "$(dirname "$target")" && sha256sum "$(basename "$target")" >"$(basename "$target").sha256")
  chmod 0444 "$target.sha256"
}
```

驗證：`$ATTEMPT` 必須是本次唯一新目錄；任何 resume 都使用新的 attempt，不重導向到舊 log。

## 2. Fresh／resume 狀態機

只允許三種入口：

1. `RUN_ROOT` 不存在：fresh bootstrap。
2. PRE 已完成、contract 尚未 freeze：重新驗 PRE 與 live exact-11 後完成 freeze；不重跑 PRE。
3. contract 已存在：只執行 frozen verify-only，然後由已通過的 pilot/checkpoint/terminal receipts 推導下一步；不 re-freeze。

若存在 open/orphan batch start 或 resource gate 而沒有對應 completed checkpoint，必須 fail-closed；同一 evidence root 不得把舊 launch gate拿來重跑未完成 batch。

### 2.1 Live exact-11 allowlist（只在 freeze 前使用）

```bash
verify_live_exact11() {
  sha256sum -c <<'SHA256_ALLOWLIST'
5f910c836fd5ebaa3d5c933cfa6f0a97e36c7e4b4c38f0b206a4343e5e2d913a  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py
e378f7186b262c3629c45eb388204e65a302acf3b9f4f95ffd0be835c47862de  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py
cf016b9a046c214bbefb6a4b2509955910710fce73d3186dce27b666d5c40fc4  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py
b28d494563bea70cbfea7b13853d99c0b3aeafed7ed4f0c0a6fff3b92c84f0a9  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py
1def0de1952d127d8d013820ac7b0eabe33d6f1a66fd8c6d47a1985b14a32f77  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
66bb175404c207ef320f213c650bb10c6d5fcf3c84cbc40b8ca25e68604da767  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_ranking.py
9d15ce2bf15af5cc2c4c690cd7718b131108fd8e3946f6a72da40487b06f1578  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m2_single_task_pilot.py
87fc3ddde1052cd32b64588aa3e483faa699474047bde52e1332a07f779e5558  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py
ac4cb4107d3b3a012ee222ab2bcfb57879197e17948ab3ebc1f0c6516c56a0d7  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m2_input_identity.py
c734c1ed2142f6e2baed0155320cdaba8925e2d6b76965e60bd9e42e1bf4f7f1  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py
47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f  docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_manifest_v3.schema.json
SHA256_ALLOWLIST
}
```

### 2.2 Fresh bootstrap：只掃 process/disk，不讀 BAM

```bash
if [[ ! -e "$RUN_ROOT" ]]; then
  verify_live_exact11
  set +e
  "$PYTHON" "$LIVE_FULL_EXTRACT" \
    --manifest "$ORIGIN_MANIFEST" \
    --output-root "$RUN_ROOT/.bootstrap_probe" \
    --workers 2 --samtools-threads 1 \
    --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 \
    --task-timeout-seconds 28800 --task-timeout-grace-seconds 300 \
    --preflight-only \
    >"$ATTEMPT/bootstrap_preflight.stdout.json" \
    2>"$ATTEMPT/bootstrap_preflight.stderr.log"
  rc=$?
  set -e
  seal_diagnostic_log "$ATTEMPT/bootstrap_preflight.stdout.json"
  seal_diagnostic_log "$ATTEMPT/bootstrap_preflight.stderr.log"
  expect_rc "$rc" 0 bootstrap_preflight
  jq -e '
    .resource_gate_pass == true and
    .active_conflict_process_count == 0 and
    .active_conflict_root_count == 0 and
    .resource_gate_preview.disk_pass == true and
    .resource_gate_preview.required_reserve_bytes == 322122547200
  ' "$ATTEMPT/bootstrap_preflight.stdout.json"
  test ! -e "$RUN_ROOT"

  mkdir -p "$CONTROL"
  "$PYTHON" "$LIVE_INPUT_VERIFY" \
    --manifest "$ORIGIN_MANIFEST" --output "$PRE" \
    >"$ATTEMPT/pre_input_identity.stdout.json" \
    2>"$ATTEMPT/pre_input_identity.stderr.log"
  seal_direct_receipt "$PRE"
fi

authenticate_sidecar "$PRE"
jq -e '
  .mode == "PRE" and .all_pass == true and
  .validation_evidence_eligible == true and
  .artifact_audit.n_artifacts == 42 and
  .artifact_audit.n_unique_resolved_files == 42 and
  .artifact_audit.n_storage_identity_v1 == 7 and
  .artifact_audit.n_full_sha256 == 35 and
  .artifact_audit.n_sampled_bam_chunks == 21 and
  .artifact_audit.n_mismatches == 0
' "$PRE"
```

目前既有 all-sSNV job仍在執行時，此步預期 exit 2；正確動作是等待它自然結束，不能 kill、不能 `--ignore-resource-gate`、不能建立 formal root繞過。

## 3. Freeze exact-11 與每次 resume 的 deep verify

```bash
if [[ ! -e "$CONTRACT" ]]; then
  verify_live_exact11
  "$PYTHON" "$LIVE_FREEZER" freeze \
    --manifest "$ORIGIN_MANIFEST" \
    --pre-input-identity-receipt "$PRE" \
    --output-contract-root "$CONTRACT" \
    >"$ATTEMPT/release_freeze.stdout.json" \
    2>"$ATTEMPT/release_freeze.stderr.log"
fi

RELEASE="$CONTRACT/m2_run_manifest.json"
FROZEN_ROOT="$CONTRACT/code_snapshot/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts"
FROZEN_MANIFEST="$CONTRACT/input_contract/canonical_manifest.json"
FROZEN_PRE="$CONTRACT/input_contract/$(basename "$PRE")"
FROZEN_EXTRACTOR="$FROZEN_ROOT/extract_lossless_read_linkage.py"
FROZEN_RANKER="$FROZEN_ROOT/build_m2_patterns_and_rank.py"
FROZEN_FULL_EXTRACT="$FROZEN_ROOT/run_full_m2_extraction.py"
FROZEN_FULL_RANK="$FROZEN_ROOT/run_full_m2_ranking.py"
FROZEN_PILOT_VERIFY="$FROZEN_ROOT/verify_m2_single_task_pilot.py"
FROZEN_INPUT_VERIFY="$FROZEN_ROOT/verify_m2_input_identity.py"
FROZEN_FREEZER="$FROZEN_ROOT/freeze_m2_release_contract.py"
FROZEN_FINAL_VERIFY="$FROZEN_ROOT/verify_full_m2_receipts.py"

authenticate_sidecar "$RELEASE"
test "$(stat -c %a "$CONTRACT")" = 555
jq -e '
  .all_pass == true and .validation_evidence_eligible == true and
  .scope.technical_datasets == 7 and .scope.biological_samples == 6 and
  .scope.tasks == 154 and .source_snapshot.n_files == 11 and
  .parameters.ranking.conditional_candidate_ranking_bootstrap_replicates == 20 and
  .parameters.ranking.conditional_candidate_ranking_bootstrap_seed == 20260716
' "$RELEASE"

release_verify() {
  local label="$1"
  local output="$ATTEMPT/release_verify_${label}.json"
  "$PYTHON" "$FROZEN_FREEZER" verify-only \
    --contract-root "$CONTRACT" --output "$output" \
    >"$ATTEMPT/release_verify_${label}.stdout.json" \
    2>"$ATTEMPT/release_verify_${label}.stderr.log"
  authenticate_sidecar "$output"
  jq -e '.all_pass == true and .validation_evidence_eligible == true and .snapshot.n_files == 11' "$output"
}

release_verify resume_entry
```

## 4. 雙 direct pilots：每個 pilot 三份 resource gates

每個 pilot 固定建立：

- `resource_gates/extraction.json`
- `resource_gates/ranking_bootstrap0.json`
- `resource_gates/ranking_bootstrap20.json`

兩個 pilots 合計 6 份。B0 與 B20 的兩份 verifier各自綁定同一 pilot 的 extraction gate與對應 ranking gate。

### 4.1 Gate capture 與 direct stage 函式

```bash
capture_pilot_gate() {
  local pilot="$1" dataset="$2" chrom="$3" stage="$4" label="$5" stage_out="$6"
  local gate="$pilot/resource_gates/$label.json"
  test ! -e "$gate" && test ! -e "$gate.sha256"
  "$PYTHON" "$FROZEN_FULL_EXTRACT" \
    --manifest "$FROZEN_MANIFEST" \
    --release-contract-manifest "$RELEASE" \
    --output-root "$stage_out" \
    --workers 2 --samtools-threads 1 \
    --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 \
    --max-new-tasks 8 \
    --task-timeout-seconds 28800 --task-timeout-grace-seconds 300 \
    --preflight-only \
    --preflight-receipt-output "$gate" \
    --preflight-pilot-stage "$stage" \
    --preflight-pilot-label "$label" \
    --preflight-pilot-dataset "$dataset" \
    --preflight-pilot-chrom "$chrom" \
    >"$ATTEMPT/${dataset}_${chrom}_${label}_gate.stdout.json" \
    2>"$ATTEMPT/${dataset}_${chrom}_${label}_gate.stderr.log"
  authenticate_sidecar "$gate"
  jq -e '
    .gate_scope == "pilot" and .checks.all_pass == true and
    .process_snapshot.process_count == 0 and
    .filesystem_snapshot.disk_pass == true and
    .filesystem_snapshot.required_reserve_bytes == 322122547200
  ' "$gate"
}

run_pilot_extraction() {
  local pilot="$1" dataset="$2" chrom="$3"
  local gate="$pilot/resource_gates/extraction.json"
  local receipt="$pilot/extraction/receipt.json"
  if [[ -e "$receipt" ]]; then
    authenticate_sidecar "$gate"
    authenticate_sidecar "$receipt"
    test -f "$pilot/extraction.time.txt" && test ! -L "$pilot/extraction.time.txt"
    jq -e '.schema_version == "1.2.0" and .all_pass == true' "$receipt"
    return
  fi
  test ! -e "$gate"
  mkdir -p "$pilot/resource_gates"
  capture_pilot_gate "$pilot" "$dataset" "$chrom" extraction extraction "$pilot/extraction"
  set +e
  /usr/bin/time -v -o "$pilot/extraction.time.txt" \
    timeout --signal=TERM --kill-after=5m 8h \
    nice -n 10 ionice -c2 -n7 \
    "$PYTHON" "$FROZEN_EXTRACTOR" \
      --manifest "$FROZEN_MANIFEST" --sample "$dataset" --chrom "$chrom" \
      --output-dir "$pilot/extraction" \
      --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 \
      --samtools-threads 1 \
      --require-existing-empty-output-dir \
      >"$ATTEMPT/${dataset}_${chrom}_extraction.stdout.json" \
      2>"$ATTEMPT/${dataset}_${chrom}_extraction.stderr.log"
  rc=$?
  set -e
  expect_rc "$rc" 0 "$dataset/$chrom direct extraction"
  seal_direct_receipt "$receipt"
  chmod 0444 "$pilot/extraction.time.txt"
  authenticate_sidecar "$receipt"
  jq -e '.schema_version == "1.2.0" and .all_pass == true' "$receipt"
}

run_pilot_ranking() {
  local pilot="$1" dataset="$2" chrom="$3" label="$4" bootstrap="$5" verify_output="$6"
  local gate="$pilot/resource_gates/$label.json"
  local receipt="$pilot/$label/receipt.json"
  if [[ ! -e "$receipt" ]]; then
    test ! -e "$gate"
    capture_pilot_gate "$pilot" "$dataset" "$chrom" ranking "$label" "$pilot/$label"
    set +e
    /usr/bin/time -v -o "$pilot/$label.time.txt" \
      timeout --signal=TERM --kill-after=5m 8h \
      nice -n 10 ionice -c2 -n7 \
      "$PYTHON" "$FROZEN_RANKER" \
        --input-dir "$pilot/extraction" --output-dir "$pilot/$label" \
        --threshold 1 --threshold 2 --threshold 3 --threshold 5 \
        --component-basis PS_HP1 --component-basis PS_HP2 \
        --hp-family 1 --hp-family 2 \
        --structural-exact-pattern-minread-grid 1,2,3,5 \
        --primary-structural-exact-pattern-minread 3 \
        --exact-k-max 12 --max-vertex-sets 256 --solver-time-limit-seconds 30 \
        --fixed-error-grid 0.005,0.01,0.02,0.05 \
        --minimum-bq-error-rate 0.000001 --maximum-bq-error-rate 0.25 \
        --conditional-candidate-ranking-bootstrap-replicates "$bootstrap" \
        --conditional-candidate-ranking-bootstrap-seed 20260716 \
        --tie-tolerance 0.000001 \
        --require-existing-empty-output-dir \
        >"$ATTEMPT/${dataset}_${chrom}_${label}.stdout.json" \
        2>"$ATTEMPT/${dataset}_${chrom}_${label}.stderr.log"
    rc=$?
    set -e
    expect_rc "$rc" 0 "$dataset/$chrom $label"
    seal_direct_receipt "$receipt"
    chmod 0444 "$pilot/$label.time.txt"
  fi
  authenticate_sidecar "$gate"
  authenticate_sidecar "$receipt"
  test -f "$pilot/$label.time.txt" && test ! -L "$pilot/$label.time.txt"
  jq -e '.schema_version == "2.0.0" and .all_pass == true' "$receipt"

  if [[ ! -e "$verify_output" ]]; then
    "$PYTHON" "$FROZEN_PILOT_VERIFY" \
      --pilot-root "$pilot" --dataset "$dataset" --chrom "$chrom" \
      --ranking-dirname "$label" \
      --expected-bootstrap-replicates "$bootstrap" \
      --extraction-resource-gate-receipt "$pilot/resource_gates/extraction.json" \
      --resource-gate-receipt "$gate" \
      --output "$verify_output" \
      >"$ATTEMPT/${dataset}_${chrom}_${label}_verify.stdout.json" \
      2>"$ATTEMPT/${dataset}_${chrom}_${label}_verify.stderr.log"
  fi
  authenticate_sidecar "$verify_output"
  jq -e '
    .all_pass == true and .release_gate.verdict == "GO" and
    .checks.extraction_and_ranking_resource_gates_authenticated == true
  ' "$verify_output"
}

validate_pilot_resume_state() {
  local pilot="$1"
  [[ -e "$pilot" ]] || return 0
  test -d "$pilot" && test ! -L "$pilot"
  local child name
  for child in "$pilot"/*; do
    [[ -e "$child" ]] || continue
    name="$(basename "$child")"
    case "$name" in
      resource_gates|extraction|extraction.time.txt|ranking_bootstrap0|ranking_bootstrap0.time.txt|ranking_bootstrap20|ranking_bootstrap20.time.txt|pilot_gate_verification_receipt.json|pilot_gate_verification_receipt.json.sha256|pilot_gate_verification_receipt.ranking_bootstrap20.json|pilot_gate_verification_receipt.ranking_bootstrap20.json.sha256) ;;
      *) printf 'NO_GO unexpected pilot top-level entry: %s\n' "$child" >&2; return 92 ;;
    esac
  done
  # Existing pilot root is a resume only after the canonical extraction stage
  # completed. Empty/preseed/open-stage roots are not accepted.
  test -f "$pilot/resource_gates/extraction.json"
  test -f "$pilot/resource_gates/extraction.json.sha256"
  test -f "$pilot/extraction/receipt.json"
  test -f "$pilot/extraction/receipt.json.sha256"
  test -f "$pilot/extraction.time.txt"
  authenticate_sidecar "$pilot/resource_gates/extraction.json"
  authenticate_sidecar "$pilot/extraction/receipt.json"
  local b0_complete=0 b0_go=0 b20_complete=0 b20_go=0 stage_present=0 verify_present=0
  stage_present=0
  for child in \
    "$pilot/resource_gates/ranking_bootstrap0.json" \
    "$pilot/resource_gates/ranking_bootstrap0.json.sha256" \
    "$pilot/ranking_bootstrap0" \
    "$pilot/ranking_bootstrap0.time.txt"; do
    [[ -e "$child" ]] && stage_present=1
  done
  if [[ "$stage_present" -eq 1 ]]; then
    test -f "$pilot/resource_gates/ranking_bootstrap0.json"
    test -f "$pilot/resource_gates/ranking_bootstrap0.json.sha256"
    test -d "$pilot/ranking_bootstrap0" && test ! -L "$pilot/ranking_bootstrap0"
    test -f "$pilot/ranking_bootstrap0/receipt.json"
    test -f "$pilot/ranking_bootstrap0/receipt.json.sha256"
    test -f "$pilot/ranking_bootstrap0.time.txt"
    authenticate_sidecar "$pilot/resource_gates/ranking_bootstrap0.json"
    authenticate_sidecar "$pilot/ranking_bootstrap0/receipt.json"
    b0_complete=1
  fi
  verify_present=0
  [[ -e "$pilot/pilot_gate_verification_receipt.json" || \
     -e "$pilot/pilot_gate_verification_receipt.json.sha256" ]] && verify_present=1
  if [[ "$verify_present" -eq 1 ]]; then
    test "$b0_complete" -eq 1
    test -f "$pilot/pilot_gate_verification_receipt.json"
    test -f "$pilot/pilot_gate_verification_receipt.json.sha256"
    authenticate_sidecar "$pilot/pilot_gate_verification_receipt.json"
    b0_go=1
  fi

  stage_present=0
  for child in \
    "$pilot/resource_gates/ranking_bootstrap20.json" \
    "$pilot/resource_gates/ranking_bootstrap20.json.sha256" \
    "$pilot/ranking_bootstrap20" \
    "$pilot/ranking_bootstrap20.time.txt"; do
    [[ -e "$child" ]] && stage_present=1
  done
  if [[ "$stage_present" -eq 1 ]]; then
    test "$b0_go" -eq 1
    test -f "$pilot/resource_gates/ranking_bootstrap20.json"
    test -f "$pilot/resource_gates/ranking_bootstrap20.json.sha256"
    test -d "$pilot/ranking_bootstrap20" && test ! -L "$pilot/ranking_bootstrap20"
    test -f "$pilot/ranking_bootstrap20/receipt.json"
    test -f "$pilot/ranking_bootstrap20/receipt.json.sha256"
    test -f "$pilot/ranking_bootstrap20.time.txt"
    authenticate_sidecar "$pilot/resource_gates/ranking_bootstrap20.json"
    authenticate_sidecar "$pilot/ranking_bootstrap20/receipt.json"
    b20_complete=1
  fi
  verify_present=0
  [[ -e "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json" || \
     -e "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json.sha256" ]] && verify_present=1
  if [[ "$verify_present" -eq 1 ]]; then
    test "$b20_complete" -eq 1
    test -f "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json"
    test -f "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json.sha256"
    authenticate_sidecar "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json"
    b20_go=1
  fi
  if [[ -d "$pilot/resource_gates" ]]; then
    for child in "$pilot/resource_gates"/*; do
      [[ -e "$child" ]] || continue
      name="$(basename "$child")"
      case "$name" in
        extraction.json|extraction.json.sha256|ranking_bootstrap0.json|ranking_bootstrap0.json.sha256|ranking_bootstrap20.json|ranking_bootstrap20.json.sha256) ;;
        *) printf 'NO_GO unexpected pilot gate entry: %s\n' "$child" >&2; return 93 ;;
      esac
    done
  fi
  # Legal monotonic prefixes are extraction; +B0; +B0 GO; +B20; +B20 GO.
  # The dependency checks above reject every skipped or preseeded later state.
}

run_direct_pilot() {
  local pilot="$1" dataset="$2" chrom="$3"
  if [[ -e "$pilot" ]]; then
    validate_pilot_resume_state "$pilot"
  fi
  run_pilot_extraction "$pilot" "$dataset" "$chrom"
  run_pilot_ranking "$pilot" "$dataset" "$chrom" ranking_bootstrap0 0 \
    "$pilot/pilot_gate_verification_receipt.json"
  run_pilot_ranking "$pilot" "$dataset" "$chrom" ranking_bootstrap20 20 \
    "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json"
}

mkdir -p "$PILOTS"
run_direct_pilot "$PILOT_DORADO" HCC1395_DORADO chr6
run_direct_pilot "$PILOT_H2009" H2009 chr2
```

若 gate 已存在但對應 stage receipt 不存在，函式會因 `test ! -e "$gate"` 失敗。這是刻意的 fail-closed：舊 launch snapshot不能授權新的重跑。

### 4.2 六 gate／四 GO／B20 資源外推

```bash
for pilot in "$PILOT_DORADO" "$PILOT_H2009"; do
  for label in extraction ranking_bootstrap0 ranking_bootstrap20; do
    authenticate_sidecar "$pilot/resource_gates/$label.json"
  done
  authenticate_sidecar "$pilot/pilot_gate_verification_receipt.json"
  authenticate_sidecar "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json"
done

if [[ ! -e "$CONTROL/pilot_resource_extrapolation.json" ]]; then
  "$PYTHON" - "$PILOT_DORADO" "$PILOT_H2009" "$OUT_BASE" \
    "$CONTROL/pilot_resource_extrapolation.json" <<'PY'
import hashlib, json, math, os, pathlib, shutil, sys
dorado, h2009, out_base, output = map(pathlib.Path, sys.argv[1:])
def tree_bytes(root):
    return sum(p.stat().st_size for p in root.rglob('*') if p.is_file() and not p.is_symlink())
dorado_extraction = tree_bytes(dorado / 'extraction')
h2009_ranking_b20 = tree_bytes(h2009 / 'ranking_bootstrap20')
point = 19.2918 * dorado_extraction + 55.4308 * h2009_ranking_b20
projected = math.ceil(1.20 * point)
free_now = shutil.disk_usage(out_base).free
reserve = 300 * 1024**3
payload = {
    'schema_name': 'intersubmod.m2_dual_pilot_resource_extrapolation',
    'schema_version': '1.0.0',
    'method': '1.20*(19.2918*D_E_DORADO_chr6 + 55.4308*D_R_H2009_chr2_B20)',
    'dorado_extraction_bytes': dorado_extraction,
    'h2009_ranking_b20_bytes': h2009_ranking_b20,
    'projected_additional_full_bytes_with_20pct_margin': projected,
    'available_bytes_after_pilots': free_now,
    'required_reserve_bytes': reserve,
    'projected_remaining_bytes': free_now - projected,
    'all_pass': free_now - projected >= reserve,
    'warning': 'operational envelope, not a statistical confidence interval; each batch has its own gate',
}
raw = (json.dumps(payload, ensure_ascii=False, indent=2) + '\n').encode()
fd = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
with os.fdopen(fd, 'wb') as handle:
    handle.write(raw); handle.flush(); os.fsync(handle.fileno()); os.fchmod(handle.fileno(), 0o444)
sidecar = output.with_name(output.name + '.sha256')
side_raw = f"{hashlib.sha256(raw).hexdigest()}  {output.name}\n".encode('ascii')
fd = os.open(sidecar, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
with os.fdopen(fd, 'wb') as handle:
    handle.write(side_raw); handle.flush(); os.fsync(handle.fileno()); os.fchmod(handle.fileno(), 0o444)
if not payload['all_pass']:
    raise SystemExit(2)
PY
fi
authenticate_sidecar "$CONTROL/pilot_resource_extrapolation.json"
jq -e '.all_pass == true and .projected_remaining_bytes >= .required_reserve_bytes' \
  "$CONTROL/pilot_resource_extrapolation.json"

release_verify post_pilots
```

## 5. Full extraction：154 tasks、session＋11 batch gates

合法累積數固定為 `8,24,40,56,72,88,104,120,136,152,154`。第1批 `max-new=8`，其餘 `16`；前10批 checkpoint success exit=3，第11批 terminal exit=0。

```bash
targets=(8 24 40 56 72 88 104 120 136 152 154)

for index in "${!targets[@]}"; do
  target="${targets[$index]}"
  batch=$((index + 1))
  max_new=16
  expected_rc=3
  [[ "$batch" -eq 1 ]] && max_new=8
  [[ "$target" -eq 154 ]] && expected_rc=0
  if [[ "$target" -eq 154 ]]; then
    receipt="$EXTRACTION/full_extraction_receipt.json"
  else
    receipt="$EXTRACTION/checkpoints/checkpoint_$(printf '%03d' "$target")_of_154.json"
  fi

  if [[ ! -e "$receipt" ]]; then
    set +e
    /usr/bin/time -v -o "$ATTEMPT/extraction_batch_$(printf '%02d' "$batch").time.txt" \
      nice -n 10 ionice -c2 -n7 \
      "$PYTHON" "$FROZEN_FULL_EXTRACT" \
        --manifest "$FROZEN_MANIFEST" \
        --release-contract-manifest "$RELEASE" \
        --output-root "$EXTRACTION" \
        --workers 2 --samtools-threads 1 \
        --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 \
        --max-new-tasks "$max_new" \
        --task-timeout-seconds 28800 --task-timeout-grace-seconds 300 \
        >"$ATTEMPT/extraction_batch_$(printf '%02d' "$batch").stdout.log" \
        2>"$ATTEMPT/extraction_batch_$(printf '%02d' "$batch").stderr.log"
    rc=$?
    set -e
    expect_rc "$rc" "$expected_rc" "full extraction batch $batch"
  fi

  authenticate_sidecar "$EXTRACTION/_orchestration/resource_gates/session.json"
  authenticate_sidecar "$EXTRACTION/_orchestration/resource_gates/batch_$(printf '%03d' "$batch").json"
  jq -e '.stage == "extraction" and .checks.all_pass == true' \
    "$EXTRACTION/_orchestration/resource_gates/batch_$(printf '%03d' "$batch").json"
  authenticate_sidecar "$receipt"
  if [[ "$target" -lt 154 ]]; then
    jq -e --argjson n "$target" \
      '.checkpoint_complete == true and .passing_tasks == $n and .remaining_tasks == (154-$n)' "$receipt"
  else
    jq -e '.all_pass == true and .passing_tasks == 154 and .remaining_tasks == 0' "$receipt"
  fi
done

test "$(find "$EXTRACTION/_orchestration/resource_gates" -maxdepth 1 -type f -name '*.json' | wc -l)" -eq 12
release_verify post_extraction
```

Resume 只跳過已通過且 sidecar 正確的 receipts。Runner在下一個 batch開始前會驗完整 session→gate→batch-start→grant→completion→checkpoint chain；open/orphan/extra gate 或 task會拒絕。

## 6. Full ranking：154 tasks、session＋11 batch gates

```bash
for index in "${!targets[@]}"; do
  target="${targets[$index]}"
  batch=$((index + 1))
  max_new=16
  expected_rc=3
  [[ "$batch" -eq 1 ]] && max_new=8
  [[ "$target" -eq 154 ]] && expected_rc=0
  if [[ "$target" -eq 154 ]]; then
    receipt="$RANKING/full_ranking_receipt.json"
  else
    receipt="$RANKING/checkpoints/checkpoint_$(printf '%03d' "$target")_of_154.json"
  fi

  if [[ ! -e "$receipt" ]]; then
    set +e
    /usr/bin/time -v -o "$ATTEMPT/ranking_batch_$(printf '%02d' "$batch").time.txt" \
      nice -n 10 ionice -c2 -n7 \
      "$PYTHON" "$FROZEN_FULL_RANK" \
        --extraction-root "$EXTRACTION" \
        --release-contract-manifest "$RELEASE" \
        --output-root "$RANKING" --workers 2 \
        --thresholds 1,2,3,5 --component-bases PS_HP1,PS_HP2 --hp-families 1,2 \
        --structural-exact-pattern-minread-grid 1,2,3,5 \
        --primary-structural-exact-pattern-minread 3 \
        --exact-k-max 12 --max-vertex-sets 256 --solver-time-limit-seconds 30 \
        --fixed-error-grid 0.005,0.01,0.02,0.05 \
        --minimum-bq-error-rate 0.000001 --maximum-bq-error-rate 0.25 \
        --conditional-candidate-ranking-bootstrap-replicates 20 \
        --conditional-candidate-ranking-bootstrap-seed 20260716 \
        --tie-tolerance 0.000001 --max-new-tasks "$max_new" \
        --task-timeout-seconds 28800 --task-timeout-grace-seconds 300 \
        >"$ATTEMPT/ranking_batch_$(printf '%02d' "$batch").stdout.log" \
        2>"$ATTEMPT/ranking_batch_$(printf '%02d' "$batch").stderr.log"
    rc=$?
    set -e
    expect_rc "$rc" "$expected_rc" "full ranking batch $batch"
  fi

  authenticate_sidecar "$RANKING/_orchestration/resource_gates/session.json"
  authenticate_sidecar "$RANKING/_orchestration/resource_gates/batch_$(printf '%03d' "$batch").json"
  jq -e '.stage == "ranking" and .checks.all_pass == true' \
    "$RANKING/_orchestration/resource_gates/batch_$(printf '%03d' "$batch").json"
  authenticate_sidecar "$receipt"
  if [[ "$target" -lt 154 ]]; then
    jq -e --argjson n "$target" \
      '.checkpoint_complete == true and .passing_tasks == $n and .remaining_tasks == (154-$n)' "$receipt"
  else
    jq -e '.all_pass == true and .scope.n_results == 154' "$receipt"
  fi
done

test "$(find "$RANKING/_orchestration/resource_gates" -maxdepth 1 -type f -name '*.json' | wc -l)" -eq 12
release_verify post_ranking
```

## 7. POST_COMPARE 與 independent full verifier

```bash
if [[ ! -e "$POST" ]]; then
  "$PYTHON" "$FROZEN_INPUT_VERIFY" \
    --manifest "$FROZEN_MANIFEST" --compare-to "$FROZEN_PRE" --output "$POST" \
    >"$ATTEMPT/post_input_identity.stdout.json" \
    2>"$ATTEMPT/post_input_identity.stderr.log"
  seal_direct_receipt "$POST"
fi
authenticate_sidecar "$POST"
jq -e '
  .mode == "POST_COMPARE" and .all_pass == true and
  .validation_evidence_eligible == true and
  .compare_to.exact_snapshot_equal == true and
  .artifact_audit.n_artifacts == 42 and .artifact_audit.n_mismatches == 0
' "$POST"

if [[ ! -e "$FINAL_VERIFY" ]]; then
  /usr/bin/time -v -o "$ATTEMPT/full_m2_independent_verification.time.txt" \
    nice -n 10 ionice -c2 -n7 \
    "$PYTHON" "$FROZEN_FINAL_VERIFY" \
      --extraction-root "$EXTRACTION" --ranking-root "$RANKING" \
      --release-contract-manifest "$RELEASE" \
      --post-input-identity-receipt "$POST" --output "$FINAL_VERIFY" \
      >"$ATTEMPT/full_m2_independent_verification.stdout.json" \
      2>"$ATTEMPT/full_m2_independent_verification.stderr.log"
fi
authenticate_sidecar "$FINAL_VERIFY"
jq -e '
  .all_pass == true and ([.checks[]] | all) and
  .extraction.orchestration.n_authenticated_resource_gates == 12 and
  .ranking.orchestration.n_authenticated_resource_gates == 12
' "$FINAL_VERIFY"
```

Final verifier必須獨立讀取 154＋154 child receipts/tables、重算 R/A/X＋BQ mixture likelihood與守恆、驗兩 stages合計24 resource gates；若欄位實際路徑不同，以 verifier schema為準修正 jq，但不可降低到只看 `all_pass`。

## 8. 最終數字 summary 與 readout

```bash
test "$(sha256sum "$LIVE_NUMERIC" | awk '{print $1}')" = 8952ccb17a0e2514621ef99110b9bc589f084398855d95e19eee627f40d6a4cd
if [[ ! -e "$NUMERIC" ]]; then
  "$PYTHON" "$LIVE_NUMERIC" \
    --extraction-root "$EXTRACTION" --ranking-root "$RANKING" \
    --final-verification-receipt "$FINAL_VERIFY" --output "$NUMERIC" \
    >"$ATTEMPT/final_numeric_summary.stdout.json" \
    2>"$ATTEMPT/final_numeric_summary.stderr.log"
fi
authenticate_sidecar "$NUMERIC"
jq -e '
  .schema_name == "intersubmod.m2_final_numeric_summary" and
  .all_pass == true and
  .scope.scope_is_canonical_full_7_dataset_chr1_22 == true and
  .scope.expected_tasks_per_stage == 154 and
  .checks.tree_vertex_partition_conserves_and_binds_authenticated_T_V == true and
  (.ranking.overall_candidate_table.tree_vertex_partition.counts.T_EQ_1_V_EQ_1 +
   .ranking.overall_candidate_table.tree_vertex_partition.counts.T_GT_1_V_EQ_1 +
   .ranking.overall_candidate_table.tree_vertex_partition.counts.T_GT_1_V_GT_1) ==
   .ranking.overall_candidate_table.tree_vertex_partition.denominator
' "$NUMERIC"

jq '{
  scope, definitions, primary_parameters,
  extraction_overall_counts: .extraction.overall_counts,
  extraction_by_basis_threshold: .extraction.overall_component_by_linkage_basis_threshold,
  ranking_input_call_funnel: .ranking.overall_input_call_funnel,
  ranking_by_HP_basis_threshold: .ranking.overall_by_HP_basis_and_bridge_threshold,
  candidate_table: .ranking.overall_candidate_table,
  by_dataset: .ranking.by_dataset,
  unsupported_or_nonidentifiable, checks, all_pass
}' "$NUMERIC" >"$ATTEMPT/final_key_numbers.readout.json"
seal_diagnostic_log "$ATTEMPT/final_key_numbers.readout.json"
```

教授版必須從此summary列出絕對數、total denominator%、relative denominator%，至少涵蓋：component k／effective-k／partial funnel／h*／T與V／三個T-V buckets／unique-tied-abstain／coarse topology／parent-edge與exact-topology uniqueness，以及 overall、7 datasets、HP1/HP2、threshold 1/2/3/5。

## 9. Presentation snapshot（3 programs＋12 evidence roles）

Presentation code SHA是第二層 allowlist；任一變動都必須建立新 snapshot，不得沿用舊 manifest。

```bash
sha256sum -c <<'PRESENTATION_SHA'
8952ccb17a0e2514621ef99110b9bc589f084398855d95e19eee627f40d6a4cd  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py
613d57a8ba956f625221cf1231625424e3b011c291858fd92a3bd6e7a8467e18  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py
cf66a8cf9b5fcd408ac2a487708617016b47bcd3551ac79881e1869e963d4e16  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py
639898a882b5d72f9c58ad95701a5afdfa22c0ad37f1aebc8748320d2efbd625  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py
PRESENTATION_SHA

CANONICAL_JSON="$REPO/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json"
FUNNEL="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/current_funnel_receipt.json"
FUNNEL_VERIFY="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/independent_verification.json"
M0="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_receipt.json"
M0_VERIFY="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json"
SYMBOLIC_PILOT="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot_v3_identity_bound/pilot_receipt.json"
METHOD_AUDIT="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md"
LITERATURE_AUDIT="$REPO/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_方法學原始文獻與主張邊界_01.md"

if [[ ! -e "$PRESENTATION" ]]; then
  "$PYTHON" "$LIVE_PRESENTATION_FREEZER" freeze \
    --output-root "$PRESENTATION" \
    --canonical-json "$CANONICAL_JSON" \
    --funnel-receipt "$FUNNEL" \
    --funnel-verification-receipt "$FUNNEL_VERIFY" \
    --m0-receipt "$M0" --m0-verification-receipt "$M0_VERIFY" \
    --pilot-receipt "$SYMBOLIC_PILOT" \
    --method-audit "$METHOD_AUDIT" --literature-audit "$LITERATURE_AUDIT" \
    --m2-extraction-receipt "$EXTRACTION/full_extraction_receipt.json" \
    --m2-ranking-receipt "$RANKING/full_ranking_receipt.json" \
    --m2-verification-receipt "$FINAL_VERIFY" \
    --final-numeric-summary "$NUMERIC" \
    >"$ATTEMPT/presentation_freeze.stdout.json" \
    2>"$ATTEMPT/presentation_freeze.stderr.log"
fi

PRESENTATION_BUILDER="$PRESENTATION/code_snapshot/build_validated_html_report.py"
PRESENTATION_QA="$PRESENTATION/code_snapshot/qa_validated_html_report.py"
PRESENTATION_NUMERIC="$PRESENTATION/code_snapshot/build_final_numeric_summary.py"

authenticate_sidecar "$PRESENTATION/presentation_snapshot_manifest.json"
"$PYTHON" "$LIVE_PRESENTATION_FREEZER" verify-only \
  --snapshot-root "$PRESENTATION" \
  >"$ATTEMPT/presentation_verify_pre_html.json" \
  2>"$ATTEMPT/presentation_verify_pre_html.stderr.log"
jq -e '.all_pass == true and .programs_verified == 3 and .evidence_inputs_verified == 12' \
  "$ATTEMPT/presentation_verify_pre_html.json"
test "$(stat -c %a "$PRESENTATION")" = 555
test "$(sha256sum "$PRESENTATION_NUMERIC" | awk '{print $1}')" = 8952ccb17a0e2514621ef99110b9bc589f084398855d95e19eee627f40d6a4cd
```

## 10. FINAL HTML 與同一 frozen QA producer 的 run＋verify-only

```bash
if [[ ! -e "$HTML" ]]; then
  "$PYTHON" "$PRESENTATION_BUILDER" \
    --canonical-json "$CANONICAL_JSON" \
    --funnel-receipt "$FUNNEL" --funnel-verification-receipt "$FUNNEL_VERIFY" \
    --m0-receipt "$M0" --m0-verification-receipt "$M0_VERIFY" \
    --pilot-receipt "$SYMBOLIC_PILOT" \
    --method-audit "$METHOD_AUDIT" --literature-audit "$LITERATURE_AUDIT" \
    --m2-extraction-receipt "$EXTRACTION/full_extraction_receipt.json" \
    --m2-ranking-receipt "$RANKING/full_ranking_receipt.json" \
    --m2-verification-receipt "$FINAL_VERIFY" \
    --m2-numeric-summary "$NUMERIC" --output "$HTML" \
    >"$ATTEMPT/html_build.stdout.json" \
    2>"$ATTEMPT/html_build.stderr.log"
  jq -e '
    .generation_pass == true and .all_pass == true and .final_ready == true and
    .artifact_status == "FINAL_VALIDATION_EVIDENCE" and (.final_issues | length) == 0
  ' "$ATTEMPT/html_build.stdout.json"
  test ! -e "$HTML.sha256"
  (set -o noclobber; cd "$(dirname "$HTML")" && \
    sha256sum "$(basename "$HTML")" >"$(basename "$HTML").sha256")
  chmod 0444 "$HTML.sha256"
fi
(cd "$(dirname "$HTML")" && sha256sum -c "$(basename "$HTML").sha256")
test "$(stat -c %a "$HTML")" = 444
test "$(stat -c %h "$HTML")" = 1

if [[ ! -e "$HTML_QA" ]]; then
  "$PYTHON" "$PRESENTATION_QA" run \
    --html "$HTML" --output-dir "$HTML_QA" --expect-status final \
    >"$ATTEMPT/html_qa_run.stdout.json" \
    2>"$ATTEMPT/html_qa_run.stderr.log"
fi

"$PYTHON" "$PRESENTATION_QA" verify-only \
  --html "$HTML" --output-dir "$HTML_QA" --expect-status final \
  >"$ATTEMPT/html_qa_verify_only.stdout.json" \
  2>"$ATTEMPT/html_qa_verify_only.stderr.log"
jq -e '.all_pass == true and .outputs_verified == 3 and .exact_layout_verified == true' \
  "$ATTEMPT/html_qa_verify_only.stdout.json"
authenticate_sidecar "$HTML_QA/browser_qa_receipt.json"
test "$(stat -c %a "$HTML_QA")" = 555
```

QA receipt綁定 HTML、同一份 frozen QA producer、desktop/mobile/PDF的 path、size、SHA、mode與nlink。用 live QA script驗 frozen receipt必須失敗，不能把 producer mismatch略過。

## 11. Post-QA deep verification 與 terminal gate

```bash
"$PYTHON" "$LIVE_PRESENTATION_FREEZER" verify-only \
  --snapshot-root "$PRESENTATION" \
  >"$ATTEMPT/presentation_verify_post_qa.json" \
  2>"$ATTEMPT/presentation_verify_post_qa.stderr.log"
jq -e '.all_pass == true and .programs_verified == 3 and .evidence_inputs_verified == 12' \
  "$ATTEMPT/presentation_verify_post_qa.json"

release_verify post_qa
authenticate_sidecar "$EXTRACTION/full_extraction_receipt.json"
authenticate_sidecar "$RANKING/full_ranking_receipt.json"
authenticate_sidecar "$POST"
authenticate_sidecar "$FINAL_VERIFY"
authenticate_sidecar "$NUMERIC"
authenticate_sidecar "$PRESENTATION/presentation_snapshot_manifest.json"
authenticate_sidecar "$HTML_QA/browser_qa_receipt.json"
(cd "$(dirname "$HTML")" && sha256sum -c "$(basename "$HTML").sha256")

jq -e '.all_pass == true and ([.checks[]] | all)' "$FINAL_VERIFY"
jq -e '.all_pass == true' "$NUMERIC"
jq -e '.all_pass == true and ([.checks[]] | all)' "$HTML_QA/browser_qa_receipt.json"
```

只有上述全部 exit 0，才可把 HTML標為最終教授／論文說明證據。

## 12. Expected exit 與 resume 邊界

| Stage | Success exit | Canonical output |
|---|---:|---|
| bootstrap preflight | 0 | attempt-exclusive stdout/stderr＋SHA；formal root仍不存在 |
| PRE | 0 | `control/pre_input_identity_receipt.json{,.sha256}` |
| exact-11 freeze/deep verify | 0 | immutable release contract＋verify receipt |
| direct extraction/B0/B20 | 0 | child receipt；每pilot三 resource gates |
| four pilot verifiers | 0 | 四份 GO receipts，綁六 gates |
| extraction batch 1–10 / 11 | 3 / 0 | checkpoint 008…152 / terminal 154 |
| ranking batch 1–10 / 11 | 3 / 0 | checkpoint 008…152 / terminal 154 |
| POST / independent verifier | 0 / 0 | authenticated POST＋final verifier |
| numeric / presentation / HTML | 0 | authenticated numeric、immutable presentation、FINAL HTML |
| QA run / same-frozen verify-only | 0 / 0 | immutable 5-file QA artifact |

可以 resume：

- PRE完成但contract尚未freeze。
- pilot extraction完成後、B0或B20尚未開始。
- 任一 full batch checkpoint已完整發布後。
- terminal data完成但numeric／presentation／HTML／QA尚未開始。

不可在同root resume：

- 已存在 pilot gate但相應 child receipt不存在。
- 已存在 open/orphan full batch start、gate、grant或部分 completion而沒有 checkpoint。
- fixed receipt存在但sidecar／mode／nlink／semantic contract不符。
- frozen source、manifest、PRE/POST、output-root inode或release identity漂移。

這些情況應保留現場、停止並建立新的版本化 formal root；不可刪除、覆寫或手工補 receipt。

## 13. Partial-read 方法在最終報告中的正確說法

1. Pattern `p∈{R,A,X}^k` 有 `u` 個 X 時，full cube中概念上有 `2^u` 個 completions；observed-ALT active universe中的有效數可能較少。
2. Production不建立 `2^u` 個獨立 tree worlds，也不建立 reads間的 completion Cartesian product。每個達 exact-pattern minread 的distinct partial pattern形成一個 group constraint `N∩G(p)≠∅`。
3. 同一 unit內所有 groups必須聯合滿足；不是「逐 completion、任一成功就停止」。
4. Solver先求全域 minimum-extra h*，再以 no-good cuts列舉所有同 h* 的distinct vertex sets。只有 `complete=true` 才表示該層列完；h*+1以上目前不進 ranking。
5. Likelihood使用所有 informative molecule patterns與BQ/error，對 X／latent state做marginalization；同一批 reads衍生的VAF不再當獨立分數重複加入。
6. 相同 vertex set的不同 parent edges在目前read-pattern likelihood下同分；因此 `V=1,T>1` 代表state set唯一但parent edge不可辨識，不能寫成唯一真實拓撲。

## 14. Reader-test 必答題

在正式執行前，fresh reader應能只靠本Runbook回答：

1. 為何有7 datasets但只有6 biological samples？
2. fresh與resume怎麼分，哪些狀態不能續跑？
3. 雙 pilots為何是6 gates、4 GO receipts？
4. full extraction/ranking各為何是12 gates？
5. 前10批為何exit 3仍是成功？
6. 如何證明最後HTML來自通過的154＋154 evidence？
7. partial X是否被逐一建立成獨立樹？
8. h*、V、T各代表什麼，何者不能當clone數？
9. 為何同reads VAF不能再加一次分？
10. FINAL QA為何必須由同一份 frozen QA producer做run與verify-only？
