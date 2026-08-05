#!/usr/bin/env bash
set -euo pipefail

ROOT="/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
RESULTS="${ROOT}/results"
REVIEWS="${ROOT}/reviews"
ARCHIVE_PARENT="${ROOT}/audit_notes/failed_formal_runs"
ARCHIVE_NAME="20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
ARCHIVE="${ARCHIVE_PARENT}/${ARCHIVE_NAME}"
STAGE="${ARCHIVE_PARENT}/.${ARCHIVE_NAME}.staging.$$"

assert_file() {
    local path="$1"
    local expected_sha="$2"
    local expected_size="$3"
    local expected_mode="$4"
    local observed_sha observed_size observed_mode

    [[ -f "${path}" && ! -L "${path}" ]]
    observed_sha="$(sha256sum -- "${path}" | awk '{print $1}')"
    observed_size="$(stat -c '%s' -- "${path}")"
    observed_mode="$(stat -c '%a' -- "${path}")"
    [[ "${observed_sha}" == "${expected_sha}" ]]
    [[ "${observed_size}" == "${expected_size}" ]]
    [[ "${observed_mode}" == "${expected_mode}" ]]
}

[[ ! -e "${ARCHIVE}" && ! -e "${STAGE}" ]]
assert_file "${RESULTS}/tumor_ref_promotion_schema_recovery_authority.v11.bundle/authority.json" \
    "2fe51e9f3a746781324dfb70156d96481d59f0c673b8623f7da28df20cf924ad" 44364 444
assert_file "${RESULTS}/tumor_ref_promotion_schema_recovery_authority.v11.bundle/authority.ed25519.sig" \
    "8e95a60f23cd75bf4b8b39490c5142b7cf3fac04873beffa924f76cdc3fee69c" 64 444
assert_file "${RESULTS}/tumor_ref_promotion_schema_recovery_authority.v11.bundle/commit.json" \
    "f4900248e51d64e961454b433558fca8bb1d7f074fd257bc03bde73bbef493c8" 803 444
assert_file "${REVIEWS}/20260723_tumor_ref_schema_recovery_mendel.v11.json" \
    "ade00ba448525616bdc2b5953d3720a44f996d763ea498f594461dcfd2514007" 2713 444
assert_file "${REVIEWS}/20260723_tumor_ref_schema_recovery_nash.v11.json" \
    "ee23983fe4fbc5336ef063dc4effaffb08e621b392c7b4fa1f490a67f469c42d" 2534 444
assert_file "${REVIEWS}/20260723_tumor_ref_schema_recovery_external_claude_opus.v11.json" \
    "340a580ebc849d3b88ddbf3a50de9742b5357c4d6dbadd1d4e310fc11a459be3" 2990 444
assert_file "${RESULTS}/tumor_ref_source_receipt_promotion_verification.recovery.v11.json" \
    "011eb6afef7e0e5d20c02e241c90e19b53f7116bfc800801bf503d116f4164c1" 14017 444
assert_file "${RESULTS}/m2v5_runner_only_gate_replay.recovery.v11.json" \
    "b3ecee416e260461c19f2130269220cf064889c9850ae54e4bcd94cbf2489d29" 63842 444
assert_file "${RESULTS}/m2v5_runner_only_gate_replay.recovery.v11.log" \
    "46f5c44d1b7974484f9d8becaa8c75b43b2aba1fddad4e2eba8caf5b830b4b1f" 15743 444

mkdir -p -- "${ARCHIVE_PARENT}"
mkdir -- "${STAGE}"
mv -- "${RESULTS}/tumor_ref_promotion_schema_recovery_authority.v11.bundle" "${STAGE}/"
mv -- "${REVIEWS}/20260723_tumor_ref_schema_recovery_mendel.v11.json" "${STAGE}/"
mv -- "${REVIEWS}/20260723_tumor_ref_schema_recovery_nash.v11.json" "${STAGE}/"
mv -- "${REVIEWS}/20260723_tumor_ref_schema_recovery_external_claude_opus.v11.json" "${STAGE}/"
mv -- "${RESULTS}/tumor_ref_source_receipt_promotion_verification.recovery.v11.json" "${STAGE}/"
mv -- "${RESULTS}/m2v5_runner_only_gate_replay.recovery.v11.json" "${STAGE}/"
mv -- "${RESULTS}/m2v5_runner_only_gate_replay.recovery.v11.log" "${STAGE}/"
mv -- "${STAGE}" "${ARCHIVE}"

printf 'ARCHIVED_V11=%s\n' "${ARCHIVE}"
find "${ARCHIVE}" -maxdepth 2 -type f -printf '%m %s %p\n' | sort
