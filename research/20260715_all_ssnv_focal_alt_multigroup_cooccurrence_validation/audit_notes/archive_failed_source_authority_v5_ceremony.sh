#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly AUTHORITY_DIR="${REPO_ROOT}/docs/provenance/source_authorities"
readonly AUTHORITY="${AUTHORITY_DIR}/20260718_all_ssnv_focal_alt_release_source_authority.v5.json"
readonly APPROVAL="${AUTHORITY_DIR}/20260718_all_ssnv_focal_alt_release_source_authority.v5.approval.json"
readonly SIGNATURE="${APPROVAL}.ed25519.sig"
readonly PENDING="${SIGNATURE}.pending"
readonly DOC_ARCHIVE="${REPO_ROOT}/docs/archive/2026/07/20260718_unsigned_source_authority_v5_signer_failure_01"
readonly KEY_DIR="/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix"
readonly KEY_ARCHIVE_PARENT="/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/archive"
readonly KEY_ARCHIVE="${KEY_ARCHIVE_PARENT}/20260718_all_ssnv_v5_summary_hotfix_failed_signer_01"
readonly PUBLIC_KEY="${KEY_DIR}/ed25519_public.pem"
readonly PRIVATE_KEY="${KEY_DIR}/ed25519_private_encrypted_unrecoverable_after_signing.pem"
readonly AUTHORITY_SHA256="99d743da27779721c75bfe015b1e0d094f69345518279fe564928e9cad1df1b9"
readonly APPROVAL_SHA256="60e15e5ce92c36d025de24199e5c71b768abcd3b9d6b421b2c480cf202095aa9"
readonly PUBLIC_KEY_SHA256="cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920"
readonly PRIVATE_KEY_SHA256="bac690113ca9931058fbb037e9992bcccbb715f64c26eb56e03ee841a8dd3265"

require_identity() {
    local path="$1"
    local expected_sha256="$2"
    local expected_mode="$3"
    [[ -f "${path}" ]]
    [[ "$(/usr/bin/sha256sum "${path}" | /usr/bin/awk '{print $1}')" == "${expected_sha256}" ]]
    [[ "$(/usr/bin/stat -c '%a' "${path}")" == "${expected_mode}" ]]
}

require_identity "${AUTHORITY}" "${AUTHORITY_SHA256}" "444"
require_identity "${APPROVAL}" "${APPROVAL_SHA256}" "444"
require_identity "${PUBLIC_KEY}" "${PUBLIC_KEY_SHA256}" "444"
require_identity "${PRIVATE_KEY}" "${PRIVATE_KEY_SHA256}" "400"
[[ ! -e "${SIGNATURE}" ]]
[[ ! -e "${PENDING}" ]]
[[ -f "${DOC_ARCHIVE}/SUMMARY.md" ]]
[[ ! -e "${DOC_ARCHIVE}/UNSIGNED_NOT_AUTHORIZED.20260718_all_ssnv_focal_alt_release_source_authority.v5.json" ]]
[[ ! -e "${DOC_ARCHIVE}/UNSIGNED_NOT_AUTHORIZED.20260718_all_ssnv_focal_alt_release_source_authority.v5.approval.json" ]]
[[ ! -e "${KEY_ARCHIVE}" ]]

/usr/bin/mkdir -p "${KEY_ARCHIVE_PARENT}"
/usr/bin/mv "${AUTHORITY}" \
    "${DOC_ARCHIVE}/UNSIGNED_NOT_AUTHORIZED.20260718_all_ssnv_focal_alt_release_source_authority.v5.json"
/usr/bin/mv "${APPROVAL}" \
    "${DOC_ARCHIVE}/UNSIGNED_NOT_AUTHORIZED.20260718_all_ssnv_focal_alt_release_source_authority.v5.approval.json"
/usr/bin/mv "${KEY_DIR}" "${KEY_ARCHIVE}"

/usr/bin/printf \
    'ARCHIVE_PASS authority_sha256=%s approval_sha256=%s public_key_sha256=%s key_archive=%s\n' \
    "${AUTHORITY_SHA256}" "${APPROVAL_SHA256}" "${PUBLIC_KEY_SHA256}" "${KEY_ARCHIVE}"
