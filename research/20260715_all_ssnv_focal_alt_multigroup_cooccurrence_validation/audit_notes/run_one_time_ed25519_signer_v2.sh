#!/usr/bin/env bash
set -u -o pipefail

if [[ "$#" -ne 3 ]]; then
    printf 'usage: %s KEY_DIR PRIVATE_FILENAME EXPECTED_TARGET\n' "$0" >&2
    exit 2
fi

readonly KEY_DIR="$1"
readonly PRIVATE_FILENAME="$2"
readonly EXPECTED_TARGET="$3"
readonly PRIVATE_KEY="${KEY_DIR}/${PRIVATE_FILENAME}"
readonly PUBLIC_KEY="${KEY_DIR}/ed25519_public.pem"
readonly OPENSSL="/usr/bin/openssl"
readonly SHA256SUM="/usr/bin/sha256sum"
readonly STAT="/usr/bin/stat"
readonly CHMOD="/usr/bin/chmod"
readonly MV="/usr/bin/mv"
readonly MKDIR="/usr/bin/mkdir"

"${MKDIR}" -p "${KEY_DIR}"
umask 077

if [[ -e "${PRIVATE_KEY}" || -e "${PUBLIC_KEY}" ]]; then
    printf 'REFUSE_EXISTING_KEY_MATERIAL key_dir=%s\n' "${KEY_DIR}" >&2
    exit 3
fi

PASS=$("${OPENSSL}" rand -hex 48) || exit 4
export PASS
"${OPENSSL}" genpkey -algorithm ED25519 -aes-256-cbc \
    -pass env:PASS -out "${PRIVATE_KEY}" || exit 5
"${OPENSSL}" pkey -in "${PRIVATE_KEY}" -passin env:PASS \
    -pubout -out "${PUBLIC_KEY}" || exit 6
"${CHMOD}" 400 "${PRIVATE_KEY}"
"${CHMOD}" 444 "${PUBLIC_KEY}"

printf 'SIGNER_READY public=%s sha256=%s expected_target=%s\n' \
    "${PUBLIC_KEY}" \
    "$("${SHA256SUM}" "${PUBLIC_KEY}" | /usr/bin/awk '{print $1}')" \
    "${EXPECTED_TARGET}"

while IFS= read -r line; do
    if [[ "${line}" != "SIGN ${EXPECTED_TARGET}" ]]; then
        printf 'WAITING_FOR_EXACT_SIGN_TARGET\n'
        continue
    fi

    readonly SIGNATURE="${EXPECTED_TARGET}.ed25519.sig"
    readonly PENDING_SIGNATURE="${SIGNATURE}.pending"
    if [[ ! -s "${EXPECTED_TARGET}" || -e "${SIGNATURE}" || -e "${PENDING_SIGNATURE}" ]]; then
        printf 'SIGN_PRECONDITION_FAILED target=%s\n' "${EXPECTED_TARGET}" >&2
        continue
    fi

    if ! "${OPENSSL}" pkeyutl -sign -rawin -inkey "${PRIVATE_KEY}" \
        -passin env:PASS -in "${EXPECTED_TARGET}" -out "${PENDING_SIGNATURE}"; then
        printf 'SIGN_OPERATION_FAILED_KEY_RETAINED\n' >&2
        continue
    fi
    if ! "${OPENSSL}" pkeyutl -verify -rawin -pubin -inkey "${PUBLIC_KEY}" \
        -sigfile "${PENDING_SIGNATURE}" -in "${EXPECTED_TARGET}"; then
        printf 'SIGN_VERIFY_FAILED_KEY_RETAINED pending=%s\n' "${PENDING_SIGNATURE}" >&2
        continue
    fi

    "${MV}" "${PENDING_SIGNATURE}" "${SIGNATURE}"
    "${CHMOD}" 444 "${SIGNATURE}"
    unset PASS
    "${CHMOD}" 000 "${PRIVATE_KEY}"
    printf 'SIGNED signature=%s sha256=%s private_mode=%s\n' \
        "${SIGNATURE}" \
        "$("${SHA256SUM}" "${SIGNATURE}" | /usr/bin/awk '{print $1}')" \
        "$("${STAT}" -c '%a' "${PRIVATE_KEY}")"
    exit 0
done

printf 'SIGNER_INPUT_CLOSED_WITHOUT_SIGNATURE\n' >&2
exit 7
