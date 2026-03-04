#!/usr/bin/env bash
set -euo pipefail

# 安全搬移 InterSubMod 輸出目錄到 /bip8_disk。
# - 不刪除資料，只做搬移（mv）
# - 預設跳過 symlink（避免移走 output/bip8_disk_output）
# - 產生 TSV manifest 便於追蹤

SOURCE_ROOT="/big8_disk/liaoyoyo2001/InterSubMod_runs/output"
DEST_ROOT="/bip8_disk/liaoyoyo2001/InterSubMod_out/output"
MANIFEST=""
DRY_RUN=false
INCLUDE_SYMLINK=false

usage() {
    cat <<'EOF'
Usage:
  migrate_output_to_bip8.sh [options]

Options:
  --source DIR            Source output root
  --dest DIR              Destination output root
  --manifest FILE         TSV manifest path (default: ./output_migration_YYYYmmdd_HHMMSS.tsv)
  --dry-run               Print actions only
  --include-symlink       Also move symlink entries (default: skip symlink)
  -h, --help              Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source) SOURCE_ROOT="$2"; shift 2 ;;
        --dest) DEST_ROOT="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --include-symlink) INCLUDE_SYMLINK=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${MANIFEST}" ]]; then
    MANIFEST="./output_migration_$(date +%Y%m%d_%H%M%S).tsv"
fi

if [[ ! -d "${SOURCE_ROOT}" ]]; then
    echo "[ERROR] Source root not found: ${SOURCE_ROOT}" >&2
    exit 1
fi

mkdir -p "${DEST_ROOT}"
mkdir -p "$(dirname "${MANIFEST}")"

echo -e "timestamp\tname\tsource_path\tdestination_path\tsize_kb_before\tsize_kb_after\tstatus\tnote" > "${MANIFEST}"

for entry in "${SOURCE_ROOT}"/*; do
    [[ -e "${entry}" ]] || continue

    name="$(basename "${entry}")"
    now="$(date '+%Y-%m-%d %H:%M:%S')"

    if [[ -L "${entry}" ]] && [[ "${INCLUDE_SYMLINK}" != true ]]; then
        echo -e "${now}\t${name}\t${entry}\t-\t0\t0\tskip\tsymlink" >> "${MANIFEST}"
        continue
    fi

    if [[ ! -d "${entry}" ]]; then
        echo -e "${now}\t${name}\t${entry}\t-\t0\t0\tskip\tnot_directory" >> "${MANIFEST}"
        continue
    fi

    size_kb_before="$(du -sk "${entry}" | awk '{print $1}')"
    target="${DEST_ROOT}/${name}"
    note=""

    if [[ -e "${target}" ]]; then
        target="${DEST_ROOT}/${name}_from_move_$(date +%Y%m%d_%H%M%S)"
        note="name_conflict_renamed"
    fi

    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] mv ${entry} ${target}"
        echo -e "${now}\t${name}\t${entry}\t${target}\t${size_kb_before}\t0\tdry_run\t${note}" >> "${MANIFEST}"
        continue
    fi

    mv "${entry}" "${target}"
    size_kb_after="$(du -sk "${target}" | awk '{print $1}')"
    echo -e "${now}\t${name}\t${entry}\t${target}\t${size_kb_before}\t${size_kb_after}\tmoved\t${note}" >> "${MANIFEST}"
done

echo "[DONE] Migration manifest: ${MANIFEST}"

