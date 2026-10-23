#!/usr/bin/env bash
# Loop filter_to_final_cells.py over every dataset in cohort_analysis/.
#
# Usage:
#   ./run_filter_to_final_cells.sh [--overwrite] [prefix ...]
#
#   With no prefixes, all datasets below are processed. Existing outputs are
#   skipped unless --overwrite is given.
#
#   DATA_DIR can be set to override the default data location.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DATA_DIR:-${SCRIPT_DIR}/../ground_truth/cohort_analysis}"

ALL_PREFIXES=(
    Substantia_Nigra
    Putamen
    Amygdala
    Middle_Temporal_Gyrus
    Middle_Frontal_Gyrus
    Anterior_Cingulate_Gyrus
    Inferior_Parietal_Lobule
    Cingulate_Gyrus
    Hippocampus
    )

SOURCE=asap-cohort.merged_cleaned_unfiltered.h5ad


EXTRA_ARGS=()
PREFIXES=()
for arg in "$@"; do
    case "$arg" in
        --overwrite) EXTRA_ARGS+=(--overwrite) ;;
        -h|--help)   sed -n '2,10p' "$0"; exit 0 ;;
        *)           PREFIXES+=("$arg") ;;
    esac
done
[[ ${#PREFIXES[@]} -eq 0 ]] && PREFIXES=("${ALL_PREFIXES[@]}")

failed=()
for prefix in "${PREFIXES[@]}"; do
    echo "=== Processing ${prefix} ==="
    if ! python "${SCRIPT_DIR}/filter_to_final_cells.py" \
            --final      "${DATA_DIR}/${prefix}.h5ad" \
            --unfiltered "${DATA_DIR}/${SOURCE}" \
            --output     "${DATA_DIR}/${prefix}.mmc.h5ad" \
            "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"; then
        echo "!!! ${prefix} failed" >&2
        failed+=("$prefix")
    fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "Failed: ${failed[*]}" >&2
    exit 1
fi
echo "All done."
