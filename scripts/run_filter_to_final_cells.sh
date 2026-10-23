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
    team-hafler       # hafler-pmdbs-sn-rnaseq-pfc
    team-lee          # lee-pmdbs-sn-rnaseq
    team-jakobsson    # jakobsson-pmdbs-sn-rnaseq
    team-scherzer     # scherzer-pmdbs-sn-rnaseq-mtg
    team-hardy        # hardy-pmdbs-sn-rnaseq
    team-sulzer       # sulzer-pmdbs-sn-rnaseq
    asap-cohort       # cohort-pmdbs-sc-rnaseq
)

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
            --final      "${DATA_DIR}/${prefix}.final.h5ad" \
            --unfiltered "${DATA_DIR}/${prefix}.merged_cleaned_unfiltered.h5ad" \
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
