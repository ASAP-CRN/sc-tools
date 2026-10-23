#!/usr/bin/env bash
# Run MapMyCells (sc_tools.mmc) over every {prefix}.mmc.h5ad in cohort_analysis/.
#
# Usage:
#   ./run_mmc.sh [--overwrite] [prefix ...]
#
#   With no prefixes, all datasets below are processed. Datasets whose
#   results CSV already exists are skipped unless --overwrite is given.
#
#   Environment overrides: DATA_DIR, RESULTS_DIR, REFERENCE, PRECOMPUTED_STATS,
#   QUERY_MARKERS (set to "" to compute markers on the fly), N_PROCESSORS,
#   MAX_GB, TMP_DIR, SC_TOOLS_SRC, PYTHON.

set -euo pipefail

SC_TOOLS_DIR="${SC_TOOLS_DIR:-/home/ergonyc/Projects/ASAP/sc-tools/src/sc_tools}"
BENCHMARK_ROOT="${BENCHMARK_ROOT:-/home/ergonyc/Projects/ASAP/benchmark}"
DATA_DIR="${DATA_DIR:-${BENCHMARK_ROOT}/ground_truth/cohort_analysis}"
ABC_ROOT="${ABC_ROOT:-/home/ergonyc/Projects/ABC/cell_type_mapper}"
ABC_DATA="${ABC_ROOT}/examples/data/abc_atlas_data"
REFERENCE="${REFERENCE:-WHB}"

PRECOMPUTED_STATS="${PRECOMPUTED_STATS:-${ABC_DATA}/precomputed_stats.siletti.training.h5}"
QUERY_MARKERS="${QUERY_MARKERS:-${ABC_DATA}/query_markers.n10.20240221800.json}"
RESULTS_DIR="${RESULTS_DIR:-${ABC_ROOT}/MMC/${REFERENCE}_RESULTS}"
N_PROCESSORS="${N_PROCESSORS:-8}"
MAX_GB="${MAX_GB:-48}"
TMP_DIR="${TMP_DIR:-${ABC_ROOT}/tmp}"

echo "TMP_DIR: ${TMP_DIR}"
echo "RESULTS_DIR: ${RESULTS_DIR}"

# sc_tools from the local sc-tools checkout (no install required)
# SC_TOOLS_SRC="${SC_TOOLS_SRC:-${SCRIPT_DIR}/../../sc-tools/src}"
# export PYTHONPATH="${SC_TOOLS_SRC}${PYTHONPATH:+:${PYTHONPATH}}"
# PYTHON="${PYTHON:-python}"

ALL_PREFIXES=(
    # team-hafler       # hafler-pmdbs-sn-rnaseq-pfc
    # team-lee          # lee-pmdbs-sn-rnaseq
    # team-jakobsson    # jakobsson-pmdbs-sn-rnaseq
    # team-scherzer     # scherzer-pmdbs-sn-rnaseq-mtg
    # team-hardy        # hardy-pmdbs-sn-rnaseq
    # team-sulzer       # sulzer-pmdbs-sn-rnaseq
    # asap-cohort       # cohort-pmdbs-sc-rnaseq
    # Substantia_Nigra
    # Putamen
    # Amygdala
    Middle_Temporal_Gyrus
    Middle_Frontal_Gyrus
    Anterior_Cingulate_Gyrus
    Inferior_Parietal_Lobule
    Cingulate_Gyrus
    Hippocampus
)

OVERWRITE=0
PREFIXES=()
for arg in "$@"; do
    case "$arg" in
        --overwrite) OVERWRITE=1 ;;
        -h|--help)   sed -n '2,13p' "$0"; exit 0 ;;
        *)           PREFIXES+=("$arg") ;;
    esac
done
[[ ${#PREFIXES[@]} -eq 0 ]] && PREFIXES=("${ALL_PREFIXES[@]}")

mkdir -p "${RESULTS_DIR}"

MMC_ARGS=(
    --mmc-precomputed-stats "${PRECOMPUTED_STATS}"
    --n-processors "${N_PROCESSORS}"
    --max-gb "${MAX_GB}"
)
[[ -n "${QUERY_MARKERS}" ]] && MMC_ARGS+=(--mmc-marker-genes "${QUERY_MARKERS}")
[[ -n "${TMP_DIR}" ]] && MMC_ARGS+=(--tmp-dir "${TMP_DIR}")

failed=()
for prefix in "${PREFIXES[@]}"; do
    input="${DATA_DIR}/${prefix}.mmc.h5ad"
    output_prefix="${RESULTS_DIR}/${prefix}.mmc.${REFERENCE}"

    if [[ ! -f "${input}" ]]; then
        echo "!!! ${prefix}: missing ${input}" >&2
        failed+=("$prefix")
        continue
    fi
    if [[ -f "${output_prefix}.results.csv" && ${OVERWRITE} -eq 0 ]]; then
        echo "Skipping ${prefix}: ${output_prefix}.results.csv exists (use --overwrite to replace)"
        continue
    fi

    echo "=== Processing ${prefix} ==="
    if ! python "${SC_TOOLS_DIR}/mmc.py" \
            --adata-input   "${input}" \
            --output-prefix "${output_prefix}" \
            "${MMC_ARGS[@]}"; then
        echo "!!! ${prefix} failed" >&2
        failed+=("$prefix")
    fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "Failed: ${failed[*]}" >&2
    exit 1
fi
echo "All done."
