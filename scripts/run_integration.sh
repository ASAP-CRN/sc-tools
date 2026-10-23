#!/usr/bin/env bash
# Normalize, attach MapMyCells phenotypes, and scVI-integrate every
# {prefix}.mmc.h5ad in cohort_analysis/ -- mirrors the cohort_analysis WDL:
#
#   1. process.py                    {prefix}.mmc.h5ad -> {prefix}.normalized.h5ad
#   2. transcriptional_phenotype.py  + MMC results.csv -> {prefix}.mmc_annotated.{REF}.h5ad
#   3. integrate_scvi.py                               -> {prefix}.scvi_integrated.{REF}.h5ad
#
# Requires run_mmc.sh to have produced {prefix}.mmc.{REF}.results.csv first.
#
# Usage:
#   ./run_integration.sh [--overwrite] [prefix ...]
#
#   With no prefixes, all datasets below are processed. Steps whose output
#   already exists are skipped unless --overwrite is given.
#
#   Environment overrides: DATA_DIR, MMC_RESULTS_DIR, OUT_DIR, REFERENCE,
#   WORKFLOW_NAME, BATCH_KEY, NORM_TARGET_SUM, N_TOP_GENES, N_COMPS,
#   SCVI_LATENT_KEY, SC_TOOLS_SRC.

set -euo pipefail

BENCHMARK_ROOT="${BENCHMARK_ROOT:-/home/ergonyc/Projects/ASAP/benchmark}"
DATA_DIR="${DATA_DIR:-${BENCHMARK_ROOT}/ground_truth/cohort_analysis}"
ABC_ROOT="${ABC_ROOT:-/home/ergonyc/Projects/ABC/cell_type_mapper}"
REFERENCE="${REFERENCE:-WHB}"
MMC_RESULTS_DIR="${MMC_RESULTS_DIR:-${ABC_ROOT}/MMC/${REFERENCE}_RESULTS}"
MMC_OUT_DIR="${MMC_OUT_DIR:-${MMC_RESULTS_DIR}}"
OUT_DIR="${OUT_DIR:-${ABC_ROOT}/MMC/}"

# Defaults from sc-rnaseq-wf workflows/main.wdl
WORKFLOW_NAME="${WORKFLOW_NAME:-pmdbs_sc_rnaseq}"
BATCH_KEY="${BATCH_KEY:-batch_id}"
NORM_TARGET_SUM="${NORM_TARGET_SUM:-10000}"
N_TOP_GENES="${N_TOP_GENES:-3000}"
N_COMPS="${N_COMPS:-30}"
SCVI_LATENT_KEY="${SCVI_LATENT_KEY:-X_scVI}"

# sc_tools from the local sc-tools checkout (scripts import sc_tools.helpers)
SC_TOOLS_SRC="${SC_TOOLS_SRC:-/home/ergonyc/Projects/ASAP/sc-tools/src}"
SC_TOOLS_DIR="${SC_TOOLS_SRC}/sc_tools"
export PYTHONPATH="${SC_TOOLS_SRC}${PYTHONPATH:+:${PYTHONPATH}}"

echo "MMC_RESULTS_DIR: ${MMC_RESULTS_DIR}"
echo "OUT_DIR: ${OUT_DIR}"

ALL_PREFIXES=(
    team-hafler       # hafler-pmdbs-sn-rnaseq-pfc
    team-lee          # lee-pmdbs-sn-rnaseq
    team-jakobsson    # jakobsson-pmdbs-sn-rnaseq
    team-scherzer     # scherzer-pmdbs-sn-rnaseq-mtg
    team-hardy        # hardy-pmdbs-sn-rnaseq
    team-sulzer       # sulzer-pmdbs-sn-rnaseq
    asap-cohort       # cohort-pmdbs-sc-rnaseq
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

OVERWRITE=0
PREFIXES=()
for arg in "$@"; do
    case "$arg" in
        --overwrite) OVERWRITE=1 ;;
        -h|--help)   sed -n '2,19p' "$0"; exit 0 ;;
        *)           PREFIXES+=("$arg") ;;
    esac
done
[[ ${#PREFIXES[@]} -eq 0 ]] && PREFIXES=("${ALL_PREFIXES[@]}")

mkdir -p "${OUT_DIR}"

# skip_step <output>: true if output exists and --overwrite was not given
skip_step() {
    [[ -e "$1" && ${OVERWRITE} -eq 0 ]] && echo "  skipping: $1 exists (use --overwrite to replace)"
}

failed=()
for prefix in "${PREFIXES[@]}"; do
    input="${DATA_DIR}/${prefix}.mmc.h5ad"
    mmc_results="${MMC_RESULTS_DIR}/${prefix}.mmc.${REFERENCE}.results.csv"

    normalized="${OUT_DIR}/${prefix}.normalized.h5ad"
    all_genes="${OUT_DIR}/${prefix}.all_genes.csv"
    hvg_genes="${OUT_DIR}/${prefix}.hvg_genes.csv"
    annotated="${OUT_DIR}/${prefix}.mmc_annotated.${REFERENCE}.h5ad"
    cell_types="${MMC_OUT_DIR}/${prefix}.mmc_results.${REFERENCE}.parquet"
    integrated="${MMC_OUT_DIR}/${prefix}.scvi_integrated.${REFERENCE}.h5ad"
    scvi_dir="${MMC_OUT_DIR}/${prefix}.${REFERENCE}_scvi_model"

    for f in "${input}" "${mmc_results}"; do
        if [[ ! -f "${f}" ]]; then
            echo "!!! ${prefix}: missing ${f}" >&2
            failed+=("$prefix")
            continue 2
        fi
    done

    echo "=== Processing ${prefix} ==="

    echo "--- [1/3] process (normalize + HVG + PCA)"
    if ! skip_step "${normalized}"; then
        if ! python "${SC_TOOLS_DIR}/process.py" \
                --adata-input      "${input}" \
                --workflow-name    "${WORKFLOW_NAME}" \
                --batch-key        "${BATCH_KEY}" \
                --norm-target-sum  "${NORM_TARGET_SUM}" \
                --n-top-genes      "${N_TOP_GENES}" \
                --n-comps          "${N_COMPS}" \
                --adata-output     "${normalized}" \
                --output-all-genes "${all_genes}" \
                --output-hvg-genes "${hvg_genes}"; then
            echo "!!! ${prefix} failed at process" >&2
            failed+=("$prefix")
            continue
        fi
    fi

    # echo "--- [2/3] transcriptional_phenotype (attach MMC labels)"
    # if ! skip_step "${annotated}"; then
    #     if ! python "${SC_TOOLS_DIR}/transcriptional_phenotype.py" \
    #             --adata-input            "${normalized}" \
    #             --workflow-name          "${WORKFLOW_NAME}" \
    #             --mmc-results            "${mmc_results}" \
    #             --adata-output           "${annotated}" \
    #             --output-cell-types-file "${cell_types}"; then
    #         echo "!!! ${prefix} failed at transcriptional_phenotype" >&2
    #         failed+=("$prefix")
    #         continue
    #     fi
    # fi

    echo "--- [3/3] integrate_scvi"
    if ! skip_step "${integrated}"; then
        if ! python "${SC_TOOLS_DIR}/integrate_scvi.py" \
                --latent-key      "${SCVI_LATENT_KEY}" \
                --batch-key       "${BATCH_KEY}" \
                --adata-input     "${normalized}" \
                --adata-output    "${integrated}" \
                --output-scvi-dir "${scvi_dir}"; then
            echo "!!! ${prefix} failed at integrate_scvi" >&2
            failed+=("$prefix")
            continue
        fi
    fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "Failed: ${failed[*]}" >&2
    exit 1
fi
echo "All done."
