#!/usr/bin/env bash
# Re-label every scVI-integrated {prefix} with each MapMyCells reference in
# ALL_REFS (sc_tools/annotate.py: MMC phenotypes -> scANVI -> Harmony -> scib):
#
#   {prefix}.scvi_integrated.{INTEGRATION_REF}.h5ad + {prefix}.{INTEGRATION_REF}_scvi_model
#     + {ref}_RESULTS/{prefix}.mmc.{ref}.results.csv
#     -> {ref}_RESULTS/{prefix}.annotated.{ref}.h5ad (+ scANVI model, cell types, scib report)
#
# Requires run_integration.sh (REFERENCE=INTEGRATION_REF) and run_mmc*.sh for each ref.
#
# Usage:
#   ./run_annotation.sh [--overwrite] [prefix ...]
#
#   With no prefixes, all datasets below are processed. Outputs that already
#   exist are skipped unless --overwrite is given.
#
#   Environment overrides: DATA_DIR, MMC_ROOT, INTEGRATION_REF, WORKFLOW_NAME,
#   BATCH_KEY, SCANVI_LATENT_KEY, SCANVI_PREDICTIONS_KEY, SC_TOOLS_SRC.

set -euo pipefail

ABC_ROOT="${ABC_ROOT:-/home/ergonyc/Projects/ABC/cell_type_mapper}"

MMC_ROOT="${MMC_ROOT:-${ABC_ROOT}/MMC}"
# scVI integration is reference-independent; reuse the one run_integration.sh made
INTEGRATION_DIR="${MMC_ROOT}"



# Defaults from sc-rnaseq-wf workflows/main.wdl
WORKFLOW_NAME="${WORKFLOW_NAME:-pmdbs_sc_rnaseq}"

BATCH_KEY="${BATCH_KEY:-batch_id}"

SCANVI_LATENT_KEY="${SCANVI_LATENT_KEY:-X_scANVI}"
SCANVI_PREDICTIONS_KEY="${SCANVI_PREDICTIONS_KEY:-C_scANVI}"

# sc_tools from the local sc-tools checkout (scripts import sc_tools.*)
SC_TOOLS_SRC="${SC_TOOLS_SRC:-/home/ergonyc/Projects/ASAP/sc-tools/src}"
SC_TOOLS_DIR="${SC_TOOLS_SRC}/sc_tools"
export PYTHONPATH="${SC_TOOLS_SRC}${PYTHONPATH:+:${PYTHONPATH}}"

echo "INTEGRATION_DIR: ${INTEGRATION_DIR}"

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
    Amygdala
    Middle_Temporal_Gyrus
    Middle_Frontal_Gyrus
    Anterior_Cingulate_Gyrus
    Inferior_Parietal_Lobule
    Cingulate_Gyrus
    Hippocampus
)

ALL_REFS=(
    WHB
    SEAAD
)

OVERWRITE=0
PREFIXES=()
for arg in "$@"; do
    case "$arg" in
        --overwrite) OVERWRITE=1 ;;
        -h|--help)   sed -n '2,18p' "$0"; exit 0 ;;
        *)           PREFIXES+=("$arg") ;;
    esac
done
[[ ${#PREFIXES[@]} -eq 0 ]] && PREFIXES=("${ALL_PREFIXES[@]}")

failed=()
for prefix in "${PREFIXES[@]}"; do
    integrated="${INTEGRATION_DIR}/${prefix}.scvi_integrated.WHB.h5ad" # the WHB is a mistake 
    scvi_dir="${INTEGRATION_DIR}/${prefix}.WHB_scvi_model" # the WHB is a mistake 

    for ref in "${ALL_REFS[@]}"; do
        mmc_out_dir="${MMC_ROOT}/${ref}_RESULTS"
        mmc_results="${mmc_out_dir}/${prefix}.mmc.${ref}.results.csv"

        annotated="${mmc_out_dir}/${prefix}.annotated.${ref}.h5ad"
        scanvi_dir="${mmc_out_dir}/${prefix}.${ref}_scanvi_model"
        cell_types="${mmc_out_dir}/${prefix}.scanvi_cell_types.${ref}.parquet"
        report_dir="${mmc_out_dir}/${prefix}.scib_report.${ref}"

        echo "=== Annotating ${prefix} with ${ref} ==="

        missing=0
        for f in "${integrated}" "${scvi_dir}" "${mmc_results}"; do
            if [[ ! -e "${f}" ]]; then
                echo "!!! ${prefix}/${ref}: missing ${f}" >&2
                missing=1
            fi
        done
        if [[ ${missing} -eq 1 ]]; then
            failed+=("${prefix}/${ref}")
            continue
        fi

        if [[ -e "${annotated}" && ${OVERWRITE} -eq 0 ]]; then
            echo "  skipping: ${annotated} exists (use --overwrite to replace)"
            continue
        fi

        mkdir -p "${mmc_out_dir}"
        if ! python "${SC_TOOLS_DIR}/annotate.py" \
                --workflow-name          "${WORKFLOW_NAME}" \
                --batch-key              "${BATCH_KEY}" \
                --latent-key             "${SCANVI_LATENT_KEY}" \
                --predictions-key        "${SCANVI_PREDICTIONS_KEY}" \
                --adata-input            "${integrated}" \
                --scvi-outputs-dir       "${scvi_dir}" \
                --mmc-results            "${mmc_results}" \
                --adata-output           "${annotated}" \
                --output-scanvi-dir      "${scanvi_dir}" \
                --output-cell-types-file "${cell_types}" \
                --output-report-dir      "${report_dir}"; then
            echo "!!! ${prefix}/${ref} failed at annotate" >&2
            failed+=("${prefix}/${ref}")
        fi
    done
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "Failed: ${failed[*]}" >&2
    exit 1
fi
echo "All done."
