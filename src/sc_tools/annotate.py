#!/usr/bin/env python3
# Re-label a scVI-integrated cohort with an alternate MMC reference:
#   attach MMC transcriptional phenotypes -> scANVI -> Harmony -> scib metrics
# NOTE: not upstream -- sync.py wipes src/sc_tools/, so keep a copy elsewhere.

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scvi

from sc_tools.add_harmony import add_harmony

# need this to use alternate reference labels
from sc_tools.transcriptional_phenotype import read_csv_results2, summarize_mmc_results
from sc_tools.artifact_metrics import get_artifact_metrics
from sc_tools.label_scanvi import label_with_scanvi


def lowercase_taxonomy(
    mmc_results: pd.DataFrame, taxonomy_levels: list[str]
) -> tuple[pd.DataFrame, list[str]]:
    """
    summarize_mmc_results expects lower-case level names, but newer references
    (e.g. SEA-AD Multiregion 2026: ["Class", "Subclass", "Supertype"]) capitalize them.
    """
    rename = {}
    for level in taxonomy_levels:
        for col in mmc_results.columns:
            if col.startswith(f"{level}_"):
                rename[col] = level.lower() + col[len(level) :]
    return mmc_results.rename(columns=rename), [lvl.lower() for lvl in taxonomy_levels]


def replace_mmc_results(adata: ad.AnnData, mmc_results: pd.DataFrame) -> ad.AnnData:
    """
    Replace any existing MMC columns in adata.obs with mmc_results (aligned on cell id).
    Cells without an MMC call are labelled "Unknown" so scANVI treats them as unlabelled.
    """
    mmc_results = mmc_results.reindex(adata.obs_names)
    n_missing = int(mmc_results["cell_type"].isna().sum())
    if n_missing:
        print(f"WARNING: {n_missing} cells have no MMC result; labelling 'Unknown'")
    mmc_results["cell_type"] = mmc_results["cell_type"].fillna("Unknown").astype(str)

    localobs = adata.obs.drop(columns=mmc_results.columns, errors="ignore")
    adata.obs = localobs.join(mmc_results)
    return adata


def main(args: argparse.Namespace):
    # Set the number of data loader workers
    # num_workers = max(1, os.cpu_count() - 1)
    num_workers = 0  # Pytorch bug unable to mmap solution https://github.com/pytorch/pytorch/issues/92134
    scvi.settings.dl_num_workers = num_workers
    print(f"Using {scvi.settings.dl_num_workers} workers")

    # 0. Load data
    adata = ad.read_h5ad(args.adata_input)  # type: ignore

    # 1. Force transcriptional phenotype from the (alternate) reference MMC results
    mmc_results, taxonomy_levels = read_csv_results2(args.mmc_results)
    mmc_results, taxonomy_levels = lowercase_taxonomy(mmc_results, taxonomy_levels)
    transformed_results = summarize_mmc_results(
        mmc_results, args.workflow_name, taxonomy_levels
    )
    adata = replace_mmc_results(adata, transformed_results)

    # 2. Load the scVI model trained on this adata
    model = scvi.model.SCVI.load(
        dir_path=args.scvi_outputs_dir,
        adata=adata,
    )

    # 3. Get scANVI model
    adata, scanvi_model = label_with_scanvi(
        adata,
        model,
        num_workers,
        args.workflow_name,
        args.latent_key,
        args.predictions_key,
    )
    scanvi_model.save(args.output_scanvi_dir, overwrite=True)

    # 4. Save the cell types to parquet
    adata.obs[[args.predictions_key]].to_parquet(
        args.output_cell_types_file, compression="gzip"
    )

    # 5. Add harmony (X_pca_harmony is one of the embeddings benchmarked below)
    add_harmony(adata, args.batch_key)

    # 6. Save the annotated adata
    adata.write_h5ad(filename=args.adata_output, compression="gzip")

    # 7. Save the artifact metrics
    report_dir = Path.cwd() / args.output_report_dir
    if not report_dir.exists():
        report_dir.mkdir(parents=True, exist_ok=True)

    report_df = get_artifact_metrics(
        adata, args.batch_key, args.predictions_key, report_dir
    )
    report_df.to_csv((report_dir / "scib_report.csv"), index=True)


def cli():
    parser = argparse.ArgumentParser(
        description="Re-label a scVI-integrated dataset with MMC results from an alternate reference via scANVI, then add Harmony and scib metrics"
    )
    parser.add_argument(
        "--workflow-name",
        type=str,
        required=True,
        help="Workflow name to detect sample source and set scANVI training parameters (e.g. 'pmdbs_sc_rnaseq', 'mouse_sc_rnaseq')",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        required=True,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--latent-key",
        type=str,
        required=True,
        help="Latent key to save the scANVI latent to (must differ from the scVI latent key)",
    )
    parser.add_argument(
        "--predictions-key",
        type=str,
        required=True,
        help="scANVI cell type predictions column name in AnnData object",
    )
    parser.add_argument(
        "--adata-input",
        type=str,
        required=True,
        help="scVI-integrated AnnData object for a dataset",
    )
    parser.add_argument(
        "--scvi-outputs-dir",
        type=str,
        required=True,
        help="Saved scVI model folder for --adata-input",
    )
    parser.add_argument(
        "--mmc-results",
        type=str,
        required=True,
        help="Path to MMC results CSV for the reference to label with",
    )
    parser.add_argument(
        "--adata-output",
        type=str,
        required=True,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--output-scanvi-dir",
        type=str,
        required=True,
        help="Output folder to save scANVI model",
    )
    parser.add_argument(
        "--output-cell-types-file",
        type=str,
        required=True,
        help="Output file to write cell types to",
    )
    parser.add_argument(
        "--output-report-dir",
        type=str,
        required=True,
        help="Output folder to save `scib` report",
    )

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    cli()
