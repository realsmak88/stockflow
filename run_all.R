# =============================================================================
# stockflow :: run_all.R  — complete, sequential pipeline
# Run from the package root:  source("run_all.R")
# =============================================================================
# 0. Install dependencies (once): source("scripts/install_packages.R")

message(">>> 01 Cleaning & QC");            source("R/prepare/01_clean_data.R")
message(">>> 02 Biological parameters");    source("R/analysis/02_biological_params.R")
message(">>> 03 Standardized CPUE");         source("R/analysis/03_cpue_standardisation.R")
message(">>> 04 Stock assessment & BRP"); source("R/analysis/04_stock_assessment.R")
message(">>> 05 Management: TAC/quotas/closures/TMD"); source("R/analysis/05_management.R")
# message(">>> 06 openMSE MSE (long)");      source("R/analysis/06_mse_openMSE.R")
message(">>> Pipeline complete — see output/tables/ and output/figures/")