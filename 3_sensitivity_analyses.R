####################################################################################################
## Author: Adrien Allorant
##
## Description:
##    This script performs a sensitivity analysis on previously fitted INLA models by refitting them 
##    with certain predictors (wealth, education) removed. It loops through a defined set of indicators, 
##    regions, and sexes. For each scenario, it attempts to:
##       1. Load the best model identified from previous runs.
##       2. Load the corresponding adjacency matrix.
##       3. Update the model formula by removing certain covariates (e.g., wealth, education).
##       4. Refit the model and save the results for comparison.
##
## Requirements:
##    - R version >= 4.0.0
##    - Packages: INLA, tidyverse, stringr
##    - The best-fit models and adjacency matrices should be previously generated and saved.
##
## Inputs:
##    - Best model fit files: "<indicator>-<tag>-best-logit-sae-fit.rds" (per region, sex)
##    - Adjacency matrices: "adjM_<indicator>_<region_with_underscores>.rds"
##
## Outputs:
##    - Sensitivity analysis fit results: "<indicator>-<tag_new>-sensitivity-<model_name>-fit.rds"
##
## Notes:
##    - Adjust directory paths and environment checks as needed.
##    - tag and tag_new differentiate between original and sensitivity analysis runs.
##
####################################################################################################

# Load the necessary packages
library(INLA)
library(tidyverse)
library(stringr)

#####################################
## 00: Setup and Configuration      ##
#####################################

if (run_rstudio) {
  # Example: using a local directory structure
  root <- "."  # or set this to your project root using here::here()
} else {
  args <- commandArgs(trailingOnly = TRUE)
  root <- args[1]
  indicator <- args[3]
}
# Define output directory
out_dir <- file.path(main_dir, 'output')

# Define regions, sexes, indicators, and tags
regions <- c("Eastern Africa", "Western Africa", "Southern Africa", "Central Africa")
sexes <- c("Men", "Women")
indicators <- c("hivstatus", "test12m", "arv")

# Original run tag and new tag for sensitivity analysis
tag <- "August2nd"
tag_new <- "November18th"

#####################################
## 01: Sensitivity Analysis Function
#####################################

# This function performs sensitivity analyses by refitting the best model after removing certain predictors.
perform_sensitivity_analysis <- function(region, sex, indicator, tag, tag_new, out_dir) {
  
  # Construct file paths for model and adjacency matrix
  model_file <- file.path(out_dir, region, sex, paste0(indicator, '-', tag, "-best-logit-sae-fit.rds"))
  adjM_file <- file.path(out_dir, paste0('adjM_', indicator, '_', str_replace_all(region, " ", "_"), '.rds'))
  
  # Check model file existence
  if (!file.exists(model_file)) {
    message("Model file does not exist:", model_file, " Skipping...")
    return(NULL)
  }

  # Load the best model
  cat("Loading best model from:", model_file, "\n")
  best_model <- readRDS(model_file)
  cat("Best model loaded successfully.\n")

  # Extract the formula and data from the best model object
  best_formula <- best_model$.args$formula
  data <- best_model$.args$data
  
  # Load the adjacency matrix if it exists
  if (!file.exists(adjM_file)) {
    message("Adjacency matrix file does not exist:", adjM_file, " Skipping...")
    return(NULL)
  }
  
  cat("Loading adjacency matrix from:", adjM_file, "\n")
  adjM <- readRDS(adjM_file)
  cat("Adjacency matrix loaded successfully.\n")

  # Create updated formulas for sensitivity tests (removing wealth or education)
  # Note: Replace 'wealthq_idx' and 'edu_idx' with the exact variables from your formula if needed.
  formula_no_wealth <- update(best_formula, . ~ . - factor(wealthq_idx))
  formula_no_edu <- update(best_formula, . ~ . - factor(edu_idx))

  # Ensure the formula environments are set to the global environment
  environment(formula_no_wealth) <- globalenv()
  environment(formula_no_edu) <- globalenv()

  # List of sensitivity models to fit
  sensitivity_models <- list(
    "No Wealth" = formula_no_wealth,
    "No Education" = formula_no_edu
  )

  sensitivity_results <- list()

  # Fit each sensitivity model
  for (model_name in names(sensitivity_models)) {
    cat("Running sensitivity analysis for:", model_name, "\n")
    formula <- sensitivity_models[[model_name]]

    # Try fitting the model
    result <- tryCatch({
      inla(
        formula = formula,
        data = data,
        control.predictor = best_model$.args$control.predictor,
        control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
      )
    }, error = function(e) {
      message("Error fitting model for ", model_name, ": ", e$message)
      return(NULL)
    })

    # If successful, save and store the results
    if (!is.null(result)) {
      output_file <- file.path(out_dir, region, sex, paste0(indicator, '-', tag_new, '-sensitivity-', model_name, '-fit.rds'))
      saveRDS(result, output_file)
      cat("Sensitivity analysis result saved for:", model_name, "at:", output_file, "\n")
      sensitivity_results[[model_name]] <- result
    } else {
      cat("Sensitivity analysis failed for:", model_name, "\n")
    }
  }

  return(sensitivity_results)
}

#####################################
## 02: Main Loop Over Indicators, Regions, Sexes
#####################################

for (indicator in indicators) {
  cat("Starting sensitivity analysis for indicator:", indicator, "\n")

  for (region in regions) {
    for (sex in sexes) {
      cat(sprintf("\nProcessing region: %s, sex: %s\n", region, sex))

      # Use tryCatch to handle errors gracefully
      tryCatch({
        results <- perform_sensitivity_analysis(
          region = region,
          sex = sex,
          indicator = indicator,
          tag = tag,
          tag_new = tag_new,
          out_dir = out_dir
        )
        if (!is.null(results)) {
          cat("Sensitivity analysis completed for", indicator, region, sex, "\n")
        }
      }, error = function(e) {
        message("An error occurred for ", indicator, " ", region, " ", sex, ": ", e$message)
      })
    }
  }
  cat("Completed sensitivity analysis for indicator:", indicator, "\n")
}
