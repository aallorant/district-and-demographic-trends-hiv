####################################################################################################
## Author: Adrien Allorant
##
## Description:
##   This script computes average predictive comparisons (APCs) initially called average marginal effects (AMEs)
##   from previously fitted hierarchical logistic models. It:
##     - Identifies the best-fitting model for each region, sex, and indicator.
##     - Samples from the model posterior to compute predicted probabilities.
##     - Calculates AMEs across various covariates and their levels.
##     - Saves and plots the results (e.g., odds ratios, marginal effects) for further analysis.
##
## Requirements:
##   - R version >= 4.0.0
##   - Packages: data.table, tidyverse, lubridate, ggplot2, INLA, arm, tile, RColorBrewer,
##               brms, PNWColors, Hmisc, gridExtra, knitr, kableExtra
##
## Inputs:
##   - Prepped data files generated in previous steps (CSV and RDS files for model fits).
##   - Sector and geography data (referenced in script).
##
## Outputs:
##   - RDATA files containing combined AME/OR results.
##   - PNG plots for OR forest plots, AME plots across time and categories.
##   - HTML tables summarizing OR results.
##
## Notes:
##   - Adjust file paths, indicator sets, and region sets as per your environment.
##   - The code includes extensive transformations and plots; customize as needed.
####################################################################################################

#####################################
## 00: Setup Environment and Paths ##
#####################################

# Remove environment clearing for reproducibility
# rm(list=ls())

library(data.table)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(INLA)
library(arm)
library(tile)
library(RColorBrewer)
library(brms)
library(PNWColors)
library(Hmisc)
library(gridExtra)
library(knitr)
library(kableExtra)

# Set main directory based on environment
if (run_rstudio) {
  # Example: using a local directory structure
  root <- "."  # or set this to your project root using here::here()
} else {
  args <- commandArgs(trailingOnly = TRUE)
  root <- args[1]
  indicator <- args[3]
}

# Define directories
post_proc_dir <- file.path(main_dir, 'prepped_data')
loc_areas_dir <- file.path(main_dir, 'geography_data')
out_dir       <- file.path(main_dir, 'output')
fig_dir       <- file.path(main_dir, 'figures')

# Create figure subdirectory for AMEs if doesn't exist
ame_dir <- file.path(fig_dir, 'AMEs')
if (!dir.exists(ame_dir)) dir.create(ame_dir)

# Define plotting colors
colors_region <- pnw_palette(name = 'Sailboat', n = 6, type = 'discrete')[c(1:3,5)]

#####################################
## 01: Configuration Variables     ##
#####################################

tag <- 'August2nd'
tag_new <- Sys.Date()

regions <- c('Eastern Africa', 'Western Africa', 'Southern Africa', 'Central Africa')
sexes <- c('Men', 'Women')
S <- 100 # number of samples from posterior
indicators <- c('hivstatus', 'arv', 'test12m')
rerun <- TRUE

##########################################
## 02: Main Computation Loop per Indicator
##########################################

for (indicator in indicators) {
  cat(sprintf("Processing indicator: %s\n", indicator))
  
  if (rerun) {
    combined_ors <- NULL
    combined_ames <- NULL
    combined_ames_across_time <- NULL
    
    # Loop over sexes and regions
    for (sex. in sexes) {
      for (region in regions) {
        cat(sprintf("  Processing %s in %s...\n", sex., region))
        
        # Identify best model by WAIC
        ic_file <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, '-logit-information-criteria.csv'))
        ic_df <- read_csv(ic_file)
        best_model <- ic_df %>% filter(waic == min(waic)) %>% pull(model)
        if(length(best_model)>1) best_model <- best_model[1]
        
        # Load data and fits
        resdf_path <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, '-best-logit-sae.csv'))
        modfit_path <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, '-best-logit-sae-fit.rds'))
        
        resdf <- read_csv(resdf_path)
        modfit <- readRDS(modfit_path)
        
        # Posterior samples from INLA fit
        samples <- inla.posterior.sample(n = S, result = modfit)
        
        # Extract relevant data fields
        resfit <- resdf[, c('iso3_idx', 'area_idx', 'period_idx', 'age_idx', 'wealthq_idx', 'edu_idx', 'res_idx','N')]
        predictors <- c('wealthq_idx', 'edu_idx', 'age_idx', 'res_idx', 'period_idx', 'iso3_idx', 'area_idx')
        
        #########################
        ## Extract U Effects   ##
        #########################
        
        # Extract samples for each predictor's random effects
        extracted_samples <- list()
        ors_list <- NULL
        
        for (var in predictors) {
          # Identify correct indices based on model and var
          var_pattern <- paste0(var)
          
          var_samples <- matrix(
            unlist(lapply(samples, function(x) {
              x$latent[grep(var_pattern, rownames(x$latent))]
            })),
            nrow = S, byrow = TRUE
          )
          
          extracted_samples[[var]] <- var_samples
          
          # Compute OR summary
          or_means <- apply(var_samples, 2, mean)
          or_cis <- t(apply(var_samples, 2, quantile, probs = c(0.025, 0.975)))
          ors_list <- bind_rows(ors_list,
                                data.frame(var = rep(var, nrow(or_cis)),
                                           level = 1:nrow(or_cis),
                                           OR = exp(or_means),
                                           CI_Lower = exp(or_cis[, 1]),
                                           CI_Upper = exp(or_cis[, 2])))
        }
        
        # Extract intercept
        intercept_index <- grep("Intercept", rownames(samples[[1]]$latent))
        intercept_samples <- unlist(lapply(samples, function(x) x$latent[intercept_index]))
        
        #################################################
        ## Average Marginal Effects (AMES) Calculation ##
        #################################################
        
        ame_df <- resfit[complete.cases(resfit[, c(predictors, 'N')]), ]
        ame_df$area_idx <- match(ame_df$area_idx, unique(ame_df$area_idx))
        
        # Function to compute weights for AME
        compute_weights <- function(df, target_variable) {
          df <- setDT(df)
          other_vars <- setdiff(names(df), c(target_variable, "area_idx", "iso3_idx","N"))
          df <- df[, n := sum(N), by = c(other_vars)]
          df <- df[, weight := n / sum(n), by = 'period_idx']
          return(df$weight)
        }
        
        ames <- NULL
        ames_across_time <- NULL
        
        # Variables to compute AME on
        ame_vars <- c('wealthq_idx', 'edu_idx', 'age_idx', 'res_idx','period_idx', 'iso3_idx', 'area_idx')
        
        # Loop over each target variable
        for (target_variable in ame_vars) {
          levels_ <- sort(unique(ame_df[[target_variable]]))
          if(length(levels_) == 1) next
          
          # Probabilities for each level
          prob_list <- list()
          
          for (level_ in levels_) {
            modified_df <- ame_df
            modified_df[[target_variable]] <- level_
            expanded_df <- distinct(modified_df)
            expanded_df <- expanded_df[!is.na(expanded_df$area_idx), ]
            
            expanded_df$weight <- compute_weights(expanded_df, target_variable)
            
            # Construct linear predictor from samples
            # For simplicity, use a generic formula. If the model differs, adjust accordingly.
            linear_predictor <- intercept_samples +
              extracted_samples[['iso3_idx']][, expanded_df$iso3_idx] +
              extracted_samples[['area_idx']][, expanded_df$area_idx] +
              extracted_samples[['age_idx']][, expanded_df$age_idx] +
              extracted_samples[['res_idx']][, expanded_df$res_idx] +
              extracted_samples[['period_idx']][, expanded_df$period_idx] +
              extracted_samples[['edu_idx']][, expanded_df$edu_idx] +
              extracted_samples[['wealthq_idx']][, expanded_df$wealthq_idx]
            
            prob_list[[level_]] <- 1 / (1 + exp(-linear_predictor))
          }
          
          # Compute differences from reference level
          if (target_variable %in% c('iso3_idx', 'area_idx')) {
            # Multi-category varying intercepts: compute mean absolute differences
            prob_matrix <- do.call(cbind, prob_list)
            diff_prob <- combn(levels_, 2, function(pair) {
              mean(abs(prob_matrix[, pair[1]] - prob_matrix[, pair[2]]))
            })
            
            summary_diff_prob <- data.frame(
              mean_diff = mean(diff_prob),
              ci_lower = quantile(diff_prob, 0.1),
              ci_upper = quantile(diff_prob, 0.9),
              var = target_variable
            )
            
            ames <- bind_rows(ames, summary_diff_prob)
            ames_across_time <- bind_rows(ames_across_time, summary_diff_prob)
            
          } else {
            # Binary or ordered categories: differences from the first category as reference
            prob_ref <- prob_list[[levels_[1]]]
            for (lev in levels_[-1]) {
              prob <- prob_list[[lev]]
              diff_prob <- prob - prob_ref
              
              expanded_df$diff_prob <- diff_prob
              summary_diff_prob <- expanded_df %>%
                group_by(period_idx) %>%
                summarize(
                  mean_diff = weighted.mean(diff_prob, weight),
                  ci_lower = wtd.quantile(diff_prob, weights = weight, probs = 0.05, normwt = TRUE),
                  ci_upper = wtd.quantile(diff_prob, weights = weight, probs = 0.95, normwt = TRUE),
                  .groups = 'drop'
                ) %>%
                mutate(var = target_variable, level = lev)
              
              ames <- bind_rows(ames, summary_diff_prob)
              
              # Across time average
              ames_across_time <- bind_rows(ames_across_time,
                                            data.frame(
                                              mean_diff = mean(diff_prob),
                                              ci_lower = quantile(diff_prob, 0.025),
                                              ci_upper = quantile(diff_prob, 0.975),
                                              var = target_variable,
                                              level = lev
                                            ))
            }
          }
        }
        
        # Combine results for this region/sex
        combined_ors <- bind_rows(combined_ors, ors_list %>% mutate(sex = sex., region = region))
        combined_ames <- bind_rows(combined_ames, ames %>% mutate(sex = sex., region = region))
        combined_ames_across_time <- bind_rows(combined_ames_across_time,
                                               ames_across_time %>% mutate(sex = sex., region = region))
        
      } # end region
    } # end sex
    
    # Save combined results
    save(combined_ames, combined_ames_across_time, combined_ors,
         file = file.path(out_dir, paste0(tag, '-', indicator, '-combined_ames_for_plot.RDATA')))
    
  } else {
    # If not rerunning, just load results
    load(file.path(out_dir, paste0(tag, '-', indicator, '-combined_ames_for_plot.RDATA')))
  }
  
  # End of indicator loop
}
