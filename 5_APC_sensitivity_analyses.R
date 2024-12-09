####################################################################################################
## Author: Adrien Allorant
##
## Description:
##    This script compares the average marginal effects (AMEs) of certain predictors (e.g., wealth or 
##    education) in a main model and sensitivity models where these predictors are omitted. 
##    The comparison helps understand how removing wealth or education from the model 
##    influences the estimated AMEs of other predictors.
##
## Requirements:
##    - R version >= 4.0.0
##    - Packages: data.table, tidyverse, lubridate, ggplot2, INLA, arm, tile, RColorBrewer,
##                brms, PNWColors, Hmisc
##
## Inputs:
##    - Results from previously fitted models ("best-logit-sae-fit.rds", "sensitivity-No Education-fit.rds",
##      "sensitivity-No Wealth-fit.rds") saved under the output directory structure.
##    - Corresponding survey data ("best-logit-sae.csv").
##
## Outputs:
##    - A combined RDATA file with AMEs from both main and sensitivity models.
##    - Example ggplot visualizations comparing AMEs across models.
##
## Notes:
##    - Adjust paths and environment detection as needed.
##    - Ensure that the main and sensitivity model files and data CSV files exist before running.
####################################################################################################

#####################################
## 00: Setup Environment and Paths ##
#####################################

# Load necessary libraries
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

# Set random seed for reproducibility
set.seed(123)
options(mc.cores = parallel::detectCores())

if (run_rstudio) {
  # Example: using a local directory structure
  root <- "."  # or set this to your project root using here::here()
} else {
  args <- commandArgs(trailingOnly = TRUE)
  root <- args[1]
  indicator <- args[3]
}

# Define data directories
post_proc_dir <- file.path(main_dir, 'prepped_data')
loc_areas_dir <- file.path(main_dir, 'geography_data')
out_dir       <- file.path(main_dir, 'output')
fig_dir       <- file.path(main_dir, 'figures')

#####################################
## 01: Configuration Variables     ##
#####################################

indicators <- c('hivstatus','arv','test12m')
sexes <- c("Men", "Women")
regions <- c("Eastern Africa", "Western Africa", "Southern Africa", "Central Africa")

S <- 100   # Number of posterior samples from INLA fits
rerun <- TRUE

# Tag used in previous scripts (adjust if needed)
tag <- 'August2nd'

#####################################
## 02: Supporting Functions        ##
#####################################

# This function calculates AMEs for a given predictor based on posterior samples from an INLA model.
# It sets the predictor to each of its levels in turn, calculates predicted probabilities, and 
# then computes differences in probabilities relative to a reference level.
calculate_ame <- function(samples, resdf, predictor) {
  # Extract intercept samples
  intercept_index <- grep("Intercept", rownames(samples[[1]]$latent))
  intercept_samples <- sapply(samples, function(x) x$latent[intercept_index])
  
  # Determine the levels of the predictor
  levels_ <- sort(unique(resdf[[predictor]]))
  ref_level <- levels_[1]  # Reference category is the first level
  
  ame_list <- list()
  S <- length(samples)
  
  # Loop over predictor levels
  for (lev in levels_) {
    # Modify dataset: set all observations to this predictor level
    modified_df <- resdf
    modified_df[[predictor]] <- lev
    modified_df <- modified_df[complete.cases(modified_df), ]
    
    N_mod <- nrow(modified_df)
    predicted_probabilities <- matrix(NA, nrow = N_mod, ncol = S)
    
    # Compute predicted probabilities for each posterior sample
    for (s in 1:S) {
      sample <- samples[[s]]
      
      # Initialize linear predictor eta. For simplicity, assume only intercept and 
      # the predictor of interest are considered. If other covariates are included, 
      # you must extend this logic.
      eta <- intercept_samples[s]
      
      # Add predictor's effect if level > reference
      if (lev > 1) {
        predictor_coef_name <- paste0(predictor, ":", lev)
        predictor_index <- which(rownames(sample$latent) == predictor_coef_name)
        
        # If direct coefficient not found, try partial matches (for interactions)
        if (length(predictor_index) == 0) {
          predictor_index <- grep(paste0(predictor, ".*:", lev), rownames(sample$latent))
          if (length(predictor_index) == 0) {
            stop(paste("Coefficient for", predictor_coef_name, "not found in sample"))
          }
        }
        eta <- eta + sample$latent[predictor_index]
      }
      
      # Convert linear predictor to probability
      p <- 1 / (1 + exp(-eta))
      predicted_probabilities[, s] <- p
    }
    
    # Mean probability per observation across all samples
    mean_prob <- rowMeans(predicted_probabilities)
    ame_list[[as.character(lev)]] <- mean_prob
  }
  
  # Compute AME as difference from reference level probabilities
  prob_ref <- ame_list[[as.character(ref_level)]]
  
  ame_results <- data.frame()
  for (lev in levels_[-1]) {
    prob <- ame_list[[as.character(lev)]]
    diff_prob <- prob - prob_ref
    ame_mean <- mean(diff_prob)
    ci <- quantile(diff_prob, probs = c(0.025, 0.975))
    ame_results <- rbind(ame_results, 
                         data.frame(level = lev, 
                                    ame = ame_mean, 
                                    ci_lower = ci[1], 
                                    ci_upper = ci[2]))
  }
  
  return(ame_results)
}


#####################################
## 03: Main Loop Over Indicators   ##
#####################################

for (indicator in indicators) {
  if (rerun) {
    combined_ames <- NULL
    
    # Loop over each sex and region
    for (sex. in sexes) {
      for (region in regions) {
        cat(sprintf("Processing %s in %s for %s...\n", sex., region, indicator))
        
        # File paths for main and sensitivity models
        main_resdf_path <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, "-best-logit-sae.csv"))
        main_modfit_path <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, "-best-logit-sae-fit.rds"))
        
        no_edu_modfit_path <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, "-sensitivity-No Education-fit.rds"))
        no_wealth_modfit_path <- file.path(out_dir, region, sex., paste0(indicator, '-', tag, "-sensitivity-No Wealth-fit.rds"))
        
        # Check file existence
        if (!(file.exists(main_modfit_path) && file.exists(no_edu_modfit_path) && file.exists(no_wealth_modfit_path) && file.exists(main_resdf_path))) {
          cat("One or more model/data files are missing. Skipping...\n")
          next
        }
        
        # Load models and data
        main_modfit <- readRDS(main_modfit_path)
        no_edu_modfit <- readRDS(no_edu_modfit_path)
        no_wealth_modfit <- readRDS(no_wealth_modfit_path)
        resdf <- read_csv(main_resdf_path)
        
        # Keep only necessary columns
        resdf <- resdf %>% 
          dplyr::select(iso3_idx, area_idx, period_idx, age_idx, wealthq_idx, edu_idx, res_idx, N)
        
        # Extract posterior samples
        main_samples <- inla.posterior.sample(n = S, result = main_modfit)
        no_edu_samples <- inla.posterior.sample(n = S, result = no_edu_modfit)
        no_wealth_samples <- inla.posterior.sample(n = S, result = no_wealth_modfit)
        
        # Compute AMEs for wealthq_idx with main model and no-edu model
        if ("wealthq_idx" %in% colnames(resdf)) {
          main_wealth_ame <- calculate_ame(main_samples, resdf, "wealthq_idx") %>%
            mutate(estimate_type = "Main", var = "wealthq_idx", sex = sex., region = region)
          
          no_edu_wealth_ame <- calculate_ame(no_edu_samples, resdf, "wealthq_idx") %>%
            mutate(estimate_type = "No Education", var = "wealthq_idx", sex = sex., region = region)
          
          wealth_ame_combined <- bind_rows(main_wealth_ame, no_edu_wealth_ame)
          combined_ames <- bind_rows(combined_ames, wealth_ame_combined)
        }
        
        # Compute AMEs for edu_idx with main model and no-wealth model
        if ("edu_idx" %in% colnames(resdf)) {
          main_edu_ame <- calculate_ame(main_samples, resdf, "edu_idx") %>%
            mutate(estimate_type = "Main", var = "edu_idx", sex = sex., region = region)
          
          no_wealth_edu_ame <- calculate_ame(no_wealth_samples, resdf, "edu_idx") %>%
            mutate(estimate_type = "No Wealth", var = "edu_idx", sex = sex., region = region)
          
          edu_ame_combined <- bind_rows(main_edu_ame, no_wealth_edu_ame)
          combined_ames <- bind_rows(combined_ames, edu_ame_combined)
        }
      }
    }
    
    # Save combined results
    save(combined_ames, file = file.path(out_dir, paste0(tag, '-combined_ames.RDATA')))
  } else {
    # If not rerunning, just load results
    load(file.path(out_dir, paste0(tag, '-combined_ames.RDATA')))
  }
}


#####################################
## 04: Visualization (Example)     ##
#####################################

# Load combined AMEs
load(file.path(out_dir, paste0(tag, '-combined_ames.RDATA')))

# Example plot for wealthq_idx
wealth_ame <- combined_ames %>% filter(var == 'wealthq_idx')

ggplot(wealth_ame, aes(x = factor(level), y = ame * 100, color = estimate_type, group = estimate_type)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100), width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ sex + region, scales = "free") +
  labs(title = 'AMEs for wealthq_idx across Main and No Education Models',
       x = 'Wealth Quintile Level',
       y = 'Average Marginal Effect (%)') +
  theme_bw() +
  theme(legend.position = "bottom")

# Example plot for edu_idx
edu_ame <- combined_ames %>% filter(var == 'edu_idx')

ggplot(edu_ame, aes(x = factor(level), y = ame * 100, color = estimate_type, group = estimate_type)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100), width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ sex + region, scales = "free") +
  labs(title = 'AMEs for edu_idx across Main and No Wealth Models',
       x = 'Education Level',
       y = 'Average Marginal Effect (%)') +
  theme_bw() +
  theme(legend.position = "bottom")