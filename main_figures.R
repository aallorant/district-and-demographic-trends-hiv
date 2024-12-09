####################################################################################################
## Author: Adrien Allorant
##
## Description:
##   This script loads combined AMEs (Average Marginal Effects) results from previously saved RDATA files
##   and creates a visualization of expected differences in probabilities for certain health indicators 
##   (HIV status, ART coverage, recent HIV testing) across various demographic, socio-economic, and 
##   geographic factors, including changes over time.
##
## Requirements:
##   - R version >= 4.0.0
##   - Packages: data.table, tidyverse, lubridate, ggplot2, INLA, arm, tile, RColorBrewer,
##               brms, PNWColors, ggtext
##
## Inputs:
##   - RDATA files containing `combined_ames_across_time` and related objects for specified indicators.
##
## Outputs:
##   - PNG plots saved in the specified `fig_dir` under 'AMEs' subdirectory.
##
## Notes:
##   - Adjust the `main_dir`, `tag`, `indicators`, and plotting parameters as needed.
####################################################################################################

#####################################
## 00: Setup Environment and Paths ##
#####################################

# Load required libraries
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
library(ggtext)

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

# Create AME directory if it doesn't exist
ame_dir <- file.path(fig_dir, 'AMEs')
if (!dir.exists(ame_dir)) dir.create(ame_dir)

#####################################
## 01: Configuration Variables     ##
#####################################

# Define colors, tag, sexes, regions, and indicators
colors_region <- pnw_palette(name = 'Sailboat', n = 6, type = 'discrete')[c(1:3,5)]
tag <- 'August2nd'
regions <- c('Eastern Africa', 'Western Africa', 'Southern Africa', 'Central Africa')
sexes <- c('Men', 'Women')
indicators <- c('hivstatus', 'arv', 'test12m')
S <- 100
rerun <- FALSE  # Set to TRUE if you need to rerun computations; currently unused here

#####################################
## 02: Label and Type Mappings     ##
#####################################

# Map variable names to broader categories
type_mapping <- c(
  'edu_idx' = 'Socio-Economic',
  'wealthq_idx' = 'Socio-Economic',
  'res_idx' = 'Demographic',
  'age_idx' = 'Demographic',
  'period_idx' = 'Time',
  'iso3_idx' = 'Geographic',
  'area_idx' = 'Geographic'
)

# Map variable names to more descriptive labels
label_mapping <- c(
  'wealthq_idx' = 'Relative Wealth Quintile',
  'age_idx' = 'Age',
  'res_idx' = 'Place of Residence',
  'edu_idx' = 'Education',
  'iso3_idx' = 'Country',
  'area_idx' = 'District'
)

#####################################
## 03: Main Loop Over Indicators    ##
#####################################

for (indicator in indicators) {
  
  # Load combined AMEs data for the current indicator
  data_file <- file.path(out_dir, paste0(tag, '-', indicator, '-combined_ames_for_plot.RDATA'))
  load(data_file)
  
  # Adjust category labeling based on indicator (arv vs others)
  margEffect <- if (indicator != 'arv') {
    combined_ames_across_time %>%
      mutate(
        significance = ifelse(abs(mean_diff * 100) > 0.1 & (ci_lower > 0 | ci_upper < 0), "Significant", "Not Significant"),
        category = case_when(
          var == 'wealthq_idx' & level == 5 ~ 'Relative Wealth\nQ5 vs Q1',
          var == 'age_idx' & level == 3 ~ 'Age\n35-49 vs 15-24',
          var == 'res_idx' & level == 2 ~ 'Residence\nUrban vs Rural',
          var == 'edu_idx' & level == 4 ~ 'Education\nTertiary vs less than Primary',
          var == 'area_idx' ~ 'District intercepts',
          var == 'iso3_idx' ~ 'Country intercepts',
          var == 'period_idx' & level == 4 ~ 'Period\n2018-2023 vs 2003-2007',
          TRUE ~ NA_character_
        )
      )
  } else {
    combined_ames_across_time %>%
      mutate(
        significance = ifelse(abs(mean_diff * 100) > 0.1 & (ci_lower > 0 | ci_upper < 0), "Significant", "Not Significant"),
        category = case_when(
          var == 'wealthq_idx' & level == 5 ~ 'Relative Wealth\nQ5 vs Q1',
          var == 'age_idx' & level == 3 ~ 'Age\n35-49 vs 15-24',
          var == 'res_idx' & level == 2 ~ 'Residence\nUrban vs Rural',
          var == 'edu_idx' & level == 4 ~ 'Education\nTertiary vs less than Primary',
          var == 'area_idx' ~ 'District intercepts',
          var == 'iso3_idx' ~ 'Country intercepts',
          var == 'period_idx' & level == 2 ~ 'Period\n2018-2023 vs 2013-2017',
          TRUE ~ NA_character_
        )
      )
  }
  
  # Filter out comparisons not of interest and apply mapping
  margEffect <- margEffect %>%
    filter(!is.na(category)) %>%
    mutate(
      estimate_type = dplyr::recode(var, !!!type_mapping),
      var_labelled = dplyr::recode(var, !!!label_mapping)
    ) %>%
    group_by(region, sex) %>%
    mutate(category = forcats::fct_reorder(category, mean_diff, .desc = FALSE)) %>%
    ungroup()
  
  # Determine indicator title
  indicator_title <- c('HIV status', 'Recent HIV testing', 'ART coverage')[which(indicator == c('hivstatus','test12m','arv'))]
  
  #####################################
  ## 04: Plotting the Results        ##
  #####################################
  
  p <- ggplot(data = margEffect, aes(x = 100 * mean_diff, y = category, shape = significance)) +
    geom_point(aes(color = estimate_type), size = 2, position = position_jitter(width = 0.5, height = 0)) +
    geom_errorbar(aes(xmin = 100 * ci_lower, xmax = 100 * ci_upper, color = estimate_type),
                  width = 0.2, position = position_jitter(width = 0.5, height = 0)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text(aes(label = paste0(round(100 * mean_diff, 1), '%'), color = estimate_type),
              vjust = -.5, size = 3.5) +
    scale_shape_manual(values = c("Significant" = 19, "Not Significant" = 1)) +
    scale_x_continuous(labels = scales::percent_format(scale = 1),
                       limits = c(min(margEffect$ci_lower * 100) - 0.5, max(margEffect$ci_upper * 100) + 0.5)) +
    scale_color_manual(values = colors_region) +
    labs(
      title = paste0('Expected differences in ', indicator_title, ' by factor.'),
      subtitle = 'Given the changes in covariates listed at the left,\nthe difference in probability is expected to be...',
      y = '', x = 'Expected difference in probability (%)'
    ) +
    guides(shape = "none", color = 'none') +
    facet_grid(sex ~ region, scales = 'free') +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = rel(1.15)),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0.5, face = "italic"),
      legend.text = element_text(size = rel(1.05)),
      title = element_text(size = rel(1.25), face = "bold"),
      axis.text.y = element_text(size = rel(0.95), hjust = 0.5, vjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = rel(0.9)),
      strip.text = element_text(size = rel(1.05), angle = 0, face = 'bold'),
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "white", color = "black")
    )
  
  # Save the plot for the current indicator
  ggsave(plot = p,
         filename = file.path(ame_dir, paste0(tag, '-', indicator, '-ame-across-time.png')),
         width = 247, height = 247, units = "mm", dpi = 300)
}
