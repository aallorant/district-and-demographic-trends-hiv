####################################################################################################
## Author: Adrien Allorant
##
## Description:
##    This script prepares and aggregates HIV-related survey data from multiple sources.
##    It reads geolocated individual-level survey datasets (previously pre-processed),
##    harmonizes and aggregates them by relevant demographic and geographic strata.
##
## Requirements:
##    - R version >= 4.0.0
##    - Packages: sf, data.table, foreign, tidyverse, Matrix, stringr
##    - Input data: 
##      * Pre-processed survey data files and geographic shapefiles stored in data directories.
##      * A merged area shapefile (`areas.geojson`) with hierarchical identifiers.
##
## Outputs:
##    - Aggregated CSV or Rdata files of HIV indicators (`test12m`, `hivstatus`, `arv`) by 
##      geographic units, year, sex, age groups, and other dimensions.
##
## Notes:
##    - The code is modularized: configuration, data loading, processing, and saving steps 
##      are clearly separated.
##    - Adjust directory paths as needed.
####################################################################################################


#####################################
## 00: Setup and Configuration      ##
#####################################

## Clear environment (optional; not always best practice in scripts)
# rm(list = ls())

## Load required packages
library(sf)
library(data.table)
library(foreign)
library(tidyverse)
library(Matrix)
library(stringr)

## Configuration Flags
run_rstudio <- TRUE   # Set to FALSE if running from command line with arguments
rerun <- FALSE         # If TRUE, will re-load and re-process raw data sources

## Indicators to process
indicators <- c('test12m', 'hivstatus', 'arv')

## Define root directory and relative paths
## Adjust these paths to match your repository structure.
if (run_rstudio) {
  # Example: using a local directory structure
  root <- "."  # or set this to your project root using here::here()
} else {
  args <- commandArgs(trailingOnly = TRUE)
  root <- args[1]
  indicator <- args[3]
}

## Relative paths (adjust as needed)
data_dir        <- file.path(root, "data")
processed_dir   <- file.path(data_dir, "Processed")
geography_dir   <- file.path(data_dir, "geography_data")
prepped_data_dir <- file.path(data_dir, "prepped_data")
code_dir        <- file.path(root, "code")
model_dir       <- file.path(code_dir, "stan_models")

## Input data files
areas_file      <- file.path(data_dir, "areas.geojson")
matched_tot_file <- file.path(data_dir, "data_matched_tot.rds")
geomatched_file  <- file.path(data_dir, "geomatched_svy_dataset.rds")
moz2022_file     <- file.path(data_dir, "moz_2022_geomatched.rds")
new_phia_file    <- file.path(data_dir, "new_phia_geomatched.rds")

## Ensure that output directories exist (optional)
dir.create(prepped_data_dir, showWarnings = FALSE, recursive = TRUE)


#####################################
## 01: Load and Pre-Process Data    ##
#####################################

## Load merged area-level data
area_merged <- setDT(read_sf(areas_file))

if (rerun) {
  # If re-running data extraction and merging:
  # 1. Load a combined survey dataset
  data_matched_tot <- setDT(readRDS(geomatched_file)) %>%
    dplyr::select(-c(aware, artself, vls, recent, cd4)) %>%
    filter(!(survey_id == 'BWA2021BAIS'))

  # Append Mozambique 2022 data
  moz2022 <- readRDS(moz2022_file)
  data_matched_tot <- data_matched_tot %>% 
    bind_rows(moz2022 %>%
      mutate(
        country = 'Mozambique',
        iso3 = 'MOZ',
        sex = as.numeric(as.factor(sex)),
        aware = NA,
        artself = NA,
        recent = NA,
        cd4 = NA,
        vls = NA,
        survey_region_id = region,
        longitude = NA, latitude = NA, geoloc_distance = NA
      ) %>%
      dplyr::select(all_of(colnames(data_matched_tot))) %>%
      st_drop_geometry()
    )

  # Append additional PHIA data
  new_phia_matched <- readRDS(new_phia_file) %>%
    rename(longitude = center_x, latitude = center_y, geoloc_area_id = area_id) %>%
    dplyr::select(-c(geometry, area_level))
  
  data_matched_tot <- bind_rows(data_matched_tot, new_phia_matched)

  # Clean up missing ISO3 codes
  data_matched_tot[is.na(iso3) & survey_id == 'MWI2020PHIA', iso3 := 'MWI']

  # Filter out rows where geoloc_area_id is missing
  data_matched_tot <- data_matched_tot %>%
    filter(!is.na(geoloc_area_id)) %>%
    mutate(area_level = as.integer(str_split(geoloc_area_id, "_", simplify = TRUE)[,2]))

  saveRDS(data_matched_tot, file = matched_tot_file)
  
} else {
  # If not re-running, load pre-saved combined dataset
  data_matched_tot <- readRDS(matched_tot_file)
}


#####################################
## 02: Harmonize Geographic Levels  ##
#####################################

## Identify minimum area_level per ISO3
min_level_by_iso3 <- data_matched_tot %>%
  group_by(iso3) %>%
  summarise(min_level = min(area_level, na.rm = TRUE))

## Distinct area-level info
distinct_data <- data_matched_tot %>%
  left_join(min_level_by_iso3, by = "iso3") %>%
  dplyr::select(iso3, geoloc_area_id, area_level, min_level) %>%
  distinct()

## Recursive function to find lowest parent area ID at a specified level
find_lowest_parent_id <- function(area_id, area_df, target_level) {
  if (is.na(area_id)) return(NA)
  
  current_info <- area_df %>% filter(area_id == !!area_id) %>% head(1)
  
  if (nrow(current_info) == 0 || current_info$area_level <= target_level) {
    return(area_id)
  } else {
    parent_id <- current_info$parent_area_id
    return(find_lowest_parent_id(parent_id, area_df, target_level))
  }
}

## Apply function to match geoloc_area_id to the lowest available level
matched_distinct_data <- distinct_data %>%
  group_by(iso3) %>%
  rowwise() %>%
  mutate(matched_geoloc_area_id = find_lowest_parent_id(geoloc_area_id, area_merged, min_level)) %>%
  ungroup()

data_matched_tot <- data_matched_tot %>%
  left_join(
    matched_distinct_data %>% dplyr::select(iso3, geoloc_area_id, matched_geoloc_area_id), 
    by = c("iso3", "geoloc_area_id")
  ) %>%
  mutate(geoloc_area_id = coalesce(matched_geoloc_area_id, geoloc_area_id)) %>%
  dplyr::select(-matched_geoloc_area_id)


#####################################
## 03: Define Age Groups and Other  ##
#####################################

# Collapse age groups
data_matched_tot[age >= 15 & age < 25, age_start := 15]
data_matched_tot[age >= 25 & age < 35, age_start := 25]
data_matched_tot[age >= 35 & age < 50, age_start := 35]
data_matched_tot[age >= 50, age_start := 50]

# Extract year from survey_id
data_matched_tot[, year := str_extract(survey_id, '[0-9]+')]
data_matched_tot[, survey_type := sapply(str_split(survey_id, "[0-9]+"), "[[", 2)]

# Normalize sex variable
data_matched_tot[sex %in% c('1','ale','male'), sex_num := 1]
data_matched_tot[sex %in% c('2','female'), sex_num := 2]
data_matched_tot[, sex := sex_num]

# Keep only non-missing age, residence type, and sex
data_matched_tot <- data_matched_tot[!is.na(age_start) & !is.na(res_type) & !is.na(sex)]

# Normalize HIV status
data_matched_tot[hivstatus %in% c('positive','1'), hivstatus_num := 1]
data_matched_tot[hivstatus %in% c('negative','0'), hivstatus_num := 0]
data_matched_tot[, hivstatus := hivstatus_num]

## Recalculate wealth quintiles by area
mean_wealth_area <- data_matched_tot[, .(mean_wealth = mean(wealths, na.rm = TRUE)),
                                     by = c('res_type','iso3','survey_id')]
data_matched_tot <- merge(data_matched_tot, mean_wealth_area, by = c('res_type','iso3','survey_id'))
data_matched_tot[, wealths_adj := wealths - mean_wealth]
data_matched_tot[, wealthq_adj := ntile(wealths_adj, 5), by = c('res_type','iso3','survey_id')]

## Recalculate weights
weights <- data_matched_tot[, .(M = sum(!is.na(indweight)),
                                Ms = sum(indweight, na.rm = TRUE)^2 / sum(indweight^2, na.rm = TRUE),
                                w = sum(indweight, na.rm = TRUE),
                                m_w = mean(indweight, na.rm = TRUE),
                                Mhiv = sum(!is.na(hivweight)),
                                Mhiv_s = sum(hivweight, na.rm = TRUE)^2 / sum(hivweight^2, na.rm = TRUE),
                                whiv = sum(hivweight, na.rm = TRUE),
                                m_whiv = mean(hivweight, na.rm = TRUE)),
                             by = c('iso3','survey_id','res_type','survey_region_id')]

data_matched_tot <- merge(data_matched_tot, weights, by = c('iso3','survey_id','res_type','survey_region_id'))
data_matched_tot[, indweight_k := (indweight/m_w)*(M/Ms)]
data_matched_tot[, hivweight_k := (hivweight/m_whiv)*(Mhiv/Mhiv_s)]


#####################################
## 04: Aggregate Indicators         ##
#####################################

## Function to aggregate a given indicator
aggregate_indicator <- function(ind_name, dt, adjusted_wealth = TRUE) {
  cat("Processing indicator: ", ind_name, "\n")
  
  # Decide which weight to use based on the indicator
  weight_var <- if (ind_name == 'test12m') 'indweight_k' else 'hivweight_k'
  
  # Choose wealth variable (adjusted or original)
  wealth_var <- if (adjusted_wealth) 'wealthq_adj' else 'wealthq'
  
  # Subset columns
  keep_vars <- c("iso3","sex","geoloc_area_id","survey_region_id","survey_id","res_type",
                 "age_start","edu", wealth_var, weight_var, ind_name)
  
  dt_sub <- dt[, ..keep_vars]
  setnames(dt_sub, weight_var, "weight_k")
  
  # Melt data for aggregation
  value_vars <- ind_name
  id_vars <- setdiff(colnames(dt_sub), value_vars)
  data_melted <- melt.data.table(dt_sub, id.vars = id_vars, measure.vars = value_vars)
  
  # Aggregate data
  data_aggregated <- data_melted[!is.na(weight_k), .(
    Y = sum(as.numeric(value)*weight_k, na.rm = TRUE),
    N = sum(as.numeric(!is.na(value))*weight_k, na.rm = TRUE)
  ), by = c('geoloc_area_id','survey_region_id','iso3','survey_id','sex','age_start','res_type','edu', wealth_var)]
  
  data_aggregated[, indicator := ind_name]
  data_aggregated <- data_aggregated[N > 0, ]
  
  # Save results
  suffix <- if (adjusted_wealth) "" else "_wealthq_original"
  save_file <- file.path(prepped_data_dir, paste0(ind_name, suffix, ".rdata"))
  save(data_aggregated, file = save_file)
  
  cat("Saved aggregated data for ", ind_name, " to ", save_file, "\n")
}

## Aggregate indicators with adjusted wealth quintiles
for (ind in indicators) {
  aggregate_indicator(ind, data_matched_tot, adjusted_wealth = TRUE)
}

## Aggregate indicators with original wealth quintiles
## Note: This assumes that 'wealthq' exists in the original dataset.
## If not, you may need to adjust code above to retain original wealthq.
for (ind in indicators) {
  aggregate_indicator(ind, data_matched_tot, adjusted_wealth = FALSE)
}
