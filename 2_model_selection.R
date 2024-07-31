rm(list=ls())

library(data.table)
library(dplyr)
library(stringr)
library(tidyverse)
library(INLA)
library(countrycode)
library(sf)
# library(rnaturalearth)
# library(rnaturalearthdata)
library(haven)
library(PNWColors)

if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}
library(R.utils)

# Determine the main directory based on the environment
if (file.exists('/cloud/project')) {
  main_dir <- '/cloud/project'
} else if (file.exists('/Users/adrienallorant')) {
  main_dir <- '/Users/adrienallorant/Documents/McGill/Research/HIV_inequality_paper'
} else if (file.exists('/Users/aallorant')) {
  main_dir <- '/Users/aallorant/Desktop/HIV_inequality_paper'
} else {
  stop("Unrecognized environment: Ensure the directory structure is correct.")
}

post_proc_dir <- paste0(main_dir, '/prepped_data')
loc_areas_dir <- paste0(main_dir, '/geography_data')
out_dir <- paste0(main_dir, '/output')
fig_dir <- paste0(main_dir, '/figures')

if(!dir.exists(out_dir)) dir.create(out_dir)
if(!dir.exists(fig_dir)) dir.create(fig_dir)

pal <- rev(pnw_palette("Sailboat", 35))

indicators <- c('test12m','hivstatus','arv')

source(paste0(main_dir,'/hiv_inequality_paper/fit_logistic.R'))

prg_dat <- setDT(readRDS(paste0(main_dir, '/data/recoded_prg_final/all_naomi.rds')))

prg_dat_wide <- dcast(prg_dat, 
                      area_level + area_id + area_name + sex + age_group_label + quarter_label ~ indicator, 
                      value.var = "mean")

tag <- 'july31st'

# for (miss in c('complete','with_missing_data')) {
  
  # for (indicator in indicators) {
    indicator <- 'hivstatus'
    cat('indicator ', indicator)
    
    # read in the data
    load(paste0(post_proc_dir,'/' ,indicator, ".rdata"))
    cat('indicator ', indicator, 'observations ', sum(data_aggregated$N))
    
    
    sector <- readRDS(paste0(main_dir,'/sector_key.rds'))
    sector$iso3 <- countrycode(sector$country, "country.name","iso3c")
    sector <- sector %>%
      rbind(data.frame(sector = 'Southern Africa', country = 'Botswana', iso3 = 'BWA'))
    
    data_aggregated <- data_aggregated %>%
      mutate(iso3 = ifelse(iso3 == 'SAB', 'ZAF', iso3),
             year = as.numeric(str_extract(survey_id, '[0-9]+')),
             period = case_when(
               year < 2008 ~ 1,
               year < 2013 ~ 2,
               year < 2018 ~ 3,
               TRUE ~ 4
             )) %>%
      left_join(sector)
    
    area_merged <- setDT(read_sf(paste0(loc_areas_dir,"/areas.geojson")))
    if('MDG' %in% unique(data_aggregated$iso3)){
      area_merged <- area_merged %>%
        bind_rows(setDT(read_sf(paste0(loc_areas_dir,"/gadm41_MDG_shp/gadm41_MDG_2.shp"))) %>%
                    dplyr::select(NAME_2, geometry) %>%
                    mutate(iso3 = 'MDG', area_level = 2, area_name = tolower(NAME_2),
                           geoloc_area_id = paste0(iso3, '_2_', match(area_name, unique(area_name))),
                           area_id = geoloc_area_id) %>%
                    dplyr::select(-NAME_2))
      
      region_to_area <- unique(data_aggregated %>%
                                 filter(iso3 == 'MDG') %>%
                                 dplyr::select(survey_region_id))
    }
    
    # Extract country code and area level
    data_aggregated  <-  data_aggregated %>%
      mutate(area_level = as.numeric(sub(".*_([0-9]+)_.*", "\\1", geoloc_area_id)),
             area_level = ifelse(is.na(area_level), 1, area_level))
    
    # Group by country and year, and find the highest area level for each country
    lowest_area_level <- data_aggregated %>%
      group_by(iso3) %>%
      summarize(lowest_area_level = min(area_level))
    
    # Filter the original data to keep only rows with the lowest area level for each country-year
    areas_filtered <- area_merged %>%
      inner_join(lowest_area_level, #%>% filter(iso3 != 'CIV'),
                 by = c("iso3","area_level" = "lowest_area_level")) %>%
      bind_rows(area_merged %>%
                  right_join(data_aggregated %>% filter(iso3 == 'CIV' & geoloc_area_id != 'CIV') %>%
                               dplyr::select(iso3, geoloc_area_id) %>% distinct(), by = c("iso3","area_id" = "geoloc_area_id"))) %>%
      #' Add an integer index for INLA
      arrange(area_sort_order) %>%
      mutate(area_idx = row_number()) %>%
      #' Add an integer index for INLA
      arrange(area_sort_order) %>%
      mutate(area_idx = row_number())

    # Step 1: Clean the Geometry
    areas_clean <- st_make_valid(areas_filtered %>% st_as_sf()) # Fix invalid geometries
    areas_clean <- st_simplify(areas_clean, preserveTopology = TRUE) # Simplify geometries
    
    rm(areas_filtered)
    # Step 2: Check for Validity
    if (all(st_is_valid(areas_clean))) {
      # Step 3: Create the Adjacency Matrix
      adjM <- spdep::poly2nb(st_as_sf(areas_clean))
    } else {
      stop("Some geometries are still invalid.")
    }
    
    adjM <- spdep::nb2mat(adjM, style = "B", zero.policy = TRUE)
    colnames(adjM) <- rownames(adjM)
    
    # for (region in c('Central Africa', 'Eastern Africa', 'Southern Africa', 'Western Africa')) {
      region <- 'Central Africa'
      print(region)
      
      dt <- data_aggregated %>%
        filter(sector == region)
      
      print(nrow(dt))
      
      
      if(!dir.exists(paste0(out_dir, '/', region))) 
        dir.create(paste0(out_dir, '/', region))
      
      # for (sex. in unique(dt$sex)) {
        sex. <- 2
        print(c('Men', 'Women')[sex.])
        if(!dir.exists(paste0(out_dir, '/', region, '/', c('Men', 'Women')[sex.]))) 
          dir.create(paste0(out_dir, '/', region, '/', c('Men', 'Women')[sex.]))
        
        tmp_dir <- paste0(out_dir, '/', region, '/', c('Men', 'Women')[sex.])
        dt_s <- dt %>%
          filter(sex == sex.) %>%
          mutate(
            Y = round(as.numeric(zap_labels(Y))),
            N = round(as.numeric(zap_labels(N)))
          )
        
        # Create a mapping for each variable based on its unique values
        iso3_map <- unique(dt_s$iso3)
        area_map <- unique(dt_s$geoloc_area_id)
        age_map <- sort(unique(dt_s$age_start))
        res_map <- unique(dt_s$res_type)
        edu_map <- sort(unique(dt_s$edu))
        wealthq_map <- sort(unique(dt_s$wealthq_adj))
        sur_map <- unique(dt_s$survey_id)
        year_map <- unique(dt_s$year)
        period_map <- sort(unique(dt_s$period))
        
        # if(region != 'Southern Africa'){
        dt_s <- dt_s %>%
          # filter(!is.na(edu), !is.na(wealthq_adj),
          #        !is.na(age_start), !is.na(res_type)) %>%
          mutate(
            iso3_idx = match(iso3, iso3_map),
            area_idx = match(geoloc_area_id, area_map),
            age_idx = match(age_start, age_map),
            res_idx = match(res_type, res_map),
            edu_idx = match(edu, edu_map),
            wealthq_idx = match(wealthq_adj, wealthq_map),
            # sur_idx = match(survey_id, sur_map),
            # period_idx = match(period, period_map),
            period_idx = match(period, period_map)
          )
        
        #' Create the scaffolding for the estimates
        df <- crossing(
          
          #' All of the different periods
          # period_idx = unique(dt_s$period_idx),
          period_idx = unique(dt_s$period_idx),
          #' Predictors
          sex = unique(dt_s$sex),
          
          res_idx = unique(dt_s$res_idx),
          
          wealthq_idx = unique(dt_s$wealthq_idx),
          
          edu_idx = unique(dt_s$edu_idx),
          
          age_idx = unique(dt_s$age_idx),
          
          #' Both the areas in the model and the aggregate country
          areas_clean %>%
            st_drop_geometry() %>%
            dplyr::select(iso3, area_id, area_name, area_level,
                          area_sort_order, center_x, center_y) %>%
            filter(iso3 %in% unique(dt_s$iso3))) #%>%
        
        #' Merge district observations into df
        df <- df %>%
          left_join(
            dt_s %>%
              dplyr::select(area_id = geoloc_area_id, survey_id, age_idx,
                            sex, age_idx,
                            edu_idx, wealthq_idx, res_idx,
                            period_idx, #year_idx,
                            Y, N),
            by = c("sex","res_idx","edu_idx","wealthq_idx",
                   "age_idx","area_id","period_idx")
          )

        df$urban <- ifelse(df$res_idx == 1, 1, 0)

        ## Period -- factor
        df$period <- scale(df$period_idx)
        
        #' Check that all the surveys are into df
        #' (The +1 is for the NA that's included in df but not ind)
        length(unique(df$survey_id)) == length(unique(dt_s$survey_id)) + 1 
        
        rm(dt_s)
        
        
        df <- df %>%
          mutate(
            iso3_idx = match(iso3, iso3_map),
            area_idx = match(area_id, area_map),
          )
        
        
        # Baseline formula
        formula_baseline <- Y ~ 1 + urban + factor(age_idx) +
          factor(edu_idx) + factor(wealthq_idx) +
            f(period_idx, model = "iid",
            constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
          f(iso3_idx, model = "iid",
            constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))

          # adding iid random effects for areas
          formula1 <- update(formula_baseline,
                               . ~ . + f(area_idx, model = "iid", constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))))
            
          # adding area-year and country-year interaction
          formula2 <- update(formula_baseline,
                               . ~ . -f(iso3_idx, model = "iid",
                                        constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
                               f(area_idx, model = "iid", constr = TRUE, 
                                       group = period_idx, control.group = list(model = "iid"),
                                       hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
                               f(iso3_idx, model = "iid", constr = TRUE, 
                                 group = period_idx, control.group = list(model = "iid"),
                                 hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))))

          # adding country-year interaction with the four demographic covariates
          formula3 <- Y ~ 1 + f(urban, constr = TRUE, 
                                group = period_idx, control.group = list(model = "iid"),
                                hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(age_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(edu_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(wealthq_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(area_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(iso3_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
            f(period_idx, model = "iid", constr = TRUE,
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))
          
          # formula 3 + country interaction on urban/rural   
          formula4 <- Y ~ 1 + f(urban, constr = TRUE, 
                                group = period_idx, control.group = list(model = "iid"),
                                replicate = iso3_idx,
                                hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(age_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(edu_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(wealthq_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(area_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))+
            f(iso3_idx, model = "iid", constr = TRUE, 
              group = period_idx, control.group = list(model = "iid"),
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
            f(period_idx, model = "iid", constr = TRUE,
              hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))
          
              # same as formula 4 but using replicate everywhere
          formula5 <- Y ~ 1 + f(period_idx,model = "iid", constr = TRUE,
                                hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
                                  f(iso3_idx, model = "iid", constr = TRUE, group = period_idx, control.group = list(model = "iid"),
                                hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
            f(urban, model = "iid", constr = TRUE, group = period_idx, control.group = list(model = "iid"),
                                            replicate = iso3_idx, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
                                  f(age_idx, model = "iid", constr = TRUE, group = period_idx, control.group = list(model = "iid"),
                                    replicate = iso3_idx, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
                                      f(edu_idx, model = "iid", constr = TRUE, group = period_idx, control.group = list(model = "iid"),
                                        replicate = iso3_idx, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
                                          f(wealthq_idx, model = "iid", constr = TRUE, group = period_idx, control.group = list(model = "iid"),
                                            replicate = iso3_idx, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01)))) +
            f(area_idx, model = "iid", constr = TRUE, group = period_idx, control.group = list(model = "iid"),
              replicate = iso3_idx, hyper = list(prec = list(prior = "pc.prec", param = c(2.5, 0.01))))
        
          # lumping all the formulas together
          
            formulas <- parse(text = paste0("list(", paste0("formula", 1:5, collapse = ", "), ")")) %>% eval()
            models <- paste0("Model ", 1:5) %>% as.list()
            
           
            library(INLA)
            library(tidyverse)

            #' Fit the models
            res <- purrr::pmap(
              list(formula = formulas, model_name = models),
              ~ logistic_model(..1, ..2)
            )

            res_df <- lapply(res, "[[", 1) %>% bind_rows()
            res_fit <- lapply(res, "[[", 2)

            #' save model objects
            saveRDS(res_fit, paste0(tmp_dir, '/',indicator,'-', tag, "-logit-sae-fits.rds"))
            saveRDS(res_df, paste0(tmp_dir,'/',indicator,'-', tag, "-logit-sae-df.rds"))

            # calculate goodness of fit metrics
            
            ic_df <- sapply(res_fit, function(fit) {
              local_dic <- na.omit(fit$dic$local.dic)
              local_waic <- na.omit(fit$waic$local.waic)
              local_cpo <- na.omit(fit$cpo$cpo)

              c("dic" = sum(local_dic),
                "dic_se" = stats::sd(local_dic) * sqrt(length(local_dic)),
                "waic" = sum(local_waic),
                "waic_se" = stats::sd(local_waic) * sqrt(length(local_waic)),
                "cpo" = sum(local_cpo),
                "cpo_se" = stats::sd(local_cpo) * sqrt(length(local_cpo)))
            }) %>%
              t() %>%
              round() %>%
              as.data.frame() %>%
              mutate(
                model = unlist(models),
                .before = dic
              )

            write_csv(ic_df, paste0(tmp_dir,'/',indicator, '-', tag,"-logit-information-criteria.csv"), na = "")

            #' Which model has the highest CPO?
            # keep model with covariates
            res_df_best <- filter(res_df, model == paste0(models[[which.max(ic_df$cpo)]]))
            res_fit_best <- res_fit[[which.max(ic_df$cpo)]]

            #' Artefacts: Best fitted model results and object
            write_csv(res_df_best, paste0(tmp_dir, '/',indicator, '-', tag,"-best-logit-sae.csv", na = ""))
            saveRDS(res_fit_best, paste0(tmp_dir, '/',indicator, '-', tag,"-best-logit-sae-fit.rds"))

            #' Artefact: Random effect variance parameter posterior means
            variance_df <- tryCatch(
              map(res_fit, function(fit)
                fit$marginals.hyperpar %>%
                  #' Calculate the expectation of the variance
                  map_df(function(x) inla.emarginal(fun = function(y) 1 / y, x)) %>%
                  #' Rename Precision to variance
                  rename_all(list(~ str_replace(., "Precision for ", "variance_")))
              ) %>%
                bind_rows() %>%
                #' Some of the models have other hyperparameters (e.g. rho)
                dplyr::select(starts_with("Variance")) %>%
                #' Sum of variance means
                mutate(total_variance = rowSums(., na.rm = TRUE)) %>%
                #' Create new columns with the percentage variance
                mutate(
                  across(
                    .cols = starts_with("Variance"),
                    .fns = list(percentage = ~ . / total_variance),
                    .names = "{fn}_{col}"
                  )
                ) %>%
                #' Add model identifier column
                mutate(
                  model = unlist(models),
                  .before = everything()
                ),
              error = function(e) {
                message("Error!")
                return(NULL)
              }
            )

            write_csv(variance_df, paste0(tmp_dir,'/',indicator, '-', tag,"variance-proportions.csv", na = ""))

      }
    }
  }
}
