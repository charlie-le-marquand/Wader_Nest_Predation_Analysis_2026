# ===================================================================================================================================================
# CONTACT: Charlotte Le Marquand (charlotte.lemarquand@york.ac.uk)
# DATE: 08/05/2026
#
# SCRIPT: 02_visualisation_tables.R
# DESCRIPTION: 
#   This script automates the extraction of posterior distributions and model selection 
#   metrics from the fitted INLA models. It produces standardized tables for 
#   manuscript reporting, including:
#     1. Posterior means, medians, and 95% Credible Intervals (CI) for fixed and random effects.
#     2. Mapping of technical variable names to human-readable manuscript labels.
#     3. Integration of DIC (Deviance Information Criterion) values to show variable importance.
#
# ===================================================================================================================================================

# 1. LOAD LIBRARIES ---------------------------------------------------------------------------------------------------------------------------------
library(tidyverse) # Core data manipulation

# 2. RUN CONTROL ------------------------------------------------------------------------------------------------------------------------------------
# Load the master control table to iterate through all species/model runs
control_table <- read_csv("model_runs.csv")

for(i in 1:nrow(control_table)) {
  
  current_run <- control_table[i, ]
  run_id  <- current_run$run_id
  out_dir <- current_run$output_dir
  
  # Detect woodland variables specific to this model run for dynamic labeling
  wood_vars_string <- current_run$woodland_covariates
  indiv_wood_vars <- unlist(strsplit(wood_vars_string, " \\+ "))
  
  message("--------------------------------------------------")
  message(paste0("EXTRACTING RESULTS FOR: ", run_id))
  message("--------------------------------------------------")
  
  # 3. LOAD MODEL OUTPUTS ---------------------------------------------------------------------------------------------------------------------------
  model_file <- readRDS(file.path(out_dir, "fit_model.rds"))
  dic_values <- read_csv(file.path(out_dir, "dic_comparisons.csv"))
  df         <- readRDS(file.path(out_dir, "script1_data.rds"))
  
  # 4. PARAMETER EXTRACTION -------------------------------------------------------------------------------------------------------------------------
  
  # Extract Posterior Summaries for Fixed Effects
  fixed_effects_df <- model_file$summary.fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    mutate(Effect_Type = "Fixed") %>%
    rename(Mean = mean, SD = sd, Lower_95_CI = '0.025quant', Upper_95_CI = '0.975quant')
  
  # Extract factor levels from the original data to map back to the model indices
  habitat_levels <- levels(as.factor(df$habitat))
  camera_levels <- levels(as.factor(df$camera))
  project_year_levels <- levels(as.factor(df$project_year))
  
  # Extract Factor Contrasts (Habitat and Camera effects)
  habitat_df <- model_file$summary.random$habitat_ef %>% as.data.frame() %>%
    mutate(Term = habitat_levels[as.integer(ID)], Effect_Type = "Factor Contrast - Habitat") %>%
    rename(Mean = mean, SD = sd, Lower_95_CI = '0.025quant', Upper_95_CI = '0.975quant')
  
  camera_df <- model_file$summary.random$camera_ef %>% as.data.frame() %>%
    mutate(Term = camera_levels[as.integer(ID)], Effect_Type = "Factor Contrast - Camera") %>%
    rename(Mean = mean, SD = sd, Lower_95_CI = '0.025quant', Upper_95_CI = '0.975quant')
  
  # Extract Grouped Random Effects (Project-Year)
  project_year_df <- model_file$summary.random$project_year_ef %>% as.data.frame() %>%
    mutate(Term = project_year_levels[as.integer(ID)], Effect_Type = "Random Effect - Project-Year") %>%
    rename(Mean = mean, SD = sd, Lower_95_CI = '0.025quant', Upper_95_CI = '0.975quant') %>%
    select(-ID)
  
  # 5. CONSOLIDATE RESULTS --------------------------------------------------------------------------------------------------------------------------
  
  final_summary_df <- bind_rows(fixed_effects_df, habitat_df, camera_df, project_year_df) %>%
    rename(median = "0.5quant") %>%
    mutate(
      Lower_95_CI = round(Lower_95_CI, 5), 
      Upper_95_CI = round(Upper_95_CI, 5),
      median = round(median, 5), 
      Term = str_trim(Term)
    )
  
  # 6. MANUSCRIPT LABEL MAPPING ---------------------------------------------------------------------------------------------------------------------
  # Map technical internal variable names to descriptive labels for publication tables.
  
  final_summary_df <- final_summary_df %>%
    mutate(
      covariate = case_when(
        Term %in% habitat_levels ~ "Habitat",
        grepl("_20[0-9]{2}", Term) ~ "Project-Year",
        Term == "Intercept" ~ "Intercept",
        Term == "dist_all_wood" ~ "Distance to All Woodland",
        Term == "dist_broad_wood" ~ "Distance to Broadleaved",
        Term == "dist_conifer_wood" ~ "Distance to Coniferous",
        Term == "area_1000_all_wood" ~ "Area of Woodland",
        Term == "edge_1000_all_wood" ~ "Woodland Perimeter Length",
        Term == "linear_feat_dist" ~ "Distance to Woody Feature",
        Term == "slope" ~ "Slope",
        Term == "wader_density" ~ "Wader Density",
        Term == "corvid_density" ~ "Corvid Density",
        Term == "Y" ~ "Camera",
        TRUE ~ "UNCLASSIFIED"
      )
    )
  
  # Apply identical mapping to the DIC comparison results for seamless joining
  dic_values <- dic_values %>%
    mutate(
      covariate = case_when(
        dropped_covariate == "None" ~ "Intercept",
        dropped_covariate == "dist_all_wood" ~ "Distance to All Woodland",
        dropped_covariate == "dist_broad_wood" ~ "Distance to Broadleaved",
        dropped_covariate == "dist_conifer_wood" ~ "Distance to Coniferous",
        dropped_covariate == "area_1000_all_wood" ~ "Area of Woodland",
        dropped_covariate == "edge_1000_all_wood" ~ "Woodland Perimeter Length",
        dropped_covariate == "linear_feat_dist" ~ "Distance to Woody Feature",
        dropped_covariate == "slope" ~ "Slope",
        dropped_covariate == "wader_density" ~ "Wader Density",
        dropped_covariate == "corvid_density" ~ "Corvid Density",
        dropped_covariate == "project_year_ef" ~ "Project-Year",
        dropped_covariate == "habitat_ef" ~ "Habitat",
        dropped_covariate == "camera_ef" ~ "Camera",
        TRUE ~ dropped_covariate
      )
    )
  
  # 7. FINAL INTEGRATION & EXPORT -------------------------------------------------------------------------------------------------------------------
  
  # Join posterior estimates with variable importance (Delta DIC)
  final_summary_df <- final_summary_df %>%
    left_join(dic_values, by = "covariate") %>% 
    mutate(
      dic = round(dic, 2),
      dic_difference = round(dic_difference, 2),
      variable_label = paste(covariate, Term, sep = " - ")
    )
  
  # Select and format columns for the final manuscript-ready CSV
  details_for_table <- final_summary_df %>% 
    mutate(variable = factor(variable_label, levels = unique(variable_label))) %>%
    dplyr::select(variable, median, Lower_95_CI, Upper_95_CI, dic_difference)
  
  # Export the consolidated summary
  write.csv(details_for_table, file = file.path(out_dir, paste0(run_id, "_model_summary_full.csv")), row.names = FALSE)
  
  # Save handover object for subsequent plotting scripts (05 & 06)
  saveRDS(list(
    final_summary_df = final_summary_df,
    habitat_levels = habitat_levels,
    camera_levels = camera_levels,
    species = current_run$species 
  ), file.path(out_dir, "table_handover.rds"))
  
  message(paste0("Summary table generated for: ", run_id))
}

print("Workflow complete: All model summaries exported.")