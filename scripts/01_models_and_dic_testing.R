# ===================================================================================================================================================
# CONTACT: Charlotte Le Marquand (charlotte.lemarquand@york.ac.uk)
# DATE: 08/05/2026

# SCRIPT: 01_models_and_dic_testing.R
# DESCRIPTION: 
#   This script performs Bayesian spatially-explicit Generalized Linear Mixed Models (GLMMs) 
#   to analyze wader nest predation. The workflow includes:
#     1. Spatial data preparation and projection (British National Grid).
#     2. Construction of a stochastic partial differential equation (SPDE) mesh for spatial autocorrelation.
#     3. Model fitting via Integrated Nested Laplace Approximation (INLA) using the 'inlabru' package.
#     4. Generation of predictive plots and spatial field visualizations.
#     5. Leave-one-out DIC comparison to evaluate variable importance.
#
# DATA NOTE: 
#   The input data has been jittered and attributes scrambled to protect sensitive nest locations.
#   Results are for reproducibility of the workflow and will vary from the manuscript results.
# ===================================================================================================================================================

# 1. LOAD LIBRARIES ---------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)      # Data wrangling and visualization
library(inlabru)        # Bayesian spatial modeling wrapper for INLA
library(sf)             # Spatial data handling (Simple Features)
library(rnaturalearth)  # Base map data for the UK
library(sp)             # Spatial objects (required for certain INLA functions)
library(INLA)           # Core Bayesian inference engine

# 2. RUN CONTROL & ITERATION ------------------------------------------------------------------------------------------------------------------------

# A control table is used to automate runs for different species and covariate combinations
control_table <- read_csv("model_runs.csv")

for (i in 1:nrow(control_table)) {
  
  current_run <- control_table[i, ]
  run_id <- current_run$run_id
  out_dir <- current_run$output_dir
  species <- current_run$species
  
  message("--------------------------------------------------")
  message(paste0("STARTING ANALYSIS: ", run_id, " (Species: ", species, ")"))
  message(paste0("Run ", i, " of ", nrow(control_table)))
  message("--------------------------------------------------")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 3. DATA PREPARATION -----------------------------------------------------------------------------------------------------------------------------
  
  df <- readRDS(current_run$data_in)
  
  # Create a grouped random effect term for Year nested within Project
  df <- df %>%
    mutate(project_year = as.factor(paste(project, year, sep = "_")))
  
  # Filter out specific groupings where zero variance in outcome prevents model convergence
  to_drop <- strsplit(current_run$levels_to_drop, ", ")[[1]]
  if (to_drop[1] != "None") {
    df <- df %>% filter(!(project_year %in% to_drop)) 
  }
  df <- df %>% droplevels() 
  
  # Define categorical variables and set reference level for habitat comparisons.
  # improved_grassland is set because this is the largest category by far.
  df$camera <- as.factor(df$camera)
  df$habitat <- as.factor(df$habitat)
  df$habitat <- relevel(df$habitat, ref = "improved_grassland")
  
  # LOGISTIC-EXPOSURE ADJUSTMENT:
  # Adjust exposure days to account for the timing of the predation event.
  df$mod_exposure_days <- (df$exposure_days - df$binary)
  model_data <- df
  
  # Apply proximity thresholds for woodland distance - if further than 10km away, exclude - too far to accurately include.
  if("dist_all_wood" %in% names(model_data)) model_data <- model_data %>% filter(dist_all_wood < 10.0001)
  if("dist_broad_wood" %in% names(model_data)) model_data <- model_data %>% filter(dist_broad_wood < 10.0001)
  if("dist_conifer_wood" %in% names(model_data)) model_data <- model_data %>% filter(dist_conifer_wood < 10.0001)
  
  saveRDS(model_data, file.path(out_dir, "processed_model_data.rds"))
  
  # 4. SPATIAL DOMAIN & MESH CONSTRUCTION -----------------------------------------------------------------------------------------------------------
  
  # Convert to Spatial Features and project to British National Grid (meters)
  points_sf <- st_as_sf(model_data, coords = c("lon", "lat"), crs = 4326)
  points_bng <- st_transform(points_sf, crs = 27700)
  
  # Define study domain: 5km buffer for high-density mesh, 100km buffer for boundary effects
  buffers <- st_buffer(points_bng, dist = 5000)
  single_domain <- st_union(buffers)
  
  uk_map_sf <- ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")
  uk_map_bng <- st_transform(uk_map_sf, crs = 27700)
  outer_buffer <- st_buffer(points_bng, dist = 100000)
  outer_domain <- st_intersection(st_union(outer_buffer), uk_map_bng)
  
  # SCALE CONVERSION FOR NUMERICAL STABILITY:
  # Coordinates are divided by 10,000 to convert units to 10km, improving INLA convergence.
  scale_factor <- 10000
  
  points_scaled <- points_bng 
  st_geometry(points_scaled) <- st_geometry(points_scaled) / scale_factor 
  st_crs(points_scaled) <- NA 
  
  single_domain_scaled <- single_domain / scale_factor 
  st_crs(single_domain_scaled) <- NA 
  
  outer_domain_scaled <- outer_domain / scale_factor 
  st_crs(outer_domain_scaled) <- NA 
  
  # Construct SPDE mesh with dual-density triangulation
  max_edge_original <- c(5000, 40000) # 5km and 40km in original meters
  max_edge_scaled <- max_edge_original / scale_factor
  
  mesh <- fmesher::fm_mesh_2d_inla(
    loc = st_coordinates(points_scaled), 
    max.edge = max_edge_scaled, 
    boundary = list(sf::as_Spatial(single_domain_scaled), 
                    sf::as_Spatial(outer_domain_scaled))
  )
  
  saveRDS(mesh, file.path(out_dir, "spatial_mesh.rds"))
  
  # 5. MODEL SPECIFICATION & FITTING ----------------------------------------------------------------------------------------------------------------
  
  bru_data_df <- as.data.frame(model_data)
  coords_bng <- st_coordinates(st_transform(st_as_sf(bru_data_df, coords = c("lon", "lat"), crs = 4326), crs = 27700))
  
  # Store scaled coordinates in the data frame
  bru_data_df$x <- coords_bng[,1] / scale_factor   
  bru_data_df$y <- coords_bng[,2] / scale_factor
  
  # Identify covariates for standardization
  wood_vars_string <- current_run$woodland_covariates
  indiv_wood_vars <- unlist(strsplit(wood_vars_string, " \\+ "))
  numerical_predictors <- c(indiv_wood_vars, "area_1000_all_wood", "edge_1000_all_wood", 
                            "linear_feat_dist", "slope", "wader_density", "corvid_density")
  
  # Standardize continuous variables (mean = 0, sd = 1) for better model estimation
  bru_data_df_scaled <- bru_data_df
  centers <- list(); scales <- list()
  for (var in numerical_predictors) {
    temp_scaled <- scale(bru_data_df[[var]])
    bru_data_df_scaled[[var]] <- as.vector(temp_scaled)
    centers[[var]] <- attr(temp_scaled, "scaled:center")
    scales[[var]] <- attr(temp_scaled, "scaled:scale")
  }
  
  # Prepare factor indices
  bru_data_df_scaled$project_year_id <- as.numeric(as.factor(bru_data_df_scaled$project_year))
  bru_data_df_scaled$habitat_id <- as.numeric(bru_data_df_scaled$habitat)
  bru_data_df_scaled$camera_id <- as.numeric(bru_data_df_scaled$camera)
  
  # Define Spatial and Random Effect Priors
  spde <- inla.spde2.pcmatern(mesh, prior.range = c(1, 0.2), prior.sigma = c(1, 0.2))
  pc_prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
  
  # Dynamic Formula Construction
  comp_formula_string <- paste0(
    "~ Intercept(1) + ", wood_vars_string, " + area_1000_all_wood + edge_1000_all_wood + ",
    "linear_feat_dist + slope + wader_density + corvid_density + ",
    "project_year_ef(project_year_id, model = 'iid', hyper = pc_prec_prior) + ",
    "habitat_ef(habitat_id, model = 'factor_contrast') + ",
    "camera_ef(camera_id, model = 'factor_contrast') + ",
    "field(main = cbind(x, y), model = spde)"
  )
  
  set.seed(42)
  fit <- bru(
    components = as.formula(comp_formula_string),
    likelihoods = bru_obs(
      formula = binary ~ .,
      family = "binomial",
      data = bru_data_df_scaled,
      Ntrials = bru_data_df_scaled$mod_exposure_days + 1
    ),
    options = list(control.compute = list(dic = TRUE))
  )
  
  saveRDS(fit, file.path(out_dir, "fitted_inla_model.rds"))
  
  # 6. POST-HOC VISUALIZATION -----------------------------------------------------------------------------------------------------------------------
  
  # Function to visualize marginalized fixed effects back-transformed to probability scale
  plot_continuous_effect <- function(variable_name) {
    
    # Create prediction frame across the observed range of the covariate
    x_range <- seq(min(df[[variable_name]], na.rm = TRUE), max(df[[variable_name]], na.rm = TRUE), length.out = 100)
    pred_df <- data.frame(
      x = mean(bru_data_df$x), y = mean(bru_data_df$y),
      habitat_id = 1, camera_id = 1, project_year_id = 1
    )
    pred_df[[variable_name]] <- x_range
    
    # Scale prediction data using center/scale from original fitting
    pred_df_scaled <- pred_df
    for (var in numerical_predictors) {
      if (var %in% names(pred_df)) {
        pred_df_scaled[[var]] <- (pred_df[[var]] - centers[[var]]) / scales[[var]]
      }
    }  
    
    # Predict daily predation probability
    predictions <- predict(fit, pred_df_scaled, formula = as.formula(paste0("~ Intercept + ", variable_name)))
    
    daily_mean_prob <- exp(predictions$mean) / (1 + exp(predictions$mean))
    daily_low_prob <- exp(as.numeric(predictions$q0.025)) / (1 + exp(as.numeric(predictions$q0.025)))
    daily_high_prob <- exp(as.numeric(predictions$q0.975)) / (1 + exp(as.numeric(predictions$q0.975)))
    
    # Transform daily probability to cumulative predation risk over the full incubation period
    inc_days <- case_when(species == "curlew" ~ 29, species == "lapwing" ~ 34, 
                          species == "oystercatcher" ~ 27, species == "redshank" ~ 24, TRUE ~ 0)
    
    predictions$mean_prob <- 1 - (1 - daily_mean_prob)^inc_days
    predictions$low_prob  <- 1 - (1 - daily_low_prob)^inc_days
    predictions$high_prob <- 1 - (1 - daily_high_prob)^inc_days
    predictions[[variable_name]] <- x_range
    
    ggplot(predictions) +
      aes_string(x = variable_name, y = "mean_prob") +
      geom_ribbon(aes(ymin = low_prob, ymax = high_prob), alpha = 0.2) +
      geom_line() +
      labs(title = paste("Effect of", variable_name), x = variable_name, y = "Predation Probability") +
      theme_classic()
  }
  
  # Export plots as high-resolution TIFFs for publication
  continuous_vars <- c(indiv_wood_vars, "area_1000_all_wood", "edge_1000_all_wood", "linear_feat_dist")
  for (var in continuous_vars) {
    p <- plot_continuous_effect(var)
    ggsave(file.path(out_dir, paste0("effect_", var, ".tiff")), plot = p, device = "tiff", dpi = 600)
  }
  
  # 7. MODEL SELECTION (DIC COMPARISONS) ------------------------------------------------------------------------------------------------------------
  # Iteratively drop each covariate to evaluate its contribution to the Deviance Information Criterion (DIC)
  
  vars_to_drop <- c(indiv_wood_vars, numerical_predictors[-c(1:length(indiv_wood_vars))], 
                    "project_year_ef", "habitat_ef", "camera_ef")
  
  dic_results <- data.frame(model = "full", dropped = "None", dic = fit$dic$dic)
  
  for (v in vars_to_drop) {
    message(paste0("Testing variable importance: ", v))
    
    # Re-build formula excluding current variable
    rem_fixed <- numerical_predictors[numerical_predictors != v]
    py_ef <- if (v == "project_year_ef") "" else "project_year_ef(project_year_id, model = 'iid', hyper = pc_prec_prior) + "
    hab_ef <- if (v == "habitat_ef") "" else "habitat_ef(habitat_id, model = 'factor_contrast') + "
    cam_ef <- if (v == "camera_ef") "" else "camera_ef(camera_id, model = 'factor_contrast') + "
    
    sub_formula <- as.formula(paste0("~ Intercept(1) + ", paste(rem_fixed, collapse = " + "), " + ",
                                     py_ef, hab_ef, cam_ef, "field(main = cbind(x, y), model = spde)"))
    
    set.seed(42)
    fit_sub <- bru(components = sub_formula, likelihoods = bru_obs(formula = binary ~ ., family = "binomial", 
                                                                   data = bru_data_df_scaled, Ntrials = bru_data_df_scaled$mod_exposure_days + 1),
                   options = list(control.compute = list(dic = TRUE)))
    
    dic_results <- rbind(dic_results, data.frame(model = paste0("minus_", v), dropped = v, dic = fit_sub$dic$dic))
    gc()
  }
  
  write_csv(dic_results, file.path(out_dir, "model_selection_dic.csv"))
  message(paste("Run ID", run_id, "completed successfully."))
} 

print("Analysis workflow complete.")