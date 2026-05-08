# ===================================================================================================================================================
# CONTACT: Charlotte Le Marquand (charlotte.lemarquand@york.ac.uk)
# DATE: 08/05/2026

# SCRIPT: 02_visualisation_prediction_plots.R
# DESCRIPTION: 
#   This script generates publication-ready visualizations from the fitted INLA models.
#   It processes the marginalized effects of woodland covariates on nest predation 
#   probability and creates multi-panel figures. Key features include:
#     1. Back-transformed probabilities (daily risk to cumulative incubation risk).
#     2. 'Rug plots' showing the distribution of raw data (quantiles).
#     3. Biological thresholds (e.g., 1km proximity markers).
#     4. High-resolution multi-panel assembly for manuscript submission.

# DATA NOTE: 
#   The input data has been jittered and attributes scrambled to protect sensitive nest locations.
#   Results are for reproducibility of the workflow and will vary from the manuscript results.

# ===================================================================================================================================================

# 1. LOAD LIBRARIES ---------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)      # Data manipulation and grammar of graphics
library(patchwork)      # Tools for combining separate ggplots into panels
library(stats)          # Statistical functions for quantile calculation

# 2. RUN CONTROL & ITERATION ------------------------------------------------------------------------------------------------------------------------

control_table <- read_csv("model_runs.csv")

for(i in 1:nrow(control_table)) {
  
  current_run <- control_table[i, ]
  run_id  <- current_run$run_id
  out_dir <- current_run$output_dir
  species <- current_run$species
  
  message("--------------------------------------------------")
  message(paste0("GENERATING PANELS: ", run_id, " (", species, ")"))
  message("--------------------------------------------------")
  
  # Load processed observation data and the base predictive plots from Script 01
  df <- readRDS(file.path(out_dir, "processed_model_data.rds"))
  all_plots <- readRDS(file.path(out_dir, "continuous_plots.rds"))
  
  # 3. DATA DENSITY PROCESSING (RUG PLOTS) ----------------------------------------------------------------------------------------------------------
  
  # To visualize the distribution of our observations, we calculate 101 quantiles 
  # for each covariate, split by the binary outcome (Predated vs. Hatched).
  vars_to_process <- names(all_plots)
  quantile_all <- list()
  
  for(v in vars_to_process) {
    # Distribution for predated nests (Binary = 1)
    quantile_all[[paste0(v, "_TRUE")]] <- quantile(df[[v]][df$binary == 1], seq(0, 1, length = 101), na.rm = TRUE)
    # Distribution for hatched/survived nests (Binary = 0)
    quantile_all[[paste0(v, "_FALSE")]] <- quantile(df[[v]][df$binary == 0], seq(0, 1, length = 101), na.rm = TRUE)
  }
  quantile_all <- as.data.frame(quantile_all)
  
  # 4. PLOT REFINEMENT & AESTHETICS -----------------------------------------------------------------------------------------------------------------
  
  edited_plots_list <- list()
  
  for(v in vars_to_process) {
    
    # Retrieve the base marginal effect plot (Line + 95% Credible Interval)
    p_base <- all_plots[[v]]
    
    # Apply manuscript-standard theme and add visual indicators
    p_edited <- p_base +
      # Bottom Rug: Distribution of predated nests
      geom_rug(data = quantile_all, aes(x = .data[[paste0(v, "_TRUE")]]), 
               inherit.aes = FALSE, sides = "b", color = "black", alpha = 0.4) + 
      # Top Rug: Distribution of surviving nests
      geom_rug(data = quantile_all, aes(x = .data[[paste0(v, "_FALSE")]]), 
               inherit.aes = FALSE, sides = "t", color = "black", alpha = 0.4) + 
      # Reference line for 50% predation probability
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "blue", alpha = 0.5) + 
      theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 18)
      )
    
    # CUSTOM LABELS AND BIOLOGICAL THRESHOLDS:
    # Adding 1km 'edge effect' threshold for distance-based covariates
    
    if(v == "dist_all_wood") {
      p_edited <- p_edited + geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
        labs(title = NULL, x = "Distance to All Woodland (km)", y = "Cumulative Predation Probability")
    }
    
    if(v == "dist_broad_wood") {
      p_edited <- p_edited + geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
        labs(title = NULL, x = "Distance to Broadleaved Woodland (km)", y = "Cumulative Predation Probability")
    }
    
    if(v == "dist_conifer_wood") {
      p_edited <- p_edited + geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
        labs(title = NULL, x = "Distance to Coniferous Woodland (km)", y = "Cumulative Predation Probability")
    }
    
    if(v == "area_1000_all_wood") {
      p_edited <- p_edited + labs(title = NULL, x = "Proportion of Woodland (1km buffer)", y = "Cumulative Predation Probability")
    }
    
    if(v == "edge_1000_all_wood") {
      p_edited <- p_edited + labs(title = NULL, x = "Woodland Perimeter (1km buffer)", y = "Cumulative Predation Probability")
    }
    
    if(v == "linear_feat_dist") {
      p_edited <- p_edited + geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
        labs(title = NULL, x = "Distance to Woody Linear Feature (km)", y = "Cumulative Predation Probability")
    }
    
    edited_plots_list[[v]] <- p_edited
    saveRDS(p_edited, file.path(out_dir, paste0(v, "_final_plot.rds")))
  }
  
  # 5. ASSEMBLE MULTI-PANEL FIGURES -----------------------------------------------------------------------------------------------------------------
  
  # Combine plots using patchwork and add alphabetical tagging (A, B, C...)
  all_plots_panel <- wrap_plots(edited_plots_list, ncol = 2) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(size = 20, face = "bold", family = "serif"),
      plot.tag.position = c(0, 1), 
      plot.tag.margin = margin(t = 2, r = 2, b = 0, l = 0)
    )
  
  # Export final figure at high resolution (600 DPI) for publication
  if(!dir.exists(file.path(out_dir, "figures"))) dir.create(file.path(out_dir, "figures"), recursive = TRUE)
  
  ggsave(file.path(out_dir, "figures", paste0(run_id, "_main_figure.jpg")), 
         plot = all_plots_panel, width = 15, height = 20, dpi = 600)
  
  message(paste("Successfully saved figure panel for:", run_id))
}

print("Workflow complete: All figures generated.")