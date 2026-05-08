# ===================================================================================================================================================
# CONTACT: Charlotte Le Marquand (charlotte.lemarquand@york.ac.uk)
# DATE: 08/05/2026

# SCRIPT: 05_visualisation_woodland_forest_plots.R
# DESCRIPTION: 
#   This script generates a specialized "Woodland-Only" version of the forest plots.
#   By filtering out control covariates (e.g., slope, density, camera effects), it 
#   provides a focused visualization of the primary research hypothesis: the 
#   relationship between woodland metrics and nest predation risk.
#
# ===================================================================================================================================================

# 1. LIBRARIES -----------------------------------------------------------------------------------------------------
library(forestplot) # Visualization of estimates and credible intervals
library(dplyr)      # Data manipulation
library(grid)       # Graphical formatting (gpar)
library(stringr)    # Label cleaning
library(readr)      # Loading control table

# 2. RUN CONTROL & ITERATION ---------------------------------------------------------------------------------------
control_table <- read_csv("model_runs.csv")

for(i in 1:nrow(control_table)) {
  
  current_run <- control_table[i, ]
  run_id  <- current_run$run_id
  out_dir <- current_run$output_dir
  species_name <- current_run$species
  
  message("--------------------------------------------------")
  message(paste0("GENERATING TARGETED WOODLAND PLOT: ", run_id))
  message("--------------------------------------------------")
  
  # 3. DATA LOADING ------------------------------------------------------------------------------------------------
  handover_path <- file.path(out_dir, "table_handover.rds")
  
  if(!file.exists(handover_path)) {
    warning(paste("Handover file missing for", run_id, "- Skipping."))
    next
  }
  
  handover <- readRDS(handover_path)
  final_summary_df <- handover$final_summary_df
  
  # 4. TARGETED FILTERING ------------------------------------------------------------------------------------------
  # Define the specific woodland metrics of interest
  woodland_covariates <- c(
    "Distance to All Woodland",
    "Distance to Broadleaved", 
    "Distance to Coniferous",
    "Area of Woodland",           
    "Woodland Perimeter Length", 
    "Distance to Woody Feature"
  )
  
  # Isolate only the woodland variables for this focused plot
  woodland_df <- final_summary_df %>%
    filter(covariate %in% woodland_covariates)
  
  if(nrow(woodland_df) == 0) {
    warning(paste("No woodland-specific data found for", run_id))
    next
  }
  
  # 5. TIDYING LABELS & ORDERING -----------------------------------------------------------------------------------
  woodland_df <- woodland_df %>%
    mutate(
      variable_display = case_when(
        covariate == "Distance to All Woodland" ~ "Distance to All Wood (km)",
        covariate == "Distance to Broadleaved" ~ "Distance to Broadleaved Wood (km)",
        covariate == "Distance to Coniferous" ~ "Distance to Coniferous Wood (km)",
        covariate == "Area of Woodland" ~ "Proportion of All Wood within 1km", 
        covariate == "Woodland Perimeter Length" ~ "Length of Woodland Perimeter within 1km",
        covariate == "Distance to Woody Feature" ~ "Distance to Woody Linear Feature (km)",
        TRUE ~ Term
      )
    )
  
  # Establish the master ordering for the forest plot y-axis
  variable_order <- c(
    "Distance to All Wood (km)",
    "Distance to Broadleaved Wood (km)",
    "Distance to Coniferous Wood (km)",
    "Proportion of All Wood within 1km", 
    "Length of Woodland Perimeter within 1km", 
    "Distance to Woody Linear Feature (km)"
  )
  
  woodland_df$variable_display <- factor(woodland_df$variable_display, 
                                         levels = variable_order, 
                                         ordered = TRUE)
  
  # 6. PLOT TABLE CONSTRUCTION -------------------------------------------------------------------------------------
  # Creating the specific data structure required by the 'forestplot' package
  processed_data <- woodland_df %>%
    filter(!is.na(variable_display)) %>%
    mutate(dic_val = as.character(dic_difference)) %>%
    group_by(Effect_Type) %>%
    do({
      group_header <- data.frame(
        variable_display = unique(.$Effect_Type), 
        median = NA_real_, Lower_95_CI = NA_real_, Upper_95_CI = NA_real_, 
        Effect_Type = "HEADER", dic_val = NA_character_,
        Original_Effect_Type = unique(.$Effect_Type) 
      )
      bind_rows(group_header, mutate(., Original_Effect_Type = Effect_Type))
    }) %>%
    ungroup()
  
  plot_data <- bind_rows(
    data.frame(variable_display = "Woodland Variable", median = NA_real_, Lower_95_CI = NA_real_, 
               Upper_95_CI = NA_real_, Effect_Type = "TITLE", dic_val = "DIC Diff",
               Original_Effect_Type = "TITLE"),
    processed_data
  ) %>%
    mutate(
      group_rank = match(Original_Effect_Type, c("TITLE", "Fixed")),
      type_rank = ifelse(Effect_Type %in% c("TITLE", "HEADER"), 1, 2),
      var_rank = match(variable_display, c("Woodland Variable", variable_order))
    ) %>%
    arrange(group_rank, type_rank, var_rank) %>%
    mutate(
      DIC_Label = case_when(
        is.na(dic_val) ~ "",
        dic_val == "DIC Diff" ~ "DIC Diff",
        suppressWarnings(!is.na(as.numeric(dic_val))) ~ sprintf("%.2f", as.numeric(dic_val)),
        TRUE ~ dic_val
      ),
      Term = as.character(variable_display) 
    )
  
  is_header <- plot_data$Effect_Type %in% c("TITLE", "HEADER")
  total_rows <- NROW(plot_data)
  
  # Aesthetic parameters: Horizontal lines for section separation
  fixed_top_lines <- list("1" = gpar(lwd = 1.5, col="black"), "2" = gpar(lwd = 1.5, col="black"))
  header_rows <- which(plot_data$Effect_Type == "HEADER")
  group_separator_lines <- setNames(
    lapply(seq_along(header_rows), function(i) gpar(lwd = 0.8, col = "grey60", lty = 2)),
    as.character(header_rows)
  )
  bottom_line <- setNames(list(gpar(lwd = 1.5, col="black")), as.character(total_rows + 1))
  all_hrzl_lines <- c(fixed_top_lines, group_separator_lines, bottom_line)
  
  # 7. GENERATE & SAVE TARGETED PLOT -------------------------------------------------------------------------------
  
  generate_wood_forest_plot <- function() {
    plot_data %>%
      forestplot(
        labeltext = c(Term, DIC_Label),
        mean = median,          
        lower = Lower_95_CI,
        upper = Upper_95_CI,
        is.summary = is_header,
        zero = 0,
        boxsize = 0.2,
        xlab = "Posterior Coefficient Estimate (Linear Predictor)",
        col = fpColors(box = "black", line = "#444444", zero = "#CCCCCC", summary = "black"),
        lineheight = "auto",
        hrzl_lines = all_hrzl_lines,
        align = c("l", "c"),
        txt_gp = fpTxtGp(
          label = gpar(cex = 0.9),
          ticks = gpar(cex = 0.8),
          xlab  = gpar(cex = 1.2),
          title = gpar(cex = 1.5)
        )
      )
  }
  
  fig_dir <- file.path(out_dir, "figures")
  file_base <- file.path(fig_dir, paste0(run_id, "_Woodland_Only_Plot"))
  
  # Export focused plot in PDF and PNG
  plot_height <- max(6, nrow(plot_data) * 0.4) 
  
  pdf(paste0(file_base, ".pdf"), width = 12, height = plot_height)
  print(generate_wood_forest_plot())
  dev.off()
  
  png(paste0(file_base, ".png"), width = 1400, height = plot_height * 100, res = 150)
  print(generate_wood_forest_plot())
  dev.off()
  
  message(paste("Success: Targeted woodland plot saved for", run_id))
}

print("Workflow complete: Specialized woodland analysis figures exported.")