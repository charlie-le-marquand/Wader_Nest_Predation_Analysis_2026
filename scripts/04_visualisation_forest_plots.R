# ===================================================================================================================================================
# CONTACT: Charlotte Le Marquand (charlotte.lemarquand@york.ac.uk)
# DATE: 08/05/2026 

# SCRIPT: 04_visualisation_forest_plots.R
# DESCRIPTION: 
#   This script generates high-resolution Forest Plots for the manuscript.
#   It visualizes the posterior estimates (coefficient means) and 95% Credible Intervals 
#   for all model parameters. Key features:
#     1. Automatic grouping of variables (Fixed Effects, Habitat, Camera, Random Effects).
#     2. Integration of variable importance metrics (Delta DIC).
#     3. Dynamic plot scaling to accommodate different numbers of factor levels per species.
#     4. Publication-quality formatting using the 'forestplot' package.
#
# ===================================================================================================================================================

# 1. LIBRARIES -----------------------------------------------------------------------------------------------------
library(forestplot) # Specialized tool for forest plot visualization
library(dplyr)      # Data manipulation
library(grid)       # Graphical parameter control (gpar)
library(stringr)    # String cleaning and formatting
library(readr)      # Fast loading of CSV control tables

# 2. RUN CONTROL & ITERATION ---------------------------------------------------------------------------------------
control_table <- read_csv("model_runs.csv")

for(i in 1:nrow(control_table)) {
  
  current_run <- control_table[i, ]
  run_id  <- current_run$run_id
  out_dir <- current_run$output_dir
  species_name <- current_run$species
  
  message("--------------------------------------------------")
  message(paste0("GENERATING FOREST PLOT: ", run_id))
  message("--------------------------------------------------")
  
  # 3. DATA PREPARATION --------------------------------------------------------------------------------------------
  
  # Load the handover object from Script 04 containing summarized model results
  handover_path <- file.path(out_dir, "table_handover.rds")
  
  if(!file.exists(handover_path)) {
    warning(paste("Handover file missing for", run_id, "- Ensure Script 04 ran successfully."))
    next
  }
  
  handover <- readRDS(handover_path)
  final_summary_df <- handover$final_summary_df
  habitat_levels   <- handover$habitat_levels
  
  # Ensure DIC values are present for labels; fill with NA if missing
  if(!"dic_difference" %in% colnames(final_summary_df)){
    final_summary_df$dic_difference <- NA
  }
  
  # 4. DESCRIPTIVE LABELING -----------------------------------------------------------------------------------------
  # Map technical model terms to human-readable labels for the plot y-axis.
  
  final_summary_df <- final_summary_df %>%
    mutate(
      variable_display = case_when(
        # Continuous Predictors
        Term == "dist_all_wood" ~ "Distance to All Wood (km)",
        Term == "dist_broad_wood" ~ "Distance to Broadleaved Wood (km)",
        Term == "dist_conifer_wood" ~ "Distance to Coniferous Wood (km)",
        Term == "area_1000_all_wood" ~ "Proportion of All Wood within 1km",
        Term == "edge_1000_all_wood" ~ "Length of Woodland Perimeter within 1km",
        Term == "linear_feat_dist" ~ "Distance to Woody Linear Feature (km)",
        Term == "slope" ~ "Slope (Degrees)",
        Term == "wader_density" ~ "Wader Density",
        Term == "corvid_density" ~ "Corvid Density",
        
        # Categorical Effects (Formatted for title case)
        Term %in% habitat_levels ~ paste0("Habitat: ", str_to_title(str_replace_all(Term, "_", " "))),
        Term == "Y" ~ "Camera - Yes",
        
        # Random Effects (Project-Year groupings)
        grepl("_20[0-9]{2}", Term) ~ paste0("Project-Year: ", str_to_title(str_replace_all(Term, "_", " "))),
        
        TRUE ~ Term
      )
    )
  
  # 5. ORDERING & GROUPING LOGIC -----------------------------------------------------------------------------------
  # We define a strict order to ensure consistency across species-specific plots.
  
  fixed_order <- c(
    "Distance to All Wood (km)", "Distance to Broadleaved Wood (km)", "Distance to Coniferous Wood (km)",
    "Proportion of All Wood within 1km", "Length of Woodland Perimeter within 1km", 
    "Distance to Woody Linear Feature (km)", "Slope (Degrees)", "Wader Density", 
    "Corvid Density", "Camera - Yes"
  )
  
  # Capture all dynamic levels for factors to include them in the ordering
  habitat_names <- final_summary_df$variable_display[grepl("Habitat: ", final_summary_df$variable_display)] %>% unique() %>% sort()
  project_year_names <- final_summary_df$variable_display[grepl("Project-Year: ", final_summary_df$variable_display)] %>% unique() %>% sort()
  
  variable_order <- c(fixed_order, habitat_names, project_year_names)
  
  final_summary_df$variable_display <- factor(final_summary_df$variable_display, 
                                              levels = variable_order, 
                                              ordered = TRUE)
  
  # Define the hierarchy of effect types
  group_order <- c(
    "Fixed", 
    "Factor Contrast - Habitat", 
    "Factor Contrast - Camera",
    "Random Effect - Project-Year Combination"
  )
  
  # CREATE PLOT TABLE STRUCTURE:
  # This section generates "Header" rows (empty rows with bold titles) for each effect category.
  processed_data <- final_summary_df %>%
    filter(!is.na(variable_display)) %>%
    mutate(dic_val = as.character(dic_difference)) %>%
    group_by(Effect_Type) %>%
    do({
      group_header <- data.frame(
        variable_display = unique(.$Effect_Type), 
        median = NA, Lower_95_CI = NA, Upper_95_CI = NA,
        Effect_Type = "HEADER", dic_val = NA,
        Original_Effect_Type = unique(.$Effect_Type) 
      )
      bind_rows(group_header, mutate(., Original_Effect_Type = Effect_Type))
    }) %>%
    ungroup()
  
  # Final table construction including the plot title row
  plot_data <- bind_rows(
    data.frame(variable_display = paste(species_name, "Variable"), median = NA, Lower_95_CI = NA, 
               Upper_95_CI = NA, Effect_Type = "TITLE", dic_val = "DIC",
               Original_Effect_Type = "TITLE"),
    processed_data
  ) %>%
    mutate(
      group_rank = match(Original_Effect_Type, c("TITLE", group_order)),
      type_rank = ifelse(Effect_Type %in% c("TITLE", "HEADER"), 1, 2),
      var_rank = match(variable_display, c(paste(species_name, "Variable"), variable_order))
    ) %>%
    arrange(group_rank, type_rank, var_rank)
  
  is_header <- plot_data$Effect_Type %in% c("TITLE", "HEADER")
  total_rows <- NROW(plot_data)
  
  # Format DIC labels for display (rounded to 2 decimal places)
  table_display <- plot_data %>%
    dplyr::select(variable_display, dic_val) %>%
    rename(Term = variable_display, 'DIC_Label' = dic_val) %>%
    mutate(
      DIC_Label = case_when(
        is.na(DIC_Label) ~ "",
        DIC_Label == "DIC" ~ "DIC",
        suppressWarnings(!is.na(as.numeric(DIC_Label))) ~ sprintf("%.2f", as.numeric(DIC_Label)),
        TRUE ~ DIC_Label
      )
    )
  
  # 6. VISUAL CUSTOMIZATION (LINES & STYLES) ------------------------------------------------------------------------
  
  fixed_top_lines <- list("1" = gpar(lwd = 1.5, col="black"), "2" = gpar(lwd = 1.5, col="black"))
  header_rows <- which(plot_data$Effect_Type == "HEADER")
  
  # Add subtle dashed separators between major categories
  group_separator_lines <- setNames(
    lapply(seq_along(header_rows), function(i) gpar(lwd = 0.8, col = "grey60", lty = 2)),
    as.character(header_rows)
  )
  
  bottom_line <- setNames(list(gpar(lwd = 1.5, col="black")), as.character(total_rows + 1))
  all_hrzl_lines <- c(fixed_top_lines, group_separator_lines, bottom_line)
  
  # FOREST PLOT FUNCTION:
  # Encapsulating the plot logic for easy export to multiple file formats
  generate_my_forest_plot <- function() {
    table_display %>%
      forestplot(
        labeltext = c(Term, DIC_Label),
        boxsize = 0.2,
        mean = plot_data$median,
        lower = plot_data$Lower_95_CI,
        upper = plot_data$Upper_95_CI,
        zero = 0,
        xlab = "Posterior Coefficient Estimate (Linear Predictor)",
        col = fpColors(box = "black", line = "#444444", zero = "#CCCCCC", summary = "black"),
        is.summary = is_header,
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
  
  # 7. EXPORT PLOTS -------------------------------------------------------------------------------------------------------
  
  fig_dir <- file.path(out_dir, "figures")
  if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  file_base <- file.path(fig_dir, paste0(run_id, "_Forest_Plot"))
  
  # Export PDF: Dynamic height adjustment ensures plots remain legible regardless of factor counts
  plot_height <- max(12, nrow(plot_data) * 0.3) 
  pdf(paste0(file_base, ".pdf"), width = 12, height = plot_height)
  print(generate_my_forest_plot())
  dev.off()
  
  # Export PNG: High resolution (150 DPI) for presentations/quick viewing
  png(paste0(file_base, ".png"), width = 1400, height = plot_height * 100, res = 150)
  print(generate_my_forest_plot())
  dev.off()
  
  message(paste("Success: Forest plot generated for", run_id))
}

print("Workflow complete: All forest plots exported.")