# Wader_Nest_Predation_Analysis_2026
Scripts and artificial data for the manuscript: Predation-mediated edge effects on wader nests in rural landscapes vary with distance to woodland and woodland type.


# ===================================================================================================================================================
# PROJECT: Woodland Proximity and Wader Nest Predation Analysis
# REPO: https://github.com/charlotte-lemarquand/wader-woodland-predation
# AUTHOR: Charlotte Le Marquand (charlotte.lemarquand@york.ac.uk)
# LAST EDIT: 08/05/2026
#
# DESCRIPTION: 
#   This repository contains the analytical pipeline and visualization scripts for 
#   investigating the effects of woodland features on wader nest predation. The analysis 
#   utilizes a Bayesian framework (INLA) to account for spatial autocorrelation and 
#   varied environmental predictors.
# ===================================================================================================================================================

# 1. SCRIPT WORKFLOW --------------------------------------------------------------------------------------------------------------------------------
# Scripts are designed to be run sequentially. Each script reads from a central 
# 'model_runs.csv' control file to automate iterations across species and data groups.

# [01_run_models_and_dic_testing.R]
#   The core analytical script. Handles data cleaning, spatial mesh construction (SPDE), 
#   and model fitting via 'inlabru'. Includes a DIC-based variable importance loop.

# [02_visualisation_prediction_plots.R]
#   Translates model outputs into marginal effect plots. Back-transforms daily 
#   risk into cumulative predation probability over the full incubation period. 
#   Includes data density rug plots and biological thresholds (e.g., 1km markers).

# [03_visualisation_tables.R]
#   Automates extraction of posterior summaries (means, medians, 95% CIs). 
#   Maps technical internal variable names to descriptive manuscript labels.

# [04_visualisation_forest_plots.R]
#   Generates high-resolution forest plots of all model parameters (Fixed, 
#   Factor Contrasts, and Random Effects) for a comprehensive summary.

# [05_visualisation_forest_plots.R]
#   A targeted visualization isolating primary woodland-related research 
#   hypotheses (proximity, area, and composition).


# 2. DATA AVAILABILITY & ETHICS ---------------------------------------------------------------------------------------------------------------------

# ARTIFICIAL DATASET:
#   The data provided in the '/data' folder is ARTIFICIAL. It mimics the structure 
#   and properties of the original dataset while protecting sensitive 
#   nest locations and landowner privacy. It allows for workflow reproduction, 
#   though coefficients will differ from published results. 

# ACCESSING PRIMARY DATA:
#   Primary data belongs to the RSPB (Royal Society for the Protection of Birds) 
#   and the GWCT (Game & Wildlife Conservation Trust).
#   
#   Due to the sensitivity of wader nesting locations:
#     - Raw data is NOT hosted in this public repository.
#     - Reasonable requests for the actual data for research purposes will 
#       be considered by the RSPB and GWCT.
#     - Inquiries should be directed to the corresponding author.


# 3. REQUIREMENTS -----------------------------------------------------------------------------------------------------------------------------------

# DEPENDENCIES:
#   The following packages are required:
#   install.packages(c("tidyverse", "sf", "patchwork", "forestplot", "readr", "inlabru))
#   
#   INLA must be installed from the R-INLA repository:
#   install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)


# 4. CITATION ---------------------------------------------------------------------------------------------------------------------------------------

# If using this code for your research, please cite:
# [Citation will be updated on publication of research, if this has not been updated please contact the corresponding author].


# ===================================================================================================================================================
# END OF README
# ===================================================================================================================================================
