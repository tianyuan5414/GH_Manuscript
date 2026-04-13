###########################################################
#
#  Configure R Environment and Workspace
#
###########################################################


#-----------------------------------------------------------
# 1. Load packages & internal functions
#-----------------------------------------------------------
# provide package list
package_list <- list(
  "assertthat",
  "tidyverse",
  "ggpubr",
  "here",
  "rlang",
  "usethis",
  "renv",
  "remotes",
  "future",
  "furrr",
  "purrr",
  "janitor",
  "quarto",
  "kableExtra",
  "languageserver",
  "analogue",
  "vegan",
  "rioja",
  "paletteer",
  "gWidgets2",
  "gWidgets2tcltk",
  "MASS",
  "devtools",
  "EMSC",
  "hyperSpec",
  "baseline",
  "pls",
  "MASS",
  "emmeans",
  "patchwork",
  "ggdist",
  "viridisLite"
)


# Load helper package `palaeochem` (preferred: installed package).
if (requireNamespace("palaeochem", quietly = TRUE)) {
  suppressPackageStartupMessages(library(palaeochem))
} else {
  palaeochem_path <- Sys.getenv("PALAEOCHEM_R", unset = NA)
  if (!is.na(palaeochem_path) && nzchar(palaeochem_path) && dir.exists(palaeochem_path)) {
    devtools::load_all(palaeochem_path, quiet = TRUE)
  } else {
    message(
      "The helper package 'palaeochem' is not installed and PALAEOCHEM_R is not set.\n",
      "Please either install it from GitHub with: remotes::install_github('PalaeoChem/palaeochem')\n",
      "or set PALAEOCHEM_R in your .Renviron to point to a local clone of the package."
    )
  }
}

# check and load R packages
check_and_load_packages(package_list)

library(showtext)
font_add(family = "arial", regular = file.path('Packages_R', 'Fonts', 'arial.ttf'))
showtext_auto()

#-----------------------------------------------------------
# 2. Load data compilation
#-----------------------------------------------------------
# load most recent data compilation

#-----------------------------------------------------------
# 3. Select datasets to process
#-----------------------------------------------------------
# Load project constants (dataset names, processing steps, colour palettes)
# — identical symbols shared with _targets.R
source(here::here("Code_R", "00_constants.r"))
source(here::here("Code_R", "Plot_PLS_Scatter.r"))
source(here::here("Code_R", "Plot_PLS_Coefficients.r"))
source(here::here("Code_R", "Plot_PLS_Loadings.r"))
source(here::here("Code_R", "Plot_PLS_Variance.r"))
source(here::here("Code_R", "Plot_PLS_Residual.r"))

#-----------------------------------------------------------
# 4. Interactive-only: ggplot theme
#-----------------------------------------------------------
ggplot2::theme_set(
  ggplot2::theme_bw()
)
