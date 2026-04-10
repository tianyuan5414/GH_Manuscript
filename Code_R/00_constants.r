###########################################################
#
#  Project-wide constants: dataset names, processing step
#  lists, and colour palettes.
#
#  No side effects: no library(), no data loading, no I/O.
#  Sourced by both 00_Config_file.r (interactive) and
#  _targets.R (pipeline) so both share identical symbols.
#
###########################################################

#-----------------------------------------------------------
# 1. Dataset name vectors
#-----------------------------------------------------------

# FTIR ASE experiments
ftir_ase_exp <- c(
  'FTIR_Greenhouse_experiment_2025'
)

#-----------------------------------------------------------
# 2. FTIR (ATR) preprocessing stepss
#-----------------------------------------------------------

#0 - No treatment
processing_steps_ftir0 <- list(
  functions = c(
  ),
  parameters = list(
  )
)


#1 - No treatment - only extract wavenumbers
processing_steps_ftir1 <- list(
  functions = c(
    get_filter_wavenumbers
  ),
  parameters = list(
    list(min_wn = 600, max_wn = 3996)
  )
)

#2 - MSC
processing_steps_ftir2 <- list(
  functions = c(
    run_msc,
    get_filter_wavenumbers,
    get_baseline_offset
  ),
  parameters = list(
    list(),
    list(min_wn = 600, max_wn = 3996),
    list()
  )
)

#3 - ALS
processing_steps_ftir3 <- list(
  functions = c(
    run_als_baseline_correction,
    get_filter_wavenumbers,
    get_baseline_offset
  ),
  parameters = list(
    list(lambda = 5, p = 0.001),
    list(min_wn = 600, max_wn = 3996),
    list()
  )
)

#4 - ALS-MSC
processing_steps_ftir4 <- list(
  functions = c(
    run_als_baseline_correction,
    get_filter_wavenumbers,
    run_msc,
    get_baseline_offset
  ),
  parameters = list(
    list(lambda = 5, p = 0.001),
    list(min_wn = 600, max_wn = 3996),
    list(),
    list()
  )
)

#5 - SG-2nd derivative
processing_steps_ftir5 <- list(
  functions = c(
    get_savitzkygolay,
    get_filter_wavenumbers
  ),
  parameters = list(
    list(width = 7, poly = 2, deriv = 2),
    list(min_wn = 600, max_wn = 3996)
  )
)

#6 - signature region - MSC
processing_steps_ftir6 <- list(
  functions = c(
    run_msc,
    get_filter_wavenumbers,
    get_baseline_offset
  ),
  parameters = list(
    list(),
    list(min_wn = 800, max_wn = 1800),
    list()
  )
)

#7 - signature region - ALS
processing_steps_ftir7 <- list(
  functions = c(
    run_als_baseline_correction,
    get_filter_wavenumbers,
    get_baseline_offset
  ),
  parameters = list(
    list(lambda = 5, p = 0.001),
    list(min_wn = 800, max_wn = 1800),
    list()
  )
)

#8 - signature region - ALS-MSC
processing_steps_ftir8 <- list(
  functions = c(
    run_als_baseline_correction,
    get_filter_wavenumbers,
    run_msc,
    get_baseline_offset
  ),
  parameters = list(
    list(lambda = 5, p = 0.001),
    list(min_wn = 800, max_wn = 1800),
    list(),
    list()
  )
)

#9 - signature region - GC-2nd derivative
processing_steps_ftir9 <- list(
  functions = c(
    get_savitzkygolay,
    get_filter_wavenumbers
  ),
  parameters = list(
    list(width = 7, poly = 2, deriv = 2),
    list(min_wn = 800, max_wn = 1800)
  )
)

#-----------------------------------------------------------
# 5. FTIR preprocessing variant library
#    Single standard variant; extend by adding entries.
#-----------------------------------------------------------

preprocessing_variants_ftir0 <- list(

  ftir_standard = list(
    steps         = processing_steps_ftir0,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir1 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir1,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir2 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir2,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir3 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir3,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir4 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir4,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir5 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir5,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir6 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir6,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir7 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir7,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir8 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir8,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

preprocessing_variants_ftir9 <- list(
  
  ftir_standard = list(
    steps         = processing_steps_ftir9,
    group_vars    = NULL,
    data_source   = "ftir",
    is_derivative = FALSE
  )
)

#-----------------------------------------------------------
# 6. Global plotting parameters
#-----------------------------------------------------------

# define color palette for sample treatments
colors_treatments <- setNames(
  viridisLite::viridis(
    5,
    option = "turbo",
    begin = 0.05,
    end = 0.9
  ),
  c(
    "subfossil",
    "Hexane+Acetone+MeOH_90C",
    "Hexane+Acetone",
    "Hexane",
    "Untreated"
  )
)