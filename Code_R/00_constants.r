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



#-----------------------------------------------------------
# 5. FTIR preprocessing variant library
#    Single standard variant; extend by adding entries.
#-----------------------------------------------------------

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