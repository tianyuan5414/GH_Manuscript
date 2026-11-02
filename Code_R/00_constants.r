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

multilevel_plsr <- function(df,
                            xcols,
                            ycol,
                            sitecol)
{
  
  X <- as.matrix(df[, xcols])
  y <- df[[ycol]]
  site <- factor(df[[sitecol]])
  
  ####################################
  # Between-site component
  ####################################
  
  X_between <- matrix(
    NA,
    nrow = nrow(X),
    ncol = ncol(X)
  )
  
  y_between <- rep(NA, length(y))
  
  for(s in levels(site)){
    
    idx <- which(site == s)
    
    Xmean <- colMeans(X[idx,,drop=FALSE])
    ymean <- mean(y[idx])
    
    X_between[idx,] <- matrix(
      rep(Xmean, length(idx)),
      nrow = length(idx),
      byrow = TRUE
    )
    
    y_between[idx] <- ymean
  }
  
  ####################################
  # Within-site component
  ####################################
  
  X_within <- X - X_between
  y_within <- y - y_between
  
  ####################################
  # Construct predictors
  ####################################
  
  colnames(X_between) <- paste0("B_", colnames(X))
  colnames(X_within) <- paste0("W_", colnames(X))
  
  Xall <- cbind(X_between, X_within)
  
  dat <- data.frame(y = y,
                    Xall)
  
  fit <- plsr(
    y ~ .,
    data = dat,
    validation = "CV",
    method = "simpls"
  )
  
  return(fit)
}
