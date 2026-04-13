#' This function create a scatter plot between fitted PLS values and GC-MS signals
#' 
  #' @param plsModel an object output from pls::plsr() function
  #' @param numComp number of component, which will determine the fitted values

plot_pls_scatter <- function(plsModel, numComp) {
  plsPlotHeight <- ggplot() + 
    geom_point(aes(x = plsModel$fitted.values[,,numComp],
                   y = plsModel$model$calibrated)) + 
    geom_smooth(aes(x = plsModel$fitted.values[,,numComp],
                    y = plsModel$model$calibrated),
                method = "lm") +
    geom_text(
      aes(x = 0.35, y = 0.7,
          label = lm_eqn2(numComp)), parse = TRUE) +
    labs(x = 'PLS fitted values (peak height)',
         y = 'GC-MS (ng/150gr)') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  return(plsPlotHeight)
}