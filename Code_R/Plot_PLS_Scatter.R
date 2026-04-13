#' This function create a scatter plot between fitted PLS values and GC-MS signals
#' 
  #' @param plsModel an object output, from pls::plsr() function
  #' @param numComp single numeric value, number of component, which will determine the fitted values

plot_pls_scatter <- function(plsModel, numComp) {
  plsPlotHeight <- ggplot() + 
    geom_point(aes(x = plsModel$fitted.values[,,numComp],
                   y = plsModel$model$Y_train)) + 
    geom_smooth(aes(x = plsModel$fitted.values[,,numComp],
                    y = plsModel$model$Y_train),
                method = "lm") +
    geom_text(
      aes(x = 0.35, y = 0.7,
          label = lm_eqn2(numComp, plsModel)), parse = TRUE) +
    labs(x = 'PLS fitted values (peak height)',
         y = 'GC-MS (ng/150gr)') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    )
  
  return(plsPlotHeight)
}