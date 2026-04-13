#' This function create a scatter plot on the residuals in the pls model
#' 
#' @param plsModel an object output from pls::plsr() function
#' @param numComp single numeric value, number of component
#' @param wavenumberRange numeric vector with length of 3. The first two values indicate the wavenumber
#' range of the plot, and the last value indicate the beak intervals on plot
plot_pls_residuals <- function(plsModel, numComp, wavenumberRange) {
  #Plot PLS residuals against GCMS signals under the selected number of
  #components
  plsPlotHeightResi <- ggplot() + 
    geom_point(aes(x = plsModel$fitted.values[,,numComp],
                   y = plsModel$residuals[,,numComp])) + 
    geom_hline(yintercept = 0) +
    labs(x = 'PLS fitted values (peak height)',
         y = 'Residuals') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    )
  
  return(plsPlotHeightResi)
}