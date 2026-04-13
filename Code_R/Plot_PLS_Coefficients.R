#' This function create a bar plot on the correlation coefficients at all wavenumber
#' input into the pls model
#' 
  #' @param plsModel an object output from pls::plsr() function
  #' @param numComp single numeric value, number of component, which will determine the fitted values
  #' @param wavenumberRange numeric vector with length of 3. The first two values indicate the wavenumber
    #' range of the plot, and the last value indicate the beak intervals on plot
  #' @param peakPos sequence of wavenumnber for the pls fitted values
plot_pls_coefficients <- function(plsModel, numComp, wavenumberRange,
                                  peakPos) {
  modelHeightCof <- ggplot() +
    geom_bar(aes(x = peakPos, y = plsModel$coefficients[, , numComp]),
             stat = 'identity') +
    geom_hline(yintercept = 0) +
    scale_x_reverse(limit = wavenumberRange,
                    breaks = seq(wavenumberRange[1],
                                 wavenumberRange[2], 
                                 by = wavenumberRange[3])) +
    labs(x = 'Wavenumber', y = paste('Variable coefficients ',
                                     numComp, ' components', sep = '', collapse = NULL)) +
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)

    )
  
  return(modelHeightCof)
}