#' This function create a bar plot on the variance explained by differernt components
#' in the pls model
#' 
  #' @param plsModel an object output from pls::plsr() function
  #' @param numComp single numeric value, number of component
  #' @param wavenumberRange numeric vector with length of 3. The first two values indicate the wavenumber
    #' range of the plot, and the last value indicate the beak intervals on plot
plot_pls_variance <- function(plsModel, numComp, wavenumberRange) {
  #Calculate variance explained by different components
  varCompMat <- matrix(nrow = numComp,
                       ncol = 2,
                       dimnames = list(c(),
                                       c('Accu',
                                         'Partial')))
  varCompMat[, 1] <- plsModel[[12]][1:numComp]
  for (loopI in 1:nrow(varCompMat)) {
    varCompMat[loopI, 2] <- sum(plsModel[[12]][1:loopI])
  }
  
  #Plot variance explained by different components
  modelHeightVar <- ggplot() +
    geom_bar(aes(x = 1:numComp, y = varCompMat[, 2] / sum(plsModel[[12]])),
             fill = 'red',
             stat = 'identity') +
    geom_bar(aes(x = 1:numComp, y = varCompMat[, 2] / sum(plsModel[[12]] )- 
                   varCompMat[, 1] / sum(plsModel[[12]])),
             fill = 'grey',
             stat = 'identity') +
    labs(x = 'Component', y = 'Variance explained') +
    scale_x_continuous(breaks = seq(1, numComp, by = 1)) +
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    )
  
  return(modelHeightVar)
}