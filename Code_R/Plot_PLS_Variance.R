#' This function create a bar plot on the variance explained by differernt components
#' in the pls model
#' 
  #' @param plsModel an object output from pls::plsr() function
  #' @param numComp single numeric value, number of component
plot_pls_variance <- function(plsModel, numComp) {
  #Calculate variance explained by different components
  varCompMat <- matrix(nrow = numComp,
                       ncol = 3,
                       dimnames = list(c(),
                                       c('Accu',
                                         'Partial',
                                         'Incremental')))
  varCompMat[, 2] <- explvar(plsModel)[1:numComp]
  for (loopI in 1:nrow(varCompMat)) {
    varCompMat[loopI,1] <- sum(explvar(plsModel)[1:loopI])
    if (loopI == 1) {
      varCompMat[loopI, 3] <- 0
    }else{
      varCompMat[loopI, 3] <- varCompMat[loopI - 1, 1]
    }
  }
  
  #Plot variance explained by different components
  modelHeightVar <- ggplot() +
    geom_bar(aes(x = 1:numComp, y = varCompMat[, 1]),
             fill = 'red',
             stat = 'identity') +
    geom_bar(aes(x = 1:numComp, y = varCompMat[, 3]),
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