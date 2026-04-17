#' This function aims at picking out the number of component in PLS model to be 
#' plotted
  #' @param plsModel an object output, from pls::plsr() function
  #' @param threshold single numeric value, indicating the critic border judging the portion of variances explained
pick_num_comp <- function(plsModel, threshold) {
  #Calculate variance explained by different components
  varCompMat <- matrix(nrow = length(explvar(plsModel)),
                       ncol = 3,
                       dimnames = list(c(),
                                       c('Accu',
                                         'Partial',
                                         'Incremental')))
  varCompMat[, 2] <- explvar(plsModel)
  for (loopI in 1:nrow(varCompMat)) {
    varCompMat[loopI,1] <- sum(explvar(plsModel)[1:loopI])
    if (loopI == 1) {
      varCompMat[loopI, 3] <- 0
    }else{
      varCompMat[loopI, 3] <- varCompMat[loopI - 1, 1]
    }
  }
  
  tempVec <- which.max(varCompMat[,1] > (varCompMat[nrow(varCompMat), 1] * threshold))
  
  return(tempVec)
}