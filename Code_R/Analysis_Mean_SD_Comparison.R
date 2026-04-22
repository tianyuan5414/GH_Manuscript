#' This function aims at calculating mean and sd of selected dataset

mean_sd_com <- function(
    dataList,
    dataListNameVec,
    varNameVec,
    creteriaName1Vec, 
    creteria1Vec, 
    creteriaName2Vec = rep(NA, length(dataList)),
    creteria2Vec = rep(NA, length(dataList)) 
  ){
  #Create a list to storage all t-test results
  tempListTtestResu <- matrix(nrow = length(dataListNameVec),
                              ncol = 3,
                              dimnames = list(c(),
                                              c('Group',
                                                'Mean',
                                                'SD')
                                              )
                              )
  tempListTtestResu[, 1] <- dataListNameVec
  
  #Loop to apply t test on every element in data list
  for (loopI in 1:length(dataList)) {
    #Check if there is a second criteria for current data frame
    if (is.na(creteriaName1Vec[loopI])) {
      tempVarIndex1 <- 1:length(dataList[[loopI]][[varNameVec[loopI]]])
    }else if (is.na(creteriaName2Vec[loopI])) {
      #Row indexes which meet the criteria given
      tempVarIndex1 <- which(dataList[[loopI]][[creteriaName1Vec[loopI]]] == creteria1Vec[loopI])
      
    } else {
      #Row indexes which meet the criteria given
      tempVarIndex1 <- which((dataList[[loopI]][[creteriaName1Vec[loopI]]] == creteria1Vec[loopI]) &
                               (dataList[[loopI]][[creteriaName2Vec[loopI]]] == creteria2Vec[loopI]))
    }
    
    tempListTtestResu[loopI, 2] <- mean(dataList[[loopI]][[varNameVec[loopI]]][tempVarIndex1],
                                        na.rm = TRUE)
    tempListTtestResu[loopI, 3] <- diff(range(dataList[[loopI]][[varNameVec[loopI]]][tempVarIndex1],
                                      na.rm = TRUE))
  }
  
  return(tempListTtestResu)
}