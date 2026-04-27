#' This function aims at applying student-t test between selected a pair of 
#' data
#' 
sig_test_stu_t <- function(dataList,
                           dataListNameVec,
                           varNameVec,
                           creteriaNmae1Vec, 
                           creteria1Vec, 
                           creteriaName2Vec = rep(NA, length(dataList)),
                           creteria2Vec = rep(NA, length(dataList)),
                           creteriaName3Vec = rep(NA, length(dataList)),
                           creteria3Vec = rep(NA, length(dataList))
                           ) {
  #Create a list to storage all t-test results
  tempListTtestResu <- vector(mode = 'list',
                              length = length(dataList))
  names(tempListTtestResu) <- dataListNameVec
  
  #Loop to apply t test on every element in data list
  for (loopI in 1:length(dataList)) {
    #Check if there is a second creteria for current data frame
    if (is.na(creteriaName2Vec[loopI])) {
      #Row indexes which meet the criteria given
      tempVarIndex1 <- which(dataList[[loopI]][[creteriaNmae1Vec[loopI]]] == creteria1Vec[[loopI]][1])
      tempVarIndex2 <- which(dataList[[loopI]][[creteriaNmae1Vec[loopI]]] == creteria1Vec[[loopI]][2])
      
      #Matrix for recording filtered data
      tempStuTMat <- matrix(nrow = max(length(tempVarIndex1), length(tempVarIndex2)),
                            ncol = 2,
                            dimnames = list(c(),
                                            creteria1Vec[[loopI]]))
    } else if (is.na(creteriaName3Vec[loopI])) {
      #Row indexes which meet the criteria given
      tempVarIndex1 <- which((dataList[[loopI]][[creteriaNmae1Vec[loopI]]] == creteria1Vec[[loopI]]) &
                               (dataList[[loopI]][[creteriaName2Vec[loopI]]] == creteria2Vec[[loopI]][1]))
      tempVarIndex2 <- which((dataList[[loopI]][[creteriaNmae1Vec[loopI]]] == creteria1Vec[[loopI]]) &
                               (dataList[[loopI]][[creteriaName2Vec[loopI]]] == creteria2Vec[[loopI]][2]))
      
      #Matrix for recording filtered data
      tempStuTMat <- matrix(nrow = max(length(tempVarIndex1), length(tempVarIndex2)),
                            ncol = 2,
                            dimnames = list(c(),
                                            creteria2Vec[[loopI]]))
    } else {
      #Row indexes which meet the criteria given
      tempVarIndex1 <- which((dataList[[loopI]][[creteriaNmae1Vec[loopI]]] == creteria1Vec[[loopI]]) &
                               (dataList[[loopI]][[creteriaName2Vec[loopI]]] == creteria2Vec[[loopI]]) &
                               (dataList[[loopI]][[creteriaName3Vec[loopI]]] == creteria3Vec[[loopI]][1]))
      tempVarIndex2 <- which((dataList[[loopI]][[creteriaNmae1Vec[loopI]]] == creteria1Vec[[loopI]]) &
                               (dataList[[loopI]][[creteriaName2Vec[loopI]]] == creteria2Vec[[loopI]]) &
                               (dataList[[loopI]][[creteriaName3Vec[loopI]]] == creteria3Vec[[loopI]][2]))
      
      #Matrix for recording filtered data
      tempStuTMat <- matrix(nrow = max(length(tempVarIndex1), length(tempVarIndex2)),
                            ncol = 2,
                            dimnames = list(c(),
                                            creteria3Vec[[loopI]]))
    }
    
    tempStuTMat[1:length(tempVarIndex1), 1] <- dataList[[loopI]][[varNameVec[loopI]]][tempVarIndex1]
    tempStuTMat[1:length(tempVarIndex2), 2] <- dataList[[loopI]][[varNameVec[loopI]]][tempVarIndex2]
    
    tempListTtestResu[[loopI]] <- t.test(tempStuTMat[, 1], tempStuTMat[, 2],
                                         var.equal = FALSE, alternative = 'two.sided')
  }
  
  #Create matrix for recording all significance test results
  tempStuTMatAll <- matrix(nrow = length(dataList),
                           ncol = 7,
                           dimnames = list(c(),
                                           c('Group',
                                             'UV_Mean',
                                             'CON_Mean',
                                             'Sig',
                                             'Mark',
                                             't',
                                             'n')))
  
  tempStuTMatAll[, 1] <- dataListNameVec
  
  for (loopI in 1:nrow(tempStuTMatAll)) {
    
    tempStuTMatAll[loopI, 2:3] <- round(tempListTtestResu[[loopI]][[5]], 4)
    
    tempStuTMatAll[loopI, 4] <- round(tempListTtestResu[[loopI]][[3]], 4)
    
    tempStuTMatAll[loopI, 5] <- sigMarkter(tempStuTMatAll[loopI, 4], c(0.1, 0.05, 0.01))
    
    tempStuTMatAll[loopI, 6] <- round(tempListTtestResu[[loopI]][[1]], 3)
    
    tempStuTMatAll[loopI, 7] <- round(tempListTtestResu[[loopI]][[2]], 0)
  }
  
  #Output the significance test results
  return(tempStuTMatAll)
}