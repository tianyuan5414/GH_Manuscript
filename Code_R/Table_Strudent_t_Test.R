#Load necessary packages--------------------------------------------------------
{
  #String edit
  library(stringr)
}

#Declare necessary functions----------------------------------------------------
{
  #Create a function to output significance markers
  #tempSig - single numeric value
  #gradient - numeric vector, significance criteria judging the significant 
  #markers for different levels, should be at length of 3, and in decreasing
  #order
  sigMarkter <- function(tempSig, gradient) {
    
    #Significant level 1, marked by *
    if ((as.numeric(tempSig) < gradient[1]) & 
        (as.numeric(tempSig) > gradient[2])) {
      
      return('*')
      
      #Significance level 2, marked by **  
    } else if ((as.numeric(tempSig) < gradient[2]) & 
               (as.numeric(tempSig) > gradient[3])) {
      
      return('**')
      
      #Significance level 3, marked by ***
    } else if (as.numeric(tempSig) < gradient[3]) {
      
      return('***')
      
    } else {
      
      #Not significant, marked by nothing
      return('')
      
    }
    
  }
}

#Read data of GCMS and FTIR PLS fitted values-----------------------------------
{
 #GCMS
  GCDataMean <- read.csv(file = file.path('Output_R',
                                          'GCDataMean.csv'),
                         row.names = 1)
  
  #FTIR
  pcaMatFTIRGCMalFull <- read.csv(file = 
                                        file.path(
                                          'Output_R',
                                          'pcaMatFTIRGCMalFull.csv'),
                                          row.names = 1)
}

#Statistic analysis on peak signals and malformation rates----------------------
{
  #GCMS
  {
    #GCMS whole
    tempVec1 <- which(GCDataMean[, 4] == 'UV')
    tempVec2 <- which(GCDataMean[, 4] == 'CON')
    
    stuTMatGCMS <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                          ncol = 2,
                          dimnames = list(c(),
                                          c('UV', 'CON')))
    
    stuTMatGCMS[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMS[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMS <- t.test(stuTMatGCMS[, 1], stuTMatGCMS[, 2],
                          var.equal = FALSE, alternative = 'greater')
    
    #GCMS whole subgroup UV
    tempVec1 <- which((GCDataMean[, 4] == 'UV') &
                        (is.na(GCDataMean[, 6])))
    tempVec2 <- which(GCDataMean[, 4] == 'UV'&
                        (!is.na(GCDataMean[, 6])))
    
    stuTMatGCMSSubUV <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                               ncol = 2,
                               dimnames = list(c(),
                                               c('Cul', 'Cut')))
    
    stuTMatGCMSSubUV[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSSubUV[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSSubUV <- t.test(stuTMatGCMSSubUV[, 1], stuTMatGCMSSubUV[, 2],
                               var.equal = FALSE, alternative = 'greater')
    
    #GCMS whole subgroup CON
    tempVec1 <- which((GCDataMean[, 4] == 'CON') &
                        (is.na(GCDataMean[, 6])))
    tempVec2 <- which(GCDataMean[, 4] == 'CON'&
                        (!is.na(GCDataMean[, 6])))
    
    stuTMatGCMSSubCON <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('Cul', 'Cut')))
    
    stuTMatGCMSSubCON[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSSubCON[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSSubCON <- t.test(stuTMatGCMSSubCON[, 1], stuTMatGCMSSubCON[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    #GCMS cul and cut subgroup
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        is.na(GCDataMean[, 6]))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        is.na(GCDataMean[, 6]))
    
    stuTMatGCMSCul <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                             ncol = 2,
                             dimnames = list(c(),
                                             c('UV', 'CON')))
    
    stuTMatGCMSCul[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSCul[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSCul <- t.test(stuTMatGCMSCul[, 1], stuTMatGCMSCul[, 2],
                             var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]))
    
    stuTMatGCMSCut <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                             ncol = 2,
                             dimnames = list(c(),
                                             c('UV', 'CON')))
    
    stuTMatGCMSCut[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSCut[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSCut <- t.test(stuTMatGCMSCut[, 1], stuTMatGCMSCut[, 2],
                             var.equal = FALSE, alternative = 'greater')
    
    #GCMS Plant subgroup
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '1'))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '1'))
    
    stuTMatGCMSPlant1 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatGCMSPlant1[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSPlant1[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSPlant1 <- t.test(stuTMatGCMSPlant1[, 1], stuTMatGCMSPlant1[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '2'))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '2'))
    
    stuTMatGCMSPlant2 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatGCMSPlant2[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSPlant2[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSPlant2 <- t.test(stuTMatGCMSPlant2[, 1], stuTMatGCMSPlant2[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '3'))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '3'))
    
    stuTMatGCMSPlant3 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatGCMSPlant3[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSPlant3[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSPlant3 <- t.test(stuTMatGCMSPlant3[, 1], stuTMatGCMSPlant3[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '4'))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '4'))
    
    stuTMatGCMSPlant4 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatGCMSPlant4[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSPlant4[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSPlant4 <- t.test(stuTMatGCMSPlant4[, 1], stuTMatGCMSPlant4[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '5'))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '5'))
    
    stuTMatGCMSPlant5 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatGCMSPlant5[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSPlant5[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSPlant5 <- t.test(stuTMatGCMSPlant5[, 1], stuTMatGCMSPlant5[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((GCDataMean[, 4] == 'UV')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '6'))
    tempVec2 <- which((GCDataMean[, 4] == 'CON')&
                        !is.na(GCDataMean[, 6]) &
                        (GCDataMean[, 5] == '6'))
    
    stuTMatGCMSPlant6 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatGCMSPlant6[1:length(tempVec1), 1] <- GCDataMean[tempVec1, 1]
    stuTMatGCMSPlant6[1:length(tempVec2), 2] <- GCDataMean[tempVec2, 1]
    
    stuTResGCMSPlant6 <- t.test(stuTMatGCMSPlant6[, 1], stuTMatGCMSPlant6[, 2],
                                var.equal = FALSE, alternative = 'greater')
  }
  
  #FTIR PLS
  {
    #FTIRPLS whole
    tempVec1 <- which(pcaMatFTIRGCMalFull[, 7] == 'UV')
    tempVec2 <- which(pcaMatFTIRGCMalFull[, 7] == 'CON')
    
    stuTMatFTIRPLS <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                             ncol = 2,
                             dimnames = list(c(),
                                             c('UV', 'CON')))
    
    stuTMatFTIRPLS[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLS[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLS <- t.test(stuTMatFTIRPLS[, 1], stuTMatFTIRPLS[, 2],
                             var.equal = FALSE, alternative = 'greater')
    
    #FTIRPLS whole subgroup UV
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV') &
                        (pcaMatFTIRGCMalFull[, 5] == 'cul'))
    tempVec2 <- which(pcaMatFTIRGCMalFull[, 7] == 'UV'&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut'))
    
    stuTMatFTIRPLSSubUV <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                  ncol = 2,
                                  dimnames = list(c(),
                                                  c('Cul', 'Cut')))
    
    stuTMatFTIRPLSSubUV[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSSubUV[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSSubUV <- t.test(stuTMatFTIRPLSSubUV[, 1], stuTMatFTIRPLSSubUV[, 2],
                                  var.equal = FALSE, alternative = 'greater')
    
    #FTIRPLS whole subgroup CON
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON') &
                        (pcaMatFTIRGCMalFull[, 5] == 'cul'))
    tempVec2 <- which(pcaMatFTIRGCMalFull[, 7] == 'CON'&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut'))
    
    stuTMatFTIRPLSSubCON <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('Cul', 'Cut')))
    
    stuTMatFTIRPLSSubCON[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSSubCON[1:length(tempVec2), 2] <-pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSSubCON <- t.test(stuTMatFTIRPLSSubCON[, 1], stuTMatFTIRPLSSubCON[, 2],
                                   var.equal = FALSE, alternative = 'greater')
    
    #FTIRPLS cul and cut subgroup
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cul'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cul'))
    
    stuTMatFTIRPLSCul <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatFTIRPLSCul[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSCul[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSCul <- t.test(stuTMatFTIRPLSCul[, 1], stuTMatFTIRPLSCul[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut'))
    
    stuTMatFTIRPLSCut <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                ncol = 2,
                                dimnames = list(c(),
                                                c('UV', 'CON')))
    
    stuTMatFTIRPLSCut[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSCut[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSCut <- t.test(stuTMatFTIRPLSCut[, 1], stuTMatFTIRPLSCut[, 2],
                                var.equal = FALSE, alternative = 'greater')
    
    #FTIRPLS Plant subgroup
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '1'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '1'))
    
    stuTMatFTIRPLSPlant1 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('UV', 'CON')))
    
    stuTMatFTIRPLSPlant1[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSPlant1[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSPlant1 <- t.test(stuTMatFTIRPLSPlant1[, 1], stuTMatFTIRPLSPlant1[, 2],
                                   var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '2'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '2'))
    
    stuTMatFTIRPLSPlant2 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('UV', 'CON')))
    
    stuTMatFTIRPLSPlant2[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSPlant2[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSPlant2 <- t.test(stuTMatFTIRPLSPlant2[, 1], stuTMatFTIRPLSPlant2[, 2],
                                   var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '3'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '3'))
    
    stuTMatFTIRPLSPlant3 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('UV', 'CON')))
    
    stuTMatFTIRPLSPlant3[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSPlant3[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSPlant3 <- t.test(stuTMatFTIRPLSPlant3[, 1], stuTMatFTIRPLSPlant3[, 2],
                                   var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '4'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '4'))
    
    stuTMatFTIRPLSPlant4 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('UV', 'CON')))
    
    stuTMatFTIRPLSPlant4[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSPlant4[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSPlant4 <- t.test(stuTMatFTIRPLSPlant4[, 1], stuTMatFTIRPLSPlant4[, 2],
                                   var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '5'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '5'))
    
    stuTMatFTIRPLSPlant5 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('UV', 'CON')))
    
    stuTMatFTIRPLSPlant5[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSPlant5[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSPlant5 <- t.test(stuTMatFTIRPLSPlant5[, 1], stuTMatFTIRPLSPlant5[, 2],
                                   var.equal = FALSE, alternative = 'greater')
    
    tempVec1 <- which((pcaMatFTIRGCMalFull[, 7] == 'UV')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '6'))
    tempVec2 <- which((pcaMatFTIRGCMalFull[, 7] == 'CON')&
                        (pcaMatFTIRGCMalFull[, 5] == 'cut') &
                        (as.character(as.numeric(pcaMatFTIRGCMalFull[, 8])) == '6'))
    
    stuTMatFTIRPLSPlant6 <- matrix(nrow = max(length(tempVec1), length(tempVec2)),
                                   ncol = 2,
                                   dimnames = list(c(),
                                                   c('UV', 'CON')))
    
    stuTMatFTIRPLSPlant6[1:length(tempVec1), 1] <- pcaMatFTIRGCMalFull[tempVec1, 11]
    stuTMatFTIRPLSPlant6[1:length(tempVec2), 2] <- pcaMatFTIRGCMalFull[tempVec2, 11]
    
    stuTResFTIRPLSPlant6 <- t.test(stuTMatFTIRPLSPlant6[, 1], stuTMatFTIRPLSPlant6[, 2],
                                   var.equal = FALSE, alternative = 'greater')
  }
  
  #Combine all ANOVA results
  {
    #Create matrix for recording all significance test results
    stuTMatAll <- matrix(nrow = 22,
                         ncol = 5,
                         dimnames = list(c(),
                                         c('Group',
                                           'UV_Mean',
                                           'CON_Mean',
                                           'Sig',
                                           'Mark')))
    
    stuTMatAll[, 1] <- c(paste(rep('GCMS', 11), c('', 'SubUV', 'SubCON',
                                                  'Cut', 'Cul',
                                                  paste('Plant', 1:6, sep = ''
                                                  )
    ),
    sep = ''),
    paste(rep('FTIRPLS', 11), c('', 'SubUV', 'SubCON',
                                'Cut', 'Cul',
                                paste('Plant', 1:6, sep = ''
                                )
    ),
    sep = '')
    )
    
    for (loopI in 1:nrow(stuTMatAll)) {
      
      stuTMatAll[loopI, 2:3] <- get(str_replace_all(paste('stuTRes', stuTMatAll[loopI, 1],
                                                          sep = '', collapse = NULL), ' ', ''))[[5]]
      
      stuTMatAll[loopI, 4] <- get(str_replace_all(paste('stuTRes', stuTMatAll[loopI, 1], 
                                                        sep = '', collapse = NULL), ' ', ''))[[3]]
      
      stuTMatAll[loopI, 5] <- sigMarkter(stuTMatAll[loopI, 4], c(0.1, 0.05, 0.01))
    }
  }
  
  #Output results
  {
    write.csv(stuTMatAll, 
              file = file.path( 'Output_R', 'stuTMatAll.csv'))
  }
  
  #Remove unused variables
  rm(loopI)
}