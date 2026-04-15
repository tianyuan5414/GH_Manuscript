#This script combines malformation rates with PLS FTIR fitted values, and GC-MS signals
  #' @param batchInfoMat matrix or data frame, describing the batch information for each sample
  #' @param PLSFittedValues numeric vector, providing PLS fitted values from FTIR spectra
  #' @param GCMSValuesMat numeric vector, providing GC-MS signals
  #' @param cultivarsData matrix or data frame, containing batch information 
    #' and malformation counting results from potted cultivars
  #' @param cuttingsData matrix or data frame of branch cuttings counting 
    #' results, containing group id in the first column, and the counting 
    #' results in the other columns
  #' @param cuttingsBatch matrix or data frame, containing batch information of
    #' branch cuttings
matrix_malformation <- function(batchInfoMat,
                                PLSFittedValues,
                                GCMSValuesMat,
                                cultivarsData,
                                cuttingsData,
                                cuttingsBatch) {
  #Clean up samples that do not have GCMS signals
  batchInfoMat <- batchInfoMat[which(!is.na(GCMSValuesMat[, 1])), ]
  GCMSValuesMat <- GCMSValuesMat[which(!is.na(GCMSValuesMat[, 1])), ]
  
  #Add batch information to FTIR data
  pcaMatFTIR <- cbind(batchInfoMat,
                      PLSFittedValues,
                      GCMSValuesMat)
  
  colnames(pcaMatFTIR)[-1:-10] <- c('PLS', 'GCMS', 'GCMS_SD')
  
  #Vector for recording row index in FTIR spectra matrix
  tempInd <- vector(length = nrow(pcaMatFTIR))
  
  #Combine batch information of cuttings and pot cultivars together
  tempMat <- rbind(cultivarsData[,1:6],
                   cuttingsBatch)
  
  #Add malformation data to FTIR data
  for (loopI in 1:nrow(pcaMatFTIR)) {
    
    #Check and skip the roof samples
    if (pcaMatFTIR[loopI, 3] == 'ROOF') {
      
      tempInd[loopI] <- NA
      
    }else{
      
      #For all pot cultivars
      if (pcaMatFTIR[loopI, 5] == 'cul') {
        
        if (any((tempMat[, 3] == batchInfoMat[loopI, 6]) &
                (tempMat[, 4] == batchInfoMat[loopI, 7]) &
                (tempMat[, 5] == batchInfoMat[loopI, 8]) &
                (!is.na(as.numeric(tempMat[, 6])))
        )
        ) {
          
          #Record the row index for FTIR spectra matrix
          tempInd[loopI] <- which(
            (tempMat[, 3] == pcaMatFTIR[loopI, 6]) &
              (tempMat[, 4] == pcaMatFTIR[loopI, 7]) &
              (tempMat[, 5] == 
                 as.character(as.numeric(pcaMatFTIR[loopI, 8]))) &
              (!is.na(as.numeric(tempMat[, 6])))
          )
          
        }else {
          
          tempInd[loopI] <- NA
          
        }
        
        #For all cuttings
      } else if (pcaMatFTIR[loopI, 5] == 'cut') {
        
        #Check if current FTIR sample has malformation counting results
        if (any((tempMat[, 3] == pcaMatFTIR[loopI, 6]) &
                (tempMat[, 4] == pcaMatFTIR[loopI, 7]) &
                (tempMat[, 5] == as.character(as.numeric(pcaMatFTIR[loopI, 8]))) &
                (tempMat[, 6] == pcaMatFTIR[loopI, 9])
        )
        ) {
          
          #Record the row index for FTIR spectra matrix
          tempInd[loopI] <- which(
            (tempMat[, 3] == pcaMatFTIR[loopI, 6]) &
              (tempMat[, 4] == pcaMatFTIR[loopI, 7]) &
              (tempMat[, 5] == 
                 as.character(as.numeric(pcaMatFTIR[loopI, 8]))) &
              (tempMat[, 6] == pcaMatFTIR[loopI, 9]))
          
        }else {
          
          tempInd[loopI] <- NA
          
        }
        
      } else {
        
        tempInd[loopI] <- NA
        
      }
    }
  }
  
  #Combine all counting results together from both cuttings and pot cultivars
  tempMat <- rbind(cultivarsCom[, -1:-5],
                   as.matrix(cuttingsData[, -1]))
  
  #Add malformation counting results to FTIR and GCMS matrix
  pcaMatFTIRGCMal <- cbind(pcaMatFTIR, tempMat[tempInd, c(1, 14)])
  
  #Adjest value type to numeric
  pcaMatFTIRGCMal[, 11:15] <- apply(pcaMatFTIRGCMal[, 11:15], c(1,2), as.numeric)
  
  #Calculate the malformation ratio in each sample, and add them to the whole
  #matrix
  pcaMatFTIRGCMalFull <- cbind(pcaMatFTIRGCMal, 
                               Mal_Ratio = 
                                 (pcaMatFTIRGCMal[,ncol(pcaMatFTIRGCMal)] - pcaMatFTIRGCMal[,ncol(pcaMatFTIRGCMal) - 1]) / 
                                 pcaMatFTIRGCMal[,ncol(pcaMatFTIRGCMal)])
  
  #Remove samples that do not have malformation counting results
  pcaMatFTIRGCMalFull <- pcaMatFTIRGCMalFull[which(!is.na(pcaMatFTIRGCMalFull[, 6])),]
  
  #Calculate the mean and sd of FTIR duplicates
  #Get unique batch information
  tempMat <- unique(pcaMatFTIRGCMalFull[, c(1:2, 4:8)])
  
  #Matrix for recording mean and sd of FTIR duplicates
  tempMat2 <- matrix(nrow = nrow(tempMat),
                     ncol = ncol(pcaMatFTIRGCMalFull) - 1,
                     dimnames = list(c(),
                                     c(colnames(pcaMatFTIRGCMalFull)[-9:-10], 
                                       paste(colnames(pcaMatFTIRGCMalFull)[11], '_SD')
                                     )
                     )
  )
  
  tempMat2[,c(1:2, 4:8)] <- as.matrix(tempMat)
  
  #Calculate the mean FTIR signals
  for (loopI in 1:nrow(tempMat2)) {
    
    #Row indexes for all duplicates of current sample
    tempVec <- which((pcaMatFTIRGCMalFull[, 5] == tempMat2[loopI, 5]) &
                       (pcaMatFTIRGCMalFull[, 6] == tempMat2[loopI, 6]) &
                       (pcaMatFTIRGCMalFull[, 7] == tempMat2[loopI, 7]) &
                       (as.numeric(pcaMatFTIRGCMalFull[, 8]) == as.numeric(tempMat2[loopI, 8]))
    )
    
    #Check whetehr there are duplicates
    if (length(tempVec) > 1) {
      
      #Calculate mean and sd of duplicates
      tempMat2[loopI, 9:14] <- apply(pcaMatFTIRGCMalFull[tempVec, 11:16],
                                     2,
                                     mean)
      tempMat2[loopI, 15] <- sd(pcaMatFTIRGCMalFull[tempVec, 11])
    }else {
      
      #Copy the signals to corresponding position
      tempMat2[loopI, 9:14] <- as.numeric(pcaMatFTIRGCMalFull[tempVec, 11:16])
      #Sd is assigned to 0 for no-duplicate sample
      tempMat2[loopI, 15] <- 0
      
    }
  }
  
  #Transfer the combined matrix to named matrix
  pcaMatFTIRGCMalFullMean <- tempMat2
  pcaMatFTIRGCMalFullMean <- as.data.frame(pcaMatFTIRGCMalFullMean)
  pcaMatFTIRGCMalFullMean[, 9:15] <- apply(pcaMatFTIRGCMalFullMean[, 9:15], 
                                           c(1,2),
                                           as.numeric)
  
  rm(tempInd,
     tempMat,
     tempMat2,
     loopI
  )
  
  #Create a matrix recording all signifcance results bettwen diffent subgroups
  cuttingsPlantsMat <- matrix(nrow = length(unique(
    pcaMatFTIRGCMalFullMean[which(pcaMatFTIRGCMalFullMean[, 5] == 'cut'),
                            8
    ]
  )
  ),
  ncol = 17,
  dimnames = list(c(),
                  c('Plant', 'UV_Mal', 
                    'CON_Mal',
                    'GCMS',
                    'PLS',
                    'GCMS_UV',
                    'GCMS_CON',
                    'PLS_UV',
                    'PLS_CON',
                    'UV_Mal_sd',
                    'CON_Mal_sd',
                    'GCMS_sd',
                    'PLS_sd',
                    'GCMS_UV_sd',
                    'GCMS_CON_sd',
                    'PLS_UV_sd',
                    'PLS_CON_sd')))
  
  #Fill batch information into the matrix
  cuttingsPlantsMat[, 1] <- sort(as.numeric(unique(
    pcaMatFTIRGCMalFullMean[which(pcaMatFTIRGCMalFullMean[, 5] == 'cut'),
                            8
    ])), decreasing = F)
  
  #Row indexes for all cuttings in the general matrix
  tempMat <-  pcaMatFTIRGCMalFullMean[
    which(pcaMatFTIRGCMalFullMean[, 5] == 'cut'),]
  
  #Loop to fill the matrix
  for (loopI in cuttingsPlantsMat[, 1]) {
    
    #Mean of malformation rates in all UV groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 2] <- mean(tempMat[
      which((tempMat[, 7] == 'UV' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 14],
      na.rm = TRUE)
    
    #Mean of malformation rates in all Control groups, and sub grouped by 
    #plants
    cuttingsPlantsMat[loopI, 3] <- mean(tempMat[
      which((tempMat[, 7] == 'CON' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 14],
      na.rm = TRUE)
    
    #Mean of GCMS signals sub grouped by plants
    cuttingsPlantsMat[loopI, 4] <- mean(tempMat[
      which((as.numeric(tempMat[, 8]) == loopI)), 10],
      na.rm = TRUE)
    
    #Mean of PLS predicted values sub grouped by plants
    cuttingsPlantsMat[loopI, 5] <- mean(tempMat[
      which((as.numeric(tempMat[, 8]) == loopI)), 9],
      na.rm = TRUE)
    
    #Mean of GCMS signals in all UV groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 6] <- mean(tempMat[
      which(((tempMat[, 7] == 'UV' ) & 
               as.numeric(tempMat[, 8]) == loopI)), 10],
      na.rm = TRUE)
    
    #Mean of GCMS signals in all Control groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 7] <- mean(tempMat[
      which((tempMat[, 7] == 'CON' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 10],
      na.rm = TRUE)
    
    #Mean of PLS predicted values in all UV groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 8] <- mean(tempMat[
      which((tempMat[, 7] == 'UV' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 9],
      na.rm = TRUE)
    
    #Mean of PLS predicted values in all Control groups, and sub grouped by 
    #plants
    cuttingsPlantsMat[loopI, 9] <- mean(tempMat[
      which((tempMat[, 7] == 'CON' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 9],
      na.rm = TRUE)
    
    #Sd of malformation rates in all UV groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 10] <- sd(tempMat[
      which((tempMat[, 7] == 'UV' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 14],
      na.rm = TRUE)
    
    #Sd of malformation rates in all Control groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 11] <- sd(tempMat[
      which((tempMat[, 7] == 'CON' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 14],
      na.rm = TRUE)
    
    #Sd of GCMS signals sub grouped by plants
    cuttingsPlantsMat[loopI, 12] <- sd(tempMat[
      which((as.numeric(tempMat[, 8]) == loopI)), 10],
      na.rm = TRUE)
    
    #Sd of PLS predicted values sub grouped by plants
    cuttingsPlantsMat[loopI, 13] <- sd(tempMat[
      which((as.numeric(tempMat[, 8]) == loopI)), 9],
      na.rm = TRUE)
    
    #Sd of GCMS signals in all UV groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 14] <- sd(tempMat[
      which(((tempMat[, 7] == 'UV' ) & 
               as.numeric(tempMat[, 8]) == loopI)), 10],
      na.rm = TRUE)
    
    #Sd of GCMS signals in all Control groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 15] <- sd(tempMat[
      which((tempMat[, 7] == 'CON' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 10],
      na.rm = TRUE)
    
    #Sd of PLS predicted values in all UV groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 16] <- sd(tempMat[
      which((tempMat[, 7] == 'UV' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 9],
      na.rm = TRUE)
    
    #Sd of PLS predicted values in all Control groups, and sub grouped by plants
    cuttingsPlantsMat[loopI, 17] <- sd(tempMat[
      which((tempMat[, 7] == 'CON' ) & 
              (as.numeric(tempMat[, 8]) == loopI)), 9],
      na.rm = TRUE)
  }
  
  #Transfer matrix into dataframe
  cuttingsPlantsMat <- as.data.frame(cuttingsPlantsMat)
  cuttingsPlantsMat[, 1] <- as.character(cuttingsPlantsMat[, 1])
  
  return(list(pcaMatFTIRGCMalFull,
              pcaMatFTIRGCMalFullMean,
              cuttingsPlantsMat))
}