#Read malformation data and combine them with FTIR/GCMS data--------------------
{
  #Read malformation data
  {
    #Read the cutting branch data, and replace all na with 0
    cuttingsRaw <- read.csv(file = 
                              file.path( 
                                'Data_R', 
                                'Malformations_Cuttings.csv'),
                            row.names = 1)
    
    #Replace all NAs with 0
    cuttingsRaw[is.na(cuttingsRaw)] <- 0
    
    #Extract the batch information
    batchInfoTemp <- strsplit(x = cuttingsRaw[,1], split = '_')
    cuttingsBatchInfo <- cbind(rep('Cuttings', length(batchInfoTemp)),
                               rep('2025', length(batchInfoTemp)),
                               matrix(unlist(batchInfoTemp), ncol = 4, byrow = TRUE))
    colnames(cuttingsBatchInfo) <- c('Sample', 'Year', 'Frame', 'Treatment',
                                     'Plant', 'Duplicate')
    
    #Delete unused variables
    rm(batchInfoTemp)
    
    #Read the cultivar branch data, and replace all na with 0
    cultivarsRaw <- read.csv(file = 
                               file.path( 
                                 'Data_R', 
                                 'Malformations_Cultivars.csv'),
                             row.names = 1)
    
    #Replace all NAs with 0
    cultivarsRaw[is.na(cultivarsRaw)] <- 0
    
    #Extract the batch information
    batchInfoTemp <- strsplit(x = cultivarsRaw[,1], split = '_')
    cultivarsBatchInfo <- cbind(rep('Cultivars', length(batchInfoTemp)),
                                matrix(unlist(batchInfoTemp), ncol = 4, byrow = TRUE))
    colnames(cultivarsBatchInfo) <- c('Sample', 'Year', 'Frame', 'Treatment',
                                      'Plant')
    
    #Delete unused variables
    rm(batchInfoTemp)
    
    #Combine counting results in duplicated samples, by summing them up
    tempVec <- unique(cultivarsRaw[, 1])
    
    #Matrix for all malformation data
    tempMat <- matrix(nrow = length(tempVec),
                      ncol = ncol(cultivarsRaw),
                      dimnames = list(c(),
                                      colnames(cultivarsRaw)))
    
    #Matrix for all batch information on malformation data
    tempMat2 <- matrix(nrow = length(tempVec),
                       ncol = ncol(cultivarsBatchInfo),
                       dimnames = list(c(),
                                       colnames(cultivarsBatchInfo)))
    
    #Fill sample label into matrix
    tempMat[, 1] <- tempVec
    
    #Loop to add duplicates together
    for (loopI in 1:nrow(tempMat)) {
      
      #Check  if there are duplicates in current sample
      if (length(which(cultivarsRaw[, 1] == tempMat[loopI, 1])) > 1) {
        
        #Sum up all counting results among duplicates
        tempMat[loopI, 2:ncol(tempMat)] <- apply(cultivarsRaw[
          which(cultivarsRaw[, 1] == tempMat[loopI, 1]), 
          2:ncol(tempMat)
        ], 2, sum)
        
      }else{
        
        #Transfer couting results into matrix
        tempMat[loopI, 2:ncol(tempMat)] <- cultivarsRaw[
          which(cultivarsRaw[, 1] == tempMat[loopI, 1]), 2:ncol(tempMat)
        ]
        
      }
      
      #Record the corresponding batch information
      tempMat2[loopI,] <- cultivarsBatchInfo[which(cultivarsRaw[, 1] == tempMat[loopI, 1])[1],]
      
    }
    
    #Combine batch information with counting results
    cultivarsCom <- cbind(tempMat2, tempMat[, -1])
    
    rm(loopI,
       tempMat,
       tempMat2,
       tempVec)
  }
  
  #Combine malformation data with FTIR and GC-MS data
  {
    #Read GCMS signals and FTIR peak data
    pcaMatFTIR <- read.csv(file = file.path('Output_R', 'PLS_Nor.csv'),
                           row.names = 1)[,c(1:10, 27,25:26)]
    
    colnames(pcaMatFTIR)[-1:-10] <- c('PLS', 'GCMS', 'GCMS_SD')
    
    #Vector for recording row index in FTIR spectra matrix
    tempInd <- vector(length = nrow(pcaMatFTIR))
    
    #Combine batch information of cuttings and pot cultivars together
    tempMat <- rbind(cultivarsCom[,1:6],
                     cuttingsBatchInfo)
    
    #Add malformation data to FTIR data
    for (loopI in 1:nrow(pcaMatFTIR)) {
      
      #Check and skip the roof samples
      if (pcaMatFTIR[loopI, 3] == 'ROOF') {
        
        tempInd[loopI] <- NA
        
      }else{
        
        #For all pot cultivars
        if (pcaMatFTIR[loopI, 5] == 'cul') {
          
          if (any((tempMat[, 3] == pcaMatFTIR[loopI, 6]) &
                  (tempMat[, 4] == pcaMatFTIR[loopI, 7]) &
                  (tempMat[, 5] == pcaMatFTIR[loopI, 8]) &
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
                     as.matrix(cuttingsRaw[, -1]))
    
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
    tempMat <- unique(pcaMatFTIRGCMalFull[, 1:8])
    
    #Matrix for recording mean and sd of FTIR duplicates
    tempMat2 <- matrix(nrow = nrow(tempMat),
                       ncol = ncol(pcaMatFTIRGCMalFull) - 1,
                       dimnames = list(c(),
                                       c(colnames(pcaMatFTIRGCMalFull)[-9:-10], 
                                         paste(colnames(pcaMatFTIRGCMalFull)[11], '_SD')
                                       )
                       )
    )
    
    tempMat2[,1:8] <- as.matrix(tempMat)
    
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
       loopI)
    
    #Create a matrix recording all significance results between different 
      #subgroups
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
    
    #Calculate malformation rates in pot cultivars
    tempVec <- (cultivarsRaw[, 15] - cultivarsRaw[, 2])/ cultivarsRaw[, 2]
    cultivarMean <- vector(length = 4)
    cultivarMean[1] <- mean(tempVec[which(cultivarsBatchInfo[, 4] == 'UV')])
    cultivarMean[2] <- sd(tempVec[which(cultivarsBatchInfo[, 4] == 'UV')])
    cultivarMean[3] <- mean(tempVec[which(cultivarsBatchInfo[, 4] == 'CON')])
    cultivarMean[4] <- sd(tempVec[which(cultivarsBatchInfo[, 4] == 'CON')])
    
    rm(loopI,
       tempMat,
       tempVec) 
  }
}

#Output combined data
{
  #Combined matrix of Malformation, GCMS and FTIR
  write.csv(pcaMatFTIRGCMalFull, file = 
              file.path('Output_R', 'pcaMatFTIRGCMalFull.csv'))
  
  #Combined matrix of Malformation, GCMS and FTIR - mean value in each sample
  write.csv(pcaMatFTIRGCMalFullMean, file = 
              file.path('Output_R', 'pcaMatFTIRGCMalFullMean.csv'))
  
  #Combined matrix of Malformation, GCMS and FTIR - mean value in subgroup of
    #donor plants
  write.csv(cuttingsPlantsMat, file = 
              file.path('Output_R', 'cuttingsPlantsMat.csv'))
  
  #Cultivar malformation rate - mean and sd
  write.csv(cultivarMean, file = 
              file.path('Output_R', 'cultivarMean.csv'))
}