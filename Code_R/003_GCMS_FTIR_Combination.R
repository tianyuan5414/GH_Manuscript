#Read GC-MS signals and combine them with FTIR peak signals---------------------
{
  #Read FTIR peak height data
  {
    pcaMatHeight <- read.csv(file = file.path('Output_R', 'pcaMatHeight.csv'),
                             row.names = 1)
  }
  
  #Read GC-MS data
  {
    tempDir2 <- paste(getwd(), '/Data_R/01_07_2025_gc_cultivar_data_2025.csv', 
                      sep = '', 
                      collapse = NULL)
    
    #GC-MS data on cuttings and pot cultivars
    GCDataRawCut <- as.matrix(read.csv(
      file.path( 'Data_R', '04_11_2025_gc_cutting_data_2025.csv'),
      header = TRUE,
      row.names = 1))
    GCDataRawCul <- as.matrix(read.csv(
      file.path( 'Data_R', '01_07_2025_gc_cultivar_data_2025.csv'),
      header = TRUE,
      row.names = 1))
    
    #Combine these two raw datasets, with batch information
    GCData <- rbind(GCDataRawCut[, c(8, 10:14)],
                    cbind(GCDataRawCul[, c(8, 10:12)], 
                          sampling_number = NA, 
                          GCDataRawCul[, 13])
    )
    
    #Exclude the last batch (u), as an unexpected crash happened before this
    #batch, making reseult from this batch abnormal
    GCData <- GCData[-which(GCData[, 6] == 'u'),]
    
    #Loop to calculate the mean values for each sample
    uniqueBatchInfo <- unique(GCData[, 2:5])
    
    #Create a matrix to fill in the mean GC-MS signals
    GCDataMean <- matrix(nrow = nrow(uniqueBatchInfo),
                         ncol = ncol(uniqueBatchInfo) + 2,
                         dimnames = list(c(),
                                         c('calibrated',
                                           'sd',
                                           colnames(uniqueBatchInfo)
                                         )
                         )
    )
    
    #Add batch number to the GCMS data matrix
    GCDataMean[, 3:ncol(GCDataMean)] <- uniqueBatchInfo
    
    #Loop to calculate the mean signals of GCMS duplicates
    for (loopI in 1:nrow(GCDataMean)) {
      
      #Check wheter current sample is cutting or pot cultivar 
      if (is.na(GCDataMean[loopI, 6])) {
        
        #Mean and sd for pot cultivar
        GCDataMean[loopI, 1] <- mean(as.numeric(
          GCData[which((GCData[, 2] == GCDataMean[loopI, 3]) &
                         (GCData[, 3] == GCDataMean[loopI, 4]) &
                         (GCData[, 4] == GCDataMean[loopI, 5]) 
          )
          , 1])
        )
        GCDataMean[loopI, 2] <- sd(as.numeric(
          GCData[which((GCData[, 2] == GCDataMean[loopI, 3]) &
                         (GCData[, 3] == GCDataMean[loopI, 4]) &
                         (GCData[, 4] == GCDataMean[loopI, 5])
          )
          , 1])
        )
        
      }else {
        
        #Mean and sd for cutting
        GCDataMean[loopI, 1] <- mean(as.numeric(
          GCData[which((GCData[, 2] == GCDataMean[loopI, 3]) &
                         (GCData[, 3] == GCDataMean[loopI, 4]) &
                         (GCData[, 4] == GCDataMean[loopI, 5]) &
                         (GCData[, 5] == GCDataMean[loopI, 6])
          )
          , 1])
        )
        GCDataMean[loopI, 2] <- sd(as.numeric(
          GCData[which((GCData[, 2] == GCDataMean[loopI, 3]) &
                         (GCData[, 3] == GCDataMean[loopI, 4]) &
                         (GCData[, 4] == GCDataMean[loopI, 5]) &
                         (GCData[, 5] == GCDataMean[loopI, 6])
          )
          , 1])
        )
        
      }
    }
    
    #Transfer the matrix back to data frame, and transfer values to number if
    #necessary
    GCDataMean <- as.data.frame(GCDataMean)
    GCDataMean[, 1] <- as.numeric(GCDataMean[, 1])
  }
  
  #Combine GC-MS data with FTIR data
  {
    #Create a vector to stroe the index for batch information matching
    tempInd <- vector('numeric', length = nrow(pcaMatHeight))
    
    #Loop to match the batch information between GC-MS and FTIR datasets
    for (loopI in 1:nrow(pcaMatHeight)) {
      
      #Check whether current sample is pot cultivar
      if (pcaMatHeight[loopI, 5] == 'cul') {
        
        #Check whether current FTIR spectrum has any GCMS signal
        if (any((GCDataMean[, 3] == pcaMatHeight[loopI, 6]) &
                (GCDataMean[, 4] == pcaMatHeight[loopI, 7]) &
                (GCDataMean[, 5] == pcaMatHeight[loopI, 8]) &
                (is.na(GCDataMean[, 6]))
        )
        ) {
          
          #Record the FTIR matrix row index for current GCMS signal
          tempInd[loopI] <- which((GCDataMean[, 3] == pcaMatHeight[loopI, 6]) &
                                    (GCDataMean[, 4] == pcaMatHeight[loopI, 7]) &
                                    (GCDataMean[, 5] == pcaMatHeight[loopI, 8]) &
                                    (is.na(GCDataMean[, 6])))
          
        }else {
          
          #Record nothing if there is no match between current FTIR specrum and
          #any GCMS signal
          tempInd[loopI] <- NA
          
        }
        
        #Check whether current sample is cutting
      } else if (pcaMatHeight[loopI, 5] == 'cut') {
        
        #Check whether current FTIR spectrum has any GCMS signal
        if (any((GCDataMean[, 3] == pcaMatHeight[loopI, 6]) &
                (GCDataMean[, 4] == pcaMatHeight[loopI, 7]) &
                (GCDataMean[, 5] == pcaMatHeight[loopI, 8]) &
                (GCDataMean[, 6] == pcaMatHeight[loopI, 9])
        )
        ) {
          
          #Record the FTIR matrix row index for current GCMS signal
          tempInd[loopI] <- which((GCDataMean[, 3] == pcaMatHeight[loopI, 6]) &
                                    (GCDataMean[, 4] == pcaMatHeight[loopI, 7]) &
                                    (GCDataMean[, 5] == pcaMatHeight[loopI, 8]) &
                                    (GCDataMean[, 6] == pcaMatHeight[loopI, 9]))
          
        }else {
          
          #Record nothing if there is no match between current FTIR spectrum and
          #any GCMS signal
          tempInd[loopI] <- NA
          
        }
        
      } else {
        
        #Record nothing if there is no match between current FTIR spectrum and
        #any GCMS signal
        tempInd[loopI] <- NA
        
      }
    }
    
    #Attach GCMS signals to  peak statistics (height or area)
    pcaMatHeightCom <- cbind(pcaMatHeight,
                          GCDataMean[tempInd, c(1,2)])

    rm(loopI,
       tempInd)
  }
}

#Output combined dataset with both FTIR peaks and GCMS signals------------------
{
  #Mean of GCMS signals for each sample
  write.csv(GCDataMean, file = file.path('Output_R', 'GCDataMean.csv'))
  
  #Combined dataset
  write.csv(pcaMatHeightCom, 
            file = file.path('Output_R', 'pcaMatHeightCom.csv'))
}