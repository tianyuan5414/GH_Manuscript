#Load necessary packages--------------------------------------------------------
{
  #Find peaks and extract peak statistics
  library(photobiology)
  
  #Paralell processing
  library(doSNOW)
}

#Declare necessary functions----------------------------------------------------
{
  #Benjamin A Bell github codes, for detecting and calculating peak area--------
  {
    # Load peak detection script
    source("https://raw.githubusercontent.com/benbell95/peak-detection/main/r/peak_detection.r")
    # Load peak area script
    source("https://raw.githubusercontent.com/benbell95/peak-detection/main/r/peak_area.r")
  }
  
  #Create a function to automatically identify all qualified peaks in given
  #spectrum, and extract the peak height and peak area under each peaks
  #tempSpec - one single spectra in numeric vector, with wavenumber information
  #saved in elements' names
  #tempMinPeakWid - one numeric value indicating the minimum peak width 
  #to be identified, should be bigger than the wavenumber interval in the 
  #spectrum input
  #tempMinPeakWidMax - one numeric value indicating the maximum peak width 
  #to be identified, should be bigger than tempMinPeakWid
  #plot - Boolean value, indicating whether a plot on current spectrum labelled
  #with identified peaks is output
  #tempTol - one numeric value defining the minimum peak singal to be 
  #identified
  autoPeakDecQuan <- function(
    tempSpec, 
    tempMinPeakWid, 
    tempMinPeakWidMax,
    plot = FALSE, 
    tempTol = 0
  ) {
    #Extract the spectra resolution
    tempResu <- floor(tempMinPeakWid / 
                        abs(as.numeric(names(tempSpec))[2] - as.numeric(names(tempSpec))[1])
    )
    
    tempResu2 <- floor(tempMinPeakWidMax / 
                         abs(as.numeric(names(tempSpec))[2] - as.numeric(names(tempSpec))[1])
    )
    
    #Get peak positions (central position) for all peaks identified in current
    #spectrum
    tempPeakPos <- as.numeric(names(tempSpec))[
      which(find_peaks(x = as.numeric(tempSpec),
                       local.threshold = tempTol,
                       span = tempResu,
                       strict = FALSE)
      )
    ]
    
    #Find the peak ranges, acting on the peak position identified before
    #Create a matrix storing the peak range information
    tempPeakRange <- matrix(nrow = length(tempPeakPos),
                            ncol = 5,
                            dimnames = list(c(),
                                            c('Peak_Pos',
                                              'Peak_L',
                                              'Peak_R',
                                              'Peak_Height',
                                              'Peak_Area')
                            )
    )
    
    #Fill in the matrix
    tempPeakRange[ ,1] <- tempPeakPos
    
    #Loop to go through all identified peaks - left
    for (loopI in 1:length(tempPeakPos)) {
      
      tempLeftSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                      tempPeakPos[loopI]
      )
      ]
      
      #Assign an indicator to terminate the loop for range searching
      tempIndcator <- 1
      
      tempIndLeftIndcator <- 1
      
      while (tempIndcator == 1) {
        
        if ((which(as.numeric(names(tempSpec)) == 
                   tempPeakPos[loopI]) - tempIndLeftIndcator < 1) | (tempIndLeftIndcator > tempResu2)) {
          
          #Assign the terminate indicator to 0, for preparing to stop the loop
          tempIndcator <- 0
          
          break
          
        } else {
          
          if (tempLeftSig > tempSpec[which(as.numeric(names(tempSpec)) == 
                                           tempPeakPos[loopI]
          ) - tempIndLeftIndcator]) {
            
            tempLeftSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                            tempPeakPos[loopI]
            ) - tempIndLeftIndcator]
            
            tempIndLeftIndcator <- tempIndLeftIndcator + 1
            
          }else{
            
            #Assign the terminate indicator to 0, for preparing to stop the loop
            tempIndcator <- 0
            
            #This indicate that the search might approach its local peak
            #boundary, then apply a boundary checking referring to the 
            #boundary tolerance index defined by the spectra resolution
            for (loopJ in 1:ceiling(tempResu / 4)) {
              
              if ((which(as.numeric(names(tempSpec)) == 
                         tempPeakPos[loopI]
              ) - tempIndLeftIndcator - loopJ) < 1) {
                
                tempRightSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                                 tempPeakPos[loopI]
                ) + tempIndLeftIndcator - loopJ + 1]
                
                tempIndcator <- 0
                
                tempIndLeftIndcator <- tempIndLeftIndcator + loopJ
                
                break
                
              }else {
                
                #Search for the local minimum for left section of the current peak
                if (tempLeftSig > tempSpec[which(as.numeric(names(tempSpec)) == 
                                                 tempPeakPos[loopI]
                ) - tempIndLeftIndcator - loopJ]) {
                  
                  tempLeftSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                                  tempPeakPos[loopI]
                  ) - tempIndLeftIndcator - loopJ]
                  
                  #Reassign the terminate indicator to 1, to continue the loop
                  tempIndcator <- 1
                  
                  tempIndLeftIndcator <- tempIndLeftIndcator + loopJ
                  
                  #Stop the loop for boundary detection
                  break
                }
              }
            }
          }
        }
        #Record the peak range inot matrix
        tempPeakRange[loopI, 2] <- names(tempSpec)[which(as.numeric(names(tempSpec)) == 
                                                           tempPeakPos[loopI]) - 
                                                     tempIndLeftIndcator + 1
        ]
      }
    }
    
    #Loop to go through all identified peaks - right
    for (loopI in 1:length(tempPeakPos)) {
      
      tempRightSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                       tempPeakPos[loopI]
      )
      ]
      
      #Assign an indicator to terminate the loop for range searching
      tempIndcator <- 1
      
      tempIndRightIndcator <- 1
      
      while (tempIndcator == 1) {
        
        if (which(as.numeric(names(tempSpec)) == 
                  tempPeakPos[loopI]) + tempIndRightIndcator > ncol(tempSpec) | (tempIndRightIndcator > tempResu2)) {
          
          #Assign the terminate indicator to 0, for preparing to stop the loop
          tempIndcator <- 0
          
          break
          
        } else {
          
          if (tempRightSig > tempSpec[which(as.numeric(names(tempSpec)) == 
                                            tempPeakPos[loopI]
          ) + tempIndRightIndcator]) {
            
            tempRightSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                             tempPeakPos[loopI]
            ) + tempIndRightIndcator]
            
            tempIndRightIndcator <- tempIndRightIndcator + 1
            
          }else{
            
            #Assign the terminate indicator to 0, for preparing to stop the loop
            tempIndcator <- 0
            
            #This indicate that the search might approach its local peak
            #boundary, then apply a boundary checking referring to the 
            #boundary tolerance index defined by the spectra resolution
            for (loopJ in 1:ceiling(tempResu / 4)) {
              
              if ((which(as.numeric(names(tempSpec)) == 
                         tempPeakPos[loopI]
              ) + tempIndRightIndcator + loopJ) > ncol(tempSpec)) {
                
                tempRightSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                                 tempPeakPos[loopI]
                ) + tempIndRightIndcator + loopJ - 1]
                
                tempIndcator <- 0
                
                tempIndRightIndcator <- tempIndRightIndcator + loopJ
                
                break
                
              }else {
                
                #Search for the local minimum for left section of the current peak
                if (tempRightSig > tempSpec[which(as.numeric(names(tempSpec)) == 
                                                  tempPeakPos[loopI]
                ) + tempIndRightIndcator + loopJ]) {
                  
                  tempRightSig <- tempSpec[which(as.numeric(names(tempSpec)) == 
                                                   tempPeakPos[loopI]
                  ) + tempIndRightIndcator + loopJ]
                  
                  #Reassign the terminate indicator to 1, to continue the loop
                  tempIndcator <- 1
                  
                  tempIndRightIndcator <- tempIndRightIndcator + loopJ
                  
                  #Stop the loop for boundary detection
                  break
                }
                
              }
            }
          }
        }
      }
      
      #Record the peak range inot matrix
      tempPeakRange[loopI, 3] <- names(tempSpec)[which(as.numeric(names(tempSpec)) == 
                                                         tempPeakPos[loopI]) + 
                                                   tempIndRightIndcator - 1
      ]
      
      tempPeakRange[loopI, 4] <- as.numeric(tempSpec)[which(as.numeric(names(tempSpec)) ==
                                                              tempPeakPos[loopI])
      ]
      
      tempPeakRange[loopI, 5] <- unlist(spd_area_bands(t(rev(rbind(tempSpec, rep(0, length(tempSpec))))),
                                                       bands = tempPeakRange[loopI, c(3,2)],
                                                       refine = FALSE)$area)[1]
    }
    
    
    #Plot the current spectrum with information extracted, if corresponding
    #variable is set as true
    if (plot == TRUE) {
      #View an example of specta read
      plot(x = as.numeric(names(tempSpec)),
           y = tempSpec,
           type = 'l',
           xlim = rev(range(as.numeric(names(tempSpec)
           )
           )
           )
      )
      
      #Add reference lines around peaks identified on the sample spectra
      abline(v=tempPeakPos, col="blue")
    }
    
    #Output the matrix
    return(tempPeakRange)
  }
  
  #Declare a function to combine all matrix in a list into one matrix
  #tempList - a list to be combined, each element in this list should be a 
  #matrix, all matrix in different elements should have the same column 
  #title, this list must have length bigger than 1
  combineList <- function(tempList) {
    
    #Create a temporary matrix, and put the fisrt list element into the matrix
    tempMat <- tempList[[1]]
    
    #Loop to combind all other list elements into the matrix above
    for (loopI in 2:length(tempList)) {
      tempMat <- rbind(tempMat, tempList[[loopI]])
    }
    
    #Output the combined matrix
    return(tempMat)
  }
  
  #Declare a function to combine some of rows and columns in all matrix in a 
  #list into one matrix
  #tempList - a list to be combined, each element in this list should be a 
  #matrix, all matrix in different elements should have the same column 
  #title, this list must have length bigger than 1
  #indexRow - a numeric vector indicating the indexes of all rows in each
  #matrix which will be combined, should have length bigger than 1
  #indexCol - optional, a numeric vector indicating the indexes of all columns
  #in each matrix which will be combined. Will be set as all column indexes
  #in each matrix if not specified
  combineListEle <- function(tempList, indexRow, indexCol) {
    
    #Check whether variable indexCol is input. Set it to all column indexes in  
    #single matrix if not specified
    if (exists('indexcol')) {
      indexcol <- indexCol
    }else{
      indexcol <- 1:ncol(tempList[[1]])
    }
    
    #Create a temporary matrix and copy the first matrix element in the list 
    #into it
    tempMat <- tempList[[1]][indexRow, indexCol]
    
    #Loop to combine all specified rows and columns in each matrix element in 
    #the list together
    for (loopI in 2:length(tempList)) {
      tempMat <- rbind(tempMat, tempList[[loopI]][indexRow, indexCol])
    }
    
    #Output the combined matrix
    return(tempMat)
  }
  
  #Create a function to find peaks with  similar positions among different peak
  #information matrix, and modified all original peak information matrix into
  #another one that deletes the singular peak positions (the peaks that is 
  #not identified in all of other spectra)
  #peakInfoList - the original list of peak information, each element in this
  #list is a matrix containing all peak statistics identified from function
  #'autoPeakDecQuan', should have length bigger than 1
  #tempPeakPos - numeric vector, indicating the peak position to be remained,
  #should have length bigger than 1
  #temppeakTol -  numeric vector, indicating the tolerance (drift 
  #distance between the theoretical peak position specified in variable
  #'tempPeakPos', peaks farer than this distance will not be considered as
  #the same peak and will be deleted at last. Should have the same length as
  #variable 'tempPeakPos'
  autoPeakResam <- function(peakInfoList, tempPeakPos, temppeakTol) {
    
    #Combine the input peak information list into one matrix
    tempMat <- combineList(peakInfoList)
    
    #Create a matrix to resample from the peak information matrix
    tempPeakInfoMat <- matrix(nrow = length(tempPeakPos),
                              ncol = 2 + ncol(tempMat),
                              dimnames = list(c(),
                                              c('Pre_Peak_Pos',
                                                'Peak_Pos_Tol',
                                                colnames(tempMat)
                                              )
                              )
    )
    
    #Fill in the matrix
    tempPeakInfoMat[, 1] <- tempPeakPos
    tempPeakInfoMat[, 2] <- temppeakTol
    
    #Create a list to store the resampled peak information
    tempList <- vector('list', length = length(peakInfoList))
    
    #Loop around all spectra
    for (loopI in 1:length(tempList)) {
      
      tempList[[loopI]] <- tempPeakInfoMat
      
      #Loop around possible peak positions to find the nearest peak identified in the
      #Spectra
      for (loopJ in 1:nrow(tempPeakInfoMat)) {
        
        #Check if there is at least one peak matched with any of peaks listed in
        #'tempPeakPos'
        if (length(which(
          (abs(
            as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
          ) == 
          min(
            abs(
              as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
            ))) & (abs(
              as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
            ) <= temppeakTol[loopJ]
            ))) >= 1){
          
          #Check if there is more than one peak matched with the same peak 
          #listed in 'tempPeakPos'
          if (length(which(
            (abs(
              as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
            ) == 
            min(
              abs(
                as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
              ))) & (abs(
                as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
              ) <= temppeakTol[loopJ]
              ))) > 1) {
            
            #if there is more than one peaks matching the same peak in 
            #'tempPeakPos'
            tempList[[loopI]][loopJ, 3:ncol(tempList[[loopI]])] <-
              apply(peakInfoSpec[[loopI]][which(
                (
                  abs(
                    as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
                  ) == 
                    min(
                      abs(
                        as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
                      )
                    )
                ) & (
                  abs(
                    as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
                  ) <= temppeakTol[loopJ]
                )
              ), ], 2, mean)
          }else {
            
            #If there is only one peak matching a peak in 'tempPeakPos'
            tempList[[loopI]][loopJ, 3:ncol(tempList[[loopI]])] <-
              peakInfoSpec[[loopI]][which(
                (
                  abs(
                    as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
                  ) == 
                    min(
                      abs(
                        as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
                      )
                    )
                ) & (
                  abs(
                    as.numeric(peakInfoSpec[[loopI]][, 1]) - tempPeakPos[loopJ]
                  ) <= temppeakTol[loopJ]
                )
              ), ]
          }
        }
      }
      
      #Adjust the value format into numeric
      tempList[[loopI]] <- apply(tempList[[loopI]], c(1,2), as.numeric)
    }
    
    #Output the modified list that deletes peaks that do not show in all spectra
    return(tempList)
  }
  
  #Check NAs in a vector, return true if the correconding element is NA
  #tempVec - numeric vector to be treated for deleting NAs, should have length
  #bigger than 1
  cleanNA <- function(tempVec) {
    
    for (loopI in 1:length(tempVec)) {
      
      if (is.na(tempVec[loopI])) {
        
        return(FALSE)
        
        break
        
      }else if(loopI == length(tempVec)) {
        
        return(TRUE)
        
      }
    }
  }
  
  #Detect the symbol of given value, and output 1 or -1 accordingly
  #tempVec - single numeric value to be treated
  symbolDec <- function(tempVec){
    
    tempVec2 <- tempVec
    
    for (loopI in 1:length(tempVec)) {
      if (tempVec[loopI] <= 0) {
        
        tempVec2[loopI] <- -1
        
      }else{
        
        tempVec2[loopI] <- 1
        
      }
    }
    
    return(tempVec2)
    
  }
} 

#Read FTIR data and peak signal extration---------------------------------------
{
  #Read raw spectra (smoothed and baseline corrected, z-score standardisation)
  {
    #Read standardised FTIR spectra
    bcSpec <- read.csv(file.path( 'Data_R', 'bcNorSpec.csv'),
                       row.names = 1,
                       header = TRUE)
    
    #Remove spectra of roof samples
    bcSpec <- bcSpec[which(bcSpec[, 3] != 'ROOF'),]
    
    #Rename the column names for the smoothed spectra
    colnames(bcSpec)[11:ncol(bcSpec)] <- seq(3976, 620, by = -4)
    
    #Convert spectra into numeric values
    bcSpec[, 11:ncol(bcSpec)] <- apply(bcSpec[, 11:ncol(bcSpec)], 
                                       c(1,2),
                                       as.numeric)
  }
  
  #Extract peak information (position and range) from raw spectra
  {
    #Create a list to store all peak information for all spectra
    peakInfoSpec <- vector('list', nrow(bcSpec))
    
    #Setup cores for parallel processing
    cl <- makeCluster(8, type="SOCK")
    
    #Regist multiple-core setup
    registerDoSNOW(cl)
    
    #Apply peak identification function on all spectra, and output the peak
    #information matrix into a list
    peakInfoSpec <- foreach(loopI = 1:length(peakInfoSpec),
                            .packages = c('photobiology')) %dopar% 
      {
        autoPeakDecQuan(bcSpec[loopI, -1:-10], 20, 100, FALSE
        )
      }
    
    #Quite the setup of multiple-core processing
    stopCluster(cl)
    
    #Clean up some of variables
    rm(cl)
    
    #Assign possible peak positions together with tolerance range for peak drifts
    peakPos <- c(3324, 2924, 2856, 1740, 
                 1636, 1556, 1516, 1440, 1376,
                 1324, 1268, 1036, 996, 832)
    peakTol <- c(200, rep(50, 13))
    
    #Rearrange the peak information into an uniform sequence
    peakInfoSpecUni <- autoPeakResam(peakInfoSpec, peakPos, peakTol)
  }
  
  #Rearrange all peak information into one matrix
  {
    #Create a matirx for all peak information - Height
    pcaMatHeight <- matrix(nrow = nrow(bcSpec),
                           ncol = 10 + length(peakPos),
                           dimnames = list(c(),
                                           c(colnames(bcSpec)[1:10], 
                                             peakPos)))
    
    #Fill the batch information into matrix for peak height and area
    pcaMatHeight[, 1:10] <- as.matrix(bcSpec[, 1:10])
    
    #Loop to fill in the matrix with individual spectrum on corresponding row
    for (loopI in 1:nrow(bcSpec)) {
      
      pcaMatHeight[loopI, c(11:ncol(pcaMatHeight))] <- peakInfoSpecUni[[loopI]][, 6]
   
    }
  }
}

#Output peak height information for all FTIR spectra
{
 #Peak height matrix
  write.csv(pcaMatHeight, file = file.path('Output_R', 'pcaMatHeight.csv'))
}