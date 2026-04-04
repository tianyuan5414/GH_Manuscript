#Load necessary packages
{
  #Read generic files
  library(here)
  
  #Read OPUS files
  library(opusreader)
  
  #Smoothing
  library(prospectr)
  
  #Baseline correction
  library(baseline)
}

#Declare necessary functions
{
  #Declare a function to automatically split the name of spectral file, 
  #to get all necessary
  #folderNames - single character string, the original file name to be split
  #to extract the batch information
  #seperator - single character string, the symbol used in the original file 
  #name for separating different types of batch information
  #defaultSettings - vector with length of 10, indicating the default batch
  #information
  namesSplit <- function(folderNames, seperator, defaultSettings) {
    
    splitedNames <- defaultSettings
    
    #Split the names
    splitedNamesRaw <- strsplit( x = folderNames, split = seperator)[[1]]
    
    #Check the information given in the folder name, 
    #and change the corresponding sections in output if necessary
    if(splitedNamesRaw[2] != 'GH'){
      
      splitedNames[c(1:8, 10)] <- splitedNamesRaw[c(1:9)]
      
    } else if(splitedNamesRaw[4] == 'cul'){
      
      splitedNames[c(1,3:8, 10)] <- splitedNamesRaw[c(1:8)]
      
    }else {
      
      splitedNames[c(1,3:10)] <- splitedNamesRaw[c(1:9)]
      
    }
    
    if(splitedNamesRaw[5]  == 'ROOF' | splitedNamesRaw[6] == 'ROOF') {
      
      splitedNames[c(1,3:5, 8,10)]  <- splitedNamesRaw[c(1,5,3,4,6,7)]
      
      splitedNames[c(6,7,9)]  <- NA
      
    }
    
    return(splitedNames)
  }
  
  #Declare a function to calculate the number of all elements in a list given
  #tempList - list to be used for calculation
  numOfEle <- function(tempList) {
    
    #Vector storing the number
    tempVec <- 0
    
    for (loopI in 1:length(tempList)) {
      tempVec <- tempVec + length(tempList[[loopI]])
    }
    
    return(tempVec)
  }
  
  #Create a function to re-sample data according to given sequence from a vector
  #tempVec - numeric vector, the FTIR spectrum signals at original wavenumbers
  #inputSeq - numeric vector, the sequence of wavenumber for input FTIR 
  #spectrum
  #outputSeq - numeric vector, the sequence of wavenumber to be output
  resamData <- function(tempVec, inputSeq, outputSeq){
    tempOutVec <- approx(x = inputSeq,
                         y = tempVec,
                         xout = outputSeq)[[2]]
    
    return(tempOutVec)
  }
  
  #Function to normalise the spectral absorbance to maximum intensity
  #x - numeric vector, the FTIR spectrum to be normalised
  zstandMax <- function(x) {
    return((x - mean(x))/ sd(x))
  }
}

#Read all OPUS files
{
  #List all the folder for different measuring days
  daysVec <- list.files(path = file.path('Data_R', 'FTIR_Raw_Spectra'))
  
  #Loop to list all .o files in each folder above
  #List for storing all names of spectra files
  namesList <- list()
  for (loopI in 1:length(daysVec)) {
    
    namesList[[length(namesList) + 1]] <- 
      list.files(path = file.path('Data_R', 'FTIR_Raw_Spectra',
                                    daysVec[loopI]), 
                 pattern = '\\.0$')
    
  }
  
  #Loop to read all spectra
  #List for storing all spectra
  rawSpecList <- list()
  
  #Matrix for storing all batch information
  rawSpecBatchInfo <- matrix(nrow = numOfEle(namesList),
                             ncol = 10,
                             dimnames = list(c(),
                                             c('Exp_Year',
                                               'Sample_Year',
                                               'Location',
                                               'Species',
                                               'Materials',
                                               'Group',
                                               'Treatment',
                                               'Sample_ID',
                                               'Sample_Dup',
                                               'FTIR_Dup')))
  
  
  #Two recorders on row index for matrix filling
  tempInd <- 0
  
  for (loopI in 1:length(daysVec)) {
    
    rawSpecList[[loopI]] <- list()
    
    for (loopJ in 1:length(namesList[[loopI]])) {
      rawSpecList[[loopI]][[loopJ]] <- 
        opus_read(file = file.path('Data_R','FTIR_Raw_Spectra', 
                                   daysVec[loopI], 
                                   namesList[[loopI]][loopJ]
        ),
        simplify = TRUE)
      
      #Extract the sample information from the file names at the same time
      tempVec <- namesSplit(substr(namesList[[loopI]][loopJ],
                                   1,
                                   nchar(namesList[[loopI]][loopJ]) - 2), '_',
                            c(2025,2025,
                              'GH', 'Pmu', 
                              'cul', 'A', 
                              'UV', '001',
                              '1', '1'))
      
      
      tempInd <- tempInd + 1
      
      rawSpecBatchInfo[tempInd,] <- tempVec
      
      rm(tempVec)
    }
  }
  
  #Create a matrix to store all raw spectra
  rawSpecMat <- matrix(nrow = numOfEle(namesList),
                       ncol = length(rawSpecList[[1]][[1]][['spec']]),
                       dimnames = list(c(),
                                       c(rawSpecList[[1]][[1]][['wavenumbers']])
                       )
  )
  
  #Loop to fill the matrix
  #Create a vector recording the current index
  tempIndRawSpec <- 0
  
  for (loopI in 1:length(rawSpecList)) {
    
    for (loopJ in 1:length(rawSpecList[[loopI]])) {
      
      tempIndRawSpec <- tempIndRawSpec + 1
      
      rawSpecMat[tempIndRawSpec,] <- rawSpecList[[loopI]][[loopJ]][['spec']]
    }
  }
  
  #Re-sample the spectra to a fixed sequence of wavenumber
  #Specify the target wavenumber sequence
  tempWaveSeq <- seq(from = 3996, to = 600, by = -4)
  
  rawSpecResamMat <- t(apply(rawSpecMat, 1, 
                             FUN = resamData, 
                             inputSeq = as.numeric(colnames(rawSpecMat)),
                             outputSeq = tempWaveSeq))
  
  colnames(rawSpecResamMat) <- tempWaveSeq
  
  #Remove redundant variables
  rm(loopI,
     loopJ,
     tempInd,
     tempIndRawSpec)
}

#Smoothing and baseline correction
{
  #Apply savitzky-Golay filter
  smooSpec <- savitzkyGolay(apply(rawSpecResamMat, c(1,2), as.numeric),
                            m = 0, p = 2, w = 11)
  
  #Apply polynomial linear baseline model
  bcSpec <- baseline.modpolyfit(smooSpec,
                                degree = 3, tol = 0.00002, 
                                rep = 100)[[2]]
  
  #Rename the colomn of baseline corrected matrix, 
  #using the correct wavenumber sequence
  colnames(bcSpec) <- tempWaveSeq[c(-1:-5, 
                                    (-length(tempWaveSeq) +4): 
                                      -length(tempWaveSeq))]
}

#Z-score standardisation or standatadise to the area between 2800 and 1800 cm-1
{
  #Specify the range of wavenumebr to be used as normalising the spectra
  tempRange <- c(2800, 1800)
  tempArea <- rowSums(bcSpec[,which(as.numeric(colnames(bcSpec)) 
                                    == tempRange[1]):
                               which(as.numeric(colnames(bcSpec)) 
                                     == tempRange[2])])
  
  #Normalise the spectra by the summed area above
  bcNorSpec <- sweep(bcSpec, 1, tempArea, FUN = "/")
  
  # #Z-score
  # bcNorSpec <- t(apply(bcSpec, 1, zstandMax))
  
  #Remove useless variables
  rm(tempRange,
     tempArea)
}

#Output datasets
{
  #Whole spectra
  {
    #Z-score
    write.csv(cbind(rawSpecBatchInfo, bcNorSpec), 
              file = here('Output_R', 'bcNorSpec.csv'))
  }
}