#This script aims at reading and ploting UV-B monitoring data in climate room experiment
#Load necessary packages
{
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  
  #Palette
  library(RColorBrewer)
  
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = file.path('Fonts', 'arial.ttf'))
  showtext_auto()
}

#Declare custom functions
{
  # Function to extract one spectrum
  read_spectrum <- function(block) {
    
    # Spectrum ID
    id_line <- grep("^IDData\\s*=", block, value = TRUE)
    id <- sub(".*=\\s*", "", id_line)
    
    # Find DATA section
    data_start <- grep("^\\[DATA\\]", block) + 1
    
    # Find END section
    data_end <- grep("^\\[END]\\s\\of\\s\\[DATA]", block)
    if(length(data_end) == 0)
      data_end <- length(block)
    
    data_lines <- block[data_start:(data_end[1]-1)]
    
    # Remove blank lines
    data_lines <- data_lines[nzchar(trimws(data_lines))]
    
    # Read table
    dat <- read.table(text = data_lines,
                      fill = TRUE,
                      stringsAsFactors = FALSE)
    
    # Remove metadata row
    dat <- dat[-1, ]
    
    # Assign column names
    colnames(dat) <- c("Wavelength", "Intensity", "Error", "Status")
    
    # Convert to numeric
    dat$Wavelength <- as.numeric(dat$Wavelength)
    
    dat$Intensity[dat$Intensity == "-NAN"] <- NA
    dat$Error[dat$Error == "-NAN"] <- NA
    
    dat$Intensity <- as.numeric(dat$Intensity)
    dat$Error <- as.numeric(dat$Error)
    
    # Keep wavelength, intensity and error
    dat %>%
      select(Wavelength, Intensity) %>%
      rename(
        !!id := Intensity
      )
  }
  
  namesSplit <- function(folderNames, seperator, defaultSettings) {
    
    splitedNames <- defaultSettings
    
    #Split the names
    splitedNamesRaw <- strsplit( x = folderNames, split = seperator)[[1]]
    
    splitedNames[1] <- substring(splitedNamesRaw[4],
                                 7,
                                 nchar(splitedNamesRaw[4]))
    
    splitedNames[2] <- substring(splitedNamesRaw[5],
                                 6,
                                 nchar(splitedNamesRaw[5]))
    
    splitedNames[3] <- splitedNamesRaw[6]
    
    splitedNames[4] <- substring(splitedNamesRaw[8],
                                 1,
                                 nchar(splitedNamesRaw[8]) - 2)
    
    splitedNames[5] <- splitedNamesRaw[1]
    
    return(splitedNames)
  }
  
  specPlot <- function(dataset, frameInd, modInd, treatmentInd) {
    
    if (frameInd == 1) {
      posLeg <- 'top'
    } else {
      posLeg <- 'none'
    }
    
    if (modInd == 1) {
      yAxisLab <- 1
    } else {
      yAxisLab <- NA
    }
    
    if (frameInd != length(frameID)) {
      xAxisLab <- NA
    }else{
      xAxisLab <- 1
    }
    
    if (!is.na(yAxisLab)) {
      tempY <- element_text(size = 10, family = 'arial')
    } else {
      tempY <- element_blank()
    }
    
    if (!is.na(xAxisLab)) {
      tempX <- element_text(size = 10, family = 'arial')
    } else {
      tempX <- element_blank()
    }
    
    tempPlot <- 
      ggplot(data = dataset[which(dataset[, 2] == frameID[frameInd] & 
                                       (dataset[, 5] == modeSetting[modInd]) &
                                    (dataset[, 3] == treatmentTy[treatmentInd])),],
             aes(x = wavelength,
                 y = dose,
                 colour = height)) +
      geom_line(linewidth = 0.5,
                alpha = 0.5) +
      geom_vline(xintercept = c(280,315),
                 linewidth = 0.5,
                 linetype = 'dashed',
                 colour = 'black',
                 alpha = 0.6) +
      scale_x_continuous(limits = c(250, 550),
                         breaks = seq(from = 250, to = 550, by = 50),
                         name = 'Wavelength (nm)') +
      scale_y_continuous(name = 'Dose (mW m-2 nm-1))') +
      scale_color_manual(values = brewer.pal(length(unique(dataset$height)), 'Dark2')) +
      labs(colour = 'Height (cm)',
           title = paste('Frame ', frameID[frameInd], ' ', modeSetting[modInd], ' mode ',
                         treatmentTy[treatmentInd], sep = '',
                         collapse = NULL)) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, family = 'arial'), 
            axis.title.y = tempY,
            axis.text.y = element_text(size = 8, family = 'arial'),
            axis.title.x = tempX,
            axis.text.x = element_text(size = 8, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = posLeg
            # axis.line.y = element_blank(),
            # axis.text.x = element_blank(), 
            # axis.line.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank()
      )
    
    return(tempPlot)
  }
  
  shoelace_area <- function(x, y) {
    # Ensure the polygon closes by appending the first vertex to the end
    x_closed <- c(x, x[1])
    y_closed <- c(y, y[1])
    
    # Number of vertices
    n <- length(x_closed)
    
    # Cross multiplication arrays
    sum1 <- sum(x_closed[1:(n-1)] * y_closed[2:n])
    sum2 <- sum(y_closed[1:(n-1)] * x_closed[2:n])
    
    # Calculate absolute area
    area <- 0.5 * abs(sum1 - sum2)
    return(area)
  }
  
  uvbDoseCal <- function(dataset, batchInfo, batchInd) {
    tempDF <- dataset[which((dataset[, 2] == batchInfo[batchInd, 1]) &
                                 (dataset[, 3] == batchInfo[batchInd, 2]) &
                                 (dataset[, 4] == batchInfo[batchInd, 3]) &
                                 (dataset[, 5] == batchInfo[batchInd, 4]) &
                                 (dataset[, 6] >= uvbRange[1]) &
                                 (dataset[, 6] <= uvbRange[2])),]
    
    tempUVBDose <- shoelace_area(tempDF$wavelength,
                                 tempDF$dose) * 3600 * 6 / 1000 / 1000
    
    return(tempUVBDose)
  }
}

#Read data
{
  #Access all folders
  folderList <- list.files(path = file.path('Data_R', 'UV_B_Monitoring'))
  
  #Get all addresses of raw data
  rawFileList <- c()
  for (loopI in 1:length(folderList)) {
    tempFileList <- list.files(path = file.path('Data_R', 'UV_B_Monitoring', folderList[loopI]))
    
    rawFileList <- c(rawFileList, file.path('Data_R', 'UV_B_Monitoring', folderList[loopI],
                                            tempFileList))
    
    rm(tempFileList,
       loopI)
  }
  
  #Read all UV-B spectra
  rawSpecList <- vector(mode = 'list', length = length(rawFileList))
  for (loopI in 1:length(rawSpecList)) {
    # Read file
    txt <- readLines(file.path(rawFileList[loopI]), warn = FALSE)
    
    # Locate the beginning of each spectrum
    spec_start <- grep("^\\[Spectrum\\]", txt)
    
    # Split file into spectrum blocks
    spec_blocks <- lapply(seq_along(spec_start), function(i) {
      
      start <- spec_start[i]
      end <- if (i < length(spec_start))
        spec_start[i + 1] - 1
      else
        length(txt)
      
      txt[start:end]
    })
    
    # Read all spectra
    spectra_list <- lapply(spec_blocks, read_spectrum)
    
    # Merge by wavelength
    rawSpecList[[loopI]] <- reduce(spectra_list, full_join, by = "Wavelength")
    
    names(rawSpecList)[[loopI]] <- substring(sub('.*Monitoring/2026_', '', rawFileList[loopI]),
                                             1,
                                             nchar(sub('.*Monitoring/2026_', '', rawFileList[loopI])) - 4)
    
    rm(txt,
       spec_start,
       spec_blocks,
       spectra_df,
       spectra_list,
       loopI)
  }
}

#Adjust the data structure of the raw spectra
{
  adjSpecList <- rawSpecList
  
  for (loopI in 1:length(rawSpecList)) {
    adjSpecList[[loopI]] <- cbind(wavelength = rawSpecList[[loopI]][, 1],
                                  dose = apply(rawSpecList[[loopI]][,-1], 1, mean),
                                  sd = apply(rawSpecList[[loopI]][,-1], 1, sd))
    
    adjSpecList[[loopI]][which(is.nan(adjSpecList[[loopI]]))] <- 0
    adjSpecList[[loopI]][which(is.na(adjSpecList[[loopI]]))] <- 0
  }
  
  rm(loopI)
  
  #Extract the batch information from the names of raw spectra list, and combine
    #all spectra together with batch information
  tempBatchInfo <- namesSplit(names(rawSpecList)[[1]],
                              '_',
                              vector(length = 5))
  
  combSpecDF <- cbind(date = tempBatchInfo[1],
                      frame = tempBatchInfo[2],
                      treatment = tempBatchInfo[3],
                      height = tempBatchInfo[4],
                      mode = tempBatchInfo[5],
                      adjSpecList[[1]])
  
  for (loopI in 2:length(adjSpecList)) {
    tempBatchInfo <- namesSplit(names(rawSpecList)[[loopI]],
                                '_',
                                vector(length = 5))
    
    tempDF <- cbind(date = tempBatchInfo[1],
                    frame = tempBatchInfo[2],
                    treatment = tempBatchInfo[3],
                    height = tempBatchInfo[4],
                    mode = tempBatchInfo[5],
                    adjSpecList[[loopI]])
    
    combSpecDF <- rbind(combSpecDF,
                        tempDF)
    
    rm(tempBatchInfo,
       tempDF)
  }
  
  rm(loopI)
  
  combSpecDF <- as.data.frame(combSpecDF)
  combSpecDF[, c(6:8)] <- apply(combSpecDF[, c(6:8)] , c(1,2), as.numeric)
}

#Plot dose spectra for each frame
{
  frameID <- c('A', 'B', 'C', 'D', 'E', 'F')
  modeSetting <- c('low', 'high')
  treatmentTy <- c('con', 'UV')
  
  specPlotList <- vector(mode = 'list',
                         length = length(frameID))
  
  for (loopI in 1:length(frameID)) {
    
    specPlotList[[loopI]] <- vector(mode = 'list',
                                    length = length(modeSetting))
    
    for (loopJ in 1:length(modeSetting)) {
      
      specPlotList[[loopI]][[loopJ]] <- vector(mode = 'list',
                                      length = length(treatmentTy))
      
      for (loopK in 1:length(treatmentTy)) {
        specPlotList[[loopI]][[loopJ]][[loopK]] <- specPlot(combSpecDF,
                                                   loopI, 
                                                   loopJ,
                                                   loopK)
      }
    }
  }
  
  rm(loopI,
     loopJ,
     loopK)
  
  #Output as pdf
  pdf(file = file.path('Figures_R', 'UV_B_Monitoring_con.pdf'),   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 16)
  
  (specPlotList[[1]][[1]][[1]] | specPlotList[[1]][[2]][[1]]) /
    (specPlotList[[2]][[1]][[1]] | specPlotList[[2]][[2]][[1]]) /
    (specPlotList[[3]][[1]][[1]] | specPlotList[[3]][[2]][[1]]) /
    (specPlotList[[4]][[1]][[1]] | specPlotList[[4]][[2]][[1]]) /
    (specPlotList[[5]][[1]][[1]] | specPlotList[[5]][[2]][[1]]) /
    (specPlotList[[6]][[1]][[1]] | specPlotList[[6]][[2]][[1]]) 
  
  dev.off()
  
  #Output as pdf
  pdf(file = file.path('Figures_R', 'UV_B_Monitoring_UV.pdf'),   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 16)
  
  (specPlotList[[1]][[1]][[2]] | specPlotList[[1]][[2]][[2]]) /
    (specPlotList[[2]][[1]][[2]] | specPlotList[[2]][[2]][[2]]) /
    (specPlotList[[3]][[1]][[2]] | specPlotList[[3]][[2]][[2]]) /
    (specPlotList[[4]][[1]][[2]] | specPlotList[[4]][[2]][[2]]) /
    (specPlotList[[5]][[1]][[2]] | specPlotList[[5]][[2]][[2]]) /
    (specPlotList[[6]][[1]][[2]] | specPlotList[[6]][[2]][[2]]) 
  
  dev.off()
}

#Calculate the UV-B doses for each frame
{
  #Get the unique settings for monitoring dataset
  doseFrame <- cbind(unique(combSpecDF[,2:5]),
                     dose = NA,
                     sd = NA)
  
  #Calculate the UV-B dose (280-315 nm)
  uvbRange <- c(280, 315)
  for (loopI in 1:nrow(doseFrame)) {
    doseFrame[loopI, 5] <- uvbDoseCal(combSpecDF,
                                      doseFrame,
                                      loopI)
  }
  
  rm(loopI)
  
  #Calculate the mean and sd when combing homogeneous data from different frames
  doseAll <- cbind(unique(doseFrame[,2:4]),
                   dose = NA,
                   sd = NA)
  
  for (loopI in 1:nrow(doseAll)) {
    doseAll[loopI, 4] <- mean(doseFrame[which((doseFrame[, 2] == doseAll[loopI, 1])&
                                                (doseFrame[, 3] == doseAll[loopI, 2])&
                                                (doseFrame[, 4] == doseAll[loopI, 3])), 5])
    doseAll[loopI, 5] <- sd(doseFrame[which((doseFrame[, 2] == doseAll[loopI, 1])&
                                                (doseFrame[, 3] == doseAll[loopI, 2])&
                                                (doseFrame[, 4] == doseAll[loopI, 3])), 5])
  }
  
  #Plot the daily doses
  dosePlotLow <- ggplot(data = doseAll[which(doseAll[, 3] == 'low'),],
                        aes(x = height,
                            y = dose,
                            group = treatment)) +
    geom_bar(aes(fill = treatment), stat = 'identity',position = position_dodge()) +
    geom_errorbar(aes(ymin = dose - sd, ymax = dose + sd), width = 0.2,
                  position = position_dodge(0.9)) +
    scale_fill_manual(values = brewer.pal(2, 'Set2'),
                      labels = c('Control', 'UV-B')) +
    scale_x_discrete(name = 'Height (cm)') +
    scale_y_continuous(name = 'Dose (KJ m-2 Day-1)',
                       limits = c(0, 45),
                       breaks = seq(0,45, by = 5)) +
    labs(fill = 'Treatment',
         title = 'Low mode') +
    theme_classic() +
    theme(plot.title = element_text(size = 16, family = 'arial'), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = 'none'
    )
 
  dosePlotHigh <- ggplot(data = doseAll[which(doseAll[, 3] == 'high'),],
                        aes(x = height,
                            y = dose,
                            group = treatment)) +
    geom_bar(aes(fill = treatment), stat = 'identity',position = position_dodge()) +
    geom_errorbar(aes(ymin = dose - sd, ymax = dose + sd), width = 0.2,
                  position = position_dodge(0.9)) +
    scale_fill_manual(values = brewer.pal(2, 'Set2'),
                      labels = c('Control', 'UV-B')) +
    scale_x_discrete(name = 'Height (cm)') +
    scale_y_continuous(name = 'Dose (KJ m-2 Day-1)',
                       limits = c(0, 45),
                       breaks = seq(0,45, by = 5)) +
    labs(fill = 'Treatment',
         title = 'High mode') +
    theme_classic() +
    theme(plot.title = element_text(size = 16, family = 'arial'), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    ) 
  
  #Output as pdf
  pdf(file = file.path('Figures_R', 'UV_B_Monitoring_Overall.pdf'),   # The directory you want to save the file in
      width = 16, # The width of the plot in inches
      height = 8)
  
  dosePlotLow | dosePlotHigh
  
  dev.off()
}