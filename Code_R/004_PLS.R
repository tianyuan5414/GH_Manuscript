#Load necessary packages--------------------------------------------------------
{
  #Plot
  library(ggplot2)
  library(patchwork)
  
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = file.path('Packages_R', 'Fonts', 'arial.ttf'))
  showtext_auto()
  
  #Palette
  library(RColorBrewer)
  
  #PLS
  library(pls)
}

#Declare necessary functions
{
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
  
  #Generate a string to show the statistics of a regression model
  lm_eqn2 <- function(compI){
    m <- lm(modelHeight$model$calibrated ~ modelHeight$fitted.values[,,compI]);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  
  #Detect the symbol of given value, and output 1 or -1 accordingly
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

#PLS analysis-------------------------------------------------------------------
{
  #Height
  {
    #Read combined dataset of both FTIR peaks and GCMS signals
    pcaMatHeight <- read.csv(file = 
                               file.path('Output_R', 'pcaMatHeightCom.csv'),
                             row.names = 1)
    
    #Adjust the colomn names so that they reflect the wavenumber position of the
      #peak
    colnames(pcaMatHeight)[11:24] <- substr(colnames(pcaMatHeight)[11:24], 
                                            2,
                                            nchar(colnames(pcaMatHeight)[11:24]))
    
    #Clean up the matrix to exclude NA values in GC-MS data
    pcaMatHeight <- pcaMatHeight[which(!is.na(pcaMatHeight[, ncol(pcaMatHeight) - 1])), ]
    pcaMatHeight <- pcaMatHeight[which(apply(pcaMatHeight[, -1:-10], 1, cleanNA)),]
    
    #PLS model
    modelHeight <- plsr(calibrated ~ ., 
                        data = as.data.frame(apply(pcaMatHeight[,c(-1:-10, -ncol(pcaMatHeight))],
                                                   c(1,2),
                                                   as.numeric)), 
                        scale=FALSE, validation="CV")
    
    #Investigate important results
    validationplot(modelHeight)
    validationplot(modelHeight, val.type="R2")
    
    #Setup the number of components to be analysed
    numCompsHeight <- selectNcomp(modelHeight, method = "randomization", plot = TRUE)
    
    #Plot PLS predictions with GC-MS
    plsPlotHeight <- ggplot() + 
      geom_point(aes(x = modelHeight$fitted.values[,,numCompsHeight],
                     y = modelHeight$model$calibrated)) + 
      geom_smooth(aes(x = modelHeight$fitted.values[,,numCompsHeight],
                      y = modelHeight$model$calibrated),
                  method = "lm") +
      geom_text(
        aes(x = 0.35, y = 0.7,
            label = lm_eqn2(numCompsHeight)), parse = TRUE) +
      labs(x = 'PLS fitted values (peak height)',
           y = 'GC-MS (ng/150gr)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
            # axis.line.y = element_blank(),
            # axis.text.x = element_blank(), 
            # axis.line.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank()
      )
    
    #Plot correlation coefficient under the selected number of components
    modelHeightCof <- ggplot() +
      geom_bar(aes(x = as.numeric(colnames(pcaMatHeight)[11:24]), 
                   y = modelHeight$coefficients[, , numCompsHeight]),
               stat = 'identity') +
      geom_hline(yintercept = 0) +
      geom_text(
        aes(x = as.numeric(colnames(pcaMatHeight)[11:24]), 
            y = modelHeight$coefficients[, ,numCompsHeight] + 
              min(abs(modelHeight$coefficients[, ,numCompsHeight])) *0.4 * 
              symbolDec(modelHeight$coefficients[, ,numCompsHeight]),
            label = as.numeric(colnames(pcaMatHeight)[11:24])), parse = TRUE, 
        check_overlap = T) +
      scale_x_reverse(limit = c(3500, 500),
                      breaks = seq(3500, 500, by = -500)) +
      labs(x = 'Wavenumber', y = paste('Variable coefficients ',
                                       numCompsHeight, ' components', sep = '', collapse = NULL)) +
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
            # axis.line.y = element_blank(),
            # axis.text.x = element_blank(), 
            # axis.line.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank()
      )
    
    #Plot variable loadings under the first components
    modelHeightLoa <- ggplot() +
      geom_bar(aes(x = as.numeric(colnames(pcaMatHeight)[11:24]), y = modelHeight$loadings[, 1]),
               stat = 'identity') +
      geom_hline(yintercept = 0) +
      geom_text(
        aes(x = as.numeric(colnames(pcaMatHeight)[11:24]), y = modelHeight$loadings[, 1] + 
              min(abs(modelHeight$loadings[, 1])) *0.4 * symbolDec(modelHeight$loadings[, 1]),
            label = as.numeric(colnames(pcaMatHeight)[11:24])), parse = TRUE, 
        check_overlap = T) +
      scale_x_reverse(limit = c(3500, 500),
                      breaks = seq(3500, 500, by = -500)) +
      labs(x = 'Wavenumber', y = paste('Variable loadings ',
                                       1, sep = '', collapse = NULL)) +
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
            # axis.line.y = element_blank(),
            # axis.text.x = element_blank(), 
            # axis.line.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank()
      )
    
    #Calculate variance explained by different components
    varCompMat <- matrix(nrow = numCompsHeight,
                         ncol = 2,
                         dimnames = list(c(),
                                         c('Accu',
                                           'Partial')))
    varCompMat[, 1] <- modelHeight[[12]][1:numCompsHeight]
    for (loopI in 1:nrow(varCompMat)) {
      varCompMat[loopI, 2] <- sum(modelHeight[[12]][1:loopI])
    }
    
    #Plot variance explained by different components
    modelHeightVar <- ggplot() +
      geom_bar(aes(x = 1:numCompsHeight, y = varCompMat[, 2] / sum(modelHeight[[12]])),
               fill = 'red',
               stat = 'identity') +
      geom_bar(aes(x = 1:numCompsHeight, y = varCompMat[, 2] / sum(modelHeight[[12]] )- 
                     varCompMat[, 1] / sum(modelHeight[[12]])),
               fill = 'grey',
               stat = 'identity') +
      labs(x = 'Component', y = 'Variance explained') +
      scale_x_continuous(breaks = seq(1, numCompsHeight, by = 1)) +
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
            # axis.line.y = element_blank(),
            # axis.text.x = element_blank(), 
            # axis.line.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank()
      )
    
    #Plot PLS residuals against GCMS signals under the selected number of
    #components
    plsPlotHeightResi <- ggplot() + 
      geom_point(aes(x = modelHeight$fitted.values[,,numCompsHeight],
                     y = modelHeight$residuals[,,numCompsHeight])) + 
      geom_hline(yintercept = 0) +
      labs(x = 'PLS fitted values (peak height)',
           y = 'Residuals') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
            # axis.line.y = element_blank(),
            # axis.text.x = element_blank(), 
            # axis.line.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank()
      )
    
    #Remove unused variables
    rm(loopI)
  }
}

#Output PLS results-------------------------------------------------------------
{
  #Output PLS model on current FTIR dataset
  pdf(file = file.path('Figures_R', 'PLS_Full.pdf'),   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 12)
  
  (plsPlotHeight | modelHeightVar) /
    (modelHeightCof | modelHeightLoa) /
    (plsPlotHeightResi | plsPlotHeightResi)
  
  dev.off()
  
  #Combine the PLS predicted values to the FTIR peak statistics
  tempMat <- cbind(pcaMatHeight, modelHeight$fitted.values[,, numCompsHeight])
  
  #Output PLS predicted values and FTIR statistics together
  write.csv(tempMat, file = file.path( 'Output_R', 'PLS_Nor.csv'))
  
  #Remove unused variables
  rm(tempMat)
}