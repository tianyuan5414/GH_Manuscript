#Load necessary packages--------------------------------------------------------
{
  #Plot
  library(ggplot2)
  library(patchwork)
  
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = file.path('Packages_R', 'Fonts', 'arial.ttf'))
  showtext_auto()
}

#Read GCMS and FTIR data
{
  pcaMatFTIRGCMalFull <- read.csv(file = file.path('Output_R',
                                                   'PLS_Nor.csv'),
                                  row.names = 1)
  
  #Rename the column
  colnames(pcaMatFTIRGCMalFull)[27] <- 'PLS_Fitted'
  
  GCDataMean <- read.csv(file = file.path('Output_R',
                                          'GCDataMean.csv'),
                         row.names = 1)
}

#Box plot on FTIR and GCMS data-------------------------------------------------
{
  #FTIR
  #PLS
  boxFTIRPLS <- ggplot(pcaMatFTIRGCMalFull, 
                       aes(y =  PLS_Fitted,
                         x = Treatment)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(aes(color = Materials), shape=16, position=position_jitter(0.2)) +
    labs(y = 'PLS predicted values') +
    scale_color_manual(values = c('cul' = 'red', 'cut' = 'blue'),
                       labels = c('Cultivars', 'Cuttings')) +
    theme_classic() +
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = c(0.9, 0.9)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  #Sub-divided by cultivars and cuttings
  #PLS
  boxSubFTIRPLS <- ggplot(pcaMatFTIRGCMalFull, 
                          aes(y = PLS_Fitted,
                            x = Treatment,
                            fill = Materials)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(aes(fill = Materials), shape=21, position=position_jitterdodge(0.2),
                color = 'black',
                stroke = 1,
                size = 2) +
    labs(y = 'PLS predicted values') +
    scale_fill_manual(values = c('cul' = 'red', 'cut' = 'blue'),
                      labels = c('Cultivars', 'Cuttings')) +
    theme_classic() +
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = c(0.9, 0.9)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  #Extract the group information between cuttings and pot cultivars
  GCDataSub <- is.na(GCDataMean[, 6])
  GCDataSub <- replace(GCDataSub, GCDataSub == FALSE,
                       c('Cuttings'))
  GCDataSub <- replace(GCDataSub, GCDataSub == TRUE,
                       c('Cultivars'))
  
  #GC-MS
  boxGCMS <- ggplot(GCDataMean, 
                    aes(y = calibrated,
                        x = Treatment)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(aes(color = GCDataSub), shape=16, position=position_jitter(0.2)) +
    scale_color_manual(values = c('Cuttings' = 'blue', 'Cultivars' = 'red')
    ) +
    labs(y = 'UV GC-MS (ng/100gr)', color = 'Materials') +
    theme_classic() +
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = c(0.9, 0.9)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  #GC-MS subgroup
  boxSubGCMS <- ggplot(GCDataMean, 
                       aes(y = calibrated,
                           x = Treatment,
                           fill = GCDataSub)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(aes(fill = GCDataSub) , shape=21, position=position_jitterdodge(0.2),
                color = 'black',
                stroke = 1,
                size = 2) +
    scale_color_manual(values = c('Cuttings' = 'blue', 'Cultivars' = 'red')
    ) +
    scale_fill_manual(values = c('Cuttings' = 'blue', 'Cultivars' = 'red'),
                      name = 'Materials') +
    labs(y = 'UV GC-MS (ng/100gr)', color = 'Materials') +
    theme_classic() +
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = c(0.9, 0.9)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
}

#Bar plot - plant divided
{
  #Bar plot on GC-MS
  barPlotPlantGCMS <- ggplot() +
    geom_boxplot(data = GCDataMean[c(1:65, 67:108),],
                 aes(x = as.character(TreeNumber),
                     y = calibrated,
                     fill = factor(Treatment, levels = c('UV', 'CON'))),
                 # stat = 'identity',
                 # position="dodge",
                 outliers = FALSE,
                 alpha = 0.6) +
    geom_boxplot(data = GCDataMean[109:133,],
                 aes(x = 'Cultivars',
                     y = calibrated,
                     fill = factor(Treatment, levels = c('UV', 'CON'))),
                 # stat = 'identity',
                 # position="dodge",
                 outliers = FALSE,
                 alpha = 0.6) +
    geom_jitter(data = GCDataMean[109:133,],
                aes(x = 'Cultivars',
                    y = calibrated,
                    fill = factor(Treatment, levels = c('UV', 'CON'))), shape=21, 
                position = position_jitterdodge(0.6),
                color = 'black',
                stroke = 1,
                size = 4,
                alpha = 0.6) +
    geom_jitter(data = GCDataMean[c(1:65, 67:108),],
                aes(x = as.character(TreeNumber),
                    y = calibrated,
                    fill = factor(Treatment, levels = c('UV', 'CON'))), shape=21, 
                position = position_jitterdodge(0.6),
                color = 'black',
                stroke = 1,
                size = 4,
                alpha = 0.6) +
    scale_y_continuous(name = 'UV GC-MS (ng/100gr)',
                       # breaks = seq(0, 0.5, by = 0.1)
    ) +
    scale_fill_manual(values = c('GCMS_UV' = 'red',
                                 'GCMS_CON' = 'blue',
                                 'UV' = 'red',
                                 'CON' = 'blue'),
                      labels = c('UV', 'CON', 'UV', 'CON'),
                      name = 'Treatment') +
    labs(x = 'Tree') +
    theme_classic() +
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = c(0.8, 0.9),
          legend.key = element_rect(fill = NA)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  #FTIR
  barPlotPlantPLS <- ggplot() +
    geom_boxplot(data = pcaMatFTIRGCMalFull[which(pcaMatFTIRGCMalFull[, 5] == 'cut'),],
                 aes(x = factor(as.numeric(Sample_ID), levels = c(1:6)),
                     y = PLS_Fitted,
                     fill = factor(Treatment, levels = c('UV', 'CON'))),
                 # stat = 'identity',
                 # position="dodge",
                 outliers = FALSE,
                 alpha = 0.6) +
    geom_boxplot(data = pcaMatFTIRGCMalFull[which(pcaMatFTIRGCMalFull[, 5] == 'cul'),],
                 aes(x = 'Cultivars',
                     y = PLS_Fitted,
                     fill = factor(Treatment, levels = c('UV', 'CON'))),
                 # stat = 'identity',
                 # position="dodge",
                 outliers = FALSE,
                 alpha = 0.6) +
    geom_jitter(data = pcaMatFTIRGCMalFull[which(pcaMatFTIRGCMalFull[, 5] == 'cut'),],
                aes(x = factor(as.numeric(Sample_ID), levels = c(1:6)),
                    y = PLS_Fitted,
                    fill = factor(Treatment, levels = c('UV', 'CON'))), shape=21, 
                position = position_jitterdodge(0.6),
                color = 'black',
                stroke = 1,
                size = 4,
                alpha = 0.6) +
    geom_jitter(data = pcaMatFTIRGCMalFull[which(pcaMatFTIRGCMalFull[, 5] == 'cul'),],
                aes(x = 'Cultivars',
                    y = PLS_Fitted,
                    fill = factor(Treatment, levels = c('UV', 'CON'))), shape=21, 
                position = position_jitterdodge(0.6),
                color = 'black',
                stroke = 1,
                size = 4,
                alpha = 0.6) +
    scale_y_continuous(name = 'PLS fitted value',
                       # breaks = seq(0, 3, by = 1)
    ) +
    scale_x_discrete(name = 'Plant') +
    scale_fill_manual(values = c('PLS_UV' = 'red',
                                 'PLS_CON' = 'blue',
                                 'UV' = 'red',
                                 'CON' = 'blue'),
                      labels = c('UV', 'CON', 'UV', 'CON'),
                      name = 'Treatment') +
    labs(x = 'Tree') +
    theme_classic() +
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 14, family = 'arial'),
          axis.text.y = element_text(size = 12, family = 'arial'),
          axis.title.x = element_text(size = 14, family = 'arial'),
          axis.text.x = element_text(size = 12, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          legend.position = c(0.8, 0.9),
          legend.background = element_rect(fill = NA)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
}

#Combine plots on GCMS signals and FTIR PLS fitted values
{
  #FTIR boxplot
  pdf(file = file.path( 'Figures_R', 'FTIR.pdf'),
      width = 12, # The width of the plot in inches
      height = 12)
  
  (boxFTIRPLS | boxSubFTIRPLS) /
    (barPlotPlantPLS | barPlotPlantPLS)
  
  dev.off()
  
  #GCMS boxplot
  pdf(file = file.path( 'Figures_R', 'GCMS.pdf'),
      width = 12, # The width of the plot in inches
      height = 12)
  
  (boxGCMS | boxSubGCMS) /
    (barPlotPlantGCMS | barPlotPlantGCMS)
  
  dev.off()
}