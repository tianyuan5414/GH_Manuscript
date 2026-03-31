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
}

#Read malformation data---------------------------------------------------------
{
  cuttingsPlantsMat <- read.csv(file = file.path('Output_R',
                                                 'cuttingsPlantsMat.csv'),
                                row.names = 1)
  
  #Transfer donar plant number into factor
  cuttingsPlantsMat$Plant <- as.character(cuttingsPlantsMat$Plant)
  
  cultivarMean <- as.vector(read.csv(file = file.path('Output_R',
                                           'cultivarMean.csv'),
                          row.names = 1))[[1]]
}

#Plot malformation figures which compares malformation rates for cutting 
  #branches samples and subgrouped by donor plants
{
  #Point plot to compare ratios between con and uv treatment, on different plants
  cuttingsPlantsPlot <-
    ggplot() +
    geom_point(data = cuttingsPlantsMat,
               aes(x = CON_Mal, y = UV_Mal, col = Plant)) +
    geom_point(aes(x = cultivarMean[3] , y = cultivarMean[1], 
                   color = 'Cultivars')) +
    geom_errorbar(data = cuttingsPlantsMat,
                  aes(x = CON_Mal, y = UV_Mal, col = Plant,
                      ymin = UV_Mal - UV_Mal_sd,
                      ymax = UV_Mal + UV_Mal_sd),
                  linewidth = 1,
                  width = 0.01) +
    geom_errorbar(data = cuttingsPlantsMat,
                  aes(x = CON_Mal, y = UV_Mal, col = Plant,
                      xmin = CON_Mal - CON_Mal_sd,
                      xmax = CON_Mal + CON_Mal_sd),
                  linewidth = 1,
                  width = 0.01,
                  orientation = 'y') +
    geom_errorbar(aes(x = cultivarMean[3] , y = cultivarMean[1], color = 'Cultivars',
                      ymin = cultivarMean[1] - cultivarMean[2],
                      ymax = cultivarMean[1] + cultivarMean[2]),
                  linewidth = 1,
                  width = 0.005) +
    geom_errorbar(aes(x = cultivarMean[3] , y = cultivarMean[1], color = 'Cultivars',
                      xmin = cultivarMean[3] - cultivarMean[4],
                      xmax = cultivarMean[3] + cultivarMean[4]),
                  linewidth = 1,
                  width = 0.005,
                  orientation = 'y') +
    geom_abline (slope=1, linetype = "dashed", color="grey") +
    scale_x_continuous(limits = c(-0.05, 0.15), name = 'CON malformation percentage',
                       breaks = seq(0, 0.15, by = 0.05)) +
    scale_y_continuous(limits = c(-0.05, 0.25), name = 'UV malformation percentage',
                       breaks = seq(0, 0.25, by = 0.05)) +
    scale_color_manual(name = 'Plant',
                       values = c('1' = brewer.pal(6, 'Set2')[1],
                                  '2' = brewer.pal(6, 'Set2')[2],
                                  '3' = brewer.pal(6, 'Set2')[3],
                                  '4' = brewer.pal(6, 'Set2')[4],
                                  '5' = brewer.pal(6, 'Set2')[5],
                                  '6' = brewer.pal(6, 'Set2')[6],
                                  'Cultivars' = 'black')) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(size = 20, family = 'arial'),
      axis.text.x = element_text(size = 16, family = 'arial'),
      
      axis.text.y = element_text(size = 20, family = 'arial'),
      axis.title.y = element_text(size = 16, family = 'arial'),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.2),
      legend.text = element_text(size = 14, family = 'arial')
    )
}

#Output figures on malformation rate
{
  #Malformation
  pdf(file = file.path( 'Figures_R', 'Malformation.pdf'),
      width = 12, # The width of the plot in inches
      height = 12)
  
  cuttingsPlantsPlot
  
  dev.off()
}