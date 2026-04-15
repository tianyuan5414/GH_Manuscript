#' This script aims at plotting box plots of either PLS fitted values or GC-MS signals.
#' The plot may be further ground by sample types (potted cultivars or branch cuttings)

plot_box <- function(dataMat,
                     variableName,
                     groupName = NA){
  if (variableName == 'PLS') {
    yLabel <- 'PLS predicted values'
  } else {
    yLabel <- 'UV GC-MS (ng/100gr)'
  }

  if (is.na(groupName)) {
    
    plot_output <- ggplot(dataMat, 
                          aes(y = dataMat[[variableName]],
                              x = Treatment)) +
      geom_boxplot(outliers = FALSE) +
      geom_jitter(aes(color = Materials), shape=16, position=position_jitter(0.2)) +
      labs(y = yLabel) +
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
      )
  } else {
    plot_output <- ggplot(dataMat, 
                            aes(y = dataMat[[variableName]],
                              x = Treatment,
                              fill = Materials)) +
      geom_boxplot(outliers = FALSE) +
      geom_jitter(aes(fill = Materials), shape=21, position=position_jitterdodge(0.2),
                  color = 'black',
                  stroke = 1,
                  size = 2) +
      labs(y = yLabel) +
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
      )
  }
  
  return(plot_output)
}