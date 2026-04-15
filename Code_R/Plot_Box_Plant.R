#' This script aims at plotting plant-groupped box plots of either PLS fitted values or GC-MS signals.
#' The plot may be further ground by sample types (potted cultivars or branch cuttings)

plot_box_plant <- function(dataMat,
                     variableName){
  if (variableName == 'PLS') {
    yLabel <- 'PLS predicted values'
    rowIndex <- which(dataMat[['Materials']] == 'cut')
  } else {
    yLabel <- 'UV GC-MS (ng/100gr)'
    rowIndex <- 1:nrow(dataMat)
  }
  
  plot_output <- ggplot() +
    geom_boxplot(data = dataMat[rowIndex,],
                 aes(x = factor(as.numeric(Sample_ID), levels = c(1:6)),
                     y = dataMat[[variableName]][rowIndex],
                     fill = factor(Treatment, levels = c('UV', 'CON'))),
                 outliers = FALSE,
                 alpha = 0.6) +
    geom_boxplot(data = dataMat[which(dataMat[['Materials']] == 'cul'),],
                 aes(x = 'Cultivars',
                     y = dataMat[[variableName]][which(dataMat[['Materials']] == 'cul')],
                     fill = factor(Treatment, levels = c('UV', 'CON'))),
                 outliers = FALSE,
                 alpha = 0.6) +
    geom_jitter(data = dataMat[which(dataMat[['Materials']] == 'cul'),],
                aes(x = 'Cultivars',
                    y = dataMat[[variableName]][which(dataMat[['Materials']] == 'cul')],
                    fill = factor(Treatment, levels = c('UV', 'CON'))), shape=21, 
                position = position_jitterdodge(0.6),
                color = 'black',
                stroke = 1,
                size = 4,
                alpha = 0.6) +
    geom_jitter(data = dataMat[rowIndex,],
                aes(x = factor(as.numeric(Sample_ID), levels = c(1:6)),
                    y = dataMat[[variableName]][rowIndex],
                    fill = factor(Treatment, levels = c('UV', 'CON'))), shape=21, 
                position = position_jitterdodge(0.6),
                color = 'black',
                stroke = 1,
                size = 4,
                alpha = 0.6) +
    scale_y_continuous(name = yLabel,
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
    )
  
  return(plot_output)
}