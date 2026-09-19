{
  #Load palaeochem package
  library(palaeochem)
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = file.path('Packages_R', 'Fonts', 'arial.ttf'))
  showtext_auto()
  
  #Packages needed in Vivian's code
  source(
    file.path(
      "Code_R", "00_Config_file.r"
    )
  )
  
  rm(package_list)
}

{
  #Read raw FTIR spectra
  rawFTIRSpec <- read.csv(file = file.path('Data_R',
                                           'Original_Spectra.csv'),
                          row.names = 1)
  colnames(rawFTIRSpec)[-1:-10] <- substring(colnames(rawFTIRSpec)[-1:-10], 2)
  
  plot_original_spec_cul <-
    rawFTIRSpec[rawFTIRSpec[, 5] == 'cul', c(-1:-2,-4:-10)] |> plot_spectra(
      group_var = 'unique_id',
      col_var = 'unique_id'
    ) +
    labs(title = 'Potted cultivars',
         x = '')
  
  plot_original_spec_cut <-
    rawFTIRSpec[rawFTIRSpec[, 5] == 'cut', c(-1:-2,-4:-10)] |> plot_spectra(
      group_var = 'unique_id',
      col_var = 'unique_id'
    ) +
    labs(title = 'Branch cuttings')
  
  #5 - smoothing + 2nd derivative
  processing_steps_ftir5 <- list(
    functions = c(
      get_savitzkygolay,
      get_filter_wavenumbers
    ),
    parameters = list(
      list(width = 7, poly = 2, deriv = 2),
      list(min_wn = 800, max_wn = 1800)
    )
  )
  
  data_preprocessedNew <- list(
    cbind(
      unique_id = rawFTIRSpec[, 3],
      materials = rawFTIRSpec[, 5],
      treatment = rawFTIRSpec[, 7],
      corrected5 = run_preprocessing_list(rawFTIRSpec[, c(-1:-2,-4:-10)], processing_steps_ftir5)
    )
  )
  
  data_preprocessedNew <- lapply(data_preprocessedNew, as.data.frame)
  for (loopI in 1:length(data_preprocessedNew)) {
    data_preprocessedNew[[loopI]][, -1:-3] <-
      apply(data_preprocessedNew[[loopI]][, -1:-3], c(1,2), as.numeric)
  }
  
  rm(processing_steps_ftir5)
}

#Convert the list of processed data into a dataframe
pcaDF2ndFTIR <- cbind.data.frame(data_preprocessedNew)
rownames(pcaDF2ndFTIR) <- 1:nrow(pcaDF2ndFTIR)

# Ordination
treat_pca <- ordr::ordinate(pcaDF2ndFTIR, cols=colnames(pcaDF2ndFTIR)[4:522], model= ~ prcomp(.), argument=c("materials", "treatment"))

treat_pca_prcomp <- prcomp(pcaDF2ndFTIR[, -1:-3], scale = FALSE)

fviz_eig(treat_pca_prcomp)

#Apply the PCA model
pcaMod2ndFTIR <- run_pca(pcaDF2ndFTIR[, -2:-3])

#Calculate the length of loading arrows, and filtering out the longest 100
pcaArrowLength <- pcaMod2ndFTIR$rotation[, 1]^2 + pcaMod2ndFTIR$rotation[, 2]^2
pcaArrowLongest40 <- names(pcaArrowLength)[order(pcaArrowLength ,decreasing=T)[1:40]]

pcaDF2ndFTIR %>%
  subset(select = pcaArrowLongest40) %>%
  as.matrix() %>%
  print() -> treat_data_arrows

treat_data_arrows <- treat_data_arrows * 0.03
colnames(treat_data_arrows) <- pcaArrowLongest40

lm(treat_data_arrows ~ get_rows(treat_pca)) %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  print() -> treat_data_arrows_plot

#PCA Plot
ordr::ggbiplot(treat_pca, sec.axes="cols", scale.factor=1) +
  # Add points and colour via treatment
  geom_rows_point(aes(color= materials, shape = treatment)) +
  # Add arrows
  geom_cols_vector(data = treat_data_arrows_plot,color="navy", size=1) +
  # Add ellipses - 95% confidence
  geom_mark_ellipse(aes(group = materials), 
                    stat = 'rows_ellipse',
                    size=0.2, level=.95,
                    label.buffer = unit(0.5, 'mm'),
                    con.cap = 0) +
  # Add arrow labels
  geom_cols_text_radiate(data = treat_data_arrows_plot, aes(label= round(as.numeric(name))), size=4, color="navy") +
  scale_color_manual(values = c('#f0d5d0', '#6c8f7e')) +
  # Minimal theme +
  theme_classic(base_size=16) +
  # Add grid
  geom_hline(yintercept=0, linetype="dashed", size=1) +
  geom_vline(xintercept=0, linetype="dashed", size=1) +
  theme(
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text=element_text(family="arial"),
    legend.text = element_text(size = 20, family="arial")
  )