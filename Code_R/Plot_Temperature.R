#This script aims at reading and plotting all temperature monitoring data in the GH 2025
#Load necessary packages
{
  #Read spreadsheet
  library(readxl)
  
  #Plot
  library(ggplot2)
  library(dplyr)
  
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = file.path('Packages_R', 'Fonts', 'arial.ttf'))
  showtext_auto()
}

#Read and pre-tread spreadsheet data
{
  #Raw spreadshee with temperature data
  rawTempHumName <- excel_sheets(path = file.path('Data_R', 'Temperature_Monitoring_Susanne.xlsx'))
  
  # Read all sheets into a list
  rawTempHum <- lapply(rawTempHumName, function(x) read_excel(path = file.path('Data_R', 'Temperature_Monitoring_Susanne.xlsx')
                                                              ,sheet = x
                                                              )
                       )
  
  #Rename the list elements and the column names for each element in the list
  names(rawTempHum) <- rawTempHumName
  for (loopI in 1:length(rawTempHum)) {
    colnames(rawTempHum[[loopI]]) <- c(
      'No',
      "Timestamp", 
      "Temperature", 
      "Humidity",
      'Dew_Point',
      'Serial_Num'
      )
  }
}

#Create a function to plot the temperature changes recorded in each matrix in the
  #raw data list
plot_temperature <- function(temp_list, frames = NULL) {
  # frames: character vector of sheet names to plot.
  #         Defaults to NULL, which plots all sheets.
  if (is.null(frames)) {
    frames <- names(temp_list)
  } else {
    missing_frames <- setdiff(frames, names(temp_list))
    if (length(missing_frames) > 0) {
      stop("The following frames were not found in temp_list: ",
           paste(missing_frames, collapse = ", "))
    }
  }
  
  # Returns a named list of ggplot objects for the selected frames
  lapply(frames, function(sheet_name) {
    df <- temp_list[[sheet_name]] |>
      dplyr::arrange(Timestamp)
    
    # Derive window size from the median sampling interval
    interval_mins <- median(diff(as.numeric(df$Timestamp))) / 60
    window_size   <- round((24 * 60) / interval_mins)
    # rollmean requires an odd window for centered alignment
    if (window_size %% 2 == 0) window_size <- window_size + 1
    
    df <- df |>
      dplyr::mutate(
        # Centered 3-day rolling mean
        Temp_Smooth = zoo::rollmean(Temperature, k = window_size,
                                    fill = NA, align = "center")
      )
    
    ggplot(df, aes(x = Timestamp, y = Temp_Smooth)) +
      geom_line(linewidth = 0.6) +
      labs(
        title = sub("Frame_([A-F])_.*", "Frame \\1", sheet_name),
        x     = "Date / Time",
        y     = "Temperature (°C)"
      ) +
      theme_bw(base_family = "arial")
  }) |>
    setNames(frames)
}

#Plot temperature changes for all 6 frames in a 3x2 grid (A-F order)
all_plots <- plot_temperature(rawTempHum)

# Compute global y-axis range across all smoothed series
y_min <- min(sapply(all_plots, function(p) min(p$data$Temp_Smooth, na.rm = TRUE)))
y_max <- max(sapply(all_plots, function(p) max(p$data$Temp_Smooth, na.rm = TRUE)))

# Apply unified y range; remove y-axis title from right-column plots (even positions)
sorted_names <- sort(names(all_plots))
final_plots <- lapply(seq_along(sorted_names), function(i) {
  p <- all_plots[[sorted_names[i]]] +
    coord_cartesian(ylim = c(y_min, y_max))
  if (i %% 2 == 0) p <- p + theme(axis.title.y = element_blank())
  if (i <= 4)      p <- p + theme(axis.title.x = element_blank())
  p
})

wrap_plots(final_plots, nrow = 3, ncol = 2)

#Output the figure as a PDF
ggsave(
  filename = file.path('Figures_R', 'Frames_Temperatures.pdf'),
  plot     = wrap_plots(final_plots, nrow = 3, ncol = 2),
  width    = 12,
  height   = 12,
  units    = "in"
)
