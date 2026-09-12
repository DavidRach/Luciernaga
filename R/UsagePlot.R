#' Internal for Dashboard, takes processed ApplicationLog from Cytek Aurora
#' instruments and generates a UserUsage plot for respective instrument. 
#' 
#' @param data The Data derrived from AppParse filtered for SitFlushes
#' @param TheInstrument Instrument designation in the Instrument column, example "5L"
#' @param returnType Default "ByHour", alternatively "ByFifteen"
#' @param desiredfill Default black, accepts specification for other color
#'  "lightgray" etc
#' @param textsizey Yaxis text size
#' @param textsizex Xaxis text size
#' 
#' @importFrom dplyr filter
#' @importFrom lubridate wday
#' @importFrom lubridate hour
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 lims
#' @importFrom lubridate floor_date
#' @importFrom dplyr slice
#' @importFrom lubridate hm
#' 
#' @return A ggplot2 object
#' 
#' @noRd
UsagePlot <- function(data, TheInstrument, returnType="ByHour",
 desiredfill="black", textsizey=10, textsizex=10){
  data <- data |> filter(Instrument %in% TheInstrument)
  
  data$DayOfWeek <- wday(data$DateTime, label = TRUE, abbr = TRUE)
  data$DayOfWeek <- factor(data$DayOfWeek,
   levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))

  if (returnType == "ByHour"){
  data$Hour <- hour(data$DateTime)

  dataByHourDay <- data |>
  group_by(DayOfWeek, Hour) |>
  summarise(count = n(), .groups = "drop")
  
  MaxCount <- dataByHourDay |> arrange(desc(count)) |>
    slice(1) |> pull(count)
  MaxCount <- MaxCount*1.05
  MaxCount <- round(MaxCount, 0)

  #TheTitle <- paste0("Aurora ", TheInstrument, " All Time Usage")

  plot <- ggplot(dataByHourDay, aes(x = Hour, y = count)) +
    geom_col(fill = desiredfill) + 
    labs(title = NULL,
        x = "Hour of Day",
        y = "Sit Flushes") +
    facet_grid(DayOfWeek ~ ., scales = "free_y") +
    theme_bw() +
    scale_x_continuous(breaks = 0:23) + lims(y=c(NA, MaxCount)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = textsizex),
    axis.text.y = element_text(size = textsizey),
    strip.text.y = element_text(angle = 0))
  } else {
    data$TimeBin <- floor_date(data$DateTime, "15 minutes")
    data$HourMinute <- format(data$TimeBin, "%H:%M")

    dataBy15MinDay <- data |>
      group_by(DayOfWeek, HourMinute) |>
      summarise(count = n(), .groups = "drop")
      
    MaxCount <- dataBy15MinDay |> arrange(desc(count)) |>
        slice(1) |> pull(count)
    MaxCount <- MaxCount*1.05
    MaxCount <- round(MaxCount, 0)

    #TheTitle <- paste0("Aurora ", TheInstrument, " All Time Usage")

    dataBy15MinDay$TimeNum <- as.numeric(hm(dataBy15MinDay$HourMinute))/3600

    plot <- ggplot(dataBy15MinDay, aes(x = TimeNum, y = count)) +
    geom_col(fill = desiredfill, width = 0.25) + 
    labs(title = NULL, x = NULL, y = NULL) +
    facet_grid(DayOfWeek ~ ., scales = "free_y") +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 23.75, by = 1),  
                       labels = function(x) sprintf("%02d:%02d", floor(x), (x-floor(x))*60)) +
    lims(y=c(NA, MaxCount)) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = textsizex),
    axis.text.y = element_text(size = textsizey),
    strip.text.y = element_text(angle = 0))
  }
  return(plot)
}