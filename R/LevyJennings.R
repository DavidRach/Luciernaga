
#' Internal for CytekQCPlots
#'
#' @param x  The passed column name to be plotted
#' @param FailedFlag Whether to show red flags when a detectors 
#' "Out-of-Range" is TRUE
#' @param xValue The x-axis column plotted by, default is "DateTime"
#' @param TheData The passed data.frame from which to retrieve data
#' @param Metadata The optional column name to be used for "comparison"
#' @param plotType Whether to look "individual" or "comparison"
#' @param EngineerVisits Passed data.frame of engineer visits for
#'  vertical lines, default NULL.
#'
#' @importFrom dplyr select
#' @importFrom tidyr starts_with
#' @importFrom stringr str_detect
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme labs
#'  theme_bw scale_color_manual scale_fill_manual scale_shape_manual
#'  scale_size_manual geom_vline
#' @importFrom lubridate ymd_hms
#'
#' @return The pdf and/the plots.
#'
#' @noRd
LevyJennings <- function(x, FailedFlag, xValue, TheData, Metadata,
                         plotType, YAxisLabel, EngineerVisits=NULL){

  yValue <- x

  # Select Equivalent Flag Column
  if (FailedFlag == TRUE){
    FlagValue <- paste0("Flag-", yValue)
    FlagColumn <- TheData |> select(starts_with(FlagValue)) %>% colnames(.)
    if (length(FlagColumn) >1){
      NewFlagColumn <- str_detect(FlagColumn, paste0("^", FlagValue, "$"))
      FlagColumn <- FlagColumn[which(NewFlagColumn)]
    } else (FlagColumn <- FlagColumn)
  }

  if (str_detect(yValue, "^UV\\d{1,2}-[A-Za-z]+(_Gain)?$")){
    mycolor <- "purple"
  } else if (str_detect(yValue, "^UV\\d{1,3}[-_]")){
    mycolor <- "purple"
  } else if (str_detect(yValue, "^V\\d{1,3}[-_]")){
    mycolor <- "violet"
  } else if (str_detect(yValue, "^V\\d{1,2}-[A-Za-z]+(_Gain)?$")){
    mycolor <- "violet"
  } else if (str_detect(yValue, "^B\\d{1,3}[-_]")){
    mycolor <- "blue"
  } else if (str_detect(yValue, "^B\\d{1,2}-[A-Za-z]+(_Gain)?$")){
    mycolor <- "blue"
  } else if (str_detect(yValue, "^Y\\d{1,3}[-_]")){
    mycolor <- "darkgreen"
  } else if (str_detect(yValue, "^YG\\d{1,3}[-_]")){
    mycolor <- "darkgreen"
  } else if (str_detect(yValue, "^YG\\d{1,2}-[A-Za-z]+(_Gain)?$")){
    mycolor <- "darkgreen"
  } else if (str_detect(yValue, "^R\\d{1,3}[-_]")){
    mycolor <- "darkred"
  } else if (str_detect(yValue, "^R\\d{1,2}-[A-Za-z]+(_Gain)?$")){
    mycolor <- "darkred"
  } else if (str_detect(yValue, "^Change_UV\\d{1,2}(-[A-Za-z% ]+)?$")){
    mycolor <- "purple"
  } else if (str_detect(yValue, "^Change_V\\d{1,2}(-[A-Za-z% ]+)?$")){
    mycolor <- "violet"
  } else if (str_detect(yValue, "^Change_B\\d{1,2}(-[A-Za-z% ]+)?$")){
    mycolor <- "blue"
  } else if (str_detect(yValue, "^Change_YG\\d{1,2}(-[A-Za-z% ]+)?$")){
    mycolor <- "darkgreen"
  } else if (str_detect(yValue, "^Change_R\\d{1,2}(-[A-Za-z% ]+)?$")){
    mycolor <- "darkred"
  } else {mycolor <- "black"}


  if (plotType == "individual"){

    if (FailedFlag == TRUE){

      if (any(grepl("FALSE", TheData[[FlagColumn]]))) {
        shape_qc <- c("FALSE" = 21, "TRUE" = 22)
        fill_qc <- c("FALSE" = mycolor, "TRUE" = "red")
        size_qc <- c("FALSE" = 1, "TRUE" = 3)

      } else {
        shape_qc <- c("False" = 21, "True" = 22)
        fill_qc <- c("False" = mycolor, "True" = "red")
        size_qc <- c("False" = 1, "True" = 3)
      }

      Plot <- ggplot(TheData, aes(x=.data[[xValue]], y = .data[[yValue]]))  +
        geom_line(color = mycolor, linewidth = 1) +  geom_point(aes(
        shape = .data[[FlagColumn]], size = .data[[FlagColumn]],
        fill = .data[[FlagColumn]])) +
        scale_shape_manual(values = shape_qc) +
        scale_fill_manual(values = fill_qc) +
        scale_size_manual(values = size_qc) +
        labs(title = yValue, x = NULL, y = YAxisLabel) + theme_bw() +
        theme(legend.position = "none")

    } else {Plot <- ggplot(TheData,
       aes(x=.data[[xValue]], y = .data[[yValue]],
            color = mycolor)) + geom_line(color = mycolor) + geom_point(
            color = mycolor) +
      labs(title = yValue, x = NULL, y = YAxisLabel) +
            theme_bw() + theme(legend.position = "none")

    }
  }

  if (plotType == "comparison"){
    if (mycolor != "black"){VariantColor <- c("black", mycolor)
    } else {VariantColor <- c("gray", mycolor)}

    Plot <- ggplot(TheData,
       aes(x=.data[[xValue]], y = .data[[yValue]], group = .data[[Metadata]],
      color=.data[[Metadata]])) + geom_line(aes(color = .data[[Metadata]])) +
      geom_point(aes(
        color = .data[[Metadata]])) +
      scale_color_manual(values = VariantColor) +
      labs(title = yValue, x = NULL, y = YAxisLabel) +
      theme(legend.position = "none") + theme_bw()
  }

  if (is.null(EngineerVisits)){
    Plot1 <- Plot
  } else {
    Plot1 <- Plot + geom_vline(xintercept = as.POSIXct(EngineerVisits),
                              color = "red", linetype = "dashed")
  }

  return(Plot1)
}