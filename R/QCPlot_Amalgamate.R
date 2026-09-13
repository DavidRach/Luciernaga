#' Internal for QC_Screen, produced the amalgamate plots for the passed inputs
#' 
#' @param x The identifying name for which data will be filtered for
#' @param data The passed data 
#' @param countcolumn Default is Count
#' @param samplecolumn Default is TheSample
#' @param linecolor Default is "red", used to show Averaged Signature
#' @param returnType Default is "plot"
#' @param titlename Default is NULL, alternatively provided sets title
#' @param legend Default is FALSE, when TRUE sets legend.position to "right"
#' @param normalize Default is FALSE
#' 
#' @importFrom dplyr filter pull select summarize across relocate bind_rows
#'  bind_cols mutate
#' @importFrom tidyr uncount pivot_longer
#' @importFrom tidyselect where everything all_of
#' @importFrom rlang sym
#' @importFrom stats median setNames
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs theme_linedraw
#'  theme_bw theme element_text element_blank
#' 
#' @return An individual ggplot object
#' 
#' @noRd
QCPlot_Amalgamate <- function(x,
                               data,
                               countcolumn = "Count",
                               samplecolumn = "TheSample",
                               linecolor = "red",
                               returnType = "plot",
                               titlename = NULL,
                               legend = FALSE,
                               normalize) {

  SmallPortion <- data |> filter(TheSample %in% x)
  if (is.null(titlename)) {
    TheTitle <- SmallPortion |> pull(TheSample)
  }
  NumericPortion <- SmallPortion |>
    uncount(weights = !!sym(countcolumn), .remove = TRUE) |>
    select(where(is.numeric)) |>
    select(-Ratio)
  AverageCalced <- NumericPortion |>
    summarize(across(everything(), median)) |>
    mutate(Cluster = "Average") |>
    relocate(Cluster, .before = 1)
  Rest <- SmallPortion |>
    select(Cluster, where(is.numeric)) |>
    uncount(weights = !!sym(countcolumn), .remove = TRUE) |>
    select(-Ratio)
  ComparisonTable <- bind_rows(Rest, AverageCalced)

  Numerics <- sapply(ComparisonTable, is.numeric)
  ComparisonTable[Numerics] <- lapply(ComparisonTable[Numerics], round,
                                       digits = 3)

  if (normalize == TRUE) {
    if (any(ComparisonTable |> select(where(is.numeric)) > 1)) {
      Metadata <- ComparisonTable |> select(!where(is.numeric))
      Numerics <- ComparisonTable |> select(where(is.numeric))
      n <- Numerics
      n[n < 0] <- 0
      A <- do.call(pmax, n)
      Normalized <- n / A
      ComparisonTable <- bind_cols(Metadata, Normalized)
    }
  }

  colnames(ComparisonTable) <- gsub("Comp-", "", colnames(ComparisonTable))

  LineCols <- ncol(ComparisonTable)
  DetectorOrder <- colnames(ComparisonTable)
  DetectorOrder <- DetectorOrder[-1]

  Melted <- ComparisonTable |>
    pivot_longer(all_of(2:LineCols), names_to = "Detector",
                 values_to = "value")

  Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)

  if (!is.null(titlename)) {
    name <- titlename
  } else {
    name <- TheTitle
  }

  if (legend == TRUE) {
    LegendPosition <- "right"
  } else {
    LegendPosition <- "none"
  }

  Melted[["Cluster"]] <- factor(
    Melted[["Cluster"]],
    levels = c("Average", setdiff(unique(Melted[["Cluster"]]), "Average"))
  )

  if (normalize == TRUE) {
    Expression <- "Normalized MFI"
  } else {
    Expression <- "MFI"
  }

  plot <- ggplot(Melted, aes(x = Detector, y = value, group = Cluster,
                              color = Cluster)) +
    geom_line(alpha = 0.5, size = 0.2) +
    scale_color_manual(values = setNames(
      ifelse(unique(Melted[["Cluster"]]) == "Average", linecolor, "gray"),
      unique(Melted[["Cluster"]]))) +
    labs(title = name, x = "Detectors", y = Expression) +
    theme_linedraw() +
    theme_bw() +
    theme(axis.title.x = element_text(face = "plain"),
          axis.title.y = element_text(face = "plain"),
          axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = LegendPosition)

  return(plot)
}