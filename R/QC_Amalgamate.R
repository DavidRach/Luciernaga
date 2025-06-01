#' Taking a Luciernaga_Slices (or similar data.frame) with an added Counts column, 
#' expands each row by the count, takes the median, and then returns either data
#' or a plot showing how different the average signature is from the individual ones
#' 
#' @param data The modified data.frame containing the sample, count and
#'  all the normalized detector columns
#' @param samplecolumn The column name correspondoing to the sample column,
#'  default is Percentiles
#' @param countcolumn The column name corresponding to the count column,
#'  default Count
#' @param returnType Whether to return plot (default) or data
#' @param titlename Default is NULL, specify a title if desired
#' @param linecolor Default is red, specify alternate color for the Average signature
#' @param legend Whether to return legend, default is TRUE
#' 
#' @importFrom tidyr uncount
#' @importFrom rlang sym
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom tidyselect where
#' @importFrom dplyr summarize
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr all_of
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme_linedraw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom stats setNames
#' 
#' @export
#' 
#' @return Either a plot or the underlying data
#' 
#' @examples
#' 
#' A <- 2+2
#' 
QC_Amalgamate <- function(data, samplecolumn="Percentiles", normalize=FALSE,
 countcolumn="Count", returnType="plot", titlename=NULL, linecolor="red", legend=TRUE){

Combined <- data |> tidyr::uncount(weights = !!sym(countcolumn), .remove = TRUE) |> 
     select(-all_of(samplecolumn)) |> summarize(across(everything(), median))  |>
     mutate(!!samplecolumn := "Average") |> relocate(all_of(samplecolumn), .before = 1)

ComparisonTable <- bind_rows(data, Combined)   
     
Numerics <- sapply(ComparisonTable, is.numeric)
ComparisonTable[Numerics] <- lapply(ComparisonTable[Numerics], round, digits = 3)
ComparisonTable <- ComparisonTable |> select(-all_of(countcolumn))
     
if (normalize == TRUE){
     if (any(ComparisonTable |> select(where(is.numeric)) > 1)){
      Metadata <- ComparisonTable |> select(!where(is.numeric))
      Numerics <- ComparisonTable |> select(where(is.numeric))
      n <- Numerics
      n[n < 0] <- 0
      A <- do.call(pmax, n)
      Normalized <- n/A
      ComparisonTable <- bind_cols(Metadata, Normalized)
     }
}
     
colnames(ComparisonTable) <- gsub("Comp-", "", colnames(ComparisonTable))

if (returnType == "data"){return(ComparisonTable)
} else{
     LineCols <- ncol(ComparisonTable)
     DetectorOrder <- colnames(ComparisonTable)
     DetectorOrder <- DetectorOrder[-1]
     
     Melted <- ComparisonTable |> 
          pivot_longer(all_of(2:LineCols),
          names_to = "Detector", values_to = "value")

     Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)
    
     if (!is.null(titlename)){name <- titlename
     } else {name <- NULL}

     if (legend == TRUE){LegendPosition <- "right"
     } else {LegendPosition <- "none"}

     Melted[[samplecolumn]] <- factor(Melted[[samplecolumn]],
          levels = c(setdiff(unique(Melted[[samplecolumn]]), "Average"), "Average"))

     if (normalize == TRUE){Expression <- "Normalized MFI"
     } else {Expression <- "MFI"}
     

     plot <- ggplot(Melted, aes(x = Detector, y = value, group = .data[[samplecolumn]],
               color = .data[[samplecolumn]])) + geom_line() +
               scale_color_manual(values = setNames(
                    ifelse(unique(Melted[[samplecolumn]]) == "Average", linecolor, "gray"),
                    unique(Melted[[samplecolumn]]))) +
               labs(title = name, x = "Detectors", y = Expression) +
               theme_linedraw() + theme_bw() + theme(axis.title.x = element_text(
               face = "plain"), axis.title.y = element_text(face = "plain"),
               axis.text.x = element_text(size = 5,
               angle = 45, hjust = 1), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
}

     return(plot)
}