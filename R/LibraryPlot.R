#' Internal, plots QC_LibraryParse data into plots
#'
#' @param x The passed data
#' @param references Plot argument, adds red reference signature
#' @param myfactor Plot argument, data column to group by for plotting. Default "Fluorophore".
#' @param namefactor Plot argument, data column name added to Plot Title.
#'
#' @importFrom dplyr pull rename select filter
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect where
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs
#'  theme_bw theme element_text
#' @importFrom rlang sym !!
#'
#' @return A ggplot2 object
#'
#' @noRd
LibraryPlot <- function(x, references=TRUE, myfactor="Fluorophore",
 namefactor="Sample"){

  Data <- x

  if (nrow(Data) == 1){Sample <- Data %>% pull(.data[[namefactor]])
  } else {Sample <- ""}

  NumberDetectors <- sum(sapply(Data, is.numeric))
  TheFluorophore <- Data |> pull(Fluorophore)
  TheFluorophore <- unique(TheFluorophore)
  TheDetectors <- colnames(Data)[sapply(Data, is.numeric)]

  Data <- Data |>
    pivot_longer(cols = where(is.numeric),
                 names_to = "Detector", values_to = "TheValue")

  Data$Detector <- factor(Data$Detector, levels=TheDetectors)


  ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors)
  ReferenceData <- ReferenceData |> rename(TheValue = "AdjustedY")
  ReferenceFluorList <- ReferenceData |> select(Fluorophore) %>%
    unique() %>% pull()

  ThePlot <- ggplot(Data, aes(x=Detector,
     y=TheValue, group=.data[[myfactor]])) + geom_line()  +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title=paste0(TheFluorophore, " ", Sample), y="Normalized") +
    theme_bw() + theme(plot.title = element_text(size = 8),
                       axis.title.y =  element_text(size=8),
                       axis.text.x = element_text(size=6,
                       angle = 70, hjust = 1))

  if (any(ReferenceFluorList == TheFluorophore)){

    if (references == TRUE){
      ReferenceData1 <- ReferenceData |>
        filter(Fluorophore %in% TheFluorophore) |>
        mutate(StandIn="Reference")
      ReferenceData1 <- ReferenceData1 %>%
        rename(!!sym(myfactor) := StandIn)
      ReferenceData1$Detector <- TheDetectors
      ThePlot <- ThePlot + geom_line(data = ReferenceData1,
         aes(x=Detector, y=TheValue, group=.data[[myfactor]]),
         color="red")
      return(ThePlot)
    }
  } else {
      return(ThePlot)
    }
}