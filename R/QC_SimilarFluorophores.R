#' Queries fluorophore and returns similar fluorophores.
#'
#' @param TheFluorophore The name of the Fluorophore compare, see QC_ReferenceLibrary
#' @param NumberDetectors Number of detectors of the instrument
#' @param NumberHits Number of most similar fluorophores by cosine
#' @param returnPlots Whether to also return signature plots, default is set FALSE
#'
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr ungroup
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select_if
#' @importFrom lsa cosine
#' @importFrom tibble rownames_to_column
#' @importFrom tidyselect starts_with
#' @importFrom dplyr arrange
#' @importFrom dplyr slice_head
#'
#' @returns A dataframe of similar fluorophores
#' @export
#'
#' @examples
#'
#' Results <- QC_SimilarFluorophores(TheFluorophore="Spark Blue 550",
#'  NumberDetectors=64, NumberHits = 10, returnPlots=TRUE)

QC_SimilarFluorophores <- function(TheFluorophore, NumberDetectors, NumberHits, returnPlots=FALSE) {

  ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)

  if (returnPlots == TRUE){ReferenceData1 <- ReferenceData}

  ReferenceData <- ReferenceData %>% select(-Instrument) %>%
    group_by(Fluorophore) %>% pivot_wider(
      names_from = Detector, values_from = AdjustedY) %>% ungroup()

  TheAvailableFluors <- ReferenceData %>% select(Fluorophore) %>% pull()
  if (!TheFluorophore %in% TheAvailableFluors) {stop("Fluorophore not found")}

  CombinedView <- ReferenceData
  Names <- CombinedView %>% select(Fluorophore) %>% pull()
  Numbers <- CombinedView %>% select_if(is.numeric)
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names

  CosineMatrix <- cosine(NumericsT)
  CosineMatrix <- round(CosineMatrix, 2)
  CosineFrame <- data.frame(CosineMatrix, check.names = FALSE)

  CosineFrame <- CosineFrame %>% select(all_of(TheFluorophore))
  TheData <- rownames_to_column(CosineFrame, var="Fluorophore")
  TheID <- TheData %>% select(all_of(TheFluorophore)) %>% colnames()

  TheHits <- TheData %>% filter(!Fluorophore %in% TheID) %>%
    arrange(desc(.data[[TheID]])) %>% slice_head(n=NumberHits)

  if (returnPlots==TRUE){
    TheseFluorophores <- TheHits %>% pull(Fluorophore)

    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=TheFluorophore, data=ReferenceData1)
    ReturnThese <- list(TheHits, ThePlot)
    return(ReturnThese)
  } else {return(TheHits)}
}


#' Internal for QC_SimilarFluorophores, returns plot of all similar fluorophores
#'
#' @param TheseFluorophores The similar fluorophores identified by cosine
#' @param TheFluorophore The one we were originally interested in
#' @param data The reference data of fluorophore signatures
#'
#' @importFrom dplyr filter
#' @importFrom dplyr rename
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_x_discrete
#'
#' @noRd
SimilarFluorPlots <- function(TheseFluorophores, TheFluorophore, data){

      These <- c(TheFluorophore, TheseFluorophores)

      TheData <- data %>% filter(Fluorophore %in% These) %>%
        rename(value=AdjustedY)

      Iterations <- TheData %>% filter(Fluorophore %in% These[[1]]) %>% nrow()

      if (is.character(TheData$Detector)) {
        MyVector <- TheData %>% filter(Fluorophore %in% These[[1]]) %>% pull(Detector)
      }

      if (is.numeric(TheData$Detector)) {
        MyVector <- 1:Iterations
      }

      if (any(TheData$value > 2)){YAxisLabel <- "Raw MFI"
      } else {YAxisLabel <- "Normalized Value"}

      TheData$Detector <- factor(TheData$Detector, levels=MyVector)
      TheData$Fluorophore <- factor(TheData$Fluorophore, levels=These)

      ThePlot <- ggplot(TheData, aes(x=Detector, y=value, group=Fluorophore, color = Fluorophore)) +
        geom_line() + theme_bw() + labs(title=paste0(TheFluorophore), x=NULL, y=YAxisLabel) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 6, angle = 45),
              panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.title.y =  element_text(size=8)) +
        scale_x_discrete(breaks = unique(TheData$Detector)[c(TRUE, rep(FALSE, 4))])

      return(ThePlot)
}
