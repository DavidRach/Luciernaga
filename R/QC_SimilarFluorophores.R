#' Queries fluorophore and returns similar fluorophores.
#'
#' @param TheFluorophore The name of the Fluorophore compare, see QC_ReferenceLibrary
#' @param NumberDetectors Number of detectors of the instrument
#' @param NumberHits Number of most similar fluorophores by cosine
#' @param returnSynonymns Returns only fluorophores > 0.98 cosine value, default FALSE
#' @param returnPlots Whether to also return signature plots, default is set FALSE
#' @param returnSynonyms Something
#' @param plotlinecolor Default NULL, otherwise if single line provide desired color
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
#' @importFrom tidyselect where
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
#' Results <- QC_SimilarFluorophores(TheFluorophore="Spark Blue 550",
#'  NumberDetectors=64, returnSynonymns=FALSE, NumberHits = 10, returnPlots=FALSE)

QC_SimilarFluorophores <- function(TheFluorophore, NumberDetectors,
   returnSynonyms=FALSE, NumberHits=10, returnPlots=FALSE, plotlinecolor=NULL) {

  ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
  #nrow(ReferenceData)
  #ReferenceData %>% pull(Fluorophore) %>% unique()
  #ReferenceData1 <- ReferenceData |> unique()
  #nrow(ReferenceData1)
  #ReferenceData1 %>% pull(Fluorophore) %>% unique()

  if (returnPlots == TRUE){ReferenceData1 <- ReferenceData}

  ReferenceData <- ReferenceData |> select(-Instrument) |>
    group_by(Fluorophore) |> pivot_wider(
      names_from = Detector, values_from = AdjustedY) |> ungroup()
  
  #RowNAs <- ReferenceData[rowSums(is.na(ReferenceData)) > 0, ]
  #nrow(RowNAs)

  #CleanNAs <- ReferenceData[rowSums(is.na(ReferenceData)) == 0, ]
  #nrow(CleanNAs)
  #ncol(CleanNAs)-1
  #View(CleanNAs)

  TheAvailableFluors <- ReferenceData |> pull(Fluorophore)
  if (!TheFluorophore %in% TheAvailableFluors) {stop("Fluorophore not found")}

  CombinedView <- ReferenceData
  Names <- CombinedView |> pull(Fluorophore)
  #SanitizedNames <- Luciernaga::NameCleanUp(Names, removestrings = c(",", "-", " ", "."))

  Numbers <- CombinedView %>% select(where(is.numeric))
  #zero_columns <- colSums(Numbers) == 0
  #zero_rows <- rowSums(Numbers) == 0
  #print(which(zero_columns))
  #print(which(zero_rows))
  
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names
  #colnames(NumericsT) <- SanitizedNames

  #zero_columns <- colSums(NumericsT) == 0
  #zero_rows <- rowSums(NumericsT) == 0
  #print(which(zero_columns))
  #print(which(zero_rows))

  CosineMatrix <- lsa::cosine(NumericsT)
  CosineMatrix <- round(CosineMatrix, 2)
  
  CosineFrame <- data.frame(CosineMatrix, check.names = FALSE)

  CosineFrame <- CosineFrame |> select(all_of(TheFluorophore))
  TheData <- rownames_to_column(CosineFrame, var="Fluorophore")
  TheID <- TheData |> select(all_of(TheFluorophore)) |> colnames()

  if (returnSynonyms == FALSE){
  TheHits <- TheData |> filter(!Fluorophore %in% TheID) |>
    arrange(desc(.data[[TheID]])) |> slice_head(n=NumberHits)
  } else {
    TheHits <- TheData |> filter(!Fluorophore %in% TheID) |>
    arrange(desc(.data[[TheID]])) |> filter(.data[[TheID]] > 0.98)
  }

  if (returnPlots==TRUE){
    TheseFluorophores <- TheHits |> pull(Fluorophore)

    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=TheFluorophore, data=ReferenceData1,
                                 plotlinecolor=plotlinecolor)
    ReturnThese <- list(TheHits, ThePlot)
    return(ReturnThese)
  } else {return(TheHits)}
  }


#' Internal for QC_SimilarFluorophores, returns plot of all similar fluorophores
#'
#' @param TheseFluorophores The similar fluorophores identified by cosine
#' @param TheFluorophore The one we were originally interested in
#' @param data The reference data of fluorophore signatures
#' @param legend Default TRUE, alternately removes plot legend
#' @param plotname Default NULL, alternately specifies the plot title
#' @param plotlinecolor Expects NULL, otherwise if single line provide desired color
#' @param unstained Default NULL, alternatively adds corresponding unstained signature
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom tidyselect where
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
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
#' @return An internal value
#'
#' @noRd
SimilarFluorPlots <- function(TheseFluorophores, TheFluorophore, data,
  legend=TRUE, plotname=FALSE, plotlinecolor, unstained=NULL){
      TheseFluorophores <- as.character(TheseFluorophores)
  
      if(!is.null(unstained)){
        These <- c(TheFluorophore, TheseFluorophores, "Unstained")

        if (any(unstained |> select(where(is.numeric)) > 1)){
          Metadata <- unstained |> select(!where(is.numeric))
          Numerics <- unstained |> select(where(is.numeric))
          n <- Numerics
          n[n < 0] <- 0
          A <- do.call(pmax, n)
          Normalized <- n/A
          unstained <- bind_cols(Metadata, Normalized)
       }

        UnstainedData <- unstained |> select(where(is.numeric)) |>
          pivot_longer(cols=everything(), names_to = "Detector", values_to="value") |>
          mutate(Instrument="Existing", Fluorophore="Unstained") |>
          relocate(Instrument, Fluorophore, .before="Detector")
      } else {These <- c(TheFluorophore, TheseFluorophores)}
      
      TheData <- data |> filter(Fluorophore %in% These) |>
        rename(value=AdjustedY)
  
      if(!is.null(unstained)){TheData <- bind_rows(TheData, UnstainedData)}
  
      TheData$Detector <- gsub("-A", "", TheData$Detector)

      Iterations <- TheData |> filter(Fluorophore %in% These[[1]]) %>% nrow()

      if (is.character(TheData$Detector)) {
        MyVector <- TheData |> filter(Fluorophore %in% These[[1]]) |> pull(Detector)
      }

      if (is.numeric(TheData$Detector)) {
        MyVector <- seq_len(Iterations)
      }

      if (any(TheData$value > 2)){YAxisLabel <- "Raw MFI"
      } else {YAxisLabel <- "Normalized Value"}

      TheData$Detector <- factor(TheData$Detector, levels=MyVector)
      TheData$Fluorophore <- factor(TheData$Fluorophore, levels=These)
  
      if (is.null(plotname)){TheTitle <- paste0(TheFluorophore)
      } else {TheTitle <- plotname}
  
      if (legend == TRUE){LegendLocation <- "right"
      } else {LegendLocation <- "none"}
  
      if (!is.null(plotlinecolor)){
        ThePlot <- ggplot(TheData, aes(x=Detector, y=value, group=Fluorophore)) +
          geom_line(color = plotlinecolor) + theme_bw() +
         labs(title=TheTitle, x=NULL, y=YAxisLabel) + geom_hline(yintercept = 1,
            linetype = "dashed", color = "red") +
         theme(plot.title = element_text(size = 8),
             legend.position = LegendLocation,
             axis.text.x = element_text(size = 6, angle = 45),
             panel.grid = element_blank(), axis.ticks.x = element_blank(),
             axis.title.y =  element_text(size=8)) +
       scale_x_discrete(breaks = unique(TheData$Detector)[c(TRUE, rep(FALSE, 4))])
      } else {
        ThePlot <- ggplot(TheData, aes(x=Detector, y=value, group=Fluorophore,
           color = Fluorophore)) + geom_line() + theme_bw() +
          labs(title=TheTitle, x=NULL, y=YAxisLabel) + geom_hline(yintercept = 1,
             linetype = "dashed", color = "red") +
          theme(plot.title = element_text(size = 8),
              legend.position = LegendLocation,
              axis.text.x = element_text(size = 6, angle = 45),
              panel.grid = element_blank(), axis.ticks.x = element_blank(),
              axis.title.y =  element_text(size=8)) +
        scale_x_discrete(breaks = unique(TheData$Detector)[c(TRUE, rep(FALSE, 4))])
        } 
  
      return(ThePlot)
}
