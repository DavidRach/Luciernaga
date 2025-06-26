#' Queries the available reference library for available fluorophores (and their naming conventions)
#'
#' @param FluorNameContains A character string pattern to match, example "APC"
#' @param NumberDetectors The number of detectors the instrument contains.
#' For Cytek Instruments 5L = 64, 4L_UV = 54, 4L_YG = 48, 3L=38, 2L_VB=30,
#' 2L_BR=22, 1L=14
#' For BD S8=78, S6="48_S", A5="48_A"
#' For Sony ID7000 7L=184, 6L_DUV="182_DUV", 5L=147, 4L=112, 3L=86
#' For ThermoFisher BigFoot 7L_488-561=55, 7L_532-594="52_7L", 6L_445="52_6L", 6L_785=51 
#' @param returnPlots Whether to return signature plot as well. Default FALSE.
#' @param plotlinecolor Default NULL, otherwise if single line provide desired color
#' @param plotname Default NULL, alternatively specify a title.
#' @param exact Default FALSE, else returns exact fluorophore name match. 
#' @param unstained Default NULL, alternatively provide corresponding unstained signature
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stringr str_detect
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr slice
#' @importFrom dplyr n
#' @importFrom stringr str_extract
#'
#' @return A dataframe column containing matching Fluorophores from your querry
#' @export
#'
#' @examples
#' QC_ReferenceLibrary(FluorNameContains = "FITC", NumberDetectors=64)
QC_ReferenceLibrary <- function(FluorNameContains, NumberDetectors,
                                returnPlots=FALSE, plotlinecolor=NULL,
                                legend=TRUE, plotname=NULL, exact=FALSE,
                                unstained=NULL){
  
  if (!length(NumberDetectors) == 1){
    ReferenceData <- map(.x=NumberDetectors, .f=Luciernaga:::InstrumentReferences) |>
      bind_rows()

    # Pre-filtering for intersecting fluorophores across instruments
    GroupLength <- length(NumberDetectors)
    Intermediate <- ReferenceData |> group_by(Instrument) |> 
      select(Instrument, Fluorophore) |> unique() |> ungroup()
    These <- Intermediate |> group_by(Fluorophore) |>
      mutate(Count = dplyr::n()) |> ungroup() |> 
      filter(Count == GroupLength) |> pull(Fluorophore) |> unique()
    ReferenceData <- ReferenceData |> filter(Fluorophore %in% These) 
  } else {ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)}

  # Identifying matches
  TheList <- ReferenceData |> select(Fluorophore) |> unique()
  rownames(TheList) <- NULL    

  if (exact == FALSE){
  if (length(FluorNameContains) == 1){
  Subset <- TheList |> filter(str_detect(Fluorophore, FluorNameContains))
  } else {Subset <- TheList |> 
    filter(str_detect(Fluorophore, paste(FluorNameContains, collapse = "|")))}
  } else {
    if (length(FluorNameContains) == 1){
      Subset <- TheList |> filter(str_detect(Fluorophore, FluorNameContains))
      } else {Subset <- TheList |> 
        dplyr::filter(Fluorophore %in% FluorNameContains)
    }
  }

  TheseFluors <- Subset |> pull(Fluorophore)
  Locations <- ReferenceData |> filter(Fluorophore %in% TheseFluors) |>
    group_by(Fluorophore) |> arrange(desc(AdjustedY)) |> slice(1) |> ungroup()
  Order <- c("UV", "V", "B", "YG", "R")
  Sequence <- Locations |> mutate(
    prefix = str_extract(Detector, "^[A-Z]+"),
    num = as.numeric(str_extract(Detector, "\\d+")),
    group_order = match(prefix, Order)
  ) |> arrange(group_order, num) |> pull(Fluorophore)

  Subset$Fluorophore <- factor(Subset$Fluorophore, levels=Sequence)
  Subset <- Subset |> arrange(Fluorophore)


  # Plotting if requested
  if (returnPlots==FALSE){
    return(Subset)
  } else {
    TheseFluorophores <- Subset |> pull(Fluorophore)
    if (!length(NumberDetectors) == 1){
      Instruments <- ReferenceData |> pull(Instrument) |> unique()
      ThePlot <- map(.x=Instruments, .f=SmallWrapper, data=ReferenceData, 
      TheseFluorophores=TheseFluorophores, unstained=unstained)  
    } else {
      ThePlot <- Luciernaga:::SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                   TheFluorophore=NULL, data=ReferenceData,
                                   plotlinecolor=plotlinecolor, legend=legend,
                                   plotname=plotname, unstained=unstained)  
    }
    ReturnThese <- list(Subset, ThePlot)
    return(ReturnThese)
  } 
}


#' Internal for QC_ReferenceLibrary
#' 
#' @importFrom dplyr filter
#' 
#' @noRd
SmallWrapper <- function(x, data, TheseFluorophores, unstained){
    Internal <- data |> filter(Instrument %in% x)
    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
      TheFluorophore=NULL, data=Internal, plotlinecolor=plotlinecolor,
      legend=legend, plotname=plotname, unstained=unstained)
  }
