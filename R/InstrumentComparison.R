#' Given either Instruments or Fluors, compares cosine values across instruments 
#' for all or selected fluorophores against the provided MainFluorophore
#' 
#' @param Instruments A list of instrument detectors for the respective instruments
#' @param MainFluorophore The desired fluorophore to compare against
#' @param returnType Whether to return "data" or a "plot"
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr select
#' @importFrom dplyr slice
#' @importFrom dplyr ungroup
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr bind_cols
#' 
#' @return A data.frame or in the future a plot object
#' 
#' @noRd
InstrumentComparison <- function(Instruments, MainFluorophore, returnType){

  Values <- map(.x=Instruments, .f=InstrumentNameCheck)
  Hmm <- unlist(Values)
  TheFluorophores <- map(.x=Hmm, .f=InstrumentReturn)
  Shared <- Reduce(intersect, TheFluorophores)
  Dataset <- map(.x=Hmm, .f=InstrumentData, fluorophores=Shared) |> bind_rows()
  #Dataset |> pull(Instrument) |> unique()
  
  if (!MainFluorophore %in% Shared){
      stop(MainFluorophore, " not found across selected instruments")
      }
  
  if(!is.null(Instruments)){
      OrderReference <-Hmm[1]
  } else {OrderReference <- 64} #Remember exception if not found 5L
  
  ReferenceData <- InstrumentReferences(NumberDetectors=OrderReference)
  TheReferenceList <- ReferenceData |> filter(Fluorophore %in% Shared)
  Shared_Rearranged <- TheReferenceList |> group_by(Fluorophore) |>
      arrange(desc(AdjustedY)) |> slice(1) |> select(Fluorophore, Detector) |>
      ungroup() |> arrange(Detector) |> pull(Fluorophore)
  
  LongReferences <- Dataset |>
       pivot_wider(names_from="Detector", values_from = "AdjustedY")
  
  instrument_match <- LongReferences |> pull(Instrument) |> unique()
  
  Forward <- map(.x=instrument_match, .f=MatchRearrange,
   thematch=Shared_Rearranged, data=LongReferences) |> bind_rows()
  # Forward |> pull(Instrument) |> unique()
  
  Comparison <- map(.x=instrument_match, .f=CosineReturn,
   MainFluorophore=MainFluorophore, data=Forward) |> bind_cols()
  
  if (returnType == "data"){
      return(Comparison)
      } else if (returnType == "plot"){
          message("Add gt or ggplot2 option here")
      } else {return(Comparison)}
  }

#' Internal for InstrumentComparison, handles name list,
#' still very preliminary
#' 
#' @param x The iterated list item to be checked
#' 
#' @return A formatted value to filter instruments with
#' 
#' @noRd
InstrumentNameCheck <- function(x){
  if (is.character(x)){
          if (grepl("\\d", x) && grepl("[A-Za-z]", x)){
              #message("Handle both")
          } else if (grepl("\\d", x)){
              #message("Converting numeric")
              TheY <- as.numeric(x)
          } else if (grepl("[A-Za-z]", x)){
              #message("Handle no letters")
          } else {stop("Instrument list item not character or numeric")}
      } else if (is.numeric(x)){
          #message("Numeric")
          TheY <- x
      } else {stop("Instrument list item not character or numeric")
      }
  return(TheY)
}

#' Internal for InstrumentComparison, pulls fluorophore list for respective instrument
#' 
#' @param x Number detectors corresponding desired instrument
#' 
#' @importFrom dplyr pull
#' 
#' @return The reference fluorophores for that instrument
#' 
#' @noRd
InstrumentReturn <- function(x){
  Data <- InstrumentReferences(NumberDetectors=x)
  Fluorophores <- Data |> dplyr::pull(Fluorophore) |> unique()
  return(Fluorophores)   
}

#' Internal for InstrumentComparison, filters for given fluorophores
#' info on a per instrument basis
#' 
#' @param x The number detectors corresponding to the instrument
#' @param fluorophores The vector of fluorophores to filter for
#' 
#' @importFrom dplyr filter
#' 
#' @return Data for the fluorophores
#' 
#' @noRd
InstrumentData <- function(x, fluorophores){
  Data <- InstrumentReferences(NumberDetectors=x)
  Data <- Data |> filter(Fluorophore %in% fluorophores)
  return(Data)
}

#' Internal for InstrumentComparison
#' 
#' @param x The instrument being filtered for
#' @param thematch The established laser ordered fluorophores
#' @param data The passed data for all instruments filtering from
#' 
#' @importFrom dplyr filter
#' 
#' @return Rearranged data.frame with laser ordered rows
#' 
#' @noRd
MatchRearrange <- function(x, thematch, data){
  Subset <- data |> filter(Instrument %in% x)
  MatchedReferences <- match(thematch, Subset$Fluorophore)
  MatchedReferences <- Subset[MatchedReferences, ]
  return(MatchedReferences)
}

#' Internal for InstrumentComparison, derrives the cosine data
#' for the fluors and instruments
#' 
#' @param x The instrument being filtered for
#' @param MainFluorophore The fluorophore the rest are being compared to
#' @param data The data.frame containing required data
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect where
#' 
#' @return A single column renamed for the instrument with
#'  comparisons as rows
#' 
#' @noRd
CosineReturn <- function(x, MainFluorophore, data){
  Subset <- data |> dplyr::filter(Instrument %in% x)
  Poised <- Subset |> select(-Instrument)
  Poised2 <- Poised %>% select(where(~ !all(is.na(.))))
  CosineReturn <- Luciernaga_Cosine(data=Poised2,
   returntype="data", rearrange = FALSE)
  CosineReturn <- round(CosineReturn, 2)
  Data <- data.frame(CosineReturn, check.names=FALSE)
  TheValues <- Data |> select(MainFluorophore)
  colnames(TheValues)[1] <- x
  return(TheValues)
}
