#' Internal for BrightnessIndexIterator (also an internal). When handed the gated flowframe,
#' proceeds to extract the correct positive and negative gate values for the respective single-color
#' and process the staining index. 
#' 
#' @param x The individual GatingHierarchy being iterated in
#' @param inverse.transform Whether to reverse the transformation 
#' before calculating the staining index. 
#' @param stringAppend Default is -A
#' 
#' @importFrom flowWorkspace pData gs_pop_get_data
#' @importFrom dplyr pull summarise where across select
#' @importFrom tidyselect all_of
#' @importFrom stats quantile
#' @importFrom stringr str_extract str_match str_remove_all
#' @importFrom flowCore exprs
#' 
StainIndexLocal <- function(x, inverse.transform, stringAppend){

     NameString <- pData(x) |> pull(name)
     TheAbsentFluorophore <- str_extract(NameString, "(?<=_No).*(?=Unmixed\\.fcs)")
     TheFluorophore <- str_match(NameString, "^\\S+ (.*?) \\(Beads\\)")[,2]

     FluorNameCheck <- paste0(TheFluorophore, stringAppend)

     if (!FluorNameCheck %in% colnames(x)){
          normalize <- function(s){str_remove_all(s, "\\s+")}
          TheFluorophoreNorm <- str_remove_all(FluorNameCheck, "\\s+")
          colnamesNorm <- str_remove_all(colnames(x), "\\s+")
          if (TheFluorophoreNorm %in% colnamesNorm){
               index <- which(colnamesNorm %in% TheFluorophoreNorm)
               TheFluorophore <- colnames(x)[index]
               TheFluorophore <- gsub(stringAppend, "", TheFluorophore)
          } else {warning("Can't match naming for ", FluorNameCheck)}
     }
     
     PositiveGate <- paste0(TheFluorophore, stringAppend, "+")
     NegativeGate <- paste0(TheFluorophore, stringAppend, "-")

     Positive <- gs_pop_get_data(x, PositiveGate, inverse.transform = inverse.transform)
     Positive <- exprs(Positive[[1]])
     Positive <- data.frame(Positive, check.names=FALSE)
     Positive <- Positive[,-grep(excludeThese, names(Positive))]
     TheMedian <- Positive |>
     summarise(across(where(is.numeric), \(x) quantile(x, probs = 0.95, na.rm = TRUE)))
     KeptMarkers <- colnames(TheMedian)
     Fluorophore <- names(TheMedian)[which.max(TheMedian)]
     Positive <- Positive |> select(all_of(Fluorophore))

     Negative <- gs_pop_get_data(x, NegativeGate, inverse.transform = inverse.transform)
     Negative <- exprs(Negative[[1]])
     Negative <- data.frame(Negative, check.names=FALSE)
     Negative <- Negative[,-grep(excludeThese, names(Negative))]
     Negative <- Negative |> select(all_of(Fluorophore))
  
     PositiveVals <- Positive |> pull(1)
     MFI_Pos <- quantile(PositiveVals, probs=0.5, na.rm =TRUE)
     MFI_Pos <- MFI_Pos[[1]]
     NegativeVals <- Negative |> pull(1)
     MFI_Neg <- quantile(NegativeVals, probs=0.5, na.rm =TRUE)
     MFI_Neg <- MFI_Neg[[1]]
     #hist(NegativeVals)
     MinNeg <- quantile(NegativeVals, probs=0.03, na.rm =TRUE)
     MinNeg <- MinNeg[[1]]
     MaxNeg <- quantile(NegativeVals, probs=0.97, na.rm =TRUE)
     MaxNeg <- MaxNeg[[1]]
     RSD <- (MaxNeg - MinNeg)/3.29
     StainingIndex <- (MFI_Pos - MFI_Neg)/(2*RSD)

     TheData <- data.frame(FMO=TheAbsentFluorophore,
                         Fluorophore=TheFluorophore,
                         StainIndex=StainingIndex)
     
     return(TheData)
}