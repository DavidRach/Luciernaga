#' Internal, processes individual files for signature
#' 
#' @importFrom flowCore keyword exprs
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr filter arrange desc pull select bind_cols
#' @importFrom tidyselect all_of
#' @importFrom stats quantile
#' 
#' @return A data.frame row of raw or normalized data
#' 
#' @noRd
FolderSignatureIterator <- function(x, sample.name, StringRemoval,
  fluorophore.name, Verbose, stats, PanelCuts, normalize, returnType){
  
  if (is.null(fluorophore.name)){
    FluorophoreName <- keyword(x, "TUBENAME")
    DefaultStrings <- c("DR_", " (Cells)")
    AbbreviatedFluorophore <- NameCleanUp(FluorophoreName,
       removestrings=DefaultStrings)
    fluorophore.name <- sub("^[^ ]+ ", "", AbbreviatedFluorophore)
  }
  
  if (length(sample.name) == 2){
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, sample.name)}

  if (!is.null(StringRemoval)){
    sampleName <- NameCleanUp(name, removestrings=StringRemoval)
  } else {
    DefaultStrings <- c("DR_", " (Cells)")
    sampleName <- NameCleanUp(name, removestrings=DefaultStrings)
  }

  if (Verbose == TRUE){
    message("After String Removal, sample.name is ", sampleName)
  }
  
  cs <- gs_pop_get_data(x, "root")
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)
  Data <- Data |> unique() #Precaution Zero Style Leftovers from Artificial
  TheColumns <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  DetectorOrder <- colnames(TheColumns)
  startingcells <- BiocGenerics::nrow(cs)[[1]]
  n <- TheColumns
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 1)
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  Normalized <- round(Normalized, 1)
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts),
   Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*0.0075
  Detectors <- PeakDetectorCounts |> filter(Counts > cutoff) |>
    arrange(desc(Counts))
  TheDetector <- Detectors[1,1]

  if(!is.null(PanelCuts)){
    PanelCuts <- PanelCuts
  } else {PanelCuts <- c(0,1)}

  LowerBound <- PanelCuts[1]
  UpperBound <- PanelCuts[2]

  if (!(LowerBound >= 0 & LowerBound <= 1)) {
    message("From should be between 0 and 1, proceeding to divide by
     100 on assumption it was a percentage")
    LowerBound <- LowerBound / 100
  }

  if (!(UpperBound >= 0 & UpperBound <= 1)) {
    message("To should be between 0 and 1, proceeding to divide by
     100 on assumption it was a percentage")
    UpperBound <- UpperBound / 100
  }

  QuantileData <- TheColumns |> select(all_of(TheDetector)) |>
    pull()
  LowerBoundMFI <- QuantileData %>% quantile(., LowerBound)
  UpperBoundMFI <- QuantileData %>% quantile(., UpperBound)

  ValuesInterest <- TheColumns |>
    filter(.data[[TheDetector]]  >= LowerBoundMFI &
      .data[[TheDetector]] <= UpperBoundMFI)

  if (returnType == "Signatures"){
  if (normalize == TRUE){
    Samples <- AveragedSignature(x=ValuesInterest, stats=stats,
      normalize = TRUE)
  } else {
    Samples <- AveragedSignature(x=ValuesInterest, stats=stats,
      normalize = FALSE)
  }

  Metadata <-data.frame(Fluorophore=fluorophore.name,
     Sample=sampleName, check.names=FALSE)
  Data <- bind_cols(Metadata, Samples)
  return(Data)
  } else {
    Dataset <- ValuesInterest |>
      mutate(Fluorophore = fluorophore.name) |>
      mutate(Sample=sampleName) |> 
      relocate(Fluorophore, Sample, .before=1)
    
    return(Dataset)}
}