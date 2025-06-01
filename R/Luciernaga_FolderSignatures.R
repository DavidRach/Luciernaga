#' Takes a folder of Luciernaga_QC signatures and returns a
#'  data.frame of the signatures contained within, for use
#'  in downstream data analysis. 
#' 
#' @param FolderPath Location of the Luciernaga_QC .fcs outputs
#' @param sample.name The keyword where the identifying sample name can be found
#' @param StringRemoval Default NULL, provide to remove items from sample.name
#' based on values found on the keyword
#' @param fluorophore.name Specify the name of the fluorophore
#' @param Verbose Default FALSE, provides info as it goes
#' @param stats Whether to use median or mean
#' @param PanelCuts Default NULL, provide a c(0.5,1) argument to specify 
#' the brightness percentiles to retrieve signature from for the individual files
#' @param normalize Default TRUE, whether to return normalized or
#' raw averaged MFI signatures
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' 
#' @export
#' 
#' @return A data.frame of the corresponding averaged signatures per file
#' 
#' @examples
#' A <- 2 + 2
Luciernaga_FolderSignatures <- function(FolderPath, sample.name,
   StringRemoval=NULL, fluorophore.name, Verbose=FALSE,
    stats="median", PanelCuts=NULL, normalize=TRUE){

  TheFCSFiles <- list.files(path=FolderPath, pattern="fcs",
   full.names=TRUE) 
  Selected_CS <- load_cytoset_from_fcs(TheFCSFiles,
     truncate_max_range = FALSE, transformation = FALSE)
  Selected_GS <- GatingSet(Selected_CS)

  Returns <- map(.x=Selected_GS, .f=FolderSignatureIterator,
    sample.name=sample.name, StringRemoval=StringRemoval,
    fluorophore.name=fluorophore.name, Verbose=Verbose,
    stats=stats, PanelCuts=PanelCuts, normalize=normalize,
    returnType="Signatures") |> bind_rows()

  return(Returns)
}

#' Internal, processes individual files for signature
#' 
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom stats quantile
#' @importFrom dplyr bind_cols
#' 
#' @return A data.frame row of raw or normalized data
#' 
#' @noRd
FolderSignatureIterator <- function(x, sample.name, StringRemoval,
  fluorophore.name, Verbose, stats, PanelCuts, normalize, returnType){
  
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

  QuantileData <- TheColumns |> select(all_of(TheDetector)) |> pull()
  LowerBoundMFI <- QuantileData %>% quantile(., LowerBound)
  UpperBoundMFI <- QuantileData %>% quantile(., UpperBound)

  ValuesInterest <- TheColumns |>
    filter(.data[[TheDetector]]  >= LowerBoundMFI & .data[[TheDetector]] <= UpperBoundMFI)

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
    Dataset <- ValuesInterest |> mutate(Fluorophore = fluorophore.name) |>
      mutate(Sample=sampleName) |> relocate(Fluorophore, Sample, .before=1)
    
    return(Dataset)}
  }