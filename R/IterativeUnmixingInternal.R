#' Internal for Luciernaga_IterativeUnmix
#'
#' @param x The name of the Iterated Signature data to filter for
#' @param Iteration The Iterated Signature data.frame product of
#'  Luciernaga_FolderSignatures
#' @param controlData The Luciernaga_SingleColors output that will
#'  be used for unmixing after swapping in the variant
#' @param samplecolumn Default "Sample", otherwise name of column in
#' Iterated data denoting the sample to filter for. 
#' @param FullStainedGS The GatingSet containing the full-stained samples
#'  you want to iteratively unmix
#' @param sample.name Keyword storing the samples name
#' @param removestrings  A list of character strings to remove from
#'  sample.name
#' @param subset What population to retrieve from the GatingSet,
#'  default is root.
#' @param returnType Default is set to "data", alternate "flowframe"
#'  and "fcs"
#' @param outpath Desired storage location
#' @param PanelPath Path to a reference panel for ordering of column
#'  markers in unmixed file
#'
#' @importFrom dplyr filter pull select bind_cols
#' @importFrom tidyselect where
#' @importFrom purrr map
#'
#' @return Passes unmixed full-stain .fcs files to the designated outpath
#'
#' @noRd
IterativeUnmixingInternal <- function(x, Iteration, samplecolumn,
   controlData, FullStainedGS, sample.name, removestrings, subset,
   returnType, outpath, PanelPath){

   RowInterest <- Iteration |> filter(.data[[samplecolumn]] %in% x)
   TheFluorophore <- RowInterest |> pull(Fluorophore)
   Assemble <- RowInterest |> select(where(is.numeric))
   Index <- which(controlData$Fluorophore %in% TheFluorophore)
   Metadata <- controlData[Index,] |> select(!where(is.numeric))
   Assembled <- bind_cols(Metadata, Assemble)
   controlData[Index,] <- Assembled

   TheAddOn <- paste0("_", TheFluorophore, "_", x, "_Unmixed")

   UnmixSuccess <- map(.x=FullStainedGS, .f=Luciernaga_Unmix,
      controlData=controlData, sample.name=sample.name, 
      addon=TheAddOn, subset=subset, removestrings=removestrings,
      outpath=outpath, PanelPath=PanelPath, Verbose=FALSE,
      returnType=returnType)
   
   return(UnmixSuccess)
}