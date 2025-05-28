#' Iterates over single color variant .fcs files, passing the 
#' generated signature to unmix full-stain sample
#'
#' @param IterativePath File path to folder containing the variant
#'  single color .fcs files
#' @param iterate_removestrings String values to remove from the 
#' variant .fcs files name
#' @param removestrings Regular string values to remove to match 
#' single color name in controlData
#' @param sample.name Keyword storage name
#' @param subset What population to pass from GatingSet, default 
#' is root.
#' @param PanelCuts A data.frame containing the cuts, 
#' see Luciernaga_Unmixing.
#' @param stats Whether to use "median" or "mean" to derrive 
#' average signatures
#' @param Verbose If TRUE returns steps
#' @param SignatureView Whether to visualize the SCs
#' @param returntype Default is set to "fcs" for iteration arguments
#' @param FullStainedGS The GatingSet containing the full-stained 
#' samples you want to see effect of unmixing on
#' @param controlData The Luciernaga_SingleColors output that will 
#' be used in combination with the iterated variant
#' @param outpath Desired storage location
#' @param PanelPath Path to a reference panel for ordering of 
#' column markers in unmixed file
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom Biobase pData
#' @importFrom purrr map
#'
#' @return Iterated versions of the full-stained sample
#' @export
#'
#' @examples NULL
Luciernaga_IterativeUnmixing <- function(IterativePath, IterativeSampleName,
   Iterativeremovestrings, fluorophore.name, PanelCuts=NULL, stats="median",
   samplecolumn="Sample", controlData, FullStainedGS, sample.name, removestrings,
   subset="root", returnType="fcs", outpath=outpath, PanelPath){
   
   if (!is.data.frame(IterativePath)){
   TheFCSFiles <- list.files(path=IterativePath, pattern="fcs", 
      full.names=TRUE)
   Dataset <- Luciernaga_FolderSignatures(
         FolderPath=IterativePath, sample.name=IterativeSampleName, 
         StringRemoval=removestrings, fluorophore.name=fluorophore.name,
         Verbose=FALSE, PanelCuts=PanelCuts, stats=stats, normalize=TRUE)
   } else {Dataset <- IterativePath}

   if (!is.data.frame(controlData)){
      controlData <- read.csv(controlData, check.names=FALSE)
   }

   if (any(controlData |> select(where(is.numeric)) > 1)){
      Metadata <- controlData |> select(!where(is.numeric))
      Numerics <- controlData |> select(where(is.numeric))
      n <- Numerics
      n[n < 0] <- 0
      A <- do.call(pmax, n)
      Normalized <- n/A
      controlData1 <- bind_cols(Metadata, Normalized)
   } else {controlData1 <- controlData}

   These <- Dataset |> pull(Sample)

   Stash <- map(.x=These, .f=IterativeUnmixingInternal,
      Iteration=Dataset, samplecolumn=samplecolumn, controlData=controlData1,
      FullStainedGS=FullStainedGS, sample.name=sample.name,
      removestrings=removestrings, subset=subset, returnType=returnType,
      outpath=outpath, PanelPath=PanelPath)
   
   return(Stash)
}

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
#' @param removestrings  A list of character strings to remove from sample.name
#' @param subset What population to retrieve from the GatingSet, default is root.
#' @param returnType Default is set to "data", alternate "flowframe" and "fcs"
#' @param outpath Desired storage location
#' @param PanelPath Path to a reference panel for ordering of column markers in unmixed file
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr bind_cols
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
   