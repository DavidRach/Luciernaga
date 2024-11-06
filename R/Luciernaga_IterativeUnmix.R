#' Iterates over single color variant .fcs files, passing the generated signature to unmix full-stain sample
#'
#' @param IterativePath File path to folder containing the variant single color .fcs files
#' @param iterate_removestrings String values to remove from the variant .fcs files name
#' @param removestrings Regular string values to remove to match single color name in controlData
#' @param sample.name Keyword storage name
#' @param subset What population to pass from GatingSet, default is root.
#' @param PanelCuts A data.frame containing the cuts, see Luciernaga_Unmixing.
#' @param stats Whether to use "median" or "mean" to derrive average signatures
#' @param Verbose If TRUE returns steps
#' @param SignatureView Whether to visualize the SCs
#' @param returntype Default is set to "data" for iteration arguments
#' @param FullStainedGS The GatingSet containing the full-stained samples you want to see effect of unmixing on
#' @param controlData The Luciernaga_SingleColors output that will be used in combination with the iterated variant
#' @param multiplier Default is 50000, amplifies the OLS output. Impacted by biexponential transformation
#' @param outpath Desired storage location
#' @param PanelPath Path to a reference panel for ordering of column markers in unmixed file
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
Luciernaga_IterativeUnmixing <- function(IterativePath, iterate_removestrings, removestrings, sample.name,
                                         subset="root", PanelCuts, stats, Verbose, SignatureView,
                                         returntype="data", FullStainedGS, controlData, multiplier,
                                         outpath, PanelPath){
  #Locate the FCS Files to Iterate Over
  TheFCSFiles <- list.files(path=IterativePath, pattern="fcs", full.names=TRUE)

  #Incorporate into a Gating Set
  Selected_CS <- load_cytoset_from_fcs(TheFCSFiles, truncate_max_range = FALSE, transformation = FALSE)
  Selected_GS <- GatingSet(Selected_CS)

  #Set up NewNames
  pd <- pData(Selected_GS)
  NewName <- NameCleanUp(name=pd$name, removestrings = iterate_removestrings)
  NewName <- data.frame(NewName)
  pd <- cbind(pd, NewName)
  pData(Selected_GS) <- pd

  #x <- Selected_GS[1]

  Stash <- map(.x=Selected_GS, .f=IterativeUnmixingInternal, sample.name=sample.name,
               removestrings=removestrings, subset=subset, PanelCuts=PanelCuts, stats=stats,
               Verbose=Verbose, SignatureView=SignatureView, returntype=returntype,
               FullStainedGS=FullStainedGS, controlData=controlData, multiplier=multiplier, outpath=outpath,
               PanelPath=PanelPath)

}

#' Internal for Luciernaga_IterativeUnmix
#'
#' @param x The iterated variant fcs file
#' @param sample.name Keyword storage name
#' @param removestrings  Regular string values to remove to match single color name in controlData
#' @param subset What population to pass from GatingSet, default is root.
#' @param PanelCuts A data.frame containing the cuts, see Luciernaga_Unmixing.
#' @param stats Whether to use "median" or "mean" to derrive average signatures
#' @param Verbose If TRUE returns steps
#' @param SignatureView Whether to visualize the SCs
#' @param returntype Default is set to "data" for iteration arguments
#' @param FullStainedGS The GatingSet containing the full-stained samples you want to see effect of unmixing on
#' @param controlData The Luciernaga_SingleColors output that will be used in combination with the iterated variant
#' @param multiplier Default is 50000, amplifies the OLS output. Impacted by biexponential transformation
#' @param outpath Desired storage location
#' @param PanelPath Path to a reference panel for ordering of column markers in unmixed file
#'
#' @importFrom dplyr anti_join
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#'
#' @return Passes unmixed full-stain .fcs files to the outpath
#'
#' @noRd
IterativeUnmixingInternal <- function(x, sample.name, removestrings, subset, PanelCuts, stats,
                                      Verbose, SignatureView, returntype, FullStainedGS, controlData,
                                      multiplier, outpath, PanelPath){

  SC_Reference <- Luciernaga_SingleColors(x=x, sample.name=sample.name, removestrings=removestrings,
                                          subset=subset, PanelCuts=PanelCuts, stats=stats,
                                          Verbose=Verbose, SignatureView=SignatureView,
                                          returntype="data")

  NotOverlapped <- controlData %>% anti_join(SC_Reference, by = c("Fluorophore", "Ligand"))
  controlData1 <- bind_rows(SC_Reference, NotOverlapped)

  AddOn <- pData(x)$NewName
  TheAddOn <- paste0("_", AddOn, "_Unmixed")

  UnmixSuccess <- map(.x=FullStainedGS, .f=Luciernaga_Unmix, controlData=controlData1, sample.name=TheSampleName,
                      addon=TheAddOn, subset=subset, removestrings=removestrings, multiplier=multiplier,
                      outpath=outpath, PanelPath=PanelPath, Verbose=FALSE)
}
