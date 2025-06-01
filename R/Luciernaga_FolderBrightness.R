#' Sister function to FolderSignature, returns a density plot of brightness
#'
#' @param FolderPath Location where the Luciernaga_QC .fcs file
#' outputs are stored
#' @param sample.name Keyword in .fcs corresponding samples name
#' @param StringRemoval Default NULL, remove these from the sample name
#' @param fluorophore.name The name of the fluorophore
#' @param Verbose Default FALSE
#' @param stats Default is median
#' @param PanelCuts Default is NULL
#' @param normalize Default is FALSE, will scale to peak detector
#' @param detector The detector corresponding to the fluorophore,
#' sets the x-axis
#' @param returnType Default data, providing plot returns a 
#' Luciernaga_Brightness plot
#' @param maxtik Default 1e6, for plot return
#' @param legend Default is right, set to none to remove.
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' 
#' @export
#' 
#' @return Either a data.frame of brightness values, or a plot
#' 
#' @examples
#' A <- 2+2
#' 
Luciernaga_FolderBrightness <- function(FolderPath, sample.name,
  StringRemoval=NULL, fluorophore.name, Verbose=FALSE,
   stats="median", PanelCuts=NULL, normalize=FALSE, detector,
  returnType = "data", maxtik=1e6, legend="right"){

  TheFCSFiles <- list.files(path=FolderPath, pattern="fcs",
   full.names=TRUE) 
  Selected_CS <- load_cytoset_from_fcs(TheFCSFiles,
     truncate_max_range = FALSE, transformation = FALSE)
  Selected_GS <- GatingSet(Selected_CS)

  Returns <- map(.x=Selected_GS, .f=FolderSignatureIterator,
    sample.name=sample.name, StringRemoval=StringRemoval,
    fluorophore.name=fluorophore.name, Verbose=Verbose,
    stats=stats, PanelCuts=PanelCuts, normalize=FALSE,
    returnType="Brightness") |> bind_rows()
  
  if (returnType == "data"){
    return(Returns)
  } else {
    BrightnessPlot <- Luciernaga_Brightness(data=Returns,
      fluorophore.name=fluorophore.name,
      fluorophore.column="Fluorophore", cluster.column="Sample",
      downsample=TRUE, subsample = NULL, detector=detector,
      reference=NULL, clearance=0.02, Scaled = TRUE, maxtik=maxtik,
      legend=legend)    
    return(BrightnessPlot)
  }
}
