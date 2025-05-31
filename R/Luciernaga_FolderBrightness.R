#' Sister function to FolderSignature, returns a density plot of brightness
Luciernaga_FolderBrightness <- function(FolderPath, sample.name,
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
    stats=stats, PanelCuts=PanelCuts, normalize=FALSE,
    returnType="Brightness") |> bind_rows()
  
  return(Returns)
}
