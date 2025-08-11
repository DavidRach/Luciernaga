#' Check the transformation settings to visualize .fcs files
#' 
#' @param x Iterated in GatingSet object
#' @param thesecolumns Choice of channels to transform
#' @param TransformationChoice Default flowjo_biexp_trans
#' @param channelRange Argument
#' @param maxValue Argument
#' @param pos Argument
#' @param neg Argument
#' @param widthBasis Argument
#' @param Multiple Argument
#' @param returnType Argument
#' @param optionalY Something
#' @param optionalX Something
#' @param optionalOutpath Something
#' @param optionalBins Something
#' @param optionalName Something
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom purrr map
#' 
#' @return A visualization to check transformation setting
#' 
#' @export
#' 
#' @examples
#' path <- system.file("extdata", package = "Luciernaga")
#' files <- list.files(path, pattern=".fcs", full.names=TRUE)
#' 
TransformationCheck <- function(x, thesecolumns=NULL, 
  TransformationChoice="flowjo_biexp", channelRange=256,
  maxValue=1000000, pos=4.5, neg=1, widthBasis=-500, Multiple=FALSE, 
  returnType="NxN", optionalY="BUV805-A", optionalX="BUV496-A",
  optionalOutpath=NULL, optionalBins=100, optionalName=NULL){

  if (class(x) == "character"){
  MyCytoSet <- load_cytoset_from_fcs(x, transformation=FALSE,
  truncate_max_range = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  }

  Markers <- colnames(MyCytoSet)
  KeptMarkers <- Markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", Markers)]

  if (!is.null(thesecolumns)){Retained <- KeptMarkers %in% thesecolumns
  } else {Retained <- KeptMarkers}

  InternalGS <- MyGatingSet

  if (Multiple==FALSE){
      Returns <- map(.x=InternalGS, .f=InternalTransformation, TransformationChoice=TransformationChoice, channelRange=channelRange,
      maxValue=maxValue, pos=pos, neg=neg, widthBasis=widthBasis,
      KeptMarkers=KeptMarkers)
      Flag <- FALSE
  } else {
      Returns <- InternalTransformation(x=InternalGS, TransformationChoice=TransformationChoice, channelRange=channelRange,
      maxValue=maxValue, pos=pos, neg=neg, widthBasis=widthBasis,
      KeptMarkers=KeptMarkers)
      Flag <- TRUE
  }

  if (returnType == "GatingSet"){
      message("Option GatingSet selected for returnType")
      if (Flag == TRUE){return(Returns)
      } else {message("Combining Individual GatingSets")}
  } else if (returnType == "GatingPlot"){
      message("Option Unity plot selected for returnType")
      if (Flag == TRUE){
        message("Returning Individual Plot")
          
        if (is.null(optionalOutpath)){
            outpath <- getwd()
        } else {outpath <- optionalOutpath}
          
        Plot <- Utility_GatingPlots(x=Returns[1],
             sample.name=c("GROUPNAME", "TUBENAME"),
             removestrings=".fcs",
             subset="root", gtFile=NULL,
             DesiredGates=NULL, outpath = outpath,
             returnType="pdf", therows=1, thecolumns=1,
             optionalY=optionalY, optionalX=optionalX,
             filename=optionalName)
        return(Plot)
      } else {message("Returning Multiple Plots")
        
        if (is.null(optionalOutpath)){
            outpath <- getwd()
        } else {outpath <- optionalOutpath}
          
         Plots <- map(.x=Returns, .f=Utility_GatingPlots,
            sample.name=c("GROUPNAME", "TUBENAME"),
            removestrings=".fcs",
            subset="root", gtFile=NULL,
            DesiredGates=NULL, outpath = outpath,
            returnType="pdf", therows=1, thecolumns=1,
            optionalY=optionalY, optionalX=optionalX,
            filename=optionalName)
         return(Plots)
      }
  } else if (returnType == "NxN"){
    if (Flag == TRUE){
      message("Option NxN plot selected, returning single plot") 
        if (is.null(optionalOutpath)){
            outpath <- getwd()
        } else {outpath <- optionalOutpath}
        Plot <- Utility_NbyNPlots(x=Returns[[1]],
            sample.name=c("GROUPNAME", "TUBENAME"),
            removestrings = ".fcs", marginsubset = "root",
            gatesubset = "root", ycolumn = optionalY,
            bins = optionalBins, clearance = 0.2, gatelines = FALSE,
            reference = NULL, outpath = outpath,
            returntype="pdf", filename=optionalName)
        return(Plot)
    } else {
    message("Option NxN plot selected, returning single plot") 
    if (is.null(optionalOutpath)){
        outpath <- getwd()
    } else {outpath <- optionalOutpath}
    Plot <- map(.x=Returns, .f=Utility_NbyNPlots,
        sample.name=c("GROUPNAME", "TUBENAME"),
        removestrings = ".fcs", marginsubset = "root",
        gatesubset = "root", ycolumn = optionalY,
        bins = optionalBins, clearance = 0.2, gatelines = FALSE,
        reference = NULL, outpath = outpath,
        returntype="pdf", filename=optionalName)
    return(Plot)
    }

  } else if (returnType == "Unity"){
      if (Flag == FALSE){stop("Rerun with Multiple set to True")}
      message("Option Unity plot selected for returnType")
  } else {message("No option selected for returnType")}
}

#' Internal for Transformation Check
#' 
#' @param x Iterated in GatingSet object
#' @param TransformationChoice Default flowjo_biexp_trans
#' @param channelRange Argument
#' @param maxValue Argument
#' @param pos Argument
#' @param neg Argument
#' @param widthBasis Argument
#' @param KeptMarkers Columns to include
#' 
#' @importFrom flowWorkspace flowjo_biexp_trans
#' @importFrom flowWorkspace transformerList
#' @importFrom flowWorkspace transform
#' 
#' @noRd
InternalTransformation <- function(x, TransformationChoice, channelRange,
  maxValue, pos, neg, widthBasis, KeptMarkers){

  if (TransformationChoice == "flowjo_biexp"){
      MyTransform <- flowjo_biexp_trans(channelRange = channelRange,
      maxValue = maxValue, pos = pos, neg = neg, widthBasis = widthBasis)
  }

  TransformList <- transformerList(KeptMarkers, MyTransform)
  UnmixedGatingSet <- transform(x, TransformList)
  return(UnmixedGatingSet)

}