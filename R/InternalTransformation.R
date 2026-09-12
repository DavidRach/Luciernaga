
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
#'  transformerList transform
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