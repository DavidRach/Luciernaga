#' Internal for StainBrightnessIndexCalculator, creates positive and
#' negative gates using openCyto, before iterating through the .fcs
#' files for stain index.
#'
#' @param x The file path to the respective FMO folder
#' @param excludeThese Used to exclude columns from transformation,
#'   default is "FSC|SSC|Time"
#' @param channelRange Default for biexponential transformation is 4096
#' @param maxValue Default for biexponential transformation is 4194304
#' @param pos Default for biexponential transformation is 5.62
#' @param neg Default for biexponential transformation is 0
#' @param widthBasis Default for biexponential transformation is -1000
#' @param inverse.transform Default is TRUE
#' @param stringAppend Default is -A
#'
#' @importFrom flowWorkspace load_cytoset_from_fcs GatingSet
#'   flowjo_biexp_trans transformerList transform
#' @importFrom data.table fread
#' @importFrom stringr str_detect
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom flowCore filterList
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#'
#' @return A data.frame with staining index for the files within
#'   the FMO folder
#'
#' @noRd
BrightnessIndexIterator <- function(x, excludeThese,
      channelRange, maxValue, pos, neg, widthBasis,
      inverse.transform, stringAppend) {

  files <- list.files(x, pattern=".fcs", full.names=TRUE)
  theCytoset <- load_cytoset_from_fcs(files,
    truncate_max_range = FALSE, transformation = FALSE)
  theGatingSet <- GatingSet(theCytoset)
  SFC_Parameters <- colnames(theGatingSet)
  FluorophoresOnly <- SFC_Parameters[!str_detect(
       SFC_Parameters, excludeThese)]
  Biexponential <- flowjo_biexp_trans(channelRange=channelRange,
    maxValue=maxValue, pos=pos, neg=neg, widthBasis=widthBasis)
  MyBiexTransform <- transformerList(FluorophoresOnly, Biexponential)
  transform(theGatingSet, MyBiexTransform)

  FileLocation <- system.file("extdata", package = "Luciernaga")
  UnmixedGates <- fread(file.path(path = FileLocation,
                                pattern = 'GatesUnmixed.csv'))

  Example <- UnmixedGates[6]
  Example[1,1] <- FluorophoresOnly[1]
  Example[1,2] <- "+/-"
  Example[1,3] <- "root"
  Example[1,4] <- FluorophoresOnly[1]

  Template <- Example

  for (Fluorophore in FluorophoresOnly[2:length(FluorophoresOnly)]) {
    Template1 <- Template
    Template1[1,1] <- Fluorophore
    Template1[1,4] <- Fluorophore
    Example <- rbind(Example, Template1)
  }

  UnmixedGating <- suppressMessages(gatingTemplate(Example))
  suppressMessages(
       gt_gating(UnmixedGating, theGatingSet)) #flowCore filterList

  # x <- theGatingSet[26]
  Data <- map(.x=theGatingSet, .f=StainIndexLocal,
        inverse.transform=inverse.transform, stringAppend=stringAppend)
  Data <- Data |> bind_rows()
  return(Data)
}