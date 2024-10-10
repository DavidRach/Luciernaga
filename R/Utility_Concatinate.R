#' Concatenate a gs based on subset and optional subsample
#'
#' @param gs A gating set object
#' @param sample.name Keyword specifying sample name
#' @param removestrings Value to be removed from sample name
#' @param subsets The gating hierarchy subset you want to include
#' @param subsample Default is set to NULL. If down-sample is desired, Total number of
#' events to subsample from each specimen
#' @param ReturnType Whether to return as a "data.frame", a "flow.frame", or export
#' to an "fcs" file to the outpath location
#' @param newName File Name for the Concatenate File
#' @param outpath Location to store the concatenated file
#' @param export Whether to export as a .fcs file.
#' @param inverse.transform Whether to reverse the GatingSet Transform on the data,
#' default is set to FALSE.
#'
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore write.FCS
#'
#' @return A concatenated data.frame, flow.frame or fcs file, with reference .csv for
#' the specimen information
#' @export
#'
#' @examples
#'
#' library(BiocGenerics)
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' Unmixed_FullStained <- FCS_Files[grep("Unmixed", FCS_Files)]
#' UnmixedFCSFiles <- Unmixed_FullStained[1:2]
#' UnmixedCytoSet <- load_cytoset_from_fcs(UnmixedFCSFiles[1:2],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' UnmixedGatingSet <- GatingSet(UnmixedCytoSet)
#' Markers <- colnames(UnmixedCytoSet)
#' KeptMarkers <- Markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", Markers)]
#' MyBiexponentialTransform <- flowjo_biexp_trans(channelRange = 256,
#'   maxValue = 1000000,pos = 4.5, neg = 0, widthBasis = -1000)
#' TransformList <- transformerList(KeptMarkers, MyBiexponentialTransform)
#' UnmixedGatingSet <- flowWorkspace::transform(UnmixedGatingSet, TransformList)
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' UnmixedGates <- fread(file.path(path = FileLocation, pattern = 'GatesUnmixed.csv'))
#' UnmixedGating <- gatingTemplate(UnmixedGates)
#' gt_gating(UnmixedGating, UnmixedGatingSet)
#'
#' removestrings <-  c("DTR_", ".fcs")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' ConcatinatedReturn <- Utility_Concatinate(gs=UnmixedGatingSet,
#'   sample.name = "GROUPNAME", removestrings=removestrings,
#'   subsets="live", subsample = 2000, ReturnType = "flow.frame",
#'   newName = "MyConcatinatedFile", outpath = StorageLocation, export = FALSE)
#'
Utility_Concatinate <- function(gs, sample.name, removestrings, subsets, subsample=NULL,
  ReturnType, newName, outpath=NULL, export=FALSE, inverse.transform = FALSE) {

gs <- gs
sample.name <- sample.name
removestrings <- removestrings
subsample <- subsample

if (!is.null(subsample)){
  if (!is.numeric(subsample)) {
    stop("Subsample argument is not numeric, enter just a number, no quotes needed")}
}

ConcatenatedFile <- map(gs, .f=Utility_Downsample, sample.name=sample.name,
  removestrings=removestrings, subsets=subsets, subsample=subsample,
  internal = TRUE, export=FALSE) %>% bind_rows

if (ReturnType == "data.frame"){return(ConcatenatedFile)
}

SpecimenNames <- data.frame(table(ConcatenatedFile$specimen))
colnames(SpecimenNames)[[1]] <- "specimen"
colnames(SpecimenNames)[[2]] <- "Count"
SpecimenNames <- SpecimenNames %>% mutate(Reference = as.numeric(factor(specimen)))
CSVName <- paste0("ReferenceDictionary", newName, ".csv")
if (is.null(outpath)) {outpath <- getwd()}
CSVDestination <- file.path(outpath, CSVName)

if (export == TRUE){write.csv(x=SpecimenNames, file=CSVDestination, row.names = FALSE)}

SpecimenNames <- SpecimenNames %>% select(-Count)
Referenced <- left_join(ConcatenatedFile, SpecimenNames, by="specimen")
# If Alternative Found Enter Stage Right
Referenced <- Referenced %>% select(-specimen)
Referenced <- Referenced %>% rename(specimen = Reference)

ff <- gs_pop_get_data(gs, subsets, inverse.transform = inverse.transform)
ConcatenatedExtra <- Referenced %>% select(specimen)
ConcatenatedMain <- Referenced %>% select(-specimen)

new_fcs <- Utility_ColAppend(ff=ff, DF=ConcatenatedMain, columnframe=ConcatenatedExtra)

# Renaming of the .fcs file in the internal parameters would occur here.
new_FCS_2 <- new_fcs
#Until Then

if (ReturnType == "flow.frame") {return(new_FCS_2)}

if (ReturnTye  == "fcs"){
  if (is.null(outpath)) {outpath <- getwd()}
  TheFileName <- paste0(newName, ".fcs")
  fileSpot <- file.path(outpath, TheFileName)
  write.FCS(new_FCS_2, filename = fileSpot, delimiter="#")
}

}

