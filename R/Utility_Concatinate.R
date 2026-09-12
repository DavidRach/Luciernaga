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
#' @param metadataCols column names from pData to append as metadata for the .fcs, default NULL
#' @param Override Testing purposes default FALSE
#' @param DataOverride Testing Purpose default NULL
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows select bind_cols
#' @importFrom tidyselect all_of starts_with
#' @importFrom utils write.csv
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom methods new
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
#'   newName = "MyConcatinatedFile", outpath = StorageLocation,
#'    export = FALSE, inverse.transform=TRUE)
#'
Utility_Concatinate <- function(gs, sample.name, removestrings, subsets, subsample=NULL,
  ReturnType, newName, outpath=NULL, export=FALSE, inverse.transform, metadataCols=NULL,
  Override = FALSE, DataOverride = NULL) {

if (Override == FALSE){
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
  internal = TRUE, export=FALSE, metadataCols=metadataCols, inverse.transform=inverse.transform) %>% bind_rows

if (ReturnType == "data.frame"){return(ConcatenatedFile)
}

} else {ConcatenatedFile <- DataOverride}

if (is.null(metadataCols)){metadataCols <- setdiff(colnames(ConcatenatedFile),colnames(gs[1]))}


NewData <- ConcatenatedFile %>% select(all_of(metadataCols))
Decisions <- DoWeConvert(NewData)

for(i in metadataCols){
ConcatenatedFile <- ExecuteNumerics(x=i, data=ConcatenatedFile, conversion=Decisions)
}

NewData <- ConcatenatedFile %>% select(all_of(metadataCols))
Decisions <- DoWeConvert(NewData)

if (!nrow(Decisions)==0){
CharResults <- map(.x=metadataCols, .f=ExecuteCharacters, data=ConcatenatedFile, conversion=Decisions)
CharResults <- Filter(Negate(is.null), CharResults)

if (length(CharResults)>1){Table <- bind_cols(CharResults)
  } else {Table <- CharResults[[1]]}

CSVName <- paste0("ReferenceDictionary", newName, ".csv")
if (is.null(outpath)) {outpath <- getwd()}
CSVDestination <- file.path(outpath, CSVName)

if (export == TRUE){write.csv(x=Table, file=CSVDestination, row.names = FALSE)}


Table <- Table %>% select(-starts_with("Count"))
Old <- Table %>% select(-starts_with("New")) %>% colnames()
#New <- Table %>% select(starts_with("New")) %>% colnames()

data <- ConcatenatedFile

for(i in Old){
  data <- ExecuteSwap(x=i, data=data, dictionary=Table)
}

} else {data <- ConcatenatedFile}

#str(data)
Referenced <- data

ff <- gs_pop_get_data(gs, subsets, inverse.transform = inverse.transform)

ConcatenatedExtra <- Referenced %>% select(all_of(metadataCols))
ConcatenatedMain <- Referenced %>% select(-all_of(metadataCols))

new_fcs <- Utility_ColAppend(ff=ff, DF=ConcatenatedMain, columnframe=ConcatenatedExtra)

new_kw <- new_fcs@description
UpdatedParameters <- parameters(new_fcs)
UpdatedExprs <- exprs(new_fcs)

#Temporary Fix to "Comp-" issue

UpdatedParameters@data$name <- gsub("Comp-", "", UpdatedParameters@data$name)
colnames(UpdatedExprs) <- gsub("Comp-", "", colnames(UpdatedExprs))

#for (i in seq_along(new_kw)) {
#  new_kw[[i]] <- gsub("Comp-", "", new_kw[[i]])
#}

new_FCS_2 <- new("flowFrame", exprs=UpdatedExprs, parameters=UpdatedParameters, description=new_kw)

AssembledName <- paste0(newName, ".fcs")

new_FCS_2@description$GUID <- AssembledName
new_FCS_2@description$`$FIL` <- AssembledName

if (is.null(outpath)) {outpath <- getwd()}

fileSpot <- file.path(outpath, AssembledName)

if (ReturnType == "fcs") {write.FCS(new_FCS_2, filename = fileSpot, delimiter="#")
} else {return(new_FCS_2)}
}