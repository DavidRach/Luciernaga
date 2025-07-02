#' Small wrapper, passes folder (or files) to Luciernaga_GroupHeatmap
#' 
#' @param FolderPath File path to a Folder (or files) to be used for the heatmap
#' @param cutoff Minimal ratio needed for showing, cells not making criteria rolled into other
#' @param removeFollowing Default NULL, string character start for removal via stringr::str_detect
#' @param legend Default is "right", provide "none" to remove.
#' @param returnType Default is "data", alternate is "count"
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace GatingSet
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr rename
#' @importFrom lubridate dmy
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' 
#' @export
#' 
#' @return A GroupHeatmap plot
#' 
#' @examples
#' A <- 2+2
Luciernaga_FolderGroupHeatmap <- function(FolderPath, cutoff=0.01, removeFollowing=NULL, 
  legend="right", returnType="data"){

  if (length(FolderPath > 1)){
  TheFCSFiles <- FolderPath
  } else {
    TheFCSFiles <- list.files(path=FolderPath, pattern="fcs",
    full.names=TRUE) 
  } 

  Selected_CS <- load_cytoset_from_fcs(TheFCSFiles,
   truncate_max_range = FALSE, transformation = FALSE)
  Selected_GS <- GatingSet(Selected_CS)

  Metadata <- map(.x=Selected_GS, .f=GatingSetMetadata) |> dplyr::bind_rows()
  Metadata$TUBENAME <- NameCleanUp(removestrings=" (Cells)", Metadata$TUBENAME)
  Metadata <- Metadata |> rename(Fluorophore=TUBENAME)

  Metadata$GUID <-  gsub(".*\\(Cells\\)_", "", Metadata$GUID)
  Metadata$GUID <- gsub(".fcs", "", Metadata$GUID)
  Metadata <- Metadata |> rename(Cluster=GUID)

  Metadata$Date <- lubridate::dmy(Metadata$Date)
  Metadata$Date <- factor(Metadata$Date)
  Metadata <- Metadata |> arrange(Date) #|> as.character(Date)

  # removeFollowing <- "UV910"

  if (!is.null(removeFollowing)){
      Metadata <- Metadata |>
           dplyr::filter(!stringr::str_detect(Metadata$Cluster, removeFollowing))
  }

  if (returnType == "data"){
  Data <- Luciernaga_GroupHeatmap(reports=Metadata,
   nameColumn = "Date", cutoff=cutoff, legend=legend,
   returntype="plot")
  return(Data)  
  } else {return(Metadata)
  }
  }

#' Internal Luciernaga_FolderGroupHeatmap, parses GatingSet metadata
#' 
#' @param x An iterated GatingSet object
#' 
#' @importFrom flowCore keyword
#' 
#' @noRd
GatingSetMetadata <- function(x){
  GUID <- keyword(x, "GUID")
  TUBENAME <- keyword(x, "TUBENAME")
  Date <- keyword(x, "$DATE")
  Count <- unname(nrow(x))[[1]]
  Data <- data.frame(GUID, TUBENAME, Date, Count)
  return(Data) 
}