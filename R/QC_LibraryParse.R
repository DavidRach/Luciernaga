#' Parses Library Reference Control .XML files and returns
#'
#' @param x An .XML file
#' @param returntype What to return "data" or "plots"
#' @param references Plot argument, adds red reference signature
#' @param myfactor Plot argument, data column to group by for plotting. Default "Fluorophore".
#' @param namefactor Plot argument, data column name added to Plot Title.
#'
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_text
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr rename
#'
#' @return TBD
#' @export
#'
#' @examples
#' Folder_Location <- system.file("extdata", package = "Luciernaga")
#' XML_Pattern <- ".XML$"
#' XML_Files <- list.files(path = Folder_Location, pattern = XML_Pattern,
#'                         full.names = TRUE, recursive = FALSE)
#' SinglePlot <- QC_LibraryParse(XML_Files[2], returntype="plots", references=FALSE)
QC_LibraryParse <- function(x, returntype, references=TRUE, myfactor="Fluorophore", namefactor="Sample"){

  doc <- read_xml(x)
  TheChildren <- xml_children(doc)
  Creator <- xml_text(TheChildren[6])
  Date <- xml_text(TheChildren[7])
  TheMetadataNode <- TheChildren[14]
  MetadataNodes <- xml_children(TheMetadataNode)
  Fluorochrome <- xml_text(MetadataNodes[3])
  Fluorochrome1 <- Fluorochrome #BackupForSubs
  Sample <- xml_text(MetadataNodes[9])
  TheNode <- TheChildren[15] #If Area
  #TheNode <- TheChildren[16] #If Height
  TheNoddles <- xml_children(TheNode)
  TheValue <- xml_text(TheNoddles)

  TheValue <- as.numeric(TheValue)
  TheValue <- data.frame(TheValue)
  TheValue2 <- TheValue %>% mutate(Detector=1:nrow(.)) %>%
    relocate(Detector, .before=TheValue)

  Assembling <- TheValue2 %>% pivot_wider(
      names_from = Detector, values_from = TheValue)
  Assembling <- cbind(Date, Fluorochrome1, Assembling)
  Assembling <- Assembling %>% rename(Fluorochrome=Fluorochrome1)
  Assembling$Date <- as.Date(Assembling$Date)

  Data <- ColumnNaming(Assembling)
  Data <- cbind(Data, Sample, Creator) %>% relocate(Sample, Creator, .after=Fluorophore)

  if (returntype == "data"){
   return(Data)
  } else if (returntype == "plots"){
    plot <- LibraryPlot(x=Data, references=references, myfactor=myfactor, namefactor=namefactor)
    return(plot)
  } else {return(Data)}
}