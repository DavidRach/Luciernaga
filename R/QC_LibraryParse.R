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



#' Internal, plots QC_LibraryParse data into plots
#'
#' @param x The passed data
#' @param references Plot argument, adds red reference signature
#' @param myfactor Plot argument, data column to group by for plotting. Default "Fluorophore".
#' @param namefactor Plot argument, data column name added to Plot Title.
#'
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect where
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom dplyr filter
#' @importFrom rlang sym !!
#'
#' @return A ggplot2 object
#'
#' @noRd
LibraryPlot <- function(x, references=TRUE, myfactor="Fluorophore", namefactor="Sample"){
  Data <- x

  if (nrow(Data) == 1){Sample <- Data %>% pull(.data[[namefactor]])
  } else {Sample <- ""}

  NumberDetectors <- sum(sapply(Data, is.numeric))
  TheFluorophore <- Data %>% pull(Fluorophore)
  TheFluorophore <- unique(TheFluorophore)
  TheDetectors <- colnames(Data)[sapply(Data, is.numeric)]

  Data <- Data %>%
    pivot_longer(cols = where(is.numeric),
                 names_to = "Detector", values_to = "TheValue")

  Data$Detector <- factor(Data$Detector, levels=TheDetectors)


  ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors)
  ReferenceData <- ReferenceData %>% rename(TheValue = "AdjustedY")
  ReferenceFluorList <- ReferenceData %>% select(Fluorophore) %>%
    unique() %>% pull()

  ThePlot <- ggplot(Data, aes(x=Detector, y=TheValue, group=.data[[myfactor]])) + geom_line()  +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title=paste0(TheFluorophore, " ", Sample), y="Normalized") +
    theme_bw() + theme(plot.title = element_text(size = 8),
                       axis.title.y =  element_text(size=8),
                       axis.text.x = element_text(size=6, angle = 70, hjust = 1))

  if (any(ReferenceFluorList == TheFluorophore)){

    if (references == TRUE){
      ReferenceData1 <- ReferenceData %>% filter(Fluorophore %in% TheFluorophore) %>%
        mutate(StandIn="Reference")
      ReferenceData1 <- ReferenceData1 %>% rename(!!sym(myfactor) := StandIn)
      ReferenceData1$Detector <- TheDetectors
      ThePlot <- ThePlot + geom_line(data = ReferenceData1, aes(
        x=Detector, y=TheValue, group=.data[[myfactor]]), color="red")
      return(ThePlot)
    }
  } else {
      return(ThePlot)
    }
}

#' Internal, thin wrapper that filters for Fluorophore and then group plots from Library Data
#'
#' @param x A iterated Fluorophore Name
#' @param data The data output from QC_Library
#' @param myfactor The desired factor for group
#' @param animate Whether to convert to ggplotly output, default FALSE
#'
#' @return A ggplot2 or a ggplotly object
#' @noRd
LibraryPlotWrapper <- function(x, data, myfactor, animate=FALSE){
  Subset <- data %>% filter(Fluorophore %in% x)
  plot <- Luciernaga:::LibraryPlot(x=Subset, myfactor=myfactor)
  if (animate == TRUE){
    plot <- plotly::ggplotly(plot)
  }
  return(plot)
}
