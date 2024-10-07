#' Parses Library Reference Control .XML files and returns
#'
#' @param x An .XML file
#' @param returntype What to return "dataframe" or "plots"
#' @param references Whether to add reference fluorophore signature in red
#'
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_text
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr rename
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom utils read.csv
#'
#' @return TBD
#' @export
#'
#' @examples
#'
#' Folder_Location <- system.file("extdata", package = "Luciernaga")
#' XML_Pattern <- ".XML$"
#' XML_Files <- list.files(path = Folder_Location, pattern = XML_Pattern,
#'                         full.names = TRUE, recursive = FALSE)
#' SinglePlot <- QC_LibraryParse(XML_Files[2], returntype="plots", references=FALSE)
QC_LibraryParse <- function(x, returntype, references=TRUE){

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

  if (returntype == "dataframe"){

    Assembling <- TheValue2 %>% pivot_wider(
      names_from = Detector, values_from = TheValue)
    Assembling <- cbind(Fluorochrome1, Sample, Creator, Date, Assembling)
    Assembling <- Assembling %>% rename(Fluorochrome=Fluorochrome1)
    Assembling$Date <- as.Date(Assembling$Date)
    return(Assembling)

  } else if (returntype == "plots"){

    if(references==TRUE){
      Folder_Location <- system.file("extdata", package = "Luciernaga")

      TotalDetectors <- nrow(TheValue2)

      if (TotalDetectors == 64){instrument <- "FiveLaser"
      TheFile <- file.path(Folder_Location, "CytekReferenceLibrary5L.csv")
      ReferenceData <- read.csv(TheFile, check.names = FALSE)
      } else if (TotalDetectors == 54){instrument <- "FourLaser"
      TheFile <- file.path(Folder_Location, "CytekReferenceLibrary4LUV.csv")
      ReferenceData <- read.csv(TheFile, check.names = FALSE)
      } else if (TotalDetectors == 38){instrument <- "ThreeLaser"
      TheFile <- file.path(Folder_Location, "CytekReferenceLibrary3L.csv")
      ReferenceData <- read.csv(TheFile, check.names = FALSE)
      } else {message("No References Found")}

      ReferenceData <- ReferenceData %>% rename(TheValue = "AdjustedY")
      ReferenceFluorList <- ReferenceData %>% select(Fluorophore) %>%
        unique() %>% pull()

      if (any(ReferenceFluorList == Fluorochrome)){ #TheRegularVersion
        ReferenceData1 <- ReferenceData %>% filter(Fluorophore %in% Fluorochrome)

        ThePlot <- ggplot(TheValue2, aes(x=Detector, y=TheValue)) +
          geom_line() + theme_bw() + labs(title=paste0(Fluorochrome, " ", Sample),
          y="Normalized") + geom_hline(yintercept = 1, linetype = "dashed",
          color = "red") + theme(plot.title = element_text(size = 8),
          axis.title.y =  element_text(size=8))

        ThePlot <- ThePlot + geom_line(data = ReferenceData1, aes(
          x=Detector, y=TheValue), color="red")

      } else {
        Fluorochrome <- gsub("AF", "Alexa Fluor", gsub(
          "efl", "eFl", gsub("Spk", "Spark", Fluorochrome)))
        Fluorochrome <- gsub(" ", "", gsub("-", "", gsub(
          ".", "", fixed=TRUE, Fluorochrome)))

        ReferenceData$Fluorophore <- gsub(" ", "", gsub(
          "-", "", gsub(".", "", fixed=TRUE, ReferenceData$Fluorophore)))
        ReferenceFluorList <- ReferenceData %>% select(Fluorophore) %>%
          unique() %>% pull()

        if (any(ReferenceFluorList == Fluorochrome)){ #TheCleanedVersion
          ReferenceData1 <- ReferenceData %>% filter(Fluorophore %in% Fluorochrome)

          ThePlot <- ggplot(TheValue2, aes(x=Detector, y=TheValue)) +
            geom_line() + theme_bw() + labs(title=paste0(
            Fluorochrome1, " ", Sample), y="Normalized") + geom_hline(yintercept = 1,
            linetype = "dashed", color = "red") + theme(plot.title = element_text(
            size = 8), axis.title.y =  element_text(size=8))

          ThePlot <- ThePlot + geom_line(data = ReferenceData1, aes(
            x=Detector, y=TheValue), color="red")

        } else {#The fallback version
          ThePlot <- ggplot(TheValue2, aes(x=Detector, y=TheValue)) + geom_line() +
          theme_bw() + labs(title=paste0(Fluorochrome1, " ", Sample),
          y="Normalized") + geom_hline(yintercept = 1, linetype = "dashed",
          color = "red") + theme(plot.title = element_text(size = 8),
          axis.title.y =  element_text(size=8))
        }
      }

    } else {#If no references specified
      ThePlot <- ggplot(TheValue2, aes(x=Detector, y=TheValue)) + geom_line() +
        theme_bw() + labs(title=paste0(Fluorochrome, " ", Sample), y="Normalized") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 8),
        axis.title.y =  element_text(size=8))}

    return(ThePlot)
  }
}
