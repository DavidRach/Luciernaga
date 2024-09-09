#' Takes data.frame of QC_LibraryParse, and iterating over individual users
#' returns plots
#'
#' @param x Mapped individual name present within the Creator Column
#' @param Data The data.frame object working from
#' @param NameAppend Name to add to end of the file
#' @param outpath Desired storage location
#' @param references Whether to include red reference signatures
#' @param thecolumns Passed to Patchwork, desired number of columns
#' @param therows Passed to Patchwork, desired number of rows
#' @param width Passed to Patchwork, page width
#' @param height Passed to Patchwork, page height
#' @param saveCSV Whether to return a .csv of underlying data, default is TRUE
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom tidyr gather
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr rename
#' @importFrom purrr map
#' @importFrom utils read.csv
#'
#' @return A pdf of the plots and maybe a csv
#' @export
#'
#' @examples NULL
QC_UserLibraries <- function(x, Data, NameAppend, outpath, references=TRUE,
                             thecolumns=3, therows=4, width=7, height=9, saveCSV = TRUE){
  TheUserData <- Data %>% filter(Creator %in% x)
  TheUser <- TheUserData %>% select(Creator) %>% pull() %>% unique()
  TheUserData <- TheUserData %>% arrange(Fluorochrome)
  columnLength <- ncol(TheUserData)
  TheGatheredData <- gather(TheUserData, key="Detector", value="value", all_of(
    5:columnLength))

  #Re-leveling the factor
  Iterations <- columnLength-4
  MyVector <- 1:Iterations
  TheGatheredData$Detector <- factor(TheGatheredData$Detector, levels=MyVector)

  #Identifying the iterator samples
  TheGatheredData <- TheGatheredData %>% mutate(TheSamples=paste0(
    Fluorochrome, "_", Sample, " ", Date)) %>%
    relocate(TheSamples, .before=Fluorochrome)
  TheSamples <- TheGatheredData %>% select(TheSamples) %>% unique() %>% pull()

  if(references==TRUE){
    TotalDetectors <- Iterations

    if (TotalDetectors == 64){instrument <- "FiveLaser"
    TheFile <- file.path(getwd(), "CytekReferenceLibrary5L.csv")
    ReferenceData <- read.csv(TheFile, check.names = FALSE)
    } else if (TotalDetectors == 54){instrument <- "FourLaser"
    TheFile <- file.path(getwd(), "CytekReferenceLibrary4LUV.csv")
    ReferenceData <- read.csv(TheFile, check.names = FALSE)
    } else if (TotalDetectors == 38){instrument <- "ThreeLaser"
    TheFile <- file.path(getwd(), "CytekReferenceLibrary3L.csv")
    ReferenceData <- read.csv(TheFile, check.names = FALSE)
    } else {message("No References Found")}

    ReferenceData <- ReferenceData %>% rename(TheValue = "AdjustedY")
    ReferenceData$Fluorophore <- gsub(" ", "", gsub("-", "", gsub(
      ".", "", fixed=TRUE, ReferenceData$Fluorophore)))
    #ReferenceFluorList <- ReferenceData %>% select(Fluorophore) %>%
    #unique() %>% pull()
  }

  if (references == TRUE){
    ThePlots <- map(.x=TheSamples, .f=QC_RefPlots, Data=TheGatheredData,
                    references=TRUE, refData=ReferenceData)
  } else {ThePlots <- map(.x=TheSamples, .f=QC_RefPlots,
                          Data=TheGatheredData, references=FALSE, refData=NULL)
  }

  fileName <- paste(TheUser, NameAppend, sep="_")
  StorageLocation <- file.path(outpath, fileName)

  Utility_Patchwork(x=ThePlots, filename=fileName, outfolder=outpath,
                    thecolumns=thecolumns, therows=therows)

  if (saveCSV == TRUE){
    CSVName <- paste0(StorageLocation, ".csv")
    write.csv(TheUserData, file=CSVName, row.names = FALSE)
  }

}



#' Internal Test
#'
#' @param x x
#' @param Data x
#' @param references x
#' @param refData x
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#'
#' @keywords internal

QC_RefPlots <- function(x, Data, references=FALSE, refData=NULL){

  if (!references == TRUE){
    TheValue <- Data %>% filter(TheSamples %in% x)
    TheFluorochrome <- TheValue %>% select(Fluorochrome) %>% unique %>% pull()
    ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=Sample)) +
      geom_line() + theme_bw() + labs(title=paste0(TheFluorochrome, " ", x),
                                      x=NULL, y="Normalized Value") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      theme(plot.title = element_text(size = 8), axis.text.x = element_text(
        size = 6, angle = 45), panel.grid = element_blank(), axis.ticks.x = element_blank(),
      axis.title.y =  element_text(size=8)) + scale_x_discrete(breaks = unique(
        TheValue$Detector)[c(TRUE, rep(FALSE, 4))])
  } else {
    TheValue <- Data %>% filter(TheSamples %in% x)
    TheFluorochrome <- TheValue %>% select(Fluorochrome) %>% unique %>% pull()

    #ReferenceRetrieval
    TheFluorochrome1 <- TheFluorochrome #Name Backup
    TheFluorochrome <- gsub("AF", "Alexa Fluor", gsub("efl", "eFl", gsub(
      "Spk", "Spark", TheFluorochrome)))
    TheFluorochrome <- gsub(" ", "", gsub("-", "", gsub(".", "", fixed=TRUE,
                                                        TheFluorochrome)))
    ReferenceFluorList <- refData %>% select(Fluorophore) %>% unique() %>% pull()

    if (any(ReferenceFluorList == TheFluorochrome)) { #TheCleanedVersion
      ReferenceData1 <- refData %>% filter(Fluorophore %in% TheFluorochrome) %>%
        rename(value=TheValue)
      Iterations <- nrow(ReferenceData1)
      MyVector <- 1:Iterations
      ReferenceData1$Detector <- factor(ReferenceData1$Detector, levels=MyVector)
      ReferenceData1 <- ReferenceData1 %>% mutate(TheSamples="Nope")

      ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=TheSamples)) +
        geom_line() + theme_bw() + labs(title=paste0(x), x=NULL, y="Normalized Value") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 6, angle = 45),
        panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.title.y =  element_text(size=8)) +
        scale_x_discrete(breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])

      ThePlot <- ThePlot + geom_line(data = ReferenceData1, aes(x=Detector,
                                                                y=value, group=TheSamples), color="red")
    } else {TheValue <- Data %>% filter(TheSamples %in% x)
    TheFluorochrome <- TheValue %>% select(Fluorochrome) %>% unique %>% pull()
    ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=Sample)) +
      geom_line() + theme_bw() + labs(title=paste0(TheFluorochrome, " ", x),
                                      x=NULL, y="Normalized Value") + geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 6, angle = 45), panel.grid = element_blank(),
      axis.ticks.x = element_blank(), axis.title.y =  element_text(size=8)) +
      scale_x_discrete(breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])
    }
  }
  return(ThePlot)
}
