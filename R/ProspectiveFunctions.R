#' Internal Test
#'
#' @param x x
#' @param Data x
#' @param NameAppend x
#' @param outpath x
#' @param references x
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
#' @noRd

UserLibraryReturn <- function(x, Data, NameAppend, outpath, references=TRUE){
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
    ThePlots <- map(.x=TheSamples, .f=InternalLibraryPlot, Data=TheGatheredData,
                    references=TRUE, refData=ReferenceData)
  } else {ThePlots <- map(.x=TheSamples, .f=InternalLibraryPlot,
                          Data=TheGatheredData, references=FALSE, refData=NULL)
  }

  fileName <- paste(TheUser, NameAppend, sep="_")
  StorageLocation <- file.path(outpath, fileName)

  Utility_Patchwork(x=ThePlots, filename=fileName, outfolder=outpath,
                    thecolumns=4, therows=4)

  CSVName <- paste0(StorageLocation, ".csv")
  write.csv(TheUserData, file=CSVName, row.names = FALSE)

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
#' @noRd

  InternalLibraryPlot <- function(x, Data, references=FALSE, refData=NULL){

    if (!references == TRUE){
      TheValue <- Data %>% filter(TheSamples %in% x)
      TheFluorochrome <- TheValue %>% select(Fluorochrome) %>% unique %>% pull()
      ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=Sample)) +
        geom_line() + theme_bw() + labs(title=paste0(TheFluorochrome, " ", x),
        x=NULL, y="Normalized Value") + geom_hline(yintercept = 1,
        linetype = "dashed", color = "red") + theme(plot.title = element_text(
        size = 8), axis.text.x = element_text(size = 6, angle = 45),
        panel.grid = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y =  element_text(size=8)) + scale_x_discrete(
        breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])
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
          geom_line() + theme_bw() + labs(title=paste0(x), x=NULL,
          y="Normalized Value") + geom_hline(yintercept = 1, linetype = "dashed",
          color = "red") + theme(plot.title = element_text(size = 8),
          axis.text.x = element_text(size = 6, angle = 45),
          panel.grid = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y =  element_text(size=8)) + scale_x_discrete(
          breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])

        ThePlot <- ThePlot + geom_line(data = ReferenceData1, aes(x=Detector,
                   y=value, group=TheSamples), color="red")
      } else {TheValue <- Data %>% filter(TheSamples %in% x)
      TheFluorochrome <- TheValue %>% select(Fluorochrome) %>% unique %>% pull()
      ThePlot <- ggplot(TheValue, aes(x=Detector, y=value, group=Sample)) +
        geom_line() + theme_bw() + labs(title=paste0(TheFluorochrome, " ", x),
        x=NULL, y="Normalized Value") + geom_hline(yintercept = 1, linetype = "dashed",
        color = "red") + theme(plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 45), panel.grid = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y =  element_text(size=8)) +
        scale_x_discrete(breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])
      }
    }
    return(ThePlot)
  }





#' Internal Test
#'
#' @param x x
#' @param data x
#' @param NumberHits x
#'
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr ungroup
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select_if
#' @importFrom lsa cosine
#' @importFrom tibble rownames_to_column
#' @importFrom tidyselect starts_with
#' @importFrom dplyr arrange
#' @importFrom dplyr slice_head
#'
#' @noRd

UnknownFluorWhoseThis <- function(x, data, NumberHits) {
  StartingData <- data %>% filter(Sample %in% x)
  if(nrow(StartingData) > 1){StartingData <- StartingData %>% slice(1)}

  Fluorophore <- StartingData %>% mutate(Fluorophore = paste0(
    "ID_", Fluorochrome, "_", Sample)) %>% pull()
  DetectorCols <- StartingData %>% select(where(is.numeric))
  WhoseThis <- cbind(Fluorophore, DetectorCols)

  TotalDetectors <- length(DetectorCols)

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

  ReferenceData <- ReferenceData %>% select(-Instrument) %>%
    group_by(Fluorophore) %>% pivot_wider(
      names_from = Detector, values_from = AdjustedY) %>% ungroup()

  #if (!nrow(WhoseThis) == 1){}
  CombinedView <- bind_rows(WhoseThis, ReferenceData)
  Names <- CombinedView %>% select(Fluorophore) %>% pull()
  Numbers <- CombinedView %>% select_if(is.numeric)
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names

  CosineMatrix <- cosine(NumericsT)
  CosineMatrix <- round(CosineMatrix, 2)
  CosineFrame <- data.frame(CosineMatrix, check.names = FALSE)

  #CosineFrame <- CosineFrame[1,] #If Want To Work With Cols
  CosineFrame <- CosineFrame %>% select(starts_with("ID_"))
  TheData <- rownames_to_column(CosineFrame, var="Fluorophore")
  TheID <- TheData %>% select(starts_with("ID_")) %>% colnames()


  TheHits <- TheData %>% filter(!Fluorophore %in% TheID) %>%
    arrange(desc(.data[[TheID]])) %>% slice_head(n=NumberHits)
  return(TheHits)
}



#' Internal Test
#'
#' @param x x
#' @param data x
#'
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom tidyselect starts_with
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr filter
#' @importFrom tidyselect where
#' @importFrom dplyr group_by
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr ungroup
#' @importFrom dplyr rename
#' @importFrom dplyr bind_rows
#' @importFrom tidyr gather
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot
#'
#' @noRd

ReturnPlots <- function(x, data){
  TheseFluorophores <- x %>% select(Fluorophore) %>% pull()
  TheOGID <- x %>% select(starts_with("ID_")) %>% colnames()

  TheNewData <- data %>% mutate(TheID = paste0(
    "ID_", Fluorochrome, "_", Sample)) %>% relocate(TheID, .before=Fluorochrome)
  TheValue <- TheNewData %>% filter(TheID %in% TheOGID)

  TotalDetectors <- TheValue %>% select(where(is.numeric)) %>% ncol()

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

  ReferenceData <- ReferenceData %>% select(-Instrument) %>%
    group_by(Fluorophore) %>% pivot_wider(
      names_from = Detector, values_from = AdjustedY) %>% ungroup()

  TheValue1 <- TheValue %>% select(where(is.numeric)) %>%
    mutate(TheID = TheValue$TheID) %>% relocate(TheID, .before=ncol(1)) %>%
    rename(Fluorophore=TheID)
  SelectReferences <- ReferenceData %>% filter(Fluorophore %in% TheseFluorophores)

  TheAssembly <- bind_rows(TheValue1, SelectReferences)
  columnLength <- ncol(TheAssembly)

  TheGatheredData <- gather(
    TheAssembly, key="Detector", value="value", all_of(2:columnLength))
  MyVector <- 1:columnLength
  TheGatheredData$Detector <- factor(TheGatheredData$Detector, levels=MyVector)

  ThePlot <- ggplot(TheGatheredData, aes(x=Detector, y=value, group=Fluorophore,
             color=Fluorophore)) + geom_line() + theme_bw() + labs(title=TheOGID,
             x=NULL, y="Normalized Value") + geom_hline(yintercept = 1,
             linetype = "dashed", color = "red") + theme(plot.title = element_text(
             size = 8), axis.text.x = element_text(size = 6, angle = 45),
             panel.grid = element_blank(), axis.ticks.x = element_blank(),
             axis.title.y =  element_text(size=8)) + scale_x_discrete(
             breaks = unique(TheValue$Detector)[c(TRUE, rep(FALSE, 4))])

  return(ThePlot)
}





