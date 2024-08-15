#' Main Luciernaga Function, normalized on individual cell level
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param subsets The Gating Hierarchy level you will be sampling at
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param removestrings A string of character values to remove from sample.name
#' @param Verbose Whether to print outputs as you go.
#' @param unmixingcontroltype Whether your inputs are "cells", "beads" or "both"
#' @param Unstained Set to True when running unstained samples that don't have Unstained in the Name
#' @param ratiopopcutoff A numeric ratio for peak detector inclusion, default is set to 0.01 all startingcells
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param experiment.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param condition Provide directly experiment name (ex. "JAN2024")
#' @param condition.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param AFOverlap A data.frame or a filepath to the CSV containing the Autofluorescence
#' overlap of individual fluorophores for exclusion
#' @param stats Whether to take "mean" or "median"
#' @param desiredAF Main Autofluorescence Detector (ex. "V7-A")
#' @param externalAF A data.frame row containing external autofluorescence to subtract from single colors.
#' @param ExportType Whether to return "fcs", "data"
#' @param SignatureReturnNow Short circuits the function and returns signature for specified autofluorescence.
#' @param outpath  Location where created .fcs and .csv files are sent
#' @param minimalfcscutoff A ratio determining number cells needed for .fcs export
#' @param Subtraction Whether for single color controls to use "Internal" or "External" autofluorescence
#' @param SCData Whether to return "subtracted" or "raw" data for single colors
#' @param NegativeType Whether to append a negative pop to .fcs file, "default", "artifical" or "samples"
#' @param Brightness Whether sum of detectors should be returned.
#' @param LocalMaximaRatio Height of peaks to proceed
#' @param SecondaryPeaks Number of Secondary Peaks, default is set to 2.
#' @param Increments Rounding parameter, default is set to 0.1
#' @param RetainedType Whether to return "raw" or "normalized" values for lineplots.
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr summarize_all
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom flowWorkspace keyword
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom flowCore exprs
#' @importFrom purrr map
#' @importFrom purrr set_names
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL
LuciernagaQC <- function(x, subsets, sample.name, removestrings, Verbose = FALSE,
                                  unmixingcontroltype="both", Unstained=FALSE, ratiopopcutoff=0.01,
                                  experiment = NULL, experiment.name = NULL, condition = NULL, condition.name = NULL,
                                  AFOverlap, stats="median", desiredAF=NULL, externalAF = NULL,
                                  ExportType, SignatureReturnNow=FALSE, outpath,
                                  minimalfcscutoff = 0.05, Subtraction = "Internal", SCData = "subtracted",
                                  NegativeType= "default", TotalNegatives = 500, Brightness=FALSE,
                                  LocalMaximaRatio=0.15, SecondaryPeaks=2, Increments, RetainedType="normalized", ...){

  name <- keyword(x, sample.name)
  Type <- Luciernaga:::Typing(name=name, unmixingcontroltype=unmixingcontroltype, Unstained=Unstained)
  #Type <- Typing(name=name, unmixingcontroltype=unmixingcontroltype, Unstained=Unstained)
  AggregateName <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings)
  #AggregateName <- NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, ...)

  Experiment <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, experiment = experiment,
                                           experiment.name = experiment.name, returnType = "experiment")
  #Experiment <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, returnType = "experiment", ...)

  Condition <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, condition=condition,
                                          condition.name = condition.name, returnType = "condition")
  #Condition <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, returnType = "condition", ...)

  # Internal Name Cleanup, to make sure fluorophores match
  InternalCleanupList <- c(".fcs", "Cells", "Beads", " ", "_", "-", ".", "(", ")")
  name <- Luciernaga:::NameCleanUp(name, InternalCleanupList)
  #name <- NameCleanUp(name, InternalCleanupList)

  # Brought back until overlap is converted from name reliant to type reliant.
  if (Unstained == TRUE) {if(!str_detect(name, "stained")){name <- paste0(name, "_Unstained")}}

  ###############
  # Exprs Setup #
  ###############

  #Retrieving the exprs data for the target population
  ff <- gs_pop_get_data(x, subsets)
  startingcells <- nrow(ff)[[1]]
  DF <- as.data.frame(exprs(ff[[1]]), check.names=FALSE)

  #For Future Row Reordering
  Backups <- DF %>% mutate(Backups = 1:nrow(DF)) %>% select(Backups)

  #For Future Column Reordering
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  OriginalColumnsVector <- colnames(DF)

  # Storing Other Parameter (Time, FSC, SSC) for Later Return
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  # Consolidating Desired Parameters
  n <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]

  # Next Step Will Remove Negative Values, to what extent is an ordinary .fcs file affected?
  if (Verbose == TRUE){
    TheTotal <- nrow(n) * ncol(n)
    BelowZero <- sum(apply(n, 2, function(x) x < 0))
    message(round(BelowZero/TheTotal,2), " of all events were negative and will be rounded to 0")
  }

  ##########################################
  # Generating Normalized by Peak Detector #
  ##########################################

  # Normalizing Individual Cells By Peak Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  #Normalized <- round(Normalized, 4) # Previously at 1
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  # Figuring out where Raw and Normalized Columns Start and End
  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsN + ColsNormalized

  # Bringing Together Raw And Normalized Data.Frames
  WorkAround <- cbind(n, Normalized)

  ####################################
  # Enumerating Peak Detector Counts #
  ####################################

  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  PeakDetectorCounts <- PeakDetectorCounts %>% arrange(desc(Counts))

  if (Type == "Cells_Unstained") {CellCutoff <- startingcells*ratiopopcutoff
                                 Detectors <- PeakDetectorCounts %>% filter(Counts > CellCutoff)
  }

  if (Type == "Cells") {CellCutoff <- startingcells*ratiopopcutoff
                        Detectors <- PeakDetectorCounts %>% filter(Counts > CellCutoff)
  }

  if (Type == "Beads_Unstained") {BeadCutoff <- startingcells/ColsN #EqualDistributionAssumption
                                  Detectors <- PeakDetectorCounts %>% filter(Counts > BeadCutoff)
  }


  if (Type == "Beads") {BeadCutoff <- startingcells*ratiopopcutoff
                        Detectors <- PeakDetectorCounts %>% filter(Counts > BeadCutoff)
  }

  if (Verbose == TRUE){
  print(Detectors)
  }

  #####################################################
  # Handling Fluorophore and Autofluorescent Overlaps #
  #####################################################
  if (is.data.frame(AFOverlap)){AFData <- AFOverlap
  } else {AFData <- read.csv(file=AFOverlap, check.names = FALSE)
  }

  AFChannels <- AFData %>% filter(Fluorophore %in% "Unstained") %>%
    pull(MainDetector) %>% str_split(",", simplify = TRUE)
  AFChannels <- AFChannels[1,]
  AFChannels <- gsub("-A", "", AFChannels)

  TheSCData <- AFData %>% filter(Fluorophore != "Unstained")
  TheSCData$Fluorophore <- gsub("-A", "", TheSCData$Fluorophore)
  TroubleChannels <- TheSCData %>% pull(Fluorophore)

  results <- map(.x=TroubleChannels, .f=TroubleChannelExclusion, TheSCData=TheSCData, MainDetector=MainDetector, AFChannels=AFChannels) %>%
    set_names(TroubleChannels)

  matching_names <- names(results)[str_detect(name, names(results))]
  if (length(matching_names) > 0) {ExclusionList <- results[[matching_names[1]]]
  Retained <- Detectors %>% filter(!Fluors %in% ExclusionList) %>% pull(Fluors)
  } else if (str_detect(name, "nstained")){Retained <- Detectors %>%
    pull(Fluors)
  } else {Retained <- Detectors %>% filter(!Fluors %in% AFChannels) %>%
    pull(Fluors)}

  if (length(Retained) == 0) {stop("There were no Retained detectors in ", name)}

  if (Verbose == TRUE){print(Retained)}


  ##################################
  # Main Autofluorescence Detector #
  ##################################

  # Inverse of retained, then top for peak

  if (!str_detect(name, "nstained")){
    Intermediate <- Detectors %>% filter(!Fluors %in% Retained)
    if (nrow(Intermediate) >0){
      TheMainAF <- Intermediate %>% slice(1) %>% pull(Fluors)
    }

  # Adding way to switch to alternate AFs.
  if (is.null(desiredAF)){This <- WorkAround %>% filter(.data[[TheMainAF]] == 1) %>% select(all_of(1:ColsN))
  } else {desiredAF <- Luciernaga:::NameCleanUp(desiredAF, c("-A"))
          TheMainAF <- desiredAF
          This <- WorkAround %>% filter(.data[[desiredAF]] == 1) %>% select(all_of(1:ColsN))}

  if (is.null(externalAF)){Samples <- AveragedSignature(x=This, stats=stats)
  } else {if(is.data.frame(externalAF)){
    Samples <- externalAF
    MaxVal <- do.call(pmax, Samples)
    TheNormed <- Samples/MaxVal
    TheCounts <- colSums(TheNormed == 1)
    ThePeakDetectorCounts <- data.frame(Fluors = names(TheCounts), Counts = TheCounts)
    rownames(ThePeakDetectorCounts) <- NULL
    TheMainAF <- ThePeakDetectorCounts %>% arrange(desc(Counts)) %>% slice(1) %>% pull(Fluors)
    TheMainAF <- gsub("-A", "", TheMainAF)
    } else {stop("externalAF needs to be a single row of a data.frame")}}

  }

  if (str_detect(name, "nstained")){
    Intermediate <- Detectors
    if (nrow(Intermediate) >0){
      TheMainAF <- Intermediate %>% slice(1) %>% pull(Fluors)
    }

    # Adding way to switch to alternate AFs.
    if (is.null(desiredAF)){This <- WorkAround %>% filter(.data[[TheMainAF]] == 1) %>% select(all_of(1:ColsN))
    } else {desiredAF <- Luciernaga:::NameCleanUp(desiredAF, c("-A"))
    TheMainAF <- desiredAF
    This <- WorkAround %>% filter(.data[[desiredAF]] == 1) %>% select(all_of(1:ColsN))}

    if (is.null(externalAF)){Samples <- AveragedSignature(x=This, stats=stats)
    } else {if(is.data.frame(externalAF)){
      Samples <- externalAF
      MaxVal <- do.call(pmax, Samples)
      TheNormed <- Samples/MaxVal
      TheCounts <- colSums(TheNormed == 1)
      ThePeakDetectorCounts <- data.frame(Fluors = names(TheCounts), Counts = TheCounts)
      rownames(ThePeakDetectorCounts) <- NULL
      TheMainAF <- ThePeakDetectorCounts %>% arrange(desc(Counts)) %>% slice(1) %>% pull(Fluors)
      TheMainAF <- gsub("-A", "", TheMainAF)
      } else {stop("externalAF needs to be a single row of a data.frame")}}

  }

  if (SignatureReturnNow == TRUE){return(Samples)}

  ######################################
  # Handling the Double Peak Detectors #
  ######################################

  WorkAround1 <- WorkAround %>% mutate(Backups = Backups$Backups) %>%
    relocate(Backups, .before = 1) #This will change the start/end count

  #StartN <- StartNormalizedMergedCol+1
  #EndN <- EndNormalizedMergedCol+1
  #DoublePeaks <- WorkAround1 %>% filter(rowSums(WorkAround1[, StartN:EndN] == 1.00) >= 2)
  #TheDoublePeaks <- DoublePeaks %>% select(Backups) %>% pull()
  #WorkAround2 <- WorkAround1 %>% filter(!Backups %in% TheDoublePeaks)

  ###############################
  # Back to Regular Programming #
  ###############################

  if (str_detect(name, "nstained")){
    # x <- Retained[3]
    RetainedDF <- map(.x= Retained, .f=Luciernaga:::UnstainedSignatures,
                      WorkAround1=WorkAround1, alternatename=AggregateName,
                      ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
                      EndNormalizedMergedCol=EndNormalizedMergedCol, Verbose=Verbose,
                      Increments=Increments) %>%
      bind_rows()

  } else {
    # x <- Retained[1]
    RetainedDF <- map(.x= Retained, .f=Luciernaga:::SingleStainSignatures,
                      WorkAround1=WorkAround1, AggregateName=AggregateName,
                      ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
                      EndNormalizedMergedCol=EndNormalizedMergedCol, Samples=Samples,
                      Increments=Increments, Subtraction=Subtraction, stats=stats,
                      TheMainAF=TheMainAF, Verbose = Verbose, SCData = SCData,
                      LocalMaximaRatio=LocalMaximaRatio, SecondaryPeaks=SecondaryPeaks) %>% bind_rows()
    }

  Reintegrated <- left_join(RetainedDF, StashedDF, by = "Backups")

  BackupsCol <- "Backups"
  NormalizedColumns <- colnames(Normalized)
  ClusterCol <- "Cluster"
  RearrangedColumns <- c(BackupsCol, OriginalColumnsVector, NormalizedColumns,
                         ClusterCol)

  OriginalStart <- length(BackupsCol) + 1
  OriginalEnd <- length(BackupsCol) + length(OriginalColumnsVector)

  Reintegrated1 <- Reintegrated %>% relocate(all_of(RearrangedColumns))
  Reintegrated1 <- Reintegrated1 %>% arrange(Backups)

  if (ExportType == "fcs"){

    BrightnessReturn <- Luciernaga:::Genesis(x=Reintegrated1, ff=ff, minimalfcscutoff = 0.05,
                        AggregateName=AggregateName, Brightness = Brightness,
                        outpath=outpath, OriginalStart = OriginalStart,
                        OriginalEnd = OriginalEnd, stats=stats, NegativeType=NegativeType,
                        TotalNegatives=TotalNegatives, Samples=Samples, ExportType=ExportType)
  }

  if (ExportType == "data.frame"){
    FinalDataFrame <- Reintegrated %>% select(-Backups)
    return(FinalDataFrame)
  }

  if (ExportType == "data"){
    ExportData <- RetainedDF %>% select(-Backups)
    TheData <- data.frame(table(ExportData$Cluster))
    TheData <- TheData%>% dplyr::arrange(desc(Freq))
    colnames(TheData)[1] <- "Cluster"
    colnames(TheData)[2] <- "Count"
    #Data

    TheExperiment <- as.character(Experiment)
    TheCondition <- as.character(Condition)

    TheData <- TheData %>% mutate(Sample=AggregateName)
    TheData <- TheData %>% mutate(Experiment=TheExperiment)
    TheData <- TheData %>% mutate(Condition=TheCondition)
    TheData <- TheData %>% relocate(Sample, Experiment, Condition, .before=Cluster)

    TheClusters <- TheData %>% pull(Cluster)

    TheSummary <- map(.x=TheClusters, .f=LuciernagaSmallReport, Data=ExportData, RetainedType=RetainedType,
        ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
        EndNormalizedMergedCol=EndNormalizedMergedCol, stats=stats) %>% bind_rows()

    FinalData <- left_join(TheData, TheSummary, by = "Cluster")
    return(FinalData)
  }

}



#' Internal for LuciernagaQC
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr rename
#' @noRd
LuciernagaSmallReport <- function(x, Data, RetainedType, ColsN, StartNormalizedMergedCol,
                             EndNormalizedMergedCol, stats){
    if (RetainedType == "raw"){Data <- Data %>% filter(Cluster %in% x) %>% select(all_of(1:ColsN))}
    if (RetainedType == "normalized"){Data <- Data %>% filter(Cluster %in% x) %>%select(all_of(
      StartNormalizedMergedCol:EndNormalizedMergedCol))}

    Averaged <- AveragedSignature(Data, stats)
    Summary <- cbind(x, Averaged) %>% rename(Cluster = x)
    return(Summary)
}

#' Internal for LuciernagaQC
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @noRd
TroubleChannelExclusion <- function(x, TheSCData, MainDetector, AFChannels){
  Internal <- TheSCData %>% filter(Fluorophore %in% x) %>% pull(MainDetector)
  Internal <- gsub("-A", "", Internal)
  Exclusion <- setdiff(AFChannels, Internal)
  return(Exclusion)
}

#' Internal for LuciernagaQC
#'
#' @importFrom dplyr summarize_all
#' @noRd
AveragedSignature <- function(x, stats){
  Signature <- x %>% summarize_all(stats)
  return(Signature)
}

#' Internal for LuciernagaQC
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr arrange
#'
#' @noRd
DetectorPeakCounts <- function(x, StartN, EndN){
x <- x %>% select(-Backups)
Normalized <- x %>% select(all_of(
  StartN:EndN))
Counts <- colSums(Normalized == 1)
PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
rownames(PeakDetectorCounts) <- NULL
PeakDetectorCounts <- PeakDetectorCounts %>% arrange(desc(Counts))
return(PeakDetectorCounts)
}

#' Internal LuciernagaQC detects Fluorophore Peak Detectors by Local Maxima
#'
#' @param theX A vector of detectors from 1:n
#' @param theY The corresponding y values corresponding to the measurements
#' of theX
#' @param therepeats Additional values to temporarily add to the edges to
#' allow for peak detection
#' @param w The span around which rolling will happen
#' @param alternatename The cleaned up name passed to the plots (internal)
#' @param Verbose Whether to print line plot outputs
#' @param ... Additional arguments passed to zoo package
#'
#' @importFrom stats loess
#' @importFrom zoo rollapply
#' @importFrom zoo zoo
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr select
#'
#' @return A value to be determined later
#' @noRd
LocalMaxima <- function(theX, theY, therepeats, w, alternatename, Verbose = FALSE, ...){

  #Adding Margins
  repeats <- therepeats*2
  NewX <- length(theX) + repeats
  NewX <- 1:NewX
  NewYmin <- theY[[1]]*0.99
  LengthY <- length(theY)
  NewYmax <- theY[[LengthY]]*0.80
  replicatedYmin <- rep(NewYmin, each = therepeats)
  replicatedYmax <- rep(NewYmax, each = therepeats)
  NewY <- c(replicatedYmin, theY, replicatedYmax)

  x <- NewX
  y <- NewY

  ##Local Maxima Functions
  n <- length(y)
  #y.smooth <- loess(y ~ x)$fitted
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  peaks <- list(x=x[i.max]-therepeats, i=i.max-therepeats,
                y.hat=y.smooth[(therepeats + 1):(length(y.smooth) - therepeats)])

  peak_points <- peaks$x

  MainData <- data.frame(x = theX, y = theY, yhat = peaks$y.hat)

  PointData <- MainData %>% filter(x %in% peak_points)

  Views <- ggplot(MainData, aes(x = x, y = y)) +
    geom_point(size = 2, color = "Gray") +
    geom_line(aes(y = yhat), linewidth = 1) +
    geom_point(data = PointData, aes(x, yhat),
               color = "Red", shape = 19, size = 2) +
    geom_segment(data = PointData, aes(x = x, xend = x, y = 0, yend = yhat),
                 color = "Red", linewidth = 1, linetype = "dashed") +
    labs(title = alternatename) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  if (Verbose == TRUE) {print(Views)}

  PointData <- PointData %>% select(-y)

  return(PointData)
}


#' Helper function splits Cluster into individual components
#'
#' @param x A data.frame containing column Cluster
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom dplyr relocate
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_with
#' @importFrom tidyr unnest_wider
#' @importFrom tidyr starts_with
#'
#' @return A value to be determined later
#' @noRd
RelativeBrightness <- function(x){
  Regular <- x %>% filter(!str_detect(Cluster, "-$"))

  if(nrow(Regular) > 0){
    Regular <- Regular %>% mutate(regular_split = str_split(
      as.character(Cluster), "-")) %>% relocate(regular_split, .after = Cluster)
    Regular <- Regular %>%  unnest_wider(regular_split, names_sep = "_") %>%
      rename_with(~ c("Detector1", "Detector2", "Detector3"), starts_with("regular"))

    Regular <- Regular %>% mutate(Detector1 = str_split(as.character(Detector1),
                                                        "_"))
    Regular <- Regular %>%  unnest_wider(Detector1, names_sep = "_") %>%
      rename_with(~ c("Detector1", "Detector1Value"), starts_with("Detector1"))

    Regular <- Regular %>% mutate(Detector2 = str_split(as.character(Detector2),
                                                        "_"))
    Regular <- Regular %>%  unnest_wider(Detector2, names_sep = "_") %>%
      rename_with(~ c("Detector2", "Detector2Value"), starts_with("Detector2"))

    Regular <- Regular %>% mutate(Detector3 = str_split(as.character(Detector3),
                                                        "_"))
    Regular <- Regular %>%  unnest_wider(Detector3, names_sep = "_") %>%
      rename_with(~ c("Detector3", "Detector3Value"), starts_with("Detector3"))
  }

  Exceptions <- x %>% filter(str_detect(Cluster, "-$"))

  if(nrow(Exceptions) > 0){

    Exceptions$Cluster <- gsub("-$", "", Exceptions$Cluster)
    Exceptions <- Exceptions %>% mutate(regular_split = str_split(as.character(
      Cluster), "-")) %>% relocate(regular_split, .after = Cluster)
    Exceptions <- Exceptions %>%  unnest_wider(regular_split, names_sep = "_") %>%
      rename_with(~ c("Detector1", "Detector2"), starts_with("regular"))

    Exceptions <- Exceptions %>% mutate(Detector1 = str_split(as.character(
      Detector1), "_"))
    Exceptions <- Exceptions %>%  unnest_wider(Detector1, names_sep = "_") %>%
      rename_with(~ c("Detector1", "Detector1Value"), starts_with("Detector1"))

    Exceptions <- Exceptions %>% mutate(Detector2 = str_split(as.character(
      Detector2), "_"))
    Exceptions <- Exceptions %>%  unnest_wider(Detector2, names_sep = "_") %>%
      rename_with(~ c("Detector2", "Detector2Value"), starts_with("Detector2"))

    Exceptions <- Exceptions %>% mutate(Detector3 = rep(NA, nrow(Exceptions))) %>%
      relocate(Detector3, .after = Detector2Value)
    Exceptions$Detector3 <- as.character(Exceptions$Detector3)

    Exceptions <- Exceptions %>% mutate(Detector3Value = rep(NA,
                                                             nrow(Exceptions))) %>% relocate(Detector3Value, .after = Detector3)
    Exceptions$Detector3Value <- as.character(Exceptions$Detector3Value)

  }

  Combined <- rbind(Regular, Exceptions)

  Combined$Detector1Value <- as.numeric(Combined$Detector1Value)
  Combined$Detector2Value <- as.numeric(Combined$Detector2Value)
  Combined$Detector3Value <- as.numeric(Combined$Detector3Value)

  Combined <- Combined %>% mutate(Brightness = rowSums(select(., Detector1Value,
                                                              Detector2Value, Detector3Value), na.rm = TRUE)) %>%
    relocate(Brightness, .after = Cluster)

  Combined <- Combined %>% mutate(Detector1Raw = rep(NA, nrow(Combined)),
                                  Detector2Raw = rep(NA, nrow(Combined)),
                                  Detector3Raw = rep(NA, nrow(Combined))) %>%
    relocate(Detector1Raw, .after = Detector1Value) %>%
    relocate(Detector2Raw, .after = Detector2Value) %>%
    relocate(Detector3Raw, .after = Detector3Value)

  Combined <- TheFill(Combined)

  return(Combined)
}



#' Internal for LuciernagaQC, creates .fcs files
#'
#' @param x The data.frame of Luciernaga data.
#' @param ff An individual cytoset object.
#' @param minimalfcscutoff A ratio indicating mininum of the total population needed to split off
#'  into own file, default is set to 0.05
#' @param AggregateName Passed final name with modifications from name
#' @param Brightness Whether to additionally return a brightness .csv to the outpath
#' @param outpath Location to export the fcs and .csv files to
#' @param OriginalStart Passed Argument indicating start column for Raw .fcs values
#' @param OrigingalEnd Passed argument indicating end column for raw .fcs values
#' @param stats Whether "median" or "mean", default is "median"
#' @param NegativeType Whether to append a negative pop. Args are "artificial", "internal" and "default"
#' @param TotalNegatives How many of the above rows to append, default is set to 500
#' @param Samples When Negative type = "Internal", the data.frame of averaged fluorescence per detector
#' @param ExportType Passed from above, set to "fcs" for fcs.file return
#'
#'
#' @importFrom flowCore parameters
#' @importFrom flowWorkspace keyword
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom utils write.csv
#'
#' @noRd
Genesis <- function(x, ff, minimalfcscutoff=0.05, AggregateName,
                    Brightness, outpath=NULL, OriginalStart, OriginalEnd,
                    stats = "median", NegativeType="default", TotalNegatives=500,
                    Samples=NULL, ExportType){

  # Replicate the Original FCS Parameters
  FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
  original_p <- parameters(FlowFrameTest)
  original_d <- keyword(FlowFrameTest)

  # Figure out what clusters to split from the file.
  x$Cluster <- factor(x$Cluster)
  ZZZ <- data.frame(table(x$Cluster))
  ZZZ <- ZZZ %>% arrange(desc(Freq))
  colnames(ZZZ)[1] <- "Cluster"
  colnames(ZZZ)[2] <- "Count"
  fcs_cutoff <- nrow(x)*minimalfcscutoff
  fcs_clusters <- ZZZ %>% filter(Count > fcs_cutoff) %>% pull(Cluster)

  Data <- x

  TheBrightness <- map(.x=fcs_clusters, .f=Luciernaga:::InternalGenesis, Data=Data, AggregateName=AggregateName,
                       outpath=outpath, OriginalStart=OriginalStart, OriginalEnd=OriginalEnd,
                       stats=stats, NegativeType=NegativeType, TotalNegatives=TotalNegatives,
                       Samples=Samples, ExportType=ExportType, parameters=original_p,
                       description=original_d) %>% bind_rows()

  message("TargetReached")

  if (Brightness == TRUE){
    RelativeBrightness <- Luciernaga:::RelativeBrightness(TheBrightness)
    CSVName <- paste0("RelativeBrightness", AggregateName, ".csv")
    CSVSpot <- file.path(outpath, CSVName)
    write.csv(RelativeBrightness, CSVSpot, row.names = FALSE)
  }
}


#' Internal for LuciernagaQC, creates .fcs files
#'
#' @param x Individual fluorescence cluster for filtering
#' @param Data The data.frame containing the many of the above
#' @param AggregateName Passed final name with modifications from name
#' @param outpath Location to export the fcs and .csv files to
#' @param OriginalStart Passed Argument indicating start column for Raw .fcs values
#' @param OriginalEnd Passed argument indicating end column for raw .fcs values
#' @param stats Whether "median" or "mean", default is "median"
#' @param NegativeType Whether to append a negative pop. Args are "artificial", "internal" and "default"
#' @param TotalNegatives How many of the above rows to append, default is set to 500
#' @param Samples When Negative type = "Internal", the data.frame of averaged fluorescence per detector
#' @param ExportType Passed from above, set to "fcs" for fcs.file return
#' @param parameters Passed parameters for .fcs creation
#' @param description Passed description for .fcs creation
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect one_of
#' @importFrom dplyr bind_cols
#' @importFrom dplyr relocate
#' @importFrom flowCore write.FCS
#'
#' @noRd
InternalGenesis <- function(x, Data, AggregateName, outpath=NULL, OriginalStart, OriginalEnd,
                            stats="median", NegativeType="default", TotalNegatives = 500,
                            Samples = NULL, ExportType, parameters, description){

  internalstrings <- c("-", "_")
  FCSname <- NameCleanUp(x, removestrings=internalstrings)
  FCSName <- paste(AggregateName, FCSname, sep = "_")
  #FCSName

  FCSSubset <- Data %>% filter(Cluster %in% x)

  # If Return Type No Add Ons
  RawFCSSubset <- FCSSubset %>% select(all_of(OriginalStart:OriginalEnd))

  HowBright <- AveragedSignature(RawFCSSubset, stats=stats)
  HowBright <- cbind(x, HowBright)
  colnames(HowBright)[1] <- "Cluster"
  #HowBright #Exported to bind_row with data.frame.
  #RawFCSSubset

  if (NegativeType == "artificial"){
    MeanFCS <- colMeans(RawFCSSubset)
    MeanFCS <- data.frame(t(MeanFCS), check.names = FALSE)
    mutateCols <- MeanFCS[,-grep("Time|FS|SC|SS|Original|W$|H$", names(MeanFCS))] %>% colnames(.)
    MeanFCS <- MeanFCS %>% mutate(across(all_of(mutateCols), ~ifelse(. >= 0, 0, .)))
    MeanFCS$Time <- round(MeanFCS$Time, 1)
    ArtificialNegative <- MeanFCS[rep(1, each = TotalNegatives),]
    rownames(ArtificialNegative) <- NULL
    FCSSubset <- rbind(RawFCSSubset, ArtificialNegative)
  }

  if (NegativeType == "samples"){
    if(!is.data.frame(Samples)){stop("Samples needs to be a single row of a data.frame
                                     for just the raw detectors")}
    SamplesCols <- colnames(Samples)
    MeanFCS <- colMeans(RawFCSSubset)
    MeanFCS <- data.frame(t(MeanFCS), check.names = FALSE)
    BackboneCols <- colnames(MeanFCS)
    Residual <- MeanFCS %>% select(-one_of(SamplesCols))
    Combined <- bind_cols(Residual, Samples)
    Combined <- Combined %>% relocate(all_of(BackboneCols))
    SampleNegative <- Combined[rep(1, each = TotalNegatives),]
    rownames(SampleNegative) <- NULL
    FCSSubset <- rbind(RawFCSSubset, SampleNegative)
  }

  if (NegativeType == "default"){
    FCSSubset <- RawFCSSubset
  }

  FCSSubset <- data.matrix(FCSSubset)
  new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=parameters, description=description)

  TheFileName <- paste(AggregateName, FCSname, sep="_")
  TheFileFCS <- paste0(TheFileName, ".fcs")
  if (is.null(outpath)) {outpath <- getwd()}
  fileSpot <- file.path(outpath, TheFileFCS)

  if (ExportType == "fcs") {write.FCS(new_fcs, filename = fileSpot, delimiter="#")}

  return(HowBright)
}


