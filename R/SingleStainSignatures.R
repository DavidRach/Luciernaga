#' Internal to LuciernagaQC
#'
#' @param x The Detector being mapped in
#' @param WorkAround1 The data.frame being used
#' @param Substraction Dictates what autofluorescence gets extracted. Options include "Internal", "Average"
#' @param ColsN Indicated end of Raw Value Columns
#' @param StartNormalizedMergedCol Indicated Start Normalized Columns
#' @param EndNormalizedMergedCol Indicated End Normalized Columns
#' @param Samples A data.frame row containing raw data of external Autofluorescence
#' @param Increments A numeric to round the normalized bins by. Default is 0.1
#' @param stats A passed param, usually "median" or "mean"
#' @param ratioSCcutoff The ratio dictating mininum size return for single colors, default is 0.01
#' @param TheMainAF A passed parameter dictating what detector column the autofluorescence is
#' @param AggregateName A passed parameter deriving from edits to sample.name and removestrings
#' @param Verbose Whether to return intermediate objects via print and plot for progress monitoring
#' @param SCData Whether to return "subtracted" or "raw" data at the end.
#' @param LocalMaximaRatio Height of peaks to proceed
#' @param SecondaryPeaks Number of Secondary Peaks, default is set to 2.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyselect where
#' @importFrom dplyr relocate
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_starts
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom dplyr left_join
#' @importFrom purrr map
#'
#' @noRd
SingleStainSignatures <- function(x, WorkAround1, ColsN, StartNormalizedMergedCol,
                                  EndNormalizedMergedCol, Samples, Increments,
                                  Subtraction = "Internal", stats, ratioSCcutoff=0.01, TheMainAF,
                                  AggregateName, Verbose = FALSE, SCData = "subtracted",
                                  LocalMaximaRatio, SecondaryPeaks){

  if (Subtraction == "Internal"){
  # Filtering the Non-Rounded Values First
  AutofluorescentSubset <- WorkAround1 %>% dplyr::filter(.data[[TheMainAF]] %in% 1.000) %>% arrange(desc(.data[[x]]))
  #nrow(AutofluorescentSubset)
  #DetectorPeakCounts(x=AutofluorescentSubset, StartN=65, EndN=128)

  SingleColorSubset <- WorkAround1 %>% dplyr::filter(.data[[x]] %in% 1.000) %>% arrange(desc(.data[[TheMainAF]]))
  #nrow(SingleColorSubset)
  #DetectorPeakCounts(x=SingleColorSubset, StartN=65, EndN=128)

  # Rounding to Desired Increments, no double dipping this way.
  AutofluorescentSubsetB <- AutofluorescentSubset %>% select(all_of(
    (StartNormalizedMergedCol+1):(EndNormalizedMergedCol+1))) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments)) %>%
    mutate(Backups = AutofluorescentSubset$Backups) %>% relocate(Backups, .before=1)

  SingleColorSubsetB <- SingleColorSubset %>% select(all_of(
    (StartNormalizedMergedCol+1):(EndNormalizedMergedCol+1))) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments)) %>%
    mutate(Backups = SingleColorSubset$Backups) %>% relocate(Backups, .before=1)

  #nrow(AutofluorescentSubsetB) #nrow(SingleColorSubsetB)

  #x <- Retained[1] #TheMainAF

  AutofluorescenceNamed <- Luciernaga:::SignatureCluster(Arg1=TheMainAF, Arg2=x, data=AutofluorescentSubsetB)
  SingleColorNamed <- Luciernaga:::SignatureCluster(Arg1=x, Arg2=TheMainAF, data=SingleColorSubsetB)
  SignatureData <- bind_rows(AutofluorescenceNamed, SingleColorNamed)
  TheClusters <- data.frame(table(SignatureData$Cluster), check.names=FALSE)
  colnames(TheClusters)[1] <- "Clusters"
  colnames(TheClusters)[2] <- "Counts"

  AF <- TheClusters %>% filter(str_starts(Clusters, TheMainAF)) %>% arrange(Clusters)
  SC <- TheClusters %>% filter(str_starts(Clusters, x)) %>% arrange(desc(Clusters))
  Total <- bind_rows(AF, SC)

  #TheOrder <- Total$Clusters
  #Total$Clusters <- as.character(Total$Clusters)
  #Total$Clusters <- factor(Total$Clusters, levels = TheOrder)

  Subtotal <- Total %>% filter(str_starts(Clusters, x))

  TheOrder <- Subtotal$Clusters
  Subtotal$Clusters <- as.character(Subtotal$Clusters)
  Subtotal$Clusters <- factor(Subtotal$Clusters, levels = TheOrder)

  InitialPlot <- ggplot(Subtotal, aes(x = Clusters, y = Counts)) + geom_bar(stat = "identity") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(title=paste(x, AggregateName, "Initial"), sep=" ")

  if (Verbose == TRUE) {print(Total)
                        InitialPlot
  }

  ################################
  # Autofluorescence Subtraction #
  ################################

  AFforSubtraction <- Total %>% filter(str_starts(Clusters, TheMainAF)) %>% arrange(desc(Counts)) %>% slice(1) %>% pull(Clusters) #%>% as.character()
  }


  if (Subtraction == "Internal"){

    AF_Choice <- SignatureData %>% filter(Cluster %in% AFforSubtraction) %>% select(Backups)
    Restored <- left_join(AF_Choice, WorkAround1, by="Backups")
    Restored <- Restored %>% select(-Backups)
    Restored <- Restored %>% select(all_of(1:ColsN))
    BackupNames2 <- colnames(Restored)

    Samples <- AveragedSignature(x=Restored, stats=stats)

    StoredBackups <- WorkAround1 %>% select(Backups)
    Data <- WorkAround1 %>% select(-Backups)
    Data <- Data %>% select(all_of(1:ColsN))
    Samples_replicated <- Samples[rep(1, each = nrow(Data)),]
    Test <- Data - Samples_replicated

    #We have just subtracted targeted AF from everything at this point.

  } else if (Subtraction == "Average"){
    #Samples #Either Original, desiredAF driven alternative, or external provided.
    StoredBackups <- WorkAround1 %>% select(Backups)
    Data <- WorkAround1 %>% select(-Backups)
    Data <- Data %>% select(all_of(1:ColsN))
    Samples_replicated <- Samples[rep(1, each = nrow(Data)),]
    Test <- Data - Samples_replicated

    #We have just subtracted averaged AF from everything at this point.
  }


  if (Verbose == TRUE){
    TheTotal <- nrow(Test) * ncol(Test)
    BelowZero <- sum(apply(Test, 2, function(x) x < 0))
    message(round(BelowZero/TheTotal,2), " of all events post AF subtraction were negative and will be rounded to 0")
  }

  Test[Test < 0] <- 0 #Check here for negatives...
  AA <- do.call(pmax, Test)
  Normalized2 <- Test/AA
  #Normalized2 <- round(Normalized2, 1)
  #num_na <- sum(is.na(Normalized2))
  Normalized2[is.na(Normalized2)] <- 0
  colnames(Normalized2) <- gsub("-A", "", colnames(Normalized2))
  Counts2 <- colSums(Normalized2 == 1.0)
  RevisedPeakDetectorCounts <- data.frame(Fluors = names(Counts2), Counts = Counts2)
  rownames(RevisedPeakDetectorCounts) <- NULL
  RevisedPeakDetectorCounts <- RevisedPeakDetectorCounts %>% arrange(desc(Counts))

  if (Verbose == TRUE){
    print(RevisedPeakDetectorCounts)
  }

  ### Its here that we rebind the subtracted data.
  ### Replacing Test with Data would retain original values for Raw.

  if (SCData == "subtracted"){WorkAround2 <- cbind(StoredBackups, Test, Normalized2)}
  if (SCData == "raw"){WorkAround2 <- cbind(StoredBackups, Data, Normalized2)}

  ################
  # Reclustering #
  ################

  SingleColorSubset2 <- WorkAround2 %>% dplyr::filter(.data[[x]] %in% 1.000) %>% arrange(desc(.data[[TheMainAF]]))
  #nrow(SingleColorSubset2)
  #DetectorPeakCounts(x=SingleColorSubset2, StartN=65, EndN=128)

  SingleColorSubset2B <- SingleColorSubset2 %>% select(all_of(
    (StartNormalizedMergedCol+1):(EndNormalizedMergedCol+1))) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments)) %>%
    mutate(Backups = SingleColorSubset2$Backups) %>% relocate(Backups, .before=1)

  SingleColorData <- Luciernaga:::SignatureCluster(Arg1=x, Arg2=TheMainAF, data=SingleColorSubset2B)
  TheClusters <- data.frame(table(SingleColorData$Cluster), check.names=FALSE)
  colnames(TheClusters)[1] <- "Clusters"
  colnames(TheClusters)[2] <- "Counts"
  Total <- TheClusters %>% filter(str_starts(Clusters, x)) %>% arrange(desc(Clusters))
  TheOrder <- Total$Clusters
  Total$Clusters <- as.character(Total$Clusters)
  Total$Clusters <- factor(Total$Clusters, levels = TheOrder)

  FinalPlot <- ggplot(Total, aes(x = Clusters, y = Counts)) + geom_bar(stat = "identity") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(title=paste(x, AggregateName, "Final"), sep=" ")

  if(Verbose == TRUE){
  InitialPlot + FinalPlot
  }

  ##########################
  # Local Maxima Iteration #
  ##########################
  Cutoff <- sum(TheClusters$Counts)*ratioSCcutoff
  MainClusters <- TheClusters %>% filter(Counts >= Cutoff) %>% select(Clusters) %>%
    pull() %>% as.character()

  SingleColorSubset <- SingleColorData %>% select(Backups, Cluster)
  FinalNormalized <- SingleColorData %>% select(-Backups, -Cluster)
  Ready <- left_join(SingleColorSubset, WorkAround2, by="Backups")
  Ready <- Ready %>% select(-all_of((StartNormalizedMergedCol+2):(EndNormalizedMergedCol+2)))
  Readied <- cbind(Ready, FinalNormalized)
  TheReady <- Readied %>% relocate(Cluster, .after="R8")

  TheDetector <- x

  #table(Ready$Cluster)
  #x <- MainClusters[1]

  AllData <- map(.x=MainClusters, .f=ClusterIteration, data=TheReady,
                 TheDetector=TheDetector,
                 StartNormalizedMergedCol=StartNormalizedMergedCol,
                 EndNormalizedMergedCol=EndNormalizedMergedCol,
                 ColsN=ColsN, AggregateName=AggregateName, Verbose=Verbose,
                 LocalMaximaRatio=LocalMaximaRatio, SecondaryPeaks=SecondaryPeaks) %>% bind_rows()

  NewTable <- data.frame(table(AllData$Cluster), check.names=FALSE)
  colnames(NewTable)[1] <- "Cluster"
  colnames(NewTable)[2] <- "Count"

  #NewTable %>% arrange(desc(Count))

  NewCutoff <- sum(NewTable$Count)*ratioSCcutoff
  GrabThese <- NewTable %>% filter(Count >= NewCutoff) %>% select(Cluster) %>%
    pull() %>% as.character()

  FinalData <- AllData %>% filter(Cluster %in% GrabThese)

  #SendThese <- list(Samples, FinalData)

  return(FinalData)
}


#' Internal for LuciernagaQC SingleStainSignatures
#'
#' @param x A cluster identity in the cluster column
#' @param data A data.frame
#' @param StartNormalizedMergedCol Indicated Start Normalized Columns
#' @param EndNormalizedMergedCol Indicated End Normalized Columns
#' @param ColsN Indicated end of Raw Value Columns
#' @param AggregateName The sample.name derrived name
#' @param Verbose Whether to return intermediate objects
#' @param LocalMaximaRatio Height of peaks to proceed
#' @param SecondaryPeaks Number of Secondary Peaks, default is set to 2.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom dplyr case_when
#' @importFrom dplyr near
#' @importFrom utils head
#'
#' @noRd
ClusterIteration <- function(x, data, TheDetector, StartNormalizedMergedCol,
                             EndNormalizedMergedCol, ColsN, AggregateName, Verbose,
                             LocalMaximaRatio = 0.15, SecondaryPeaks = 2){

  subset <- data %>% filter(Cluster %in% x)
  StashedIDs <- subset %>% select(Backups)
  TheNormalized <- subset %>% select(-Backups) %>% select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol))
  MyRawData <- subset %>% select(-Backups) %>% select(all_of(1:ColsN))

  #Preparation for Local Maxima
  Conversion <- data.frame(t(TheNormalized), check.names = FALSE)
  Conversion <- cbind(Detectors = rownames(Conversion), Conversion)
  rownames(Conversion) <- NULL

  #Preparing Detector Stand Ins for left_join
  Decoys <- Conversion %>% select(Detectors)
  Decoys <- Decoys %>% mutate(TheDetector = 1:nrow(Decoys)) %>% relocate(
    TheDetector, .before = Detectors)

  #Deriving an average y-vector for local maxima
  Conversion <- Conversion %>% mutate(TheSums = rowSums(.[2:ncol(.)],
                                                        na.rm = TRUE) /(ncol(Conversion) - 1)) %>% relocate(
                                                          TheSums, .after = Detectors)
  Conversion$Detectors <- 1:nrow(Conversion)
  LocalX <- Conversion$Detectors
  LocalY <- Conversion$TheSums

  #I made it export, now just need to rebuild, then remove extra :
  alternatename <- AggregateName

  PointData <- Luciernaga:::LocalMaxima(theX = LocalX, theY = LocalY, therepeats = 3,
                                                w = 3, span = 0.11, alternatename = alternatename,
                                                Verbose = Verbose)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  Newest2 <- PointData %>% filter(TheHeight > LocalMaximaRatio) %>% arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)

  if (any(These %in% TheDetector)) {These <- These[These != TheDetector]}

  if(length(These) == 0){if (Verbose == TRUE) {message("Solitary Peak")}
  } else if (length(These) > SecondaryPeaks) {
    if (Verbose == TRUE) {message("More than ", SecondaryPeaks+1, " peaks. Abbreviated.")}
    These <- head(These, SecondaryPeaks)
  }

  MyData <- cbind(StashedIDs, MyRawData, TheNormalized)
  MyData$Cluster <- paste(TheDetector, "10-", sep = "_")

  if (length(These) > 3){stop("Only currently set up to handle up to 4 fluorescence peaks per fluorophore")
  } else if (length(These) == 3){second <- These[[1]]
  third <- These[[2]]
  fourth <- These[[3]]
  } else if (length(These) == 2){second <- These[[1]]
  third <- These[[2]]
  } else if (length(These) == 1){second <- These[[1]]
  } else if (length(These) == 0){
    if (Verbose == TRUE) {message("No second peak")}}

  if (length(These) >= 1){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[second]], 0.0) ~ paste0(MyData$Cluster, second, "_00-"),
    near(MyData[[second]], 0.1) ~ paste0(MyData$Cluster, second, "_01-"),
    near(MyData[[second]], 0.2) ~ paste0(MyData$Cluster, second, "_02-"),
    near(MyData[[second]], 0.3) ~ paste0(MyData$Cluster, second, "_03-"),
    near(MyData[[second]], 0.4) ~ paste0(MyData$Cluster, second, "_04-"),
    near(MyData[[second]], 0.5) ~ paste0(MyData$Cluster, second, "_05-"),
    near(MyData[[second]], 0.6) ~ paste0(MyData$Cluster, second, "_06-"),
    near(MyData[[second]], 0.7) ~ paste0(MyData$Cluster, second, "_07-"),
    near(MyData[[second]], 0.8) ~ paste0(MyData$Cluster, second, "_08-"),
    near(MyData[[second]], 0.9) ~ paste0(MyData$Cluster, second, "_09-"),
    near(MyData[[second]], 1.0) ~ paste0(MyData$Cluster, second, "_10-")))
  }

  if(length(These) >= 2){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[third]], 0.0) ~ paste0(MyData$Cluster, third, "_00"),
    near(MyData[[third]], 0.1) ~ paste0(MyData$Cluster, third, "_01"),
    near(MyData[[third]], 0.2) ~ paste0(MyData$Cluster, third, "_02"),
    near(MyData[[third]], 0.3) ~ paste0(MyData$Cluster, third, "_03"),
    near(MyData[[third]], 0.4) ~ paste0(MyData$Cluster, third, "_04"),
    near(MyData[[third]], 0.5) ~ paste0(MyData$Cluster, third, "_05"),
    near(MyData[[third]], 0.6) ~ paste0(MyData$Cluster, third, "_06"),
    near(MyData[[third]], 0.7) ~ paste0(MyData$Cluster, third, "_07"),
    near(MyData[[third]], 0.8) ~ paste0(MyData$Cluster, third, "_08"),
    near(MyData[[third]], 0.9) ~ paste0(MyData$Cluster, third, "_09"),
    near(MyData[[third]], 1.0) ~ paste0(MyData$Cluster, third, "_10")))
  }

  if(length(These) >= 3){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[fourth]], 0.0) ~ paste0(MyData$Cluster, fourth, "_00"),
    near(MyData[[fourth]], 0.1) ~ paste0(MyData$Cluster, fourth, "_01"),
    near(MyData[[fourth]], 0.2) ~ paste0(MyData$Cluster, fourth, "_02"),
    near(MyData[[fourth]], 0.3) ~ paste0(MyData$Cluster, fourth, "_03"),
    near(MyData[[fourth]], 0.4) ~ paste0(MyData$Cluster, fourth, "_04"),
    near(MyData[[fourth]], 0.5) ~ paste0(MyData$Cluster, fourth, "_05"),
    near(MyData[[fourth]], 0.6) ~ paste0(MyData$Cluster, fourth, "_06"),
    near(MyData[[fourth]], 0.7) ~ paste0(MyData$Cluster, fourth, "_07"),
    near(MyData[[fourth]], 0.8) ~ paste0(MyData$Cluster, fourth, "_08"),
    near(MyData[[fourth]], 0.9) ~ paste0(MyData$Cluster, fourth, "_09"),
    near(MyData[[fourth]], 1.0) ~ paste0(MyData$Cluster, fourth, "_10")))
  }

  return(MyData)


}


#' Internal for LuciernagaQC SingleStainSignatures
#'
#' @param Arg1 The first Detector Column
#' @param Arg2 The second Detector Column
#' @param data A data.frame
#'
#' @noRd
SignatureCluster <- function(Arg1, Arg2, data){
  data <- data %>% mutate(Cluster = case_when(
    near(data[[Arg1]], 0.0) ~ paste(data$Cluster, Arg1, "_00-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.1) ~ paste(data$Cluster, Arg1, "_01-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.2) ~ paste(data$Cluster, Arg1, "_02-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.3) ~ paste(data$Cluster, Arg1, "_03-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.4) ~ paste(data$Cluster, Arg1, "_04-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.5) ~ paste(data$Cluster, Arg1, "_05-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.6) ~ paste(data$Cluster, Arg1, "_06-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.7) ~ paste(data$Cluster, Arg1, "_07-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.8) ~ paste(data$Cluster, Arg1, "_08-", sep = "", collapse = NULL),
    near(data[[Arg1]], 0.9) ~ paste(data$Cluster, Arg1, "_09-", sep = "", collapse = NULL),
    near(data[[Arg1]], 1.0) ~ paste(data$Cluster, Arg1, "_10-", sep = "", collapse = NULL)))

  Second <- data %>% mutate(Cluster = case_when(
    near(data[[Arg2]], 0.0) ~ paste(data$Cluster, Arg2, "_00-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.1) ~ paste(data$Cluster, Arg2, "_01-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.2) ~ paste(data$Cluster, Arg2, "_02-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.3) ~ paste(data$Cluster, Arg2, "_03-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.4) ~ paste(data$Cluster, Arg2, "_04-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.5) ~ paste(data$Cluster, Arg2, "_05-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.6) ~ paste(data$Cluster, Arg2, "_06-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.7) ~ paste(data$Cluster, Arg2, "_07-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.8) ~ paste(data$Cluster, Arg2, "_08-", sep = "", collapse = NULL),
    near(data[[Arg2]], 0.9) ~ paste(data$Cluster, Arg2, "_09-", sep = "", collapse = NULL),
    near(data[[Arg2]], 1.0) ~ paste(data$Cluster, Arg2, "_10-", sep = "", collapse = NULL)))

  return(Second)
}

