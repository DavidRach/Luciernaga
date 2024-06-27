#' Main Luciernaga Function, normalized on individual cell level
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param subsets The Gating Hierarchy level you will be sampling at
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param group.name Keyword variable which groups are stored (ex. "GROUPNAME")
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param experiment.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param stats Whether to take "mean" or "median"
#' @param Kept Whether "Raw" or "Normalized" values are retained in the
#' Luciernaga object.
#' @param external An external autofluorescence to subtract from single colors.
#' @param sourcelocation Location where .fcs creation file is stored
#' @param outpath  Location where created .fcs and .csv files are sent
#' @param artificial Whether an artificial 0 population should be added for
#'  a background autofluorescence stand in.
#' @param fcsexport Whether to export .fcs files, TRUE or FALSE
#' @param mainAF Main Autofluorescence Detector (ex. "V7-A")
#' @param AFOverlap Name of data.frame containing the Autofluorescence
#' overlap of individual fluorophores for exclusion
#' @param Beads  Whether the sample is Beads.
#' @param Brightness Whether sum of detectors should be returned.
#' @param Unstained Whether the sample is Unstained.
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize_all
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom dplyr left_join
#' @importFrom dplyr case_when
#' @importFrom dplyr rename
#' @importFrom BiocGenerics nrow
#' @importFrom flowWorkspace keyword
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom purrr map
#' @importFrom purrr set_names
#'
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL

#x <- gs[2]
#subsets = "lymph"
#sample.name = "GUID"
#removestrings = c("DR_", " (Cells)")
#mainAF = "V7-A"
#AFOverlap = AFOverlap
#stats = "median"
#Unstained = FALSE
#Beads = FALSE
#Verbose = TRUE
#external = NULL
#fcsexport = TRUE
#sourcelocation = "Genesis.R"
#outpath = MainOutPath
#artificial = TRUE
#Brightness = TRUE

#group.name = "GROUPNAME"
#experiment = NULL
#experiment.name = "$DATE"
#Kept = "Normalized"
# Kept, outpath, artificial, Brightness

Utility_SingleColorQC <- function(x, subsets, sample.name, removestrings, experiment = NULL, experiment.name = NULL,
                                  mainAF, AFOverlap, stats, Unstained=FALSE, Beads=FALSE, Verbose = FALSE,
                                  external = NULL, fcsexport, sourcelocation, outpath, artificial, Brightness=FALSE){

  ##############
  # Name Setup #
  ##############

  # Retrieving .fcs file name
  name <- keyword(x, sample.name)
  name <- NameCleanUp(name, ".fcs")

  # Predicting Reference Control Type

  if (Beads != TRUE) {
    if(str_detect(name, "(Cells)")){Type <- "Cells"
      } else if(str_detect(name, "(Beads)")){
        stop("Bead reference controls are not currently supported. Sorry! -David")
          } else {Type <- "Unknown"}
  } else {
    if(str_detect(name, "(Cells)")){Type <- "Cells"
      } else if(str_detect(name, "(Beads)")){Type <- "Beads"
        } else {Type <- "Unknown"}
  }

  name <- NameCleanUp(name, removestrings)

  # Specifying Unstained Work Around.
  if (Unstained == TRUE) {name <- paste0(name, "_Unstained")}

  # Cleaning up Cells in Case Not Specified
  name <- NameCleanUp(name, " (Cells)")

  if(exists("experiment")) {experiment <- experiment
  } else if (exists("experiment.name")) {experiment <- keyword(x, experiment.name)
  } else {experiment <- NULL}

  if (!is.null(experiment)){
    AggregateName <- paste0(name, experiment)
  } else {AggregateName <- name}

  # Removing Non-Characters and Spaces.
  ExtraSpacers <- c(" ", "_", "-", ".", "(", ")")
  AggregateName <- NameCleanUp(AggregateName, removestrings=ExtraSpacers)

  #if (!is.null(external)){external1 <- external} else {external1 <- NULL}

  ###############
  # Exprs Setup #
  ###############

  #Retrieving the exprs data for the target population
  ff <- gs_pop_get_data(x, subsets)
  startingcells <- nrow(ff)[[1]]
  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  #For Future Column Reordering
  OriginalColumnsVector <- colnames(DF)
  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  #For Future Row Reordering
  Backups <- DF %>% mutate(Backups = 1:nrow(DF)) %>% select(Backups)

  #Stashing Away Time FSC SSC For Later Use
  StashedDF <- DF[,grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  #Consolidating Columns Going Forward
  CleanedDF <- DF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  BackupNames <- colnames(CleanedDF)
  n <- CleanedDF

  #################
  # Peak Detector #
  #################

  # Enumerating Negative Values Ratio
  if (Verbose == TRUE){
  TotalCells <- nrow(n) * ncol(n)
  BelowZero <- sum(apply(n, 2, function(x) x < 0))
  message(round(BelowZero/TotalCells,2), " of all events were negative and rounded to 0")
  }

  #Normalizing By Peak Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 1)
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  #Towards Generalizing
  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsNormalized*2

  #Bringing Together Raw And Normalized and deriving Autofluorescence Negative
  WorkAround <- cbind(n, Normalized)

  ####################
  # Autofluorescence #
  ####################

  #Retrieving Main Auto fluorescent Channels signature
  if (Unstained == FALSE){MainAF <- mainAF
  MainAF <- gsub("-A", "", MainAF)
  } else{ ### NOTE ###Code this exception in later.
    MainAF <- mainAF
    MainAF <- gsub("-A", "", MainAF)
  }

  # Grabbing Main AF Peak and the Associated Raw Values
  This <- WorkAround %>% filter(.data[[MainAF]] == 1) %>%
    select(all_of(1:ColsN))

  # Deriving Middle Autofluorescence Measurement
  if(stats == "mean"){Samples <- This %>% summarize_all(mean)
  } else if (stats == "median"){Samples <- This %>%
    summarize_all(median) #%>% select(-Backups)
  } else(stop("Please specify stats parameter mean or median"))

  ################
  # Fluorophores #
  ################

  #Deriving Peak Detector Counts and Detectors of Interest
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*0.0075 #Possible Parameter Add Here
  Detectors <- PeakDetectorCounts %>% filter(Counts > cutoff) %>%
    arrange(desc(Counts))

  ################################################################
  # Bringing in known AF detectors, and overlapping fluorophores #
  ################################################################

  AFData <- AFOverlap
  AFChannels <- AFData %>% filter(Fluorophore %in% "Unstained") %>%
    pull(MainDetector) %>% str_split(",", simplify = TRUE)
  AFChannels <- AFChannels[1,]
  AFChannels <- gsub("-A", "", AFChannels)

  SCData <- AFData %>% filter(Fluorophore != "Unstained")
  SCData$Fluorophore <- gsub("-A", "", SCData$Fluorophore)
  TroubleChannels <- SCData %>% pull(Fluorophore)

  # Handling the Exceptions
  TroubleChannelExclusion <- function(x, SCData, MainDetector, AFChannels){
    Internal <- SCData %>% filter(Fluorophore %in% x) %>% pull(MainDetector)
    Internal <- gsub("-A", "", Internal)
    Exclusion <- setdiff(AFChannels, Internal)
    return(Exclusion)
  }

  results <- map(.x=TroubleChannels, .f=TroubleChannelExclusion, SCData=SCData, MainDetector=MainDetector, AFChannels=AFChannels) %>%
    set_names(TroubleChannels)

  matching_names <- names(results)[str_detect(name, names(results))]
  if (length(matching_names) > 0) {ExclusionList <- results[[matching_names[1]]]
  Retained <- Detectors %>% filter(!Fluors %in% ExclusionList) %>% pull(Fluors)
  } else if (str_detect(name, "Unstained")){Retained <- Detectors %>%
    pull(Fluors)
  } else {Retained <- Detectors %>% filter(!Fluors %in% AFChannels) %>%
    pull(Fluors)}

  if (length(Retained) == 0) {
    stop("There were no Retained detectors in ", name)
  }

  TheMainAF <- Detectors %>% filter(!Fluors %in% Retained) %>% slice(1) %>% pull(Fluors)

  if(Beads == TRUE){Retained <- Retained[[1]]}

  #####################################
  # We resume our regular programming #
  #####################################

  WorkAround1 <- WorkAround %>% mutate(Backups = Backups$Backups) %>%
    relocate(Backups, .before = 1) #This will change the start/end count

  if (str_detect(name, "Unstained")){

    RetainedDF <- map(.x= Retained, .f=UnstainedSignatures,
                      WorkAround1=WorkAround1, alternatename=AggregateName,
                      ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
                      EndNormalizedMergedCol=EndNormalizedMergedCol) %>% bind_rows()

  } else {

    RetainedDF <- map(.x= Retained, .f=SingleStainSignatures,
                      WorkAround1=WorkAround1, alternatename=AggregateName,
                      ColsN=ColsN, StartNormalizedMergedCol=StartNormalizedMergedCol,
                      EndNormalizedMergedCol=EndNormalizedMergedCol, Samples=Samples, name=name,
                      results=results, stats=stats) %>% bind_rows()
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

  if (fcsexport == TRUE){source(sourcelocation, local = TRUE)}

  return(Reintegrated1)
}

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


ModernSingleStainSignatures <- function(x, WorkAround1, alternatename, ColsN, StartNormalizedMergedCol,
                                        EndNormalizedMergedCol, Samples, name, results, Increments=0.1,
                                        Subtraction = "Internal", stats, RatioCutoff=0.01){

  WorkAround1b <- WorkAround1 %>% select(all_of(
    (StartNormalizedMergedCol+1):(EndNormalizedMergedCol+1))) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments)) %>%
    mutate(Backups = WorkAround1$Backups) %>% relocate(Backups, .before=1)

  #x <- Retained[1]
  #TheMainAF

  # The Arguments would be mapping elements in Retained, therefore Retained below would be swapped out for the x argument.
  AutofluorescentSubset <- WorkAround1b %>% dplyr::filter(.data[[TheMainAF]] == 1.000) %>% arrange(desc(.data[[x]]))
  SingleColorSubset <- WorkAround1b %>% dplyr::filter(.data[[x]] == 1.000) %>% arrange(desc(.data[[TheMainAF]])) %>%
    dplyr::filter(!.data[[TheMainAF]] == 1.00)  #To avoid double dipping

  AutofluorescenceNamed <- SignatureCluster(Arg1=TheMainAF, Arg2=x, data=AutofluorescentSubset)
  SingleColorNamed <- SignatureCluster(Arg1=x, Arg2=TheMainAF, data=SingleColorSubset)
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
  InitialPlot

  ################################
  # Autofluorescence Subtraction #
  ################################

  AFforSubtraction <- Total %>% filter(str_starts(Clusters, TheMainAF)) %>% arrange(desc(Counts)) %>% slice(1) %>% pull(Clusters) #%>% as.character()

  if (Subtraction == "Internal"){
    AF_Choice <- SignatureData %>% filter(Cluster == AFforSubtraction) %>% select(Backups)
    Restored <- left_join(AF_Choice, WorkAround1, by="Backups")
    Restored <- Restored %>% select(-Backups)
    Restored <- Restored %>% select(all_of(1:ColsN))
    BackupNames2 <- colnames(Restored)

    if(stats == "mean"){Samples <- Restored %>% summarize_all(mean)
    } else if (stats == "median"){Samples <- Restored %>%
      summarize_all(median)
    } else(stop("Please specify stats parameter mean or median"))

    StoredBackups <- WorkAround1 %>% select(Backups)
    Data <- WorkAround1 %>% select(-Backups)
    Data <- Data %>% select(all_of(1:ColsN))
    Samples_replicated <- Samples[rep(1, each = nrow(Data)),]
    Test <- Data - Samples_replicated

    Test[Test < 0] <- 0 #Check here for negatives...
    AA <- do.call(pmax, Test)
    Normalized2 <- Test/AA
    Normalized2 <- round(Normalized2, 1)
    #num_na <- sum(is.na(Normalized2))
    Normalized2[is.na(Normalized2)] <- 0
    colnames(Normalized2) <- gsub("-A", "", colnames(Normalized2))
    Counts2 <- colSums(Normalized2 == 1.0)
    #Counts2
    WorkAround2 <- cbind(StoredBackups, Test, Normalized2)

    ################
    # Reclustering #
    ################

    WorkAround2b <- WorkAround2 %>% select(all_of(
      (StartNormalizedMergedCol+1):(EndNormalizedMergedCol+1))) %>%
      mutate(across(where(is.numeric), ~ ceiling(. / Increments) * Increments)) %>%
      mutate(Backups = WorkAround2$Backups) %>% relocate(Backups, .before=1)

    SingleColorSubset <- WorkAround2b %>% dplyr::filter(.data[[x]] == 1.000) %>% arrange(desc(.data[[TheMainAF]])) #%>%
      #dplyr::filter(!.data[[TheMainAF]] == 1.00)  #Not Double Dipping Anymore Technically
    SingleColorData <- SignatureCluster(Arg1=x, Arg2=TheMainAF, data=SingleColorSubset)
    TheClusters <- data.frame(table(SingleColorData$Cluster), check.names=FALSE)
    colnames(TheClusters)[1] <- "Clusters"
    colnames(TheClusters)[2] <- "Counts"
    Total <- TheClusters %>% filter(str_starts(Clusters, x)) %>% arrange(desc(Clusters))
    TheOrder <- Total$Clusters
    Total$Clusters <- as.character(Total$Clusters)
    Total$Clusters <- factor(Total$Clusters, levels = TheOrder)

    FinalPlot <- ggplot(Total, aes(x = Clusters, y = Counts)) + geom_bar(stat = "identity") + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(title=paste(x, AggregateName, "Final"), sep=" ")
    #InitialPlot
    InitialPlot + FinalPlot

    ##########################
    # Local Maxima Iteration #
    ##########################
    Cutoff <- sum(TheClusters$Counts)*RatioCutoff
    MainClusters <- TheClusters %>% filter(Counts >= Cutoff) %>% select(Clusters) %>%
      pull() %>% as.character()

    SingleColorSubset <- SingleColorData %>% select(Backups, Cluster)
    Ready <- left_join(SingleColorSubset, WorkAround2, by="Backups") %>% relocate(Cluster, .after="R8")

    ClusterIteration <- function(x, data){
      subset <- data %>% filter(Cluster %in% x)






    }

    #x <- MainClusters[1]
    map(.x=MainClusters, .f=ClusterIteration, data=Ready)



    AF_Choice <- left_join(SingleColorData) %>% filter(Cluster == AFforSubtraction) %>% select(Backups)


  } else if (Subtraction == "Average"){

    stop("Sorry, still working on this. -David ")

  } else if (Subtraction == "External"){

    stop("Sorry, still working on this. -David ")
  }
}



  MyData <- WorkAround2 %>% select(all_of(
    StartNormalizedMergedCol:EndNormalizedMergedCol)) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / 0.2) * 0.2))
  MyRawData <- WorkAround2 %>% select(all_of(1:ColsN))

  #Preparation for Local Maxima
  Conversion <- data.frame(t(MyData), check.names = FALSE)
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
  alternatename <- alternatename
  PointData <- Luciernaga::Utility_LocalMaxima(theX = LocalX, theY = LocalY,
                                               therepeats = 3, w = 3, span = 0.11, alternatename = alternatename)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  Newest2 <- PointData %>% filter(TheHeight > 0.15) %>% arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)
  if (any(These %in% x)) {These <- These[These != x]}
  if (length(These) > 2) {These <- head(These, 2)} #If we wanted to institute a number of peaks argument, it would be here.

  MyData <- cbind(StashedIDs, MyRawData, MyData)
  MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

  if (length(These) == 2){second <- These[[1]]
  third <- These[[2]]
  } else if(length(These) == 1){second <- These[[1]]}

  if(exists("second")){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[second]], 0.0) ~ paste(MyData$Cluster, second, "_00-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.2) ~ paste(MyData$Cluster, second, "_02-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.4) ~ paste(MyData$Cluster, second, "_04-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.6) ~ paste(MyData$Cluster, second, "_06-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.8) ~ paste(MyData$Cluster, second, "_08-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 1.0) ~ paste(MyData$Cluster, second, "_10-",
                                        sep = "", collapse = NULL)))
  }

  if(exists("third")){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[third]], 0.0) ~ paste(MyData$Cluster, third, "_00",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.2) ~ paste(MyData$Cluster, third, "_02",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.4) ~ paste(MyData$Cluster, third, "_04",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.6) ~ paste(MyData$Cluster, third, "_06",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.8) ~ paste(MyData$Cluster, third, "_08",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 1.0) ~ paste(MyData$Cluster, third, "_10",
                                       sep = "", collapse = NULL)))
  }

  return(MyData)
}




UnstainedSignatures <- function(x, WorkAround1, alternatename, ColsN, StartNormalizedMergedCol, EndNormalizedMergedCol){
  MySubset <- WorkAround1 %>% dplyr::filter(.data[[x]] == 1.000)
  StashedIDs <- MySubset %>% select(Backups)
  MySubset <- MySubset %>% select(-Backups)
  DetectorName <- x
  MyData <- MySubset %>% select(all_of(
    StartNormalizedMergedCol:EndNormalizedMergedCol)) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / 0.2) * 0.2))
  MyRawData <- MySubset %>% select(all_of(1:ColsN))

  #Preparation for Local Maxima
  Conversion <- data.frame(t(MyData), check.names = FALSE)
  Conversion <- cbind(Detectors = rownames(Conversion), Conversion)
  rownames(Conversion) <- NULL

  #Preparing Detector Stand Ins for left_join
  Decoys <- Conversion %>% select(Detectors)
  Decoys <- Decoys %>% mutate(TheDetector = 1:nrow(Decoys)) %>%
    relocate(TheDetector, .before = Detectors)

  #Deriving an average y-vector for local maxima
  Conversion <- Conversion %>% mutate(TheSums = rowSums(.[2:ncol(.)],
                                                        na.rm = TRUE) /(ncol(Conversion) - 1)) %>% relocate(TheSums,
                                                                                                            .after = Detectors)
  Conversion$Detectors <- 1:nrow(Conversion)
  LocalX <- Conversion$Detectors
  LocalY <- Conversion$TheSums

  #I made it export, now just need to rebuild, then remove extra :
  alternatename <- alternatename
  PointData <- Luciernaga::Utility_LocalMaxima(theX = LocalX, theY = LocalY,
                                                therepeats = 3, w = 3, span = 0.11, alternatename = alternatename)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  Newest2 <- PointData %>% filter(TheHeight > 0.15) %>% arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)
  if (any(These %in% x)) {These <- These[These != x]}
  if (length(These) > 2) {These <- head(These, 2)} #If we wanted to institute a number of peaks argument, it would be here.

  MyData <- cbind(StashedIDs, MyRawData, MyData)
  MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

  if (length(These) == 2){second <- These[[1]]
  third <- These[[2]]
  } else if(length(These) == 1){second <- These[[1]]}

  if(exists("second")){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[second]], 0.0) ~ paste(MyData$Cluster, second, "_00-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.2) ~ paste(MyData$Cluster, second, "_02-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.4) ~ paste(MyData$Cluster, second, "_04-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.6) ~ paste(MyData$Cluster, second, "_06-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.8) ~ paste(MyData$Cluster, second, "_08-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 1.0) ~ paste(MyData$Cluster, second, "_10-",
                                        sep = "", collapse = NULL)))
  }

  if(exists("third")){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[third]], 0.0) ~ paste(MyData$Cluster, third, "_00",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.2) ~ paste(MyData$Cluster, third, "_02",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.4) ~ paste(MyData$Cluster, third, "_04",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.6) ~ paste(MyData$Cluster, third, "_06",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.8) ~ paste(MyData$Cluster, third, "_08",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 1.0) ~ paste(MyData$Cluster, third, "_10",
                                       sep = "", collapse = NULL)))
  }

  return(MyData)
}



SingleStainSignatures <- function(x, WorkAround1, alternatename, ColsN, StartNormalizedMergedCol, EndNormalizedMergedCol,
                                  Samples, name, results){

  MySubset <- WorkAround1 %>% dplyr::filter(.data[[x]] == 1.000)
  StashedIDs <- MySubset %>% select(Backups)
  MySubset <- MySubset %>% select(-Backups)
  MySubset <- MySubset %>% select(all_of(1:ColsN)) #Not Specified Internally
  BackupNames2 <- colnames(MySubset)
  DetectorName <- x

  #if (is.null(external1)) {
    Data <- MySubset
    Samples_replicated <- Samples[rep(1, each = nrow(Data)),]
    Test <- Data[, 1:ColsN] - Samples_replicated[, 1:ColsN]
    Test[Test < 0] <- 0

    AA <- do.call(pmax, Test)
    Normalized2 <- Test/AA
    Normalized2 <- round(Normalized2, 1)
    colnames(Normalized2) <- gsub("-A", "", colnames(Normalized2))

    Counts2 <- colSums(Normalized2 == 1)
    Captured <- round(Counts2[x]/sum(Counts2), 2)
    message(paste0(x, " retained ", Captured, " of the Variance"))
    WorkAround2 <- cbind(MySubset, Normalized2)

    if (any(str_detect(name, names(results)))){
      WorkAround3 <- WorkAround2 %>% mutate(Backups = StashedIDs$Backups) %>%
        relocate(Backups, .before = 1)
      WorkAroundInt <- WorkAround3 %>% dplyr::filter(.data[[x]] == 1.000)
      StashedIDs <- WorkAroundInt %>% select(Backups)
      WorkAround2 <- WorkAroundInt %>% select(-Backups)
    }
  #} else {
  #  Data <- MySubset
  #  Samples_replicated <- external1[rep(1, each = nrow(Data)),]
  #  Test <- Data[, 1:ColsN] - Samples_replicated[, 1:ColsN]
  #  Test[Test < 0] <- 0
  #
  #  AA <- do.call(pmax, Test)
  #  Normalized2 <- Test/AA
  #  Normalized2 <- round(Normalized2, 1)
  #  colnames(Normalized2) <- gsub("-A", "", colnames(Normalized2))

  #  Counts2 <- colSums(Normalized2 == 1)
  #  Captured <- round(Counts2[x]/sum(Counts2), 2)
  #  message(paste0(x, " retained ", Captured, " of the Variance"))

    #Bringing Together Raw And Subtracted Normalized
  #  WorkAround2 <- cbind(MySubset, Normalized2)

   # if (any(str_detect(name, names(results)))){
   #   WorkAround3 <- WorkAround2 %>% mutate(Backups = StashedIDs$Backups) %>%
   #     relocate(Backups, .before = 1)
   #   WorkAroundInt <- WorkAround3 %>% dplyr::filter(.data[[x]] == 1.000)
   #   StashedIDs <- WorkAroundInt %>% select(Backups)
   #   WorkAround2 <- WorkAroundInt %>% select(-Backups)
   # }
  #}

  MyData <- WorkAround2 %>% select(all_of(
    StartNormalizedMergedCol:EndNormalizedMergedCol)) %>%
    mutate(across(where(is.numeric), ~ ceiling(. / 0.2) * 0.2))
  MyRawData <- WorkAround2 %>% select(all_of(1:ColsN))

  #Preparation for Local Maxima
  Conversion <- data.frame(t(MyData), check.names = FALSE)
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
  alternatename <- alternatename
  PointData <- Luciernaga::Utility_LocalMaxima(theX = LocalX, theY = LocalY,
                                               therepeats = 3, w = 3, span = 0.11, alternatename = alternatename)

  colnames(PointData)[1] <- "TheDetector"
  colnames(PointData)[2] <- "TheHeight"

  Newest2 <- PointData %>% filter(TheHeight > 0.15) %>% arrange(desc(TheHeight))
  Assembled <- left_join(Newest2, Decoys, by = "TheDetector")
  if(nrow(Assembled) == 0){stop("Failed at Assembled, no local maxima greater than 0.15")}
  These <- Assembled %>% pull(Detectors)
  if (any(These %in% x)) {These <- These[These != x]}
  if (length(These) > 2) {These <- head(These, 2)} #If we wanted to institute a number of peaks argument, it would be here.

  MyData <- cbind(StashedIDs, MyRawData, MyData)
  MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

  if (length(These) == 2){second <- These[[1]]
  third <- These[[2]]
  } else if(length(These) == 1){second <- These[[1]]}

  if(exists("second")){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[second]], 0.0) ~ paste(MyData$Cluster, second, "_00-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.2) ~ paste(MyData$Cluster, second, "_02-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.4) ~ paste(MyData$Cluster, second, "_04-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.6) ~ paste(MyData$Cluster, second, "_06-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 0.8) ~ paste(MyData$Cluster, second, "_08-",
                                        sep = "", collapse = NULL),
    near(MyData[[second]], 1.0) ~ paste(MyData$Cluster, second, "_10-",
                                        sep = "", collapse = NULL)))
  }

  if(exists("third")){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[third]], 0.0) ~ paste(MyData$Cluster, third, "_00",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.2) ~ paste(MyData$Cluster, third, "_02",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.4) ~ paste(MyData$Cluster, third, "_04",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.6) ~ paste(MyData$Cluster, third, "_06",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 0.8) ~ paste(MyData$Cluster, third, "_08",
                                       sep = "", collapse = NULL),
    near(MyData[[third]], 1.0) ~ paste(MyData$Cluster, third, "_10",
                                       sep = "", collapse = NULL)))
  }

  return(MyData)
}
