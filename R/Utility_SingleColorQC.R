Utility_SingleColorQC <- function(x, subsets, sample.name, group.name, experiment = NULL, experiment.name = NULL, stats, Kept, external, sourcelocation, outpath, artificial, fcsexport, mainAF, AFOverlap, Beads, Brightness, Unstained){

  gc()
  x <- x

  # Retrieving .fcs file name
  name <- keyword(x, sample.name)
  name <- gsub(".fcs", "", name)

  if (Unstained == TRUE) {name <- paste0(name, " Unstained")}

  #Retrieving group name #InfantID in this case
  group <- keyword(x, group.name)

  # Guessing Type
  if(str_detect(name, "(Cells)")){Type <- "Cells"} else if(str_detect(name, "(Beads)")){Type <- "Beads"} else {Type <- "NULL"}

  #Additional Name Cleanup #Using Cytek Reference Control Output Names
  name <- gsub(" (Cells)", "", fixed = TRUE, gsub(" (Beads)", "", fixed = TRUE,name))

  #Setting Up an Alternate Name #Removing All Separators
  alternate.name <- name #Before Removing Spaces et al
  alternate.name <- gsub(" ", "", gsub("(", "", fixed = TRUE, gsub(")", "", fixed = TRUE, gsub("_", "", fixed = TRUE, gsub("-", "", fixed = TRUE, alternate.name)))))

  #Retrieving Experiment Info #Switched to an exist statement.
  if(!is.null(experiment)){Experiment <- experiment
  } else {experiment <- keyword(x, experiment.name)
  experiment <- gsub("-", "", fixed = TRUE, experiment) #Some Cleanup
  Experiment <- experiment}
  #suppressWarnings(rm(experiment)) #Being Used Somewhere

  #Retrieving the exprs data for the target population
  ff <- gs_pop_get_data(x, subsets)
  startingcells <- nrow(ff)[[1]]
  df <- flowCore::exprs(ff[[1]])
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

  #Retrieving Main Auto fluorescent Channels signature
  MainAF <- mainAF
  MainAF <- gsub("-A", "", MainAF)

  This <- WorkAround %>% dplyr::filter(.data[[MainAF]] == 1) %>% select(all_of(1:ColsN))

  if(stats == "mean"){Samples <- This %>% summarize_all(mean) #%>% select(-Backups)
  } else if (stats == "median"){Samples <- This %>% summarize_all(median) #%>% select(-Backups)
  } else(print("NA"))

  #Deriving Peak Detector Counts and Detectors of Interest
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*0.0075
  Detectors <- PeakDetectorCounts %>% filter(Counts > cutoff) %>% arrange(desc(Counts))

  ###################################################################
  # This is the very much cell and machine specific filter criteria #
  ###################################################################
  AFData <- AFOverlap
  AFChannels <- AFData %>% filter(Fluorophore %in% "Unstained") %>% pull(MainDetector) %>% str_split(",", simplify = TRUE)
  AFChannels <- AFChannels[1,]
  AFChannels <- gsub("-A", "", AFChannels)

  SCData <- AFData %>% filter(Fluorophore != "Unstained")
  SCData$Fluorophore <- gsub("-A", "", SCData$Fluorophore)
  TroubleChannels <- SCData %>% pull(Fluorophore)

  results <- list()
  for(w in TroubleChannels){Internal <- SCData %>% filter(Fluorophore %in% w) %>% pull(MainDetector)
  Internal <- gsub("-A", "", Internal)
  Exclusion <- setdiff(AFChannels, Internal)
  results[[w]] <- Exclusion
  }

  matching_names <- names(results)[str_detect(name, names(results))]
  if (length(matching_names) > 0) {ExclusionList <- results[[matching_names[1]]]
  Retained <- Detectors %>% filter(!Fluors %in% ExclusionList) %>% pull(Fluors)
  } else if (str_detect(name, "Unstained")){Retained <- Detectors %>% pull(Fluors)
  } else {Retained <- Detectors %>% filter(!Fluors %in% AFChannels) %>% pull(Fluors)}

  #Retained


  if (length(Retained) == 0) {
    stop("There were no Retained detectors in ", alternate.name)
  }

  if(Beads == TRUE){Retained <- Retained[[1]]}

  #####################################
  # We resume our regular programming #
  #####################################

  #Processing Unstaineds, no removal required

  if (str_detect(name, "Unstained")){
    #Retained
    RetainedDF <- data.frame()
    #i <- Retained[[1]]
    #i
    for(i in Retained){
      # Filter Normalized Data and Bin It
      WorkAround1 <- WorkAround %>% mutate(Backups = Backups$Backups) %>% relocate(Backups, .before = 1)
      MySubset <- WorkAround1 %>% dplyr::filter(.data[[i]] == 1.000)
      StashedIDs <- MySubset %>% select(Backups)
      MySubset <- MySubset %>% select(-Backups)
      DetectorName <- i
      MyData <- MySubset %>% select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol)) %>% mutate(across(where(is.numeric), ~ ceiling(. / 0.2) * 0.2))
      MyRawData <- MySubset %>% select(all_of(1:ColsN))

      #Preparation for Local Maxima
      Conversion <- data.frame(t(MyData), check.names = FALSE)
      Conversion <- cbind(Detectors = rownames(Conversion), Conversion)
      rownames(Conversion) <- NULL

      #Preparing Detector Stand Ins for left_join
      Decoys <- Conversion %>% select(Detectors)
      Decoys <- Decoys %>% mutate(TheDetector = 1:nrow(Decoys)) %>% relocate(TheDetector, .before = Detectors)

      #Deriving an average y-vector for local maxima
      Conversion <- Conversion %>% mutate(TheSums = rowSums(.[2:ncol(.)], na.rm = TRUE) /(ncol(Conversion) - 1)) %>% relocate(TheSums, .after = Detectors)
      Conversion$Detectors <- 1:nrow(Conversion)
      LocalX <- Conversion$Detectors
      LocalY <- Conversion$TheSums

      #Local Maxima
      Newest <- data.frame()

      PointData <- Utility_LocalMaxima(theX = LocalX, theY = LocalY, therepeats = 3, w = 3, span = 0.11, alternatename = alternate.name)

      for(l in PointData$x){Infinite <- PointData %>% filter(x %in% l)
      TheDetector <- Infinite %>% pull(x)
      TheHeight <- Infinite %>% pull(yhat)
      Measurement <- cbind(TheDetector, TheHeight)
      Newest <- rbind(Newest, Measurement)}
      Newest2 <- Newest %>% filter(TheHeight > 0.15) %>% arrange(desc(TheHeight))

      #Detector Cascade
      Assembled <- left_join(Newest2, Decoys, by = "TheDetector")

      #suppressWarnings(rm(These))
      These <- Assembled %>% pull(Detectors)

      if (any(These %in% i)) {These <- These[These != i]}
      if (length(These) > 2) {These <- head(These, 2)}

      #suppressWarnings(rm(second))

      #suppressWarnings(rm(third))

      ####################
      # Processing Start #
      ####################

      MyData <- cbind(StashedIDs, MyRawData, MyData)

      MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

      tryCatch({second <- These[[1]]}, error = function(e) {cat("Error occurred:", conditionMessage(e), "\n")})

      if(exists("second")){MyData <- MyData %>% mutate(Cluster = case_when(
        near(MyData[[second]], 0.0) ~ paste(MyData$Cluster, second, "_00-", sep = "", collapse = NULL),
        near(MyData[[second]], 0.2) ~ paste(MyData$Cluster, second, "_02-", sep = "", collapse = NULL),
        near(MyData[[second]], 0.4) ~ paste(MyData$Cluster, second, "_04-", sep = "", collapse = NULL),
        near(MyData[[second]], 0.6) ~ paste(MyData$Cluster, second, "_06-", sep = "", collapse = NULL),
        near(MyData[[second]], 0.8) ~ paste(MyData$Cluster, second, "_08-", sep = "", collapse = NULL),
        near(MyData[[second]], 1.0) ~ paste(MyData$Cluster, second, "_10-", sep = "", collapse = NULL)))
      }

      tryCatch({third <- These[[2]]}, error = function(e) {cat("Error occurred:", conditionMessage(e), "\n")})

      if(exists("third")){MyData <- MyData %>% mutate(Cluster = case_when(
        near(MyData[[third]], 0.0) ~ paste(MyData$Cluster, third, "_00", sep = "", collapse = NULL),
        near(MyData[[third]], 0.2) ~ paste(MyData$Cluster, third, "_02", sep = "", collapse = NULL),
        near(MyData[[third]], 0.4) ~ paste(MyData$Cluster, third, "_04", sep = "", collapse = NULL),
        near(MyData[[third]], 0.6) ~ paste(MyData$Cluster, third, "_06", sep = "", collapse = NULL),
        near(MyData[[third]], 0.8) ~ paste(MyData$Cluster, third, "_08", sep = "", collapse = NULL),
        near(MyData[[third]], 1.0) ~ paste(MyData$Cluster, third, "_10", sep = "", collapse = NULL)))
      }

      RetainedDF <- rbind(RetainedDF, MyData)

      suppressWarnings(rm(These))
      suppressWarnings(rm(second))
      suppressWarnings(rm(third))
    }

  } else {
    RetainedDF <- data.frame()
    #Retained
    #i <- Retained[[4]]
    #i
    for(i in Retained){
      # Filter Normalized Data and Bin It
      WorkAround1 <- WorkAround %>% mutate(Backups = Backups$Backups) %>% relocate(Backups, .before = 1)
      MySubset <- WorkAround1 %>% dplyr::filter(.data[[i]] == 1.000)
      StashedIDs <- MySubset %>% select(Backups)
      MySubset <- MySubset %>% select(-Backups)
      MySubset <- MySubset %>% select(all_of(1:ColsN))
      BackupNames2 <- colnames(MySubset)
      DetectorName <- i

      if (is.null(external)) {
        Data <- MySubset
        Samples_replicated <- Samples[rep(1, each = nrow(Data)),]
        Test <- Data[, 1:ColsN] - Samples_replicated[, 1:ColsN]
        Test[Test < 0] <- 0

        AA <- do.call(pmax, Test)
        Normalized2 <- Test/AA
        Normalized2 <- round(Normalized2, 1)
        colnames(Normalized2) <- gsub("-A", "", colnames(Normalized2))

        Counts2 <- colSums(Normalized2 == 1)
        Captured <- round(Counts2[i]/sum(Counts2), 2)
        print(paste0(i, " retained ", Captured, " of the Variance"))
        WorkAround2 <- cbind(MySubset, Normalized2)

        if (any(str_detect(name, names(results)))){
          WorkAround3 <- WorkAround2 %>% mutate(Backups = StashedIDs$Backups) %>% relocate(Backups, .before = 1)
          WorkAroundInt <- WorkAround3 %>% dplyr::filter(.data[[i]] == 1.000)
          StashedIDs <- WorkAroundInt %>% select(Backups)
          WorkAround2 <- WorkAroundInt %>% select(-Backups)
        }
        #View(WorkAround2)
      } else {
        Data <- MySubset
        Samples_replicated <- external[rep(1, each = nrow(Data)),]
        Test <- Data[, 1:ColsN] - Samples_replicated[, 1:ColsN]
        Test[Test < 0] <- 0

        AA <- do.call(pmax, Test)
        Normalized2 <- Test/AA
        Normalized2 <- round(Normalized2, 1)
        colnames(Normalized2) <- gsub("-A", "", colnames(Normalized2))

        Counts2 <- colSums(Normalized2 == 1)
        Captured <- round(Counts2[i]/sum(Counts2), 2)
        print(paste0(i, " retained ", Captured, " of the Variance"))

        #Bringing Together Raw And Subtracted Normalized
        WorkAround2 <- cbind(MySubset, Normalized2)

        if (any(str_detect(name, names(results)))){
          WorkAround3 <- WorkAround2 %>% mutate(Backups = StashedIDs$Backups) %>% relocate(Backups, .before = 1)
          WorkAroundInt <- WorkAround3 %>% dplyr::filter(.data[[i]] == 1.000)
          StashedIDs <- WorkAroundInt %>% select(Backups)
          WorkAround2 <- WorkAroundInt %>% select(-Backups)
        }

        #View(WorkAround2)
      }


      MyData <- WorkAround2 %>% select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol)) %>% mutate(across(where(is.numeric), ~ ceiling(. / 0.2) * 0.2))
      MyRawData <- WorkAround2 %>% select(all_of(1:ColsN))

      #Preparation for Local Maxima
      Conversion <- data.frame(t(MyData), check.names = FALSE)
      Conversion <- cbind(Detectors = rownames(Conversion), Conversion)
      rownames(Conversion) <- NULL

      #Preparing Detector Stand Ins for left_join
      Decoys <- Conversion %>% select(Detectors)
      Decoys <- Decoys %>% mutate(TheDetector = 1:nrow(Decoys)) %>% relocate(TheDetector, .before = Detectors)

      #Deriving an average y-vector for local maxima
      Conversion <- Conversion %>% mutate(TheSums = rowSums(.[2:ncol(.)], na.rm = TRUE) /(ncol(Conversion) - 1)) %>% relocate(TheSums, .after = Detectors)
      Conversion$Detectors <- 1:nrow(Conversion)
      LocalX <- Conversion$Detectors
      LocalY <- Conversion$TheSums

      #Local Maxima
      Newest <- data.frame()

      PointData <- Utility_LocalMaxima(theX = LocalX, theY = LocalY, therepeats = 3, w = 3, span = 0.11, alternatename = alternate.name)

      for(l in PointData$x){Infinite <- PointData %>% filter(x %in% l)
      TheDetector <- Infinite %>% pull(x)
      TheHeight <- Infinite %>% pull(yhat)
      Measurement <- cbind(TheDetector, TheHeight)
      Newest <- rbind(Newest, Measurement)}
      Newest2 <- Newest %>% filter(TheHeight > 0.15) %>% arrange(desc(TheHeight))

      #Detector Cascade
      Assembled <- left_join(Newest2, Decoys, by = "TheDetector")

      tryCatch({These <- Assembled %>% pull(Detectors)})

      if (length(These) == 1){
        if(any(These %in% i)){These <- These
        #second <- NULL
        #third <- NULL
        }
      } else if (length(These) == 0){print("No These Present")
        rm(These)
      } else {
        if (any(These %in% i)){These <- These[These != i]}
        if (length(These) > 2) {These <- head(These, 2)}
      }

      ####################
      # Processing Start #
      ####################

      MyData <- cbind(StashedIDs, MyRawData, MyData)

      MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

      tryCatch({seconded <- These[[1]]}, error = function(e) {cat("Error occurred:", conditionMessage(e), "\n")})

      if(exists("seconded")){MyData <- MyData %>% mutate(Cluster = case_when(
        near(MyData[[seconded]], 0.0) ~ paste(MyData$Cluster, seconded, "_00-", sep = "", collapse = NULL),
        near(MyData[[seconded]], 0.2) ~ paste(MyData$Cluster, seconded, "_02-", sep = "", collapse = NULL),
        near(MyData[[seconded]], 0.4) ~ paste(MyData$Cluster, seconded, "_04-", sep = "", collapse = NULL),
        near(MyData[[seconded]], 0.6) ~ paste(MyData$Cluster, seconded, "_06-", sep = "", collapse = NULL),
        near(MyData[[seconded]], 0.8) ~ paste(MyData$Cluster, seconded, "_08-", sep = "", collapse = NULL),
        near(MyData[[seconded]], 1.0) ~ paste(MyData$Cluster, seconded, "_10-", sep = "", collapse = NULL)))
      }

      tryCatch({thirded <- These[[2]]}, error = function(e) {cat("Error occurred:", conditionMessage(e), "\n")})

      if(exists("thirded")){MyData <- MyData %>% mutate(Cluster = case_when(
        near(MyData[[thirded]], 0.0) ~ paste(MyData$Cluster, thirded, "_00", sep = "", collapse = NULL),
        near(MyData[[thirded]], 0.2) ~ paste(MyData$Cluster, thirded, "_02", sep = "", collapse = NULL),
        near(MyData[[thirded]], 0.4) ~ paste(MyData$Cluster, thirded, "_04", sep = "", collapse = NULL),
        near(MyData[[thirded]], 0.6) ~ paste(MyData$Cluster, thirded, "_06", sep = "", collapse = NULL),
        near(MyData[[thirded]], 0.8) ~ paste(MyData$Cluster, thirded, "_08", sep = "", collapse = NULL),
        near(MyData[[thirded]], 1.0) ~ paste(MyData$Cluster, thirded, "_10", sep = "", collapse = NULL)))
      }

      RetainedDF <- rbind(RetainedDF, MyData)
      suppressWarnings(rm(These))
      suppressWarnings(rm(seconded))
      suppressWarnings(rm(thirded))
    }

  }

  #PerCPFragments <- data.frame(table(RetainedDF$Cluster))
  #PerCPFragments %>% arrange(desc(Freq))

  #View(RetainedDF)
  #View(OriginalColumnsIndex)
  Reintegrated <- left_join(RetainedDF, StashedDF, by = "Backups")

  BackupsCol <- "Backups"
  #OriginalColumnsVector
  NormalizedColumns <- colnames(Normalized)
  ClusterCol <- "Cluster"
  RearrangedColumns <- c(BackupsCol, OriginalColumnsVector, NormalizedColumns, ClusterCol)

  OriginalStart <- length(BackupsCol) + 1
  OriginalEnd <- length(BackupsCol) + length(OriginalColumnsVector)

  Reintegrated1 <- Reintegrated %>% relocate(all_of(RearrangedColumns))

  if (fcsexport == TRUE){source(sourcelocation, local = TRUE)}

  #Deriving Cluster Count Table
  My.Data2 <- RetainedDF %>% select(-Backups)
  My.Data2$Cluster <- factor(My.Data2$Cluster)
  AAA <- data.frame(table(My.Data2$Cluster))
  AAA <- AAA %>% dplyr::arrange(desc(Freq))
  colnames(AAA)[1] <- "Cluster"
  colnames(AAA)[2] <- "Count"
  Final <- AAA
  #View(Final)

  #Adding Metadata
  Final$sample <- name
  Final$group <- group
  Final$type <- Type
  Final$experiment <- Experiment
  Final <- Final %>% relocate(sample, group, type, experiment, .before = Cluster)
  #View(Final)

  #Extracting LinePlot Data
  TheOutputs <- Final %>% pull(Cluster)
  TheOutputsDF <- data.frame()
  #k <- "V5_10-UV7_00-"
  for(k in TheOutputs){

    if(Kept == "Raw"){RetainedSubset <- RetainedDF %>% filter(Cluster %in% k) %>% select(-Backups) %>% select(all_of(1:ColsN))
    } else if(Kept == "Normalized"){RetainedSubset <- RetainedDF %>% filter(Cluster %in% k) %>% select(-Backups) %>% select(all_of(StartNormalizedMergedCol:EndNormalizedMergedCol))}

    if(stats == "mean"){RetainedSamples <- RetainedSubset %>% summarize_all(mean)
    } else if (stats == "median"){RetainedSamples <- RetainedSubset %>% summarize_all(median)
    } else(print("NA"))

    Cluster <- k
    RetainedSamples <- cbind(k, RetainedSamples) %>% rename(Cluster = k)
    TheOutputsDF <- rbind(TheOutputsDF, RetainedSamples)
  }

  #Left Joining the Count with the LinePlot Data
  TheOutputsDF$Cluster <- factor(TheOutputsDF$Cluster)
  Final2 <- left_join(Final, TheOutputsDF, by = "Cluster")
  #View(Final2)
  return(Final2)
}
