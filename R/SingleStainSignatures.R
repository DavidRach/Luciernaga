#' Internal to Utility_SingleColorQC
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate
#' @importFrom tidyselect where
#' @importFrom dplyr relocate
#'
#'
#' @importFrom dplyr filter
#'
#' @importFrom dplyr summarize_all
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#'
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
#' @noRd

SingleStainSignatures <- function(x, WorkAround1, alternatename, ColsN, StartNormalizedMergedCol,
                                  EndNormalizedMergedCol, Samples, name, results, Increments=0.1,
                                  Subtraction = "Internal", stats, ratioSCcutoff=0.01, TheMainAF,
                                  AggregateName, Verbose = FALSE, SCData = "subtracted"){

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
                 ColsN=ColsN, AggregateName=AggregateName, Verbose=Verbose) %>% bind_rows()

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



