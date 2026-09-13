#' Main Luciernaga Function, normalized on individual cell level
#'
#' @param x A Gating Set object (ex. gs or gs[[1]])
#' @param subsets The Gating Hierarchy level you will be sampling at
#' @param sample.name Keyword variable which samples are stored (ex. "GUID")
#' @param removestrings A string of character values to remove from sample.name
#' @param Verbose Whether to print outputs as you go.
#' @param unmixingcontroltype Whether your inputs are "cells", "beads" or "both"
#' @param Unstained Set to True when running unstained samples that don't have
#' Unstained in the Name
#' @param ratiopopcutoff A numeric ratio for peak detector inclusion, default is
#' set to 0.01 all startingcells
#' @param experiment Provide directly experiment name (ex. "JAN2024")
#' @param experiment.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param condition Provide directly experiment name (ex. "JAN2024")
#' @param condition.name Keyword variable which experiment information
#' is stored (ex. "TUBENAME")
#' @param AFOverlap A data.frame or a filepath to the CSV containing the
#' Autofluorescence
#' overlap of individual fluorophores for exclusion
#' @param stats Whether to take "mean" or "median"
#' @param desiredAF Main Autofluorescence Detector (ex. "V7-A")
#' @param ExportType Whether to return "fcs", "data"
#' @param SignatureReturnNow Short circuits the function and returns signature
#' for specified autofluorescence.
#' @param outpath  Location where created .fcs and .csv files are sent
#' @param minimalfcscutoff A ratio determining number cells needed for .fcs export
#' @param Subtraction Whether for single color controls to use "Internal", "Internal_General" or
#' "External" autofluorescence
#' @param SCData Whether to return "subtracted" or "raw" data for single colors
#' @param NegativeType Whether to append a negative pop to .fcs file, "default",
#' "artifical" or "samples"
#' @param Brightness Whether sum of detectors should be returned.
#' @param LocalMaximaRatio Height of peaks to proceed
#' @param SecondaryPeaks Number of Secondary Peaks, default is set to 2.
#' @param Increments Rounding parameter, default is set to 0.1
#' @param RetainedType Whether to return "raw" or "normalized" values for lineplots.
#' @param BeadAF A passed data.frame row containing the reference for bead unstained.
#' @param BeadMainAF The detector that corresponds to the "main" bead AF, albeit dim.
#' @param CellAF A passed data.frame row containing the reference for cell unstained.
#' @param CellMainAF The detector that corresponds to the "main" cell AF.
#' @param TotalNegatives When setting NegativeType to "artificial" or "sample", how many events to add.
#' @param inverse.transform Passed to flowWorkspace, if data has been transformed and wish to return to
#' raw values, set to TRUE.
#' @param Consolidate Default NULL, alternative Cluster character string to partly match to consolidate fcs_export
#'
#' @importFrom flowCore keyword exprs
#' @importFrom stringr str_detect str_split
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr mutate select arrange desc bind_rows pull filter
#'  slice relocate left_join
#' @importFrom purrr map compact set_names
#' @importFrom utils read.csv
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @return Additional information to be added
#' @export
#'
#' @examples
#'
#' library(flowCore)
#' library(flowWorkspace)
#' library(openCyto)
#' library(data.table)
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#'
#' File_Location <- system.file("extdata", package = "Luciernaga")
#' FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
#'   full.names = TRUE)
#' CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
#' CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
#' MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
#' StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
#'
#' FileLocation <- system.file("extdata", package = "Luciernaga")
#' pattern = "AutofluorescentOverlaps.csv"
#' AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)
#'
#' SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC, subsets="lymphocytes",
#'  removestrings=removestrings, sample.name="GUID", unmixingcontroltype = "cells",
#'  Unstained = FALSE, ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
#'  stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
#'  outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
#'  experiment = "FirstExperiment", condition = "ILTPanel", Subtraction = "Internal",
#'  CellAF=TheCellAF, SCData="subtracted", NegativeType="default") |> bind_rows()
#'
Luciernaga_QC <- function(x,
                           subsets,
                           sample.name = "TUBENAME",
                           removestrings = NULL,
                           Verbose = FALSE,
                           experiment = NULL,
                           experiment.name = NULL,
                           condition = NULL,
                           condition.name = NULL,
                           AFOverlap,
                           unmixingcontroltype = "both",
                           Unstained = FALSE,
                           ratiopopcutoff = 0.01,
                           stats = "median",
                           Subtraction = "Internal",
                           desiredAF = NULL,
                           BeadAF = NULL,
                           BeadMainAF = NULL,
                           CellAF = NULL,
                           CellMainAF = NULL,
                           SignatureReturnNow = FALSE,
                           Increments = 0.1,
                           LocalMaximaRatio = 0.15,
                           SecondaryPeaks = 2,
                           Brightness = FALSE,
                           RetainedType = "raw",
                           ExportType = "data",
                           minimalfcscutoff = 0.05,
                           SCData = "subtracted",
                           NegativeType = "default",
                           TotalNegatives = 500,
                           outpath = NULL,
                           inverse.transform = TRUE,
                           Consolidate = NULL) {

  ###################
  # Metadata Module #
  ###################

  if (length(sample.name) == 2) {
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep = "_")
  } else {
    name <- keyword(x, sample.name)
  }

  Type <- Luciernaga:::Typing(name = name,
                               unmixingcontroltype = unmixingcontroltype,
                               Unstained = Unstained)
  AggregateName <- Luciernaga:::NameForSample(x = x, sample.name = sample.name,
                                                removestrings = removestrings)

  if (is.null(experiment) && is.null(experiment.name) &&
        SignatureReturnNow == FALSE) {
    message("Both experiment and experiment.name are set to NULL,
            consider adding one or the other.")
  }

  Experiment <- Luciernaga:::NameForSample(x = x, sample.name = sample.name,
                                             removestrings = removestrings,
                                             experiment = experiment,
                                             experiment.name = experiment.name,
                                             returnType = "experiment")

  if (is.null(condition) && is.null(condition.name) &&
        SignatureReturnNow == FALSE) {
    message("Both condition and condition.name are set to NULL,
            consider adding one or the other.")
  }

  Condition <- Luciernaga:::NameForSample(x = x, sample.name = sample.name,
                                            removestrings = removestrings,
                                            condition = condition,
                                            condition.name = condition.name,
                                            returnType = "condition")

  InternalCleanupList <- c(".fcs", "Cells", "Beads", " ", "_", "-", ".",
                            "(", ")")
  name <- Luciernaga:::NameCleanUp(name, InternalCleanupList)

  if (Unstained == TRUE) {
    if (!str_detect(name, "stained")) {
      name <- paste0(name, "_Unstained")
    }
    if (!str_detect(Type, "stained")) {
      Type <- paste0(Type, "_Unstained")
    }
  }

  ###############
  # Exprs Setup #
  ###############

  ff <- gs_pop_get_data(x, subsets, inverse.transform = inverse.transform)
  startingcells <- nrow(ff)[[1]]
  DF <- as.data.frame(exprs(ff[[1]]), check.names = FALSE)

  Backups <- DF |>
    mutate(Backups = 1:nrow(DF)) |>
    select(Backups)

  OriginalColumns <- colnames(DF)
  OriginalColumns <- data.frame(OriginalColumns)
  OriginalColumnsIndex <- OriginalColumns %>%
    mutate(IndexLocation = 1:nrow(.)) # TODO: `.` refers to lhs (magrittr-only)

  OriginalColumnsVector <- colnames(DF)

  StashedDF <- DF[, grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]
  StashedDF <- cbind(Backups, StashedDF)

  n <- DF[, -grep("Time|FS|SC|SS|Original|W$|H$", names(DF))]

  if (Verbose == TRUE) {
    TheTotal <- nrow(n) * ncol(n)
    BelowZero <- sum(apply(n, 2, function(x) x < 0))
    message(round(BelowZero / TheTotal, 2),
            " of all events were negative and will be rounded to 0")
  }

  #################################################
  # Normalizing Individual Cells By Peak Detector #
  #################################################

  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n / A
  colnames(Normalized) <- gsub("-A", "", colnames(Normalized))

  ColsN <- ncol(n)
  ColsNormalized <- ncol(Normalized)
  StartNormalizedMergedCol <- ColsN + 1
  EndNormalizedMergedCol <- ColsN + ColsNormalized

  WorkAround <- cbind(n, Normalized)

  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  PeakDetectorCounts <- PeakDetectorCounts |> arrange(desc(Counts))

  ##############################
  # Determining Peak Detectors #
  ##############################

  if (Type == "Cells" | Type == "Cells_Unstained" |
        Type == "Unknown_Unstained") {
    CellCutoff <- startingcells * ratiopopcutoff
    Detectors <- PeakDetectorCounts |> filter(Counts > CellCutoff)
  }

  if (Type == "Beads" | Type == "Beads_Unstained") {
    filtering <- (startingcells / ncol(n)) * 1
    TheCandidates <- PeakDetectorCounts |>
      filter(Counts > filtering) |>
      pull(Fluors)
    TheMedians <- map(.x = TheCandidates, .f = BeadDetectors,
                       data = WorkAround) |>
      bind_rows()
    if (Verbose == TRUE) {
      message("Returning Peak Bead Detector Medians")
      # YNW(TheMedians)
      TheMedians
    }
    PeakCutoffVal <- mean(TheMedians$TheMedian)
    TheRemnant <- TheMedians |>
      filter(TheMedian > PeakCutoffVal) |>
      pull(TheDetector)
    Detectors <- PeakDetectorCounts |> filter(Fluors %in% TheRemnant)
  }

  if (Verbose == TRUE) {
    # YNW(Detectors)
    Detectors
  }

  if (Type == "Beads_Unstained") {
    ReferenceUnstained <- AveragedSignature(n, stats)
    BeadPlot <- ReferenceUnstained |>
      mutate(Sample = "UnstainedBeads") |>
      relocate(Sample, .before = 1)
    Plot <- QC_ViewSignature(x = "UnstainedBeads", data = BeadPlot,
                              Normalize = TRUE)

    if (Verbose == TRUE) {
      message("Returning designated stats values for Beads_Unstained, please return
    to LuciernagaQC as BeadsAF for subtraction for single color unmixing controls")

      # YNW(Plot)
      Plot
    }
    return(ReferenceUnstained)
  }

  ###########################################################
  # Handling Fluorophore/AF overlap for cell controls cells #
  ###########################################################

  if (str_detect(Type, "ells")) {

    if (is.data.frame(AFOverlap)) {
      AFData <- AFOverlap
    } else {
      AFData <- read.csv(file = AFOverlap, check.names = FALSE)
    }

    AFChannels <- AFData |>
      filter(Fluorophore %in% "Unstained") |>
      pull(MainDetector) |>
      str_split(",", simplify = TRUE)
    AFChannels <- AFChannels[1, ]
    AFChannels <- gsub("-A", "", AFChannels)

    TheSCData <- AFData |> filter(Fluorophore != "Unstained")
    TheSCData$Fluorophore <- gsub("-A", "", TheSCData$Fluorophore)
    TroubleChannels <- TheSCData |> pull(Fluorophore)

    results <- map(.x = TroubleChannels,
                   .f = Luciernaga:::TroubleChannelExclusion,
                   TheSCData = TheSCData, MainDetector = MainDetector,
                   AFChannels = AFChannels) |>
      set_names(TroubleChannels)

    OverlapFlag <- NULL
    matching_names <- names(results)[str_detect(name, names(results))]
    if (length(matching_names) > 0) {
      OverlapFlag <- "Yep"
      ExclusionList <- results[[matching_names[1]]]
      Retained <- Detectors |>
        filter(!Fluors %in% ExclusionList) |>
        pull(Fluors)
    } else if (str_detect(name, "nstained")) {
      Retained <- Detectors |> pull(Fluors)
    } else {
      Retained <- Detectors |>
        filter(!Fluors %in% AFChannels) |>
        pull(Fluors)
    }

  } else {
    Retained <- Detectors |> pull(Fluors)
  }

  if (length(Retained) == 0) {
    stop("There were no Retained detectors in ", name)
  }

  if (Verbose == TRUE) {
    # YNW(Retained)
    Retained
  }

  ##############################################
  # Determining Main Autofluorescence Detector #
  ##############################################

  # Inverse of retained, then top for peak

  if (Type == "Cells") {
    Intermediate <- Detectors |> filter(!Fluors %in% Retained)
  }
  if (Type == "Cells_Unstained" | Type == "Unknown_Unstained") {
    Intermediate <- Detectors
  }

  if (Type == "Cells" | Type == "Cells_Unstained" |
        Type == "Unknown_Unstained") {
    InternalOverride <- FALSE

    if (nrow(Intermediate) > 0) {
      TheMainAF <- Intermediate |> slice(1) |> pull(Fluors)
    } else {

      if (Type == "Cells") {
        message("Only a single peak detector present.")
      }

      if (Type == "Cells" && Subtraction == "Internal" &&
            !is.null(CellAF)) {
        InternalOverride <- TRUE
        OverlapFlag <- "Yep"
      }

      if (Type == "Cells" && Subtraction == "Internal" &&
            is.null(CellAF)) {
        stop("Only one detector present and no external cell autofluorescence signature was provided
           for subtraction. Please provision the CellAF argument.")
      }
    }

    if (Subtraction == "Internal" | Subtraction == "Internal_General") {
      if (InternalOverride != TRUE) {

        if (!is.null(desiredAF)) {
          TheMainAF <- NameCleanUp(desiredAF, c("-A"))
        }

        This <- WorkAround |>
          filter(.data[[TheMainAF]] == 1) |>
          select(all_of(1:ColsN))

        Samples <- AveragedSignature(x = This, stats = stats)

      } else {
        Subtraction <- "External"
      }
    }

    if (Subtraction == "External") {

      if (is.data.frame(CellAF)) {
        Samples <- CellAF
        MaxVal <- do.call(pmax, Samples)
        TheNormed <- Samples / MaxVal
        TheCounts <- colSums(TheNormed == 1)
        ThePeakDetectorCounts <- data.frame(Fluors = names(TheCounts),
                                             Counts = TheCounts)
        rownames(ThePeakDetectorCounts) <- NULL
        TheMainAF <- ThePeakDetectorCounts |>
          arrange(desc(Counts)) |>
          slice(1) |>
          pull(Fluors)
        TheMainAF <- gsub("-A", "", TheMainAF)

        if (nrow(Samples) > 1) {
          Samples <- AveragedSignature(x = CellAF, stats = stats)
        }
      } else {
        stop("CellAF needs to be a data.frame object with a single row. Use SignatureReturnNow = TRUE
                 on a Luciernaga_QC() Unstained Sample, or AveragedSignature() on exprs data to generate
                 the correct format.")
      }
    }

  }

  if (Type == "Beads") {

    if (SignatureReturnNow == TRUE) {

      if (nrow(Detectors) > 1) {
        TheMainAF <- Intermediate |> slice(1) |> pull(Fluors)
      } else {
        TheMainAF <- Detectors |> pull(Fluors)
      }

      if (!is.null(desiredAF)) {
        TheMainAF <- NameCleanUp(desiredAF, c("-A"))
      }

      This <- WorkAround |>
        filter(.data[[TheMainAF]] == 1) |>
        select(all_of(1:ColsN))

      Samples <- AveragedSignature(x = This, stats = stats)

    }
  }

  if (SignatureReturnNow == TRUE) {
    ThePlot <- Samples |>
      mutate(Sample = name) |>
      relocate(Sample, .before = 1)
    Plot <- QC_ViewSignature(x = name, data = ThePlot, Normalize = TRUE)

    if (Verbose == TRUE) {
      # YNW(Detectors)
      # YNW(Plot)
      Detectors
      Plot
    }

    return(Samples)
  }

  ##################################################
  # Sending off for unmixing control type handling #
  ##################################################

  WorkAround1 <- WorkAround |>
    mutate(Backups = Backups$Backups) |>
    relocate(Backups, .before = 1) # This will change the start/end count

  if (Type == "Cells_Unstained" | Type == "Unknown_Unstained") {
    # x <- Retained[3]
    RetainedDF <- map(.x = Retained, .f = UnstainedSignatures,
                       WorkAround1 = WorkAround1,
                       alternatename = AggregateName,
                       ColsN = ColsN,
                       StartNormalizedMergedCol = StartNormalizedMergedCol,
                       EndNormalizedMergedCol = EndNormalizedMergedCol,
                       Verbose = Verbose, Increments = Increments,
                       SecondaryPeaks = SecondaryPeaks) |>
      bind_rows()

  }

  if (Type == "Cells") {
    if (Subtraction == "Internal_General" | Subtraction == "External") {
      Subtraction <- "Average"
    }

    if (!is.null(OverlapFlag)) {
      Subtraction <- "Average"
      if (Verbose == TRUE) {
        message(name, " used Average Subtraction.")
      }
    }

    # x <- Retained[4]
    RetainedDF <- map(.x = Retained, .f = SingleStainSignatures,
                       WorkAround1 = WorkAround1,
                       AggregateName = AggregateName,
                       ColsN = ColsN,
                       StartNormalizedMergedCol = StartNormalizedMergedCol,
                       EndNormalizedMergedCol = EndNormalizedMergedCol,
                       Samples = Samples, Increments = Increments,
                       Subtraction = Subtraction, stats = stats,
                       TheMainAF = TheMainAF, Verbose = Verbose,
                       SCData = SCData, LocalMaximaRatio = LocalMaximaRatio,
                       SecondaryPeaks = SecondaryPeaks)

    cleaned_results <- compact(RetainedDF)
    RetainedDF <- bind_rows(cleaned_results)

  }

  if (Type == "Beads") {

    if (!is.null(BeadAF)) {
      if (is.data.frame(BeadAF)) {
        Samples <- BeadAF
        if (nrow(Samples) > 1) {
          Samples <- AveragedSignature(x = CellAF, stats = stats)
        }
      } else {
        stop("Please provide a bead unstained reference signature in  a dataframe
      format to the BeadAF argument.")
      }
    } else {
      stop("Please provide a bead unstained reference signature in  a dataframe
      format to the BeadAF argument.")
    }

    if (!is.null(BeadMainAF)) {
      TheMainAF <- gsub("-A", "", BeadMainAF)
    } else {
      stop("Please provide a placeholder detector to the BeadMainAF argument")
    }

    RetainedDF <- map(.x = Retained, .f = SingleStainSignatures,
                       WorkAround1 = WorkAround1,
                       AggregateName = AggregateName,
                       ColsN = ColsN,
                       StartNormalizedMergedCol = StartNormalizedMergedCol,
                       EndNormalizedMergedCol = EndNormalizedMergedCol,
                       Samples = Samples, Increments = Increments,
                       Subtraction = "Average", stats = stats,
                       TheMainAF = TheMainAF, Verbose = Verbose,
                       SCData = SCData, LocalMaximaRatio = LocalMaximaRatio,
                       SecondaryPeaks = SecondaryPeaks) |>
      bind_rows()
  }

  ###################################################
  # Sending the returns to their finals destination #
  ###################################################

  Reintegrated <- left_join(RetainedDF, StashedDF, by = "Backups")

  BackupsCol <- "Backups"
  NormalizedColumns <- colnames(Normalized)
  ClusterCol <- "Cluster"
  RearrangedColumns <- c(BackupsCol, OriginalColumnsVector, NormalizedColumns,
                          ClusterCol)

  OriginalStart <- length(BackupsCol) + 1
  OriginalEnd <- length(BackupsCol) + length(OriginalColumnsVector)

  Reintegrated1 <- Reintegrated |> relocate(all_of(RearrangedColumns))
  Reintegrated1 <- Reintegrated1 |> arrange(Backups)

  if (ExportType == "fcs") {

    BrightnessReturn <- Genesis(x = Reintegrated1, ff = ff,
                                 minimalfcscutoff = minimalfcscutoff,
                                 AggregateName = AggregateName,
                                 Brightness = Brightness, outpath = outpath,
                                 OriginalStart = OriginalStart,
                                 OriginalEnd = OriginalEnd, stats = stats,
                                 NegativeType = NegativeType,
                                 TotalNegatives = TotalNegatives,
                                 Samples = Samples, ExportType = ExportType,
                                 Consolidate = Consolidate)
    ExportType <- "data"
  }

  if (ExportType == "data.frame") {
    FinalDataFrame <- Reintegrated |> select(-Backups)
    return(FinalDataFrame)
  }

  if (ExportType == "data") {
    ExportData <- RetainedDF |> select(-Backups)
    TheData <- data.frame(table(ExportData$Cluster))
    TheData <- TheData |> arrange(desc(Freq))
    colnames(TheData)[1] <- "Cluster"
    colnames(TheData)[2] <- "Count"
    # Data

    TheExperiment <- as.character(Experiment)
    TheCondition <- as.character(Condition)

    TheData <- TheData |> mutate(Sample = AggregateName)
    TheData <- TheData |> mutate(Experiment = TheExperiment)
    TheData <- TheData |> mutate(Condition = TheCondition)
    TheData <- TheData |>
      relocate(Sample, Experiment, Condition, .before = Cluster)

    TheClusters <- TheData |> pull(Cluster)

    TheSummary <- map(.x = TheClusters, .f = LuciernagaSmallReport,
                       Data = ExportData, RetainedType = RetainedType,
                       ColsN = ColsN,
                       StartNormalizedMergedCol = StartNormalizedMergedCol,
                       EndNormalizedMergedCol = EndNormalizedMergedCol,
                       stats = stats) |>
      bind_rows()

    FinalData <- left_join(TheData, TheSummary, by = "Cluster")
    return(FinalData)
  }

}