#' Calculates the single color control matrix for a given quantile cutoffs and a given
#' statistic
#'
#' @param x A Gating Set Object
#' @param sample.name The Keyword for which the Fluorophore Name is stored
#' @param removestrings Values to remove from the name
#' @param subset A desired Gating Hierarchy level of cells to filter in
#' @param PanelCuts A .csv or dataframe containing columns Fluorophore, From and To
#' Fluorophore name should match sample.name style
#' @param stats Whether to use "mean" or "median"
#' @param SignatureView Whether to also return a normalized signature plot.
#' @param Verbose Provides debugging for removestrings
#' @param returntype Allows to modify default "data" to instead return the "plots"
#'
#' @importFrom flowCore keyword
#' @importFrom stringr str_detect
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom BiocGenerics nrow
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom utils read.csv
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' @importFrom stats quantile
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#'
#' @return A tibble row for each flow object containing the summarized data.
#' @export
#'
#' @examples NULL
Luciernaga_SingleColors <- function(x, sample.name, removestrings, subset, PanelCuts,
                                    stats, SignatureView, Verbose=FALSE, returntype="data"){

  if (length(sample.name) == 2){
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    second <- keyword(x, second)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, sample.name)}

  name <- NameCleanUp(name, removestrings=c("(Cells)", "(Beads)", "Reference Group_"))
  name <- gsub("\\s+$", "", name)
  name <- NameCleanUp(name, removestrings=removestrings)

  if (Verbose == TRUE){message("After removestrings cleanup the name is ", name)}

  if (!str_detect(name, "nstained")){
  name <- sub(" ", "_", name)
  name <- sub("(.*)_(.*)", "\\2_\\1", name) #Swaps order
  TheFluorophore <- strsplit(name, "_")[[1]]
  TheFluorophores <- TheFluorophore[1]
  TheLigand <- TheFluorophore[2]
  } else {
    TheFluorophores <- "Unstained"
    TheLigand <- name
  }

  if (Verbose == TRUE){
    message("The Fluorophore is ", TheFluorophores, " and the ligand is ", TheLigand)
    }

  cs <- gs_pop_get_data(x, subset)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)
  Data <- Data %>% unique() #Precaution Artificial Leftovers

  TheColumns <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  detector_order <- colnames(TheColumns)
  startingcells <- BiocGenerics::nrow(cs)[[1]]

  # Triangulating on Detector
  n <- TheColumns
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  Normalized <- round(Normalized, 1)
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts), Counts = Counts)
  rownames(PeakDetectorCounts) <- NULL
  cutoff <- startingcells*0.0075
  Detectors <- PeakDetectorCounts %>% filter(Counts > cutoff) %>%
    arrange(desc(Counts))
  TheDetector <- Detectors[1,1]

  # Bringing in the Panel and Cut Quantiles
  if (!is.data.frame(PanelCuts)) {
    PanelCuts <- read.csv(PanelCuts, check.names = FALSE)
    } else {PanelCuts <- PanelCuts}

  PanelCuts$Fluorophore <- gsub("-A", "", PanelCuts$Fluorophore)

  TheInfo <- PanelCuts %>% dplyr::filter(Fluorophore %in% TheFluorophores)

  if(base::nrow(TheInfo) > 1){stop("More than two rows retrieved when selecting Fluorophore")}

  LowerBound <- TheInfo %>% select(From) %>% pull()
  UpperBound <- TheInfo %>% select(To) %>% pull()

  if (!(LowerBound >= 0 & LowerBound <= 1)) {
    message("From should be between 0 and 1, proceeding to divide by 100 on assumption it was a percentage")
    LowerBound <- LowerBound / 100
  }

  if (!(UpperBound >= 0 & UpperBound <= 1)) {
    message("To should be between 0 and 1, proceeding to divide by 100 on assumption it was a percentage")
    UpperBound <- UpperBound / 100
  }

  QuantileData <- TheColumns %>% select(all_of(TheDetector)) %>% pull()
  LowerBoundMFI <- QuantileData %>% quantile(., LowerBound)
  UpperBoundMFI <- QuantileData %>% quantile(., UpperBound)


  ValuesInterest <- TheColumns %>% filter(
    .data[[TheDetector]]  >= LowerBoundMFI & .data[[TheDetector]] <= UpperBoundMFI)

  Samples <- Luciernaga:::AveragedSignature(x=ValuesInterest, stats=stats)


  Data <- cbind(TheFluorophores, TheLigand, Samples) %>% rename(Fluorophore = TheFluorophores) %>%
    rename(Ligand = TheLigand)

  if (SignatureView == TRUE){

    NData <- Samples
    NumberDetectors <- ncol(NData)
    ReferenceData <- ReferenceCall(NumberDetectors)
    ReferenceData <- ReferenceData %>% select(-Instrument) %>% rename(value=AdjustedY) %>%
      dplyr::filter(Fluorophore %in%TheFluorophores)

    NData[NData < 0] <- 0
    A <- do.call(pmax, NData)
    Normalized <- NData/A
    Normalized <- round(Normalized, 1)

    thename <- paste(TheFluorophores, TheLigand, sep = " ")
    ThePlotData <- Normalized %>% mutate(Fluorophore = thename) %>% relocate(Fluorophore, .before=1)
    LineCols <- ncol(ThePlotData)

    Melted <- ThePlotData %>% pivot_longer(cols=2:LineCols, names_to = "Detector", values_to = "value")
    Melted$Detector <- factor(Melted$Detector, levels = detector_order)

    Melted <- Melted %>% mutate(Type = "Provided")

    if (!nrow(ReferenceData) == 0){

    ReferenceData <- ReferenceData %>% mutate(Type = "Reference")
    NewDetectorNames <- Melted %>% pull(Detector)
    ReferenceData <- ReferenceData %>% mutate(Detector=NewDetectorNames)
    ReferenceData$value <- round(ReferenceData$value, 2)

    Melted <- rbind(Melted, ReferenceData)
    }

    plot <- ggplot(Melted, aes(x = Detector, y = value, group = Type,
            color = Type)) + geom_line(size=1) + theme_bw() +
      labs(title = thename, x = "Detectors", y = "Normalized") +
      theme(axis.title.x = element_text(face = "plain"), axis.title.y = element_text(face = "plain"),
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.position="none") + scale_color_manual(values = c("Provided" = "blue", "Reference" = "red"))
  }

  if (SignatureView == TRUE & returntype == "plots"){
    return(plot)
  } else if (SignatureView == TRUE & returntype == "data"){
    print(plot)
    return(Data)
  } else {return(Data)}

}



#' Internal to quickly bring up reference data for comparison plotting
#'
#' @param NumberDetectors Passed Number of Columns to decipher instrument
#' @importFrom utils read.csv
#'
#' @return The corresponding reference data.frame
#'
#' @noRd
ReferenceCall <- function(NumberDetectors){
  if (NumberDetectors == 64){instrument <- "FiveLaser"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary5L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 54){instrument <- "FourLaser"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary4LUV.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 38){instrument <- "ThreeLaser"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary3L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else {message("No References Found")}

  return(ReferenceData)
}

