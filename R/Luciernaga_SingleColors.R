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
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore exprs
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise_all
#' @importFrom dplyr rename
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot
#' @importFrom utils read.csv
#'
#' @return A tibble row for each flow object containing the summarized data.
#' @export
#'
#' @examples NULL

Luciernaga_SingleColors <- function(x, sample.name, removestrings, subset, PanelCuts,
                                    stats, SignatureView, Verbose=FALSE){
  name <- keyword(x, sample.name)
  name <- NameCleanUp(name, removestrings=c("(Cells)", "(Beads)"))
  name <- gsub("\\s+$", "", name)
  name <- NameCleanUp(name, removestrings=removestrings)

  if (Verbose == TRUE){message("After removestrings cleanup the name is ", name)}

  name <- sub(" ", "_", name)
  name <- sub("(.*)_(.*)", "\\2_\\1", name) #Swaps order
  TheFluorophore <- strsplit(name, "_")[[1]]
  TheFluorophores <- TheFluorophore[1]
  TheLigand <- TheFluorophore[2]

  if (Verbose == TRUE){message("The Fluorophore is ", TheFluorophores, " and the ligand is ", TheLigand)}

  cs <- gs_pop_get_data(x, subset)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)
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

  if (!(LowerBound > 0 & LowerBound < 1)) {
    message("From should be between 0 and 1, proceeding to divide by 100 on assumption it was a percentage")
    LowerBound <- LowerBound / 100
  }

  if (!(UpperBound > 0 & UpperBound < 1)) {
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
    NData[NData < 0] <- 0
    A <- do.call(pmax, NData)
    Normalized <- NData/A
    Normalized <- round(Normalized, 1)
    ThePlotData <- cbind(name, Normalized) %>% rename(Fluorophore = name)

    LineCols <- ncol(ThePlotData)

    Melted <- gather(ThePlotData, key = "Detector", value = "value", all_of(
      2:LineCols)) #Gather is my New Best Friend

    Melted$Detector <- factor(Melted$Detector, levels = detector_order)
    Melted$Fluorophore <- factor(Melted$Fluorophore)
    Melted1 <- data.frame(Melted)

    plot <- ggplot(Melted1, aes(x = Detector, y = value, group = Fluorophore,
                                color = Fluorophore)) + geom_line() +
      scale_color_hue(direction = 1) +
      labs(title = TheDetector, x = "Detectors", y = "Normalized Values") +
      theme_linedraw() + theme_bw() +
      theme(axis.title.x = element_text(face = "plain"),
            axis.title.y = element_text(face = "plain",
                                        margin = margin(r = -120)),
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
      )

    print(plot)

  }

  return(Data)
}
