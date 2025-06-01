#' Takes a Lucierna SingleColorQC .fcs, splits into percentiles, and plots the results.
#'
#' @param x A GatingSet object
#' @param subset The desired gating hierarchy level to look at the data
#' @param sample.name Keyword for which sample name is stored
#' @param removestrings A list of values to remove from the name
#' @param stats Whether to use "mean" or "median"
#' @param returntype Whether to return a "raw" or "normalized" value lineplot.
#' @param probsratio Ratio increments to break quantiles into, default is set to 0.1.
#' @param output Whether to return "plot" or "data"
#' @param desiredAF Peak detector(ex. "V7-A") want to filter cells by before slicing, argument
#' only used to override the main peak detector when a .fcs file has more than a single peak
#' detector, default is set to NULL
#' @param legend Returns the legend, default is TRUE. 
#' @param droplowest Removes lowest bin (percentile 0), default is TRUE
#' @param titlename Default NULL, otherwise provide an alternate title. 
#' @param returncutplot Default FALSE, if true returns a histogram plot with locations
#'  where percentile slice occured.
#' @param titleplot Default NULL, sets the returncutplot title
#'
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore exprs
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom stats quantile
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom tidyr nest
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @importFrom dplyr ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_color_hue
#' @importFrom ggplot2 theme_linedraw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 annotation_logticks
#' @importFrom scales trans_format 
#' @importFrom scales math_format 
#' 
#' @return Either ggplots or the summarized data.frame object preceding
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
#'
#' APC_Example <- CellSingleColors[grep("CD16_", CellSingleColors)]
#'
#' MyCytoSet <- load_cytoset_from_fcs(APC_Example,
#'   truncate_max_range = FALSE,transformation = FALSE)
#' MyGatingSet <- GatingSet(MyCytoSet)
#' MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
#' MyGatingTemplate <- gatingTemplate(MyGates)
#' gt_gating(MyGatingTemplate, MyGatingSet)
#' removestrings <- c(".fcs")
#'
#' NormalizedSlices <- Luciernaga_LinearSlices(x=MyGatingSet[1],
#'  subset="lymphocytes",sample.name="GUID", removestrings=removestrings,
#'  stats="median", returntype="normalized", probsratio=0.1, output="plot",
#'  desiredAF="R1-A")
#'
Luciernaga_LinearSlices <- function(x, subset, sample.name, removestrings, stats,
                                    returntype, probsratio=0.1, output, desiredAF=NULL,
                                    legend=TRUE, droplowest=TRUE, titlename=NULL,
                                    returncutplot=FALSE, titleplot =NULL){
  name <- keyword(x, sample.name)
  name <- NameCleanUp(name, removestrings)

  cs <- gs_pop_get_data(x, subset)
  startingcells <- nrow(cs)[[1]]
  Data <- data.frame(exprs(cs[[1]]), check.names = FALSE)
  n <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  DetectorOrder <- colnames(n)

  # Triangulating on Detector
  n[n < 0] <- 0
  A <- do.call(pmax, n)
  Normalized <- n/A
  #Normalized <- round(Normalized, 1)
  na_counts <- colSums(is.na(Normalized))
  Normalized[is.na(Normalized)] <- 0
  Counts <- colSums(Normalized == 1)
  PeakDetectorCounts <- data.frame(Fluors = names(Counts),
                                   Counts = Counts) |> arrange(desc(Counts))
  rownames(PeakDetectorCounts) <- NULL
  Detectors <- PeakDetectorCounts |> filter(Counts > 0)

  if (nrow(Detectors) == 0){stop("No Detectors Filtered!")}

  if (nrow(Detectors) > 1){
    MultiDetector <- TRUE
    colnames(Normalized) <- gsub("-A", "", colnames(Normalized))
    data <- cbind(n, Normalized)
  } else {MultiDetector <- FALSE
          data <- n}

  if(nrow(Detectors) > 1 && is.null(desiredAF)){
    message("Luciernaga_LinearSlices is only meant to work on LuciernagaQC output .fcs files,
            your file ", name," contained two peak detectors, only the first was selected")
            Detectors <- Detectors |> slice(1)
  } else if (nrow(Detectors) > 1 && !is.null(desiredAF)){
    #desiredAF <- as.character(desiredAF)
    #Detectors$Fluors <- as.character(Detectors$Fluors)
    Detectors <- Detectors |> dplyr::filter(Fluors %in% desiredAF)
  }

  if (nrow(Detectors) == 0){stop("No Detectors Filtered at point 2! for sample ", name)}

  TheDetector <- Detectors |> pull(Fluors)

  if (MultiDetector == TRUE){
    NormDetector <- gsub("-A", "", TheDetector)
    data <- data %>% dplyr::filter(.data[[NormDetector]] == 1) # Bug is here?
    if (nrow(data) == 0){stop("We lost data at the filter step")}
    Normalized <- data |> select(!matches("-A"))
    colnames(Normalized) <- paste0(colnames(Normalized), "-A")
    data <- data |> select(matches("-A")) # Adding to Restore
  }

  #Assigning by Percentiles
  #probsratio <- 0.1
  probslabel <- probsratio*100

  percentiles <- quantile(data[[TheDetector]], probs = seq(0, 1, by = probsratio),
                          na.rm = TRUE)

  ByPercentiles <- data %>% mutate(Percentiles = cut(.data[[TheDetector]],
                      breaks = c(-Inf, percentiles),
                      labels = seq(0, 100, by = probslabel), include.lowest = TRUE))
  
  if (returncutplot == TRUE){
    breaks <- percentiles |> unname()
    Rectangles <- data.frame(xmin = breaks[-length(breaks)],
    xmax = breaks[-1], fill_group = factor(1:(length(breaks)-1))) |> 
      mutate(ymin = 0, ymax = Inf)

    custom_breaks <- c(1, 1e3, 1e4, 1e5, 1e6)

    ColorPallette <- c("#FFDDDD", "#DDFFDD", "#DDDDFF", "#FFFFDD", 
    "#FFDDFF", "#DDFFFF", "#FFCCCC", "#CCFFCC", "#CCCCFF", "#FFFFCC")

    if (!is.null(titleplot)){Title <- paste("Percentile Bins for ", TheDetector)
    } else {Title <- titleplot}

    plot <- ggplot(data, aes(x = .data[[TheDetector]])) + 
      geom_density(fill = "steelblue", alpha = 0.7) +
      geom_rect(data = Rectangles, aes(xmin = xmin,
        xmax = xmax, ymin = ymin, ymax = ymax, fill=fill_group),
        alpha = 0.65, inherit.aes = FALSE) +
      scale_fill_manual(values = ColorPallette, guide = "none") +
      scale_x_log10(limits = c(1, 1e6), breaks = custom_breaks,
      labels = trans_format("log10", math_format(10^.x))) + 
      labs(title = Title,
      x = TheDetector, y = "Density") + theme_bw() +
      annotation_logticks(sides = "b") 

    return(plot)
  }

  if (returntype == "normalized"){
    ByPercentiles <- Normalized %>% mutate(Percentiles = ByPercentiles$Percentiles)
  }

  Samples <- ByPercentiles %>% group_by(Percentiles) %>%
    nest(data = where(is.numeric)) %>%
    mutate(Averaged = map(data, .f=AveragedSignature, stats=stats)) %>%
    select(Percentiles, Averaged) %>% unnest(Averaged) %>% ungroup()

  if (droplowest == TRUE){
    Samples <- Samples |> filter(!Percentiles %in% "0")
  }

  if (output == "data"){
    return(Samples)
  }

  if (returntype == "normalized"){Expression <- "Normalized Values"
  } else {Expression <- "Raw MFI"}

  if (output == "plot"){
  LineCols <- ncol(Samples)
  Melted <- Samples %>%
    pivot_longer(all_of(2:LineCols),names_to = "Detector", values_to = "value")

  Melted$Detector <- factor(Melted$Detector, levels = DetectorOrder)
    
  if (!is.null(titlename)){name <- titlename}

  if (legend == TRUE){
  plot <- ggplot(Melted, aes(x = Detector, y = value, group = Percentiles,
          color = Percentiles)) + geom_line() + scale_color_hue(direction = 1) +
          labs(title = name, x = "Detectors", y = Expression) +
          theme_linedraw() + theme_bw() + theme(axis.title.x = element_text(
          face = "plain"), axis.title.y = element_text(face = "plain"),
          axis.text.x = element_text(size = 5,
          angle = 45, hjust = 1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  } else {
    plot <- ggplot(Melted, aes(x = Detector, y = value, group = Percentiles,
      color = Percentiles)) + geom_line() + scale_color_hue(direction = 1) +
      labs(title = name, x = "Detectors", y = Expression) +
      theme_linedraw() + theme_bw() + theme(legend.position="none", 
      axis.title.x = element_text(face = "plain"),
      axis.title.y = element_text(face = "plain"),
      axis.text.x = element_text(size = 5,
      angle = 45, hjust = 1), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  }
    

  return(plot)
  }
}
