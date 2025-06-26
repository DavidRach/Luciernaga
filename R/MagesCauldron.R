#' Playing around with Peter Mage's Hotspot calculation, extending the preprint concepts. 
#' and the .
#' 
#' @param panelfluors A vector fluorophores in your panel (names matching those in QC_ReferenceLibrary)
#' @param unstained Luciernaga_QC ReturnSignature output containing just detector columns
#' @param returnType Default is plot, alternate data
#' @param swapname Name of the Fluorophore to replace with experimental signature,
#'  see QC_ReferenceLibrary for exact formatting
#' @param swapvalue Just the detector columns for the swapname Fluorophore
#' @param outpath Default NULL, file.path to store the savePlot outputs
#' @param filename Default HotspotMatrix, specifies name to save file as
#' @param device Desired storage format, default is "png"
#' @param width Desired height for saved plot, default is 15
#' @param height Desired height for saved plot, default  is 15
#' 
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom MASS ginv
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggsave
#' 
#' @return Either a plot or the underlying matrix
#' 
#' @noRd
MagesCauldron <- function(panelfluors, unstained, returnType="plot", savePlot=FALSE,
   outpath=NULL, filename=NULL, device="png", width=15, height=15, swapname=NULL,
   swapvalue=NULL){
  
  DetectorLength <- ncol(unstained)

  if (any(unstained > 1)){
    unstained[unstained < 0] <- 0
    A <- do.call(pmax, unstained)
    unstained <- unstained/A
  }

  Vaiya <- InstrumentReferences(NumberDetectors = DetectorLength)
  TheUnstained <- unstained |> mutate(Fluorophore="Unstained") |>
       relocate(Fluorophore, .before=1)
  TheseFluorophores <- panelfluors
  Data <- Vaiya |> filter(Fluorophore %in% TheseFluorophores)
  Data <- Data |> select(-Instrument)
  Data <- Data |> pivot_wider(names_from="Detector",
   values_from="AdjustedY")
  colnames(TheUnstained) <- gsub("-A", "", gsub("-H", "", colnames(TheUnstained)))
  Data <- bind_rows(Data, TheUnstained)
  TheseFluorophores <- c(TheseFluorophores, paste0("Unstained", seq_len(nrow(unstained))))
  Data$Fluorophore <- factor(Data$Fluorophore, levels=TheseFluorophores)
  Data <- Data |> arrange(desc(Fluorophore))
  Data <- Data |> arrange(Fluorophore)

  if (!is.null(swapname)){
    Replacement <- swapvalue |> mutate("Fluorophore"=swapname) |>
      relocate(Fluorophore, .before=1)
    Index <- which(Data$Fluorophore == swapname)
    Data[Index,] <- Replacement
  }

  Similarity <- Luciernaga_Cosine(Data, returntype="data", rearrange=FALSE)

  PseudoInverse <- ginv(Similarity)
  Absolute <- abs(PseudoInverse)
  TheSQRT <- sqrt(Absolute)
  row.names(TheSQRT) <- TheseFluorophores
  Hotspots <- round(TheSQRT, 2)

  if (returnType == "data"){
       return(Hotspots)
  }

  final_matrix <- Hotspots
  lower_tri <- final_matrix
  lower_tri[upper.tri(final_matrix)] <- NA  

  Longer <- as.data.frame(lower_tri) |> rownames_to_column("row") |> 
    pivot_longer(cols = -row, names_to = "col_num", names_prefix = "V",
    values_to = "value") |> mutate(row = factor(row, levels=rev(TheseFluorophores)),
    col = as.numeric(col_num), col_label = factor(
      TheseFluorophores[col], levels = TheseFluorophores)) |>
    filter(!is.na(value)) |> arrange(row)

  Plot <- ggplot(Longer, aes(x = col_label, y = row, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(value, 2)), 
  color = ifelse(Longer$value > max(Longer$value)/2, "white", "black"),
            size = 3, fontface = "bold") +
  scale_fill_gradient(low = "white",high = "red2",na.value = "white",
  name="Hotspot") +
  scale_x_discrete(position = "bottom") + labs(x = NULL, y = NULL) +
  coord_fixed() + theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),legend.position = "right",
  axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,),
  axis.text.y = element_text(hjust = 1))

  if (savePlot == TRUE){
    if (is.null(outpath)){outpath <- getwd()}
    if (is.null(filename)){filename <- paste0("HotspotMatrix.")}
    FinalPath <- file.path(outpath, filename)
       ggsave(filename = FinalPath, plot = Plot, device = device,
       width = width, height = height, units = "in",
       dpi = 300, bg = "white")

  } else {return(Plot)}
  }

#' Internal, visualizes changes between hotspot plots
#' 
#' @param HotspotBrighter Brighter hotspot data.frame
#' @param HotspotDimmer Dimmer hotspot data.frame
#' 
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggsave
#' 
#' @return The residual plot of the difference between the two hotspots
#' 
#' @noRd
TailOfANewt <- function(HotspotBrighter, HotspotDimmer, savePlot=FALSE,
  outpath=NULL, filename=NULL, device="png", width=15, height=15){
  
  TheseFluorophores <- HotspotBrighter |> row.names()
  
  residuals <- HotspotBrighter[, -1] - HotspotDimmer[, -1]
  final_matrix <- residuals
  lower_tri <- residuals
  lower_tri[upper.tri(final_matrix)] <- NA  

  Longer <- as.data.frame(lower_tri) |> rownames_to_column("row") |> 
    pivot_longer(cols = -row, names_to = "col_num", names_prefix = "V",
    values_to = "value") |> mutate(row = factor(row, levels=rev(TheseFluorophores)),
    col = as.numeric(col_num), col_label = factor(
      TheseFluorophores[col], levels = TheseFluorophores)) |>
    filter(!is.na(value)) |> arrange(row)

  Plot <- ggplot(Longer, aes(x = col_label, y = row, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(value, 2)), 
  color = ifelse(Longer$value > max(Longer$value)/2, "white", "black"),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue2", high = "red2", mid="white", 
    midpoint=0, na.value = "white", name="Hotspot",
   limits = c(-max(abs(Longer$value)), max(abs(Longer$value)))) +
  scale_x_discrete(position = "bottom") + labs(x = NULL, y = NULL) +
  coord_fixed() + theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),legend.position = "right",
  axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,),
  axis.text.y = element_text(hjust = 1))

  if (savePlot == TRUE){
    if (is.null(outpath)){outpath <- getwd()}
    if (is.null(filename)){filename <- paste0("DifferenceInHotspotMatrix.")}

    FinalPath <- file.path(outpath, filename)
    
    ggsave(filename = FinalPath, plot = Plot, device = device,
       width = width, height = height, units = "in",
       dpi = 300, bg = "white")

  } else {return(Plot)}
}