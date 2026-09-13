#' Playing around with Peter Mage's Hotspot calculation, extending the
#' preprint concepts and the .
#'
#' @param panelfluors A vector fluorophores in your panel (names matching
#'   those in QC_ReferenceLibrary)
#' @param unstained Luciernaga_QC ReturnSignature output containing just
#'   detector columns
#' @param returnType Default is plot, alternate data
#' @param savePlot TBD
#' @param swapname Name of the Fluorophore to replace with experimental
#'   signature, see QC_ReferenceLibrary for exact formatting
#' @param swapvalue Just the detector columns for the swapname Fluorophore
#' @param outpath Default NULL, file.path to store the savePlot outputs
#' @param filename Default HotspotMatrix, specifies name to save file as
#' @param device Desired storage format, default is "png"
#' @param width Desired height for saved plot, default is 15
#' @param height Desired height for saved plot, default is 15
#' @param NumberDetectors Default NULL, used when unstained is NULL
#'
#' @importFrom dplyr mutate relocate filter select bind_rows arrange desc
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom MASS ginv
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#'   scale_x_discrete labs coord_fixed theme_minimal theme element_blank
#'   element_text ggsave
#'
#' @return Either a plot or the underlying matrix
#'
#' @noRd
MagesCauldron <- function(panelfluors, unstained=NULL, returnType="plot",
 savePlot=FALSE, outpath=NULL, filename=NULL, device="png", width=15,
 height=15, swapname=NULL, swapvalue=NULL, NumberDetectors=NULL) {

  if (!is.null(unstained)) {
    DetectorLength <- ncol(unstained)

    if (any(unstained > 1)) {
      unstained[unstained < 0] <- 0
      A <- do.call(pmax, unstained)
      unstained <- unstained/A
    }

    TheUnstained <- unstained |> mutate(Fluorophore="Unstained") |>
      relocate(Fluorophore, .before=1)
  } else {
    if (is.null(NumberDetectors)) {
      stop("When not providing unstained, provide NumberDetectors argument")
    }
    DetectorLength <- NumberDetectors
  }

  Vaiya <- InstrumentReferences(NumberDetectors = DetectorLength)
  TheseFluorophores <- panelfluors
  Data <- Vaiya |> filter(Fluorophore %in% TheseFluorophores)
  Data <- Data |> select(-Instrument)
  Data <- Data |> pivot_wider(names_from="Detector", values_from="AdjustedY")

  if (!is.null(unstained)) {
    colnames(TheUnstained) <- gsub("-A", "",
      gsub("-H", "", colnames(TheUnstained)))
    Data <- bind_rows(Data, TheUnstained)
    TheseFluorophores <- c(TheseFluorophores,
      paste0("Unstained", seq_len(nrow(unstained))))
  }

  Data$Fluorophore <- factor(Data$Fluorophore, levels=TheseFluorophores)
  Data <- Data |> arrange(desc(Fluorophore))
  Data <- Data |> arrange(Fluorophore)

  if (!is.null(swapname)) {
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

  final_matrix <- Hotspots
  lower_tri <- final_matrix
  lower_tri[upper.tri(final_matrix)] <- NA

  if (returnType == "data") {
    ReturnTri <- data.frame(lower_tri)
    colnames(ReturnTri) <- row.names(ReturnTri)
    ReturnTri <- ReturnTri |>
      rownames_to_column("Fluorophore")
    return(ReturnTri)
  }

  Longer <- as.data.frame(lower_tri) |> rownames_to_column("row") |>
    pivot_longer(cols = -row, names_to = "col_num", names_prefix = "V",
    values_to = "value") |>
    mutate(row = factor(row, levels=rev(TheseFluorophores)),
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
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_text(hjust = 1))

  if (savePlot == TRUE) {
    if (is.null(outpath)) {
      outpath <- getwd()
    }
    if (is.null(filename)) {
      filename <- paste0("HotspotMatrix.")
    }
    FinalPath <- file.path(outpath, filename)
    ggsave(filename = FinalPath, plot = Plot, device = device,
      width = width, height = height, units = "in",
      dpi = 300, bg = "white")
  } else {
    return(Plot)
  }
}