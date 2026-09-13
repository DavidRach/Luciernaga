#' Internal for QC_FlowSightPDF, parses the second page,
#'  returns a data.frame
#'
#' @param x The second page of text parsed from the QC report
#'
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr bind_cols
#'
#' @noRd
FlowSightPageTwo <- function(x) {
  lines <- strsplit(x, "\n")[[1]]
  lines <- lines[nzchar(trimws(lines))]

  SideScatterCalibrationLine <- grep("Side Scatter Calibration", lines)
  if (length(SideScatterCalibrationLine) == 1) {
    SideScatterCalibrationLine_DF <- TwoLineSandwhich(
      x = SideScatterCalibrationLine, data = lines)
    #Ignoring second (duplicated???) power listing for now
  }

  Power405nm <- grep("405nm Laser Power Test", lines)
  if (length(Power405nm) == 1) {
    Power405nm_DF <- TwoLineSandwhich(x = Power405nm, data = lines)
  }

  Power488nm <- grep("488nm Laser Power Test", lines)
  if (length(Power488nm) == 1) {
    Power488nm_DF <- TwoLineSandwhich(x = Power488nm, data = lines)
  }

  Power561nm <- grep("561nm Laser Power Test", lines)
  if (length(Power561nm) == 1) {
    Power561nm_DF <- TwoLineSandwhich(x = Power561nm, data = lines)
  }

  Power642nm <- grep("642nm Laser Power Test", lines)
  if (length(Power642nm) == 1) {
    Power642nm_DF <- TwoLineSandwhich(x = Power642nm, data = lines)
  }

  Power785nm <- grep("785nm Laser Power Test", lines)
  if (length(Power785nm) == 1) {
    Power785nm_DF <- TwoLineSandwhich(x = Power785nm, data = lines)
  }

  BrightfieldAlignmentLine <- grep("Brightfield Alignment Test", lines)
  if (length(BrightfieldAlignmentLine) == 1) {
    BrightfieldAlignmentLine <- lines[BrightfieldAlignmentLine]
    BrightfieldAlignment_DF <- TwoPartSplits(BrightfieldAlignmentLine)
  }

  BrightfieldUniformityLine <- grep("Brightfield Uniformity Test", lines)
  CameraNoiseLine <- grep("Camera Noise Test", lines)

  if (length(BrightfieldUniformityLine) == 1) {
    BrightfieldLines <- lines[
      BrightfieldUniformityLine:(CameraNoiseLine - 1)]
    BridgeData <- BridgeSplits(BrightfieldLines)
    Brightfield_DF <- BridgeData[1]
  }

  CameraNoiseLine <- grep("Camera Noise Test", lines)
  if (length(CameraNoiseLine) == 1) {
    CameraNoiseLineOne <- lines[CameraNoiseLine]
    CameraNoise_DF <- TwoPartSplits(CameraNoiseLineOne)
    CameraNoiseLineTwo <- lines[CameraNoiseLine + 1]
    Values <- as.numeric(unlist(strsplit(
      sub(".*Results:\\s*", "", CameraNoiseLineTwo), "\\s+")))
    CameraData <- data.frame(Results = Values)
  }

  AxialStabilityLine <- grep("Flow Core Axial Stability Test", lines)
  if (length(AxialStabilityLine) == 1) {
    AxialStabilityLine <- lines[AxialStabilityLine]
    AxialStability_DF <- SandwhichSplits(AxialStabilityLine)
  }

  LateralStabilityLine <- grep("Flow Core Lateral Stability Test", lines)
  if (length(LateralStabilityLine) == 1) {
    LateralStabilityLine <- lines[LateralStabilityLine]
    LateralStability_DF <- SandwhichSplits(LateralStabilityLine)
  }

  FlowCorePositionLines <- grep("Flow Core Position Test", lines)
  if (length(FlowCorePositionLines) == 1) {
    FlowCorePositionLines <- lines[
      FlowCorePositionLines:(FlowCorePositionLines + 3)]
    BridgeData <- BridgeSplits(FlowCorePositionLines)
    FlowCorePosition_DF <- BridgeData[[1]]
    BridgeDataset <- BridgeData[[2]]
    BridgeDataset$Name <- paste0(colnames(
      FlowCorePosition_DF)[1], " ", BridgeDataset$Name)
    BridgeDataset <- BridgeDataset |>
      tidyr::pivot_wider(names_from = Name, values_from = Value)
    FlowCorePosition_DF <- bind_cols(FlowCorePosition_DF, BridgeDataset)
  }

  FocusPercentageLine <- grep("Focus Percentage Test", lines)
  if (length(FocusPercentageLine) == 1) {
    FocusPercentageLine <- lines[FocusPercentageLine]
    FocusPercentage_DF <- SandwhichSplits(FocusPercentageLine)
  }

  FocusUniformityLine <- grep("Focus Uniformity Test", lines)
  if (length(FocusUniformityLine) == 1) {
    FocusUniformityLine <- lines[FocusUniformityLine]
    FocusUniformity_DF <- SandwhichSplits(FocusUniformityLine)
  }

  ImageQualityLine <- grep("Image Quality Test", lines)
  if (length(ImageQualityLine) == 1) {
    ImageQualityLine <- lines[ImageQualityLine]
    ImageQuality_DF <- SandwhichSplits(ImageQualityLine)
  }

  SecondPageDF <- bind_cols(SideScatterCalibrationLine_DF,
    Power405nm_DF, Power488nm_DF, Power561nm_DF, Power642nm_DF,
    Power785nm_DF, BrightfieldAlignment_DF, Brightfield_DF,
    CameraNoise_DF, AxialStability_DF, LateralStability_DF,
    FlowCorePosition_DF, FocusPercentage_DF, FocusUniformity_DF,
    ImageQuality_DF)

  return(SecondPageDF)
}