#' Helper function for QC PDF conversions. Tidyes two chunk items into single data.frame cell
#' 
#' @param x A string line read in from the pdf
#' 
#' @noRd
TwoPartSplits <- function(x){
    x <- gsub("mum:", "mum:  ", x)
    parts <- strsplit(trimws(x), "\\s{2,}")[[1]]
    keys <- trimws(gsub(":$", "", parts[seq(1, length(parts), by = 2)]))
    vals <- trimws(parts[seq(2, length(parts), by = 2)])
    if (length(keys) > length(vals)) {
    keys <- keys[seq_len(length(vals))]}
    data <- as.data.frame(as.list(setNames(vals, keys)), check.names=FALSE)
    return(data)
}

#' Helper function for QC PDF conversions. Combines the bread ends, then handles the middle.
#' 
#' @param x A string line read in from the pdf
#' 
#' @noRd
SandwhichSplits <- function(x){
    parts <- strsplit(trimws(x), "\\s{2,}")[[1]]
    First <- paste0(parts[1], " ", parts[2])
    Second <- paste0(parts[1], ": ", parts[3])
    parts <- c(First, Second)
    data <- as.data.frame(
        setNames(
            list(
            trimws(sub(".*:","", parts[1])),
            trimws(sub(".*:","", parts[2]))
            ),
            trimws(sub(":.*","", parts))
        ), stringsAsFactors = FALSE, check.names=FALSE)
    return(data)
    }

#' Helper function for QC PDF conversions. Combines the bread ends,
#'  then handles the middle, then tackles the middle continuation 
#'  on the subsequent row.
#' 
#' @param x A string line index from the pdf
#' @param data The lines from the pdf for filtering on
#' 
#' @importFrom dplyr bind_cols
#' 
#' @noRd
TwoLineSandwhich <- function(x, data){
    FirstLine <- data[x]
    SecondLine <- data[x+1]

    parts <- strsplit(trimws(FirstLine), "\\s{2,}")[[1]]
    First <- paste0(parts[1], " ", parts[2])
    Second <- paste0(parts[1], ": ", parts[3])
    parts <- c(First, Second)
    FirstData <- as.data.frame(
        setNames(
            list(
            trimws(sub(".*:","", parts[1])),
            trimws(sub(".*:","", parts[2]))
            ),
            trimws(sub(":.*","", parts))
        ), stringsAsFactors = FALSE, check.names=FALSE)

    SecondLine <- gsub("mum:", "mum:  ", SecondLine)
    SecondLine <- gsub("ity:", "ity:  ", SecondLine)
    SecondLine <- gsub("wer:", "wer:  ", SecondLine)
    parts <- strsplit(trimws(SecondLine), "\\s{2,}")[[1]]
    keys <- trimws(gsub(":$", "", parts[seq(1, length(parts), by = 2)]))
    vals <- trimws(parts[seq(2, length(parts), by = 2)])
    SecondData <- as.data.frame(as.list(setNames(vals, keys)), check.names=FALSE)
    colnames(SecondData) <- paste0(colnames(FirstData[2]), " ", colnames(SecondData))
    TheData <- bind_cols(FirstData, SecondData)
    return(TheData)
}


#' Helper function for QC PDF conversions. Combines the bread ends,
#'  then separates out the middle as it's own separate data.frame,
#'  returning both as components of a list. 
#' 
#' @param data The lines from the pdf for the respective section
#' 
#' @noRd
BridgeSplits <- function(data){
    Total <- length(data)
    parts <- strsplit(trimws(data), "\\s{2,}")[[1]]
    Second <- paste0(parts[1], ": ", parts[3])
    Bridge <- as.data.frame(t(setNames(strsplit(Second, ":\\s*")[[1]][2],
                               strsplit(Second, ":\\s*")[[1]][1])), check.names=FALSE)
    TheParts <- strsplit(trimws(data), "\\s{2,}")[2:6]
    TheVector <- c(parts[2], unlist(TheParts))

    BridgeData <- do.call(rbind, lapply(TheVector, function(s) {
    parts <- strsplit(s, ":\\s*")[[1]]
    data.frame(
        Name  = parts[1],
        Value = as.numeric(parts[2]),
        stringsAsFactors = FALSE
    )
    }))

    ReturnVals <- list(Bridge, BridgeData)
    return(ReturnVals)
}

#' Parses the two-page Amnis FlowSight QC report, returning tidy dataframe.
#'
#' @param x A file.path to the respective QC report pdf
#' 
#' @importFrom pdftools pdf_text
#' @importFrom dplyr bind_cols
#' 
#' @export
#' 
#' @examples A <- 2+2
QC_FlowSightPDF <- function(x){
  text <- pdftools::pdf_text(x)
  NumberPages <- length(text)
  PageOne <- FlowSightPageOne(x=text[1])
  PageTwo <- FlowSightPageTwo(x=text[2])
  Completed <- bind_cols(PageOne, PageTwo)
  return(Completed)
}


#' Internal for QC_FlowSightPDF, parses the second page, returns a data.frame
#' 
#' @param x The second page of text parsed from the QC report
#'
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr bind_cols
#' 
#' @noRd
FlowSightPageTwo <- function(x){
  lines <- strsplit(x, "\n")[[1]]
lines <- lines[nzchar(trimws(lines))]

SideScatterCalibrationLine <- grep("Side Scatter Calibration", lines)
if (length(SideScatterCalibrationLine) == 1){
    SideScatterCalibrationLine_DF <- TwoLineSandwhich(x=SideScatterCalibrationLine, data=lines)
    #Ignoring second (duplicated???) power listing for now
}

Power405nm <- grep("405nm Laser Power Test", lines)
if (length(Power405nm) == 1){
    Power405nm_DF <- TwoLineSandwhich(x=Power405nm, data=lines)
}

Power488nm <- grep("488nm Laser Power Test", lines)
if (length(Power488nm) == 1){
    Power488nm_DF <- TwoLineSandwhich(x=Power488nm, data=lines)
}

Power561nm <- grep("561nm Laser Power Test", lines)
if (length(Power561nm) == 1){
    Power561nm_DF <- TwoLineSandwhich(x=Power561nm, data=lines)
}

Power642nm <- grep("642nm Laser Power Test", lines)
if (length(Power642nm) == 1){
    Power642nm_DF <- TwoLineSandwhich(x=Power642nm, data=lines)
}

Power785nm <- grep("785nm Laser Power Test", lines)
if (length(Power785nm) == 1){
    Power785nm_DF <- TwoLineSandwhich(x=Power785nm, data=lines)
}

BrightfieldAlignmentLine <- grep("Brightfield Alignment Test", lines)
if (length(BrightfieldAlignmentLine) == 1){
    BrightfieldAlignmentLine <-lines[BrightfieldAlignmentLine]
    BrightfieldAlignment_DF <- TwoPartSplits(BrightfieldAlignmentLine)
}

BrightfieldUniformityLine <- grep("Brightfield Uniformity Test", lines)
CameraNoiseLine <- grep("Camera Noise Test", lines)

if (length(BrightfieldUniformityLine) == 1){
    BrightfieldLines <- lines[BrightfieldUniformityLine:(CameraNoiseLine-1)]
    BridgeData <- BridgeSplits(BrightfieldLines)
    Brightfield_DF <- BridgeData[1]
}

CameraNoiseLine <- grep("Camera Noise Test", lines)
if (length(CameraNoiseLine) == 1){
    CameraNoiseLineOne <- lines[CameraNoiseLine]
    CameraNoise_DF <- TwoPartSplits(CameraNoiseLineOne)
    CameraNoiseLineTwo <- lines[CameraNoiseLine+1]
    Values <- as.numeric(unlist(strsplit(
        sub(".*Results:\\s*", "", CameraNoiseLineTwo), "\\s+")))
    CameraData <- data.frame(Results = Values)
}

AxialStabilityLine <- grep("Flow Core Axial Stability Test", lines)
if (length(AxialStabilityLine) == 1){
    AxialStabilityLine <-lines[AxialStabilityLine]
    AxialStability_DF <- SandwhichSplits(AxialStabilityLine)
}

LateralStabilityLine <- grep("Flow Core Lateral Stability Test", lines)
if (length(LateralStabilityLine) == 1){
    LateralStabilityLine <-lines[LateralStabilityLine]
    LateralStability_DF <- SandwhichSplits(LateralStabilityLine)
}

FlowCorePositionLines <- grep("Flow Core Position Test", lines)
if (length(FlowCorePositionLines) == 1){
    FlowCorePositionLines <- lines[FlowCorePositionLines:(FlowCorePositionLines+3)]
    BridgeData <- BridgeSplits(FlowCorePositionLines)
    FlowCorePosition_DF <- BridgeData[[1]]
    BridgeDataset <- BridgeData[[2]]
    BridgeDataset$Name <- paste0(colnames(FlowCorePosition_DF)[1], " ", BridgeDataset$Name)
    BridgeDataset <- BridgeDataset |> tidyr::pivot_wider(names_from=Name, values_from=Value)
    FlowCorePosition_DF <- bind_cols(FlowCorePosition_DF, BridgeDataset)
}

FocusPercentageLine <- grep("Focus Percentage Test", lines)
if (length(FocusPercentageLine) == 1){
    FocusPercentageLine <-lines[FocusPercentageLine]
    FocusPercentage_DF <- SandwhichSplits(FocusPercentageLine)
}

FocusUniformityLine <- grep("Focus Uniformity Test", lines)
if (length(FocusUniformityLine) == 1){
    FocusUniformityLine <-lines[FocusUniformityLine]
    FocusUniformity_DF <- SandwhichSplits(FocusUniformityLine)
}

ImageQualityLine <- grep("Image Quality Test", lines)
if (length(ImageQualityLine) == 1){
    ImageQualityLine <-lines[ImageQualityLine]
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

#' Internal for QC_FlowSightPDF, parses the first page, returns a data.frame
#' 
#' @param x The first page of text parsed from the QC report
#'
#' @importFrom stringr str_extract
#' @importFrom lubridate mdy
#' @importFrom dplyr bind_cols
#' 
#' @noRd
FlowSightPageOne <- function(x){
lines <- strsplit(x, "\n")[[1]]
lines <- lines[nzchar(trimws(lines))]

DateTimeLine <- grep("FlowSight", lines)
if (length(DateTimeLine) == 1){
   DateTimeLine <-lines[DateTimeLine]
   DateTime_str <- str_extract(x, "[A-Za-z]+,\\s+[A-Za-z]+\\s+\\d{2},\\s+\\d{4}")
   DateTime <- lubridate::mdy(DateTime_str)
   DateTime <- data.frame(DateTime=DateTime, PDF=DateTime_str)
}

FocusAdjustLine <- grep("Focus Adjustor Calibration", lines)
if (length(FocusAdjustLine) == 1){
    FocusAdjustLine  <-lines[FocusAdjustLine]
    FocusAdjust_DF <- TwoPartSplits(FocusAdjustLine)
}

AutosamplerLine <- grep("Autosampler Nest Calibration", lines)
if (length(AutosamplerLine) == 1){
    AutosamplerLine  <-lines[AutosamplerLine]
    Autosampler_DF <- SandwhichSplits(AutosamplerLine)
}

FrameLine <- grep("Frame Offset Calibration", lines)
if (length(FrameLine) == 1){
    FrameLine <-lines[FrameLine]
    Frame_DF <- SandwhichSplits(FrameLine)
}

CoreStageLine <- grep("Core Stage Position Calibration", lines)
if (length(CoreStageLine) == 1){
    CoreStageLine <-lines[CoreStageLine]
    CoreStage_DF <- SandwhichSplits(CoreStageLine)
}

SpatialLine <- grep("Spatial Offsets Calibration", lines)
if (length(SpatialLine) == 1){
   SpatialLineOne <-lines[SpatialLine]
   SpatialLine_DF <- TwoPartSplits(SpatialLineOne)
   
   SpatialXLine <-lines[SpatialLine+1]
   SpatialYLine <-lines[SpatialLine+2]

   SpatialX <- as.numeric(unlist(strsplit(
    trimws(sub("X Offsets:", "", SpatialXLine)), "\\s+")))
   SpatialY <- as.numeric(unlist(strsplit(
    trimws(sub("Y Offsets:", "", SpatialYLine)), "\\s+")))

    SpatialData <- data.frame(
    SpatialXOffset = SpatialX,
    SpatialYOffset = SpatialY,
    stringsAsFactors = FALSE
    )
}

DarkLine <- grep("Dark Current Calibration", lines)
if (length(DarkLine) == 1){
    DarkLineOne <-lines[DarkLine]
    DarkLine_DF <- TwoPartSplits(DarkLineOne)
   
   DarkLineTwo <-lines[DarkLine+1]
    DarkLineTwo <- TwoPartSplits(DarkLineTwo)
    colnames(DarkLineTwo) <- paste0("Dark Current ", colnames(DarkLineTwo))
   DarkLineThree <-lines[DarkLine+2]
    DarkLineThree <- TwoPartSplits(DarkLineThree)
    colnames(DarkLineThree) <- paste0("Dark Current ", colnames(DarkLineThree))

    DarkLineFour <-lines[DarkLine+3]
    DarkLineFour <- as.numeric(unlist(strsplit(
        trimws(sub("Channel Means:", "", DarkLineFour)), "\\s+")))
    DarkLineData <- data.frame(DarkCurrentChannelMeans = DarkLineFour)

    DarkLine_DF <- bind_cols(DarkLine_DF, DarkLineTwo, DarkLineThree)
}

BrightfieldLine <- grep("Brightfield XTalk Coefficient Calibration", lines)
Calibration405Line <- grep("405nm Horizontal Laser Calibration", lines)

if (length(BrightfieldLine) == 1){
    BrightfieldLineOne <-lines[BrightfieldLine]
    BrightfieldLineOne_DF <- TwoPartSplits(BrightfieldLineOne)

    BrightfieldLines <- grep("Cross Talk Coefficients", lines, value = TRUE)
    chan <- sub(".*Chan([0-9]+):.*", "\\1", BrightfieldLines)
    vals <- lapply(BrightfieldLines, function(x) as.numeric(
        unlist(regmatches(x, gregexpr("[0-9\\.]+", x)))))
    BF_CrossTalk <- as.data.frame(do.call(rbind, vals))
    rownames(BF_CrossTalk) <- paste0("Chan", chan)
    colnames(BF_CrossTalk) <- paste0("ToChan", 1:ncol(BF_CrossTalk))  
}

Calibration405Line <- grep("405nm Horizontal Laser Calibration", lines)
if (length(Calibration405Line) == 1){
    Calibration405Line <-lines[Calibration405Line]
    Calibration405_DF <- SandwhichSplits(Calibration405Line)
}

Calibration488Line <- grep("488nm Horizontal Laser Calibration", lines)
if (length(Calibration488Line) == 1){
    Calibration488Line <-lines[Calibration488Line]
    Calibration488_DF <- SandwhichSplits(Calibration488Line)
}

Calibration561Line <- grep("561nm Horizontal Laser Calibration", lines)
if (length(Calibration561Line) == 1){
    Calibration561Line <-lines[Calibration561Line]
    Calibration561_DF <- SandwhichSplits(Calibration561Line)
}

Calibration642Line <- grep("642nm Horizontal Laser Calibration", lines)
if (length(Calibration642Line) == 1){
    Calibration642Line <-lines[Calibration642Line]
    Calibration642_DF <- SandwhichSplits(Calibration642Line)
}

Calibration785Line <- grep("785nm Horizontal Laser Calibration", lines)
if (length(Calibration785Line) == 1){
    Calibration785Line <-lines[Calibration785Line]
    Calibration785_DF <- SandwhichSplits(Calibration785Line)
}

RetroLine <- grep("Retro Calibration", lines)
if (length(RetroLine) == 1){
    RetroLine <-lines[RetroLine]
    Retro_DF <- SandwhichSplits(RetroLine)
}

PageOneData <- bind_cols(DateTime, FocusAdjust_DF, Autosampler_DF, Frame_DF,
   CoreStage_DF, SpatialLine_DF, DarkLine_DF, BrightfieldLineOne_DF,
   Calibration405_DF, Calibration488_DF, Calibration561_DF,
    Calibration642_DF, Calibration785_DF, Retro_DF)

return(PageOneData)
}
