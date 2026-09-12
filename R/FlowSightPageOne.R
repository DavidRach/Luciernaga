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