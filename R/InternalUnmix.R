
#' Internal for Luciernaga_Unmix
#'
#' @param cs The passed original cytoset object to extract info from
#' @param StashedDF The stored Time FSC etc
#' @param TheData The Unmixed Data
#' @param Ligands The ligand names
#'
#' @importFrom dplyr filter slice 
#' @importFrom stringr str_detect
#' @importFrom tidyr unnest
#' @importFrom tidyselect where
#' @importFrom flowCore parameters<-
#'
#' @return An internal value
#'
#' @noRd
InternalUnmix <- function(cs, StashedDF, TheData, Ligands){
  fr <- cs[[1, returnType = "flowFrame"]]
  ParamData <- fr@parameters@data
  FluorData <- ParamData |> filter(!str_detect(name, "Time|FSC|SSC")) |>
    slice(1)
  ParamData <- ParamData |> filter(str_detect(name, "Time|FSC|SSC"))

  NewColStart <- ncol(StashedDF)
  AllCols <- ncol(TheData)

  cols <- as.matrix(TheData)
  ncol <- ncol(cols)
  cn <- colnames(cols)

  new_pid <- 1
  new_pid <- seq(new_pid, length.out = ncol)
  new_pid <- paste0("$P", new_pid)

  SecondCN <- cn[NewColStart:AllCols]

  new_pd <- do.call(rbind, lapply(SecondCN, function(i){
    vec <- cols[,i]
    rg <- range(vec)
    data.frame(name = i, desc = NA, range = FluorData$range,
               minRange = FluorData$minRange, maxRange = FluorData$maxRange)
  }))

  new_pd$desc <- Ligands
  new_pd <- rbind(ParamData, new_pd)
  rownames(new_pd) <- new_pid

  new_kw <- fr@description

  NameParams <- new_kw[grepl("^\\$P\\d+N\\d*", names(new_kw))]
  VoltageParams <- new_kw[grepl("^\\$P\\d+V\\d*", names(new_kw))]
  VoltageParams <- c(NA, VoltageParams)
  DisplayParams <- new_kw[grepl("^\\P\\d+DISPLAY\\d*", names(new_kw))]
  TypeParams <- new_kw[grepl("^\\$P\\d+TYPE\\d*", names(new_kw))]

  DescriptionData <- cbind(NameParams, VoltageParams,
     DisplayParams, TypeParams)
  DescriptionData <- as.data.frame(DescriptionData)
  DescriptionData <- DescriptionData |>
    filter(str_detect(NameParams, "Time|FSC|SSC|B1-A"))
  DescriptionData <- DescriptionData %>% unnest(cols = where(is.list))
  #DescriptionData

  Test <- new_kw[!grepl("^\\$P\\d+", names(new_kw))]
  Test <- Test[!grepl("^\\P\\d+", names(Test))]
  Test <- Test[!grepl("^\\$FLOWRATE", names(Test))]
  Test <- Test[!grepl("^\\$CYTOLIB", names(Test))]

  OGLength <- length(Test)

  # Sausage Getting Made, forgive my For-loop
  for (i in new_pid){
    NoDollar <- gsub("$", "", fixed=TRUE, i)
    Test[paste0(i,"B")] <- "32"
    Test[paste0(i,"E")] <- "0,0"
    Test[paste0(i,"N")] <- new_pd[[i,1]]

    TheName <- new_pd[[i, 1]]
    if (!str_detect(TheName, "Time|FSC|SSC")) {Test[paste0(i,"V")] <- "0"
    } else {
      if (str_detect(TheName, "FSC|SSC")){
        Voltage <- DescriptionData |>
          filter(NameParams %in% TheName) |>
          pull(VoltageParams)
        Test[paste0(i,"V")] <- Voltage
      }
    }

    Test[paste0(i,"R")] <- new_pd[[i,5]]

    if (!str_detect(TheName, "FSC|SSC")) {
      Test[paste0(NoDollar,"DISPLAY")] <- "LOG"
    } else {Test[paste0(NoDollar,"DISPLAY")] <- "LIN"}

    if (str_detect(TheName, "Time")) {Test[paste0(i,"TYPE")] <- "Time"
    } else if (str_detect(TheName, "FSC")){
      Test[paste0(i,"TYPE")] <- "Forward_Scatter"
    } else if (str_detect(TheName, "SSC")){
      Test[paste0(i,"TYPE")] <- "Side_Scatter"
    } else {Test[paste0(i,"TYPE")] <- "Unmixed_Fluorescence"}
  }

  index <- which(names(Test) == "$CYTSN")
  StartLength <- OGLength+1
  FinalLength <- length(Test)

  Subset <- Test[StartLength:FinalLength]
  Residual <- Test[-(StartLength:FinalLength)]
  new_kw <- append(Residual, Subset, after = index)

  TheSpilloverNames <- cn[!grepl("Time|FSC|SSC", cn)]
  MatrixSize <- length(TheSpilloverNames)
  NewMatrix <- matrix(0, nrow = MatrixSize, ncol = MatrixSize, byrow = TRUE)
  diag(NewMatrix) <- 1
  colnames(NewMatrix) <- TheSpilloverNames

  new_kw$`$SPILLOVER` <- NewMatrix
  new_kw$`CREATOR` <- "Luciernaga 0.99.1"

  # Adding back to extract again?
  fr@exprs <- cols
  pData(parameters(fr)) <- new_pd

  UpdatedParameters <- parameters(fr)

  new_fcs <- new("flowFrame", exprs=cols, parameters=UpdatedParameters,
                 description=new_kw)

  return(new_fcs)
}