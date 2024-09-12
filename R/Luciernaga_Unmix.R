#' OLS unmixing of a Gating Set object
#'
#' @param x A Gating Set object
#' @param controlData The matrix of single color controls generated from Luciernaga
#' @param sample.name The keyword containing the fcs file name
#' @param addon Additional addon to append to the new .fcs file name
#' @param removestrings A list of values to remove from name
#' @param subset A gating hierarchy level to sort cells at, expression values retrieved
#' from these
#' @param multiplier A number to scale the OLS coefficients by
#' @param outpath The return folder for the .fcs files
#' @param Verbose For troubleshooting name after removestrings
#' @param PanelPath Location to a panel.csv containing correct order of fluorophores
#' @param Experimental Debug for unmixing de no flowframe
#'
#' @importFrom flowCore keyword
#' @importFrom dplyr pull
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom utils read.csv
#' @importFrom dplyr arrange
#' @importFrom stats lsfit
#' @importFrom flowWorkspace realize_view
#' @importFrom flowWorkspace cf_append_cols
#' @importFrom flowCore write.FCS
#'
#' @return A new .fcs file with the new columns appended
#' @export
#'
#' @examples NULL

Luciernaga_Unmix <- function(x, controlData, sample.name, addon, removestrings,
                             subset, multiplier, outpath, Verbose, PanelPath, Experimental=FALSE){

  if (length(sample.name) == 2){
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(x, first)
    first <- first #%>% pull(.)
    second <- keyword(x, second)
    second <- second #%>% pull(.)
    name <- paste(first, second, sep="_")
  } else {name <- keyword(x, sample.name)}

  name <- NameCleanUp(name, removestrings=removestrings)
  if (Verbose == TRUE){message("After removestrings, name is ", name)}

  cs <- gs_pop_get_data(x, subset)
  Data <- exprs(cs[[1]])
  Data <- data.frame(Data, check.names = FALSE)

  #For Future Column Reordering
  OriginalColumnsVector <- colnames(Data)
  OriginalColumns <- colnames(Data)
  OriginalColumns <- data.frame(OriginalColumns, check.names = FALSE)
  OriginalColumnsIndex <- OriginalColumns %>% mutate(IndexLocation = 1:nrow(.))

  #For Future Row Reordering
  Backups <- Data %>% mutate(Backups = 1:nrow(Data)) %>% select(Backups)

  #Stashing Away Time FSC SSC For Later Use
  StashedDF <- Data[,grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  StashedDF <- cbind(Backups, StashedDF)

  #Consolidating Columns Going Forward
  TheSampleData <- Data[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Data))]
  BackupNames <- colnames(TheSampleData)

  # if(ReorderCOlumns == TRUE){}

  if (!is.data.frame(PanelPath)){Panel <- read.csv(PanelPath, check.names=FALSE)
  } else {Panel <- PanelPath}

  CorrectColumnOrder <- Panel %>% pull(Fluorophore)
  CorrectColumnOrder <- gsub("-A$", "", CorrectColumnOrder)

  NewControlData <- controlData %>% arrange(match(Fluorophore, CorrectColumnOrder))


  Newest <- NewControlData %>% pull(Fluorophore)

  if (!identical(CorrectColumnOrder, Newest)){
    print(Newest)
    stop("Column Reordering Failed, printed order output for troubleshooting:")
  }


  # Ordering the Single Color Control Columns to Match the Samples
  TheControlData <- NewControlData[names(TheSampleData)] #Sans Fluor and Ligand
  NewNames <- NewControlData %>% select(Fluorophore)
  NewNames$Fluorophore <- paste0(NewNames$Fluorophore, "-A") #Restored
  NewNames <- NewNames %>% pull(Fluorophore)

  Ligands <- NewControlData %>% pull(Ligand)

  #OLS
  LeastSquares <- lsfit(x = t(TheControlData), y = t(TheSampleData), intercept = FALSE)
  UnmixedData <- t(LeastSquares$coefficients)
  UnmixedData2 <- UnmixedData*multiplier

  colnames(UnmixedData2) <- NewNames

  TheData <- cbind(StashedDF, UnmixedData2)
  TheData <- TheData %>% select(-Backups)
  rownames(TheData) <- NULL
  #TheData

  # Retrieving Raw Data Info
  fr <- cs[[1, returnType = "flowFrame"]]

  # Finding Values That Stay Put
  ParamData <- fr@parameters@data
  #ParamData$name <- as.character(ParamData$name)
  FluorData <- ParamData %>% dplyr::filter(!str_detect(name, "Time|FSC|SSC")) %>% slice(1)
  ParamData <- ParamData %>% dplyr::filter(str_detect(name, "Time|FSC|SSC"))

  # Writing a New FlowFrame
  NewColStart <- ncol(StashedDF) #Because Backup Still Present
  AllCols <- ncol(TheData)

  cols <- as.matrix(TheData) # Matrix Form Data
  ncol <- ncol(cols) # Number Columns Data
  cn <- colnames(cols) # Column Names Data

  new_pid <- 1
  new_pid <- seq(new_pid, length.out = ncol)
  new_pid <- paste0("$P", new_pid)

  #Updated $P1N to Time and Scatter
  #rownames(ParamData) <- new_pid[1:ncol(StashedDF)-1]

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

  #new_pd #Parameters (with variance to range measurements)

  new_kw <- fr@description

  NameParams <- new_kw[grepl("^\\$P\\d+N\\d*", names(new_kw))]
  VoltageParams <- new_kw[grepl("^\\$P\\d+V\\d*", names(new_kw))]
  VoltageParams <- c(NA, VoltageParams)
  DisplayParams <- new_kw[grepl("^\\P\\d+DISPLAY\\d*", names(new_kw))]
  TypeParams <- new_kw[grepl("^\\$P\\d+TYPE\\d*", names(new_kw))]

  DescriptionData <- cbind(NameParams, VoltageParams, DisplayParams, TypeParams)
  DescriptionData <- as.data.frame(DescriptionData)
  DescriptionData <- DescriptionData %>%
    dplyr::filter(str_detect(NameParams, "Time|FSC|SSC|B1-A"))
  DescriptionData <- DescriptionData %>% unnest(cols = where(is.list))
  #DescriptionData

  Test <- new_kw[!grepl("^\\$P\\d+", names(new_kw))]
  Test <- Test[!grepl("^\\P\\d+", names(Test))]
  Test <- Test[!grepl("^\\$FLOWRATE", names(Test))]
  Test <- Test[!grepl("^\\$CYTOLIB", names(Test))]

  OGLength <- length(Test)

  #i <- new_pid[1]
  #new_pd[[i, 1]]

  # Where the Sausage Gets Made
  for (i in new_pid){
    NoDollar <- gsub("$", "", fixed=TRUE, i)
    Test[paste0(i,"B")] <- "32"              #Bits?
    Test[paste0(i,"E")] <- "0,0"             #Zero
    Test[paste0(i,"N")] <- new_pd[[i,1]]     #Name

    TheName <- new_pd[[i, 1]]
    if (!str_detect(TheName, "Time|FSC|SSC")) {Test[paste0(i,"V")] <- "0"
    } else {
      if (str_detect(TheName, "FSC|SSC")){
        Voltage <- DescriptionData %>% dplyr::filter(NameParams %in% TheName) %>% pull(VoltageParams)
        Test[paste0(i,"V")] <- Voltage
      }
    }

    Test[paste0(i,"R")] <- new_pd[[i,5]]

    if (!str_detect(TheName, "FSC|SSC")) {Test[paste0(NoDollar,"DISPLAY")] <- "LOG"
    } else {Test[paste0(NoDollar,"DISPLAY")] <- "LIN"}

    if (str_detect(TheName, "Time")) {Test[paste0(i,"TYPE")] <- "Time"
    } else if (str_detect(TheName, "FSC")){Test[paste0(i,"TYPE")] <- "Forward_Scatter"
    } else if (str_detect(TheName, "SSC")){Test[paste0(i,"TYPE")] <- "Side_Scatter"
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
  new_kw$`CREATOR` <- "Luciernaga 0.0.1"

  # Adding back to extract again?
  fr@exprs <- cols
  pData(parameters(fr)) <- new_pd

  UpdatedParameters <- parameters(fr)

  new_fcs <- new("flowFrame", exprs=cols, parameters=UpdatedParameters,
                 description=new_kw)

  if (!is.null(addon)){name <- paste0(name, addon)
  }

  AssembledName <- paste0(name, ".fcs")

  if (is.null(outpath)) {outpath <- getwd()}

  fileSpot <- file.path(outpath, AssembledName)

  if (export == TRUE) {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
  } else {return(new_fcs)}
}
