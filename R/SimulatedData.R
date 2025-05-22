#' Produces a simulated .fcs file according to your specificiations.
#'
#' @param populations A data.frame containing
#' @param abundance A data.frame containing
#' @param totalevents The number of desired total events
#' @param targets A data.frame containing Fluorophores and MFI
#' @param distribution A data.frame denoting 1, 2, 3 for respective pops
#' @param NumberDetectors Aurora number of detectors, used for reference signatures
#' @param name Desired file name
#' @param addon Desired add on for filename
#' @param outpath Desired storage location for fcs file
#' @param flowWorkspace load_cytoset_from_fcs
#' @param returntype Default fcs, anything else returns a flowframe to Renviron
#'
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr relocate
#' @importFrom flowCore write.FCS
#'
#' @return An .fcs or a flow frame containing the simulated data.
#'
#' @noRd
SimulatedData <- function(populations, abundance, totalevents,
  targets, distribution, NumberDetectors, name, addon, outpath, returntype="fcs"){

TheFluorophores <- targets %>% pull(Fluorophores)
ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
Data <- ReferenceData %>% select(-Instrument) %>% filter(Fluorophore %in% TheFluorophores)
ScaledData <- map(.x=TheFluorophores, .f=FluorScaling, data=Data, targets=targets) %>% bind_rows()

IntAbund <- abundance %>% mutate(Total=Ratio*totalevents) %>% select(-Ratio)
ToAssemble <- left_join(IntAbund, populations, by="Pops")

Pops <- ToAssemble %>% pull(Pops)

ThePopulations <- map(.x=Pops, .f=DataSimulation,
 ToAssemble=ToAssemble, ScaledData=ScaledData, distribution=distribution)

Dataset <- do.call(rbind, ThePopulations)

Dataset <- data.frame(Dataset)
Dataset <- Dataset %>% mutate(DateTime="Standin", Fluorophore="Standin") %>%
  relocate(DateTime, Fluorophore, .before=1)
Dataset <- Luciernaga:::ColumnNaming(Dataset)
Dataset <- Dataset %>% select(-Fluorophore, -DateTime)
colnames(Dataset) <- paste0(colnames(Dataset), "-A")

Noise <- HouseParty(x=Dataset)
NoisyDataset <- Dataset+Noise

NoisyDataset <- ScatterParty(x=NoisyDataset)

ParameterParam <- ParameterPrep(x=NoisyDataset)
DescriptionParam <- DescriptionGenesis(x=NoisyDataset)

NoisyMatrix <- as.matrix(NoisyDataset)

path <- system.file("extdata", package = "Luciernaga")
files <- list.files(path=path, pattern="CD4_BUV805.*Cells", full.names=TRUE)
CytoSet <- load_cytoset_from_fcs(files, truncate_max_range = FALSE, transformation = FALSE)
fr <- CytoSet[[1, returnType = "flowFrame"]]
Parameter <- fr@parameters
Parameter@data <- ParameterParam

new_fcs <- new("flowFrame", exprs=NoisyMatrix, parameters=Parameter,
               description=DescriptionParam)

if (!is.null(addon)){name <- paste0(name, addon)}

AssembledName <- paste0(name, ".fcs")

new_fcs@description$GUID <- AssembledName
new_fcs@description$`$FIL` <- AssembledName

if (is.null(outpath)) {outpath <- getwd()}

fileSpot <- file.path(outpath, AssembledName)

if (returntype == "fcs") {write.FCS(new_fcs, filename = fileSpot, delimiter="#")
} else {return(new_fcs)}

}


#' Internal for Simulated Data, returns vector of noise
#'
#' @param x The data.frame to get dimensions from
#'
#' @importFrom stats rnorm
#'
#' @return A vector of noise values
#' @noRd
NoiseGenerator <- function(x){
  rnorm(nrow(x), mean=0, sd=6)
}

#' Internal for Simulated Data, generates noise data.frame
#'
#' @param x The data.frame we wish to generate noise for
#'
#' @importFrom dplyr bind_cols
#'
#' @return An equivalent data.frame of noise values
#' @noRd
HouseParty <- function(x){
  TheNames <- colnames(x)
  TheParticipants <- list()
  for (i in seq_along(x)){
    TheParticipants[[i]] <- NoiseGenerator(x)
  }
  names(TheParticipants) <- TheNames
  TheParticipants <- bind_cols(TheParticipants)
  return(TheParticipants)
}

#' Internal for Simulated Data, generates scatter data and adds to data.frame
#'
#' @param x The existing data.frame of just detectors
#'
#' @importFrom stats rnorm
#' @importFrom stringr str_detect
#' @importFrom dplyr relocate
#'
#' @return The updated data.frame now including time and scatter parameters.
#' @noRd
ScatterParty <- function(x){

  values <- 0:500000
  Time <- sample(values, size=nrow(x), replace=FALSE)

  SSCW <- data.frame("SSC-W" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)
  SSCH <- data.frame("SSC-H" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)
  SSCA <- data.frame("SSC-A" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)

  FSCW <- data.frame("FSC-W" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)
  FSCH <- data.frame("FSC-H" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)
  FSCA <- data.frame("FSC-A" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)

  SSCBW <- data.frame("SSC-B-W" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)
  SSCBH <- data.frame("SSC-B-H" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)
  SSCBA <- data.frame("SSC-B-A" = rnorm(nrow(x), mean=1000000, sd=125000), check.names=FALSE)

  TheDetectors <- colnames(x)

  Data <- cbind(Time, x)

  if (any(str_detect(TheDetectors, "UV"))){
    Data <- cbind(Data, SSCW, SSCH, SSCA) %>%
      relocate("SSC-W", "SSC-H", "SSC-A", .before="V1-A")
  } else {Data <- cbind(Data, SSCW, SSCH, SSCA) %>%
    relocate("SSC-W", "SSC-H", "SSC-A", .after="Time")
  }

  Data <- cbind(Data, FSCW, FSCH, FSCA, SSCBW, SSCBH, SSCBA) %>%
    relocate("FSC-W", "FSC-H", "FSC-A", "SSC-B-W", "SSC-B-H", "SSC-B-A", .before="B1-A")

  return(Data)
}

#' Internal for SimulatedData, adjust reference signatures to desired MFI targets
#'
#' @param x Iterated in Fluorophore
#' @param data Iterated in Fluorophore Signature Data
#' @param targets The data.frame containing Fluorophore MFI targets
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#'
#' @return MFI-scaled Fluorophore Signature
#'
#' @noRd
FluorScaling <- function(x, data, targets){
  IntData <- data %>% filter(Fluorophore %in% x)
  Multiplier <- targets %>% filter(Fluorophores %in% x) %>% pull(MFI)
  IntData <- IntData %>% mutate(AdjustedY = AdjustedY*Multiplier)
  return(IntData)
}

#' Internal for SimulatedData, generates total variant events and multiplies by signature
#'
#' @param x Iterated in Fluorophore
#' @param LocalNumber Iterated in number of events (derrived from pop and abundance)
#' @param ScaledData The Scaled signature data
#' @param distribution The desired distribution for the population (1, 2, 3)
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom stats rlnorm
#'
#' @return A matrix of detector data specific for that fluorophore by respective events
#'
#' @noRd
FluorophoreVariance <- function(x, LocalNumber, ScaledData, distribution){
LocalData <- ScaledData %>% filter(Fluorophore %in% x) %>% pull(AdjustedY)
LocalDist <- distribution %>% filter(Markers %in% x) %>% pull(Distribution)

if (LocalDist == 1){
mu <- 0 # Mean
sigma <- 0.05 #SD
Values <- rlnorm(LocalNumber, meanlog = mu, sdlog=sigma)
#mean(Values)
#hist(Values, main = "Log-normal Distribution", xlab = "Value", breaks = 50)
#LocalData

ResultMatrix <- sapply(Values, function(value) LocalData * value)
ResultMatrix <- t(ResultMatrix)

#Results <- data.frame(ResultMatrix)
#A <- do.call(pmax, Results)
#Normalized <- Results/A
#View(Normalized)

return(ResultMatrix)
}
}


#' Internal for Simulated Data, coordinates pop level and adds up across fluorophore contributions
#'
#' @param x Iterated Population
#' @param ToAssemble Combined Population data
#' @param ScaledData The scaled signatures
#' @param distribution The data.frame with desired distribution parameters.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom purrr reduce
#'
#' @return The matrices of raw values for respective population
#'
#' @noRd
DataSimulation <- function(x, ToAssemble, ScaledData, distribution){
  IntData <- ToAssemble %>% filter(Pops %in% x) %>% select(-Pops)
  LocalNumber <- IntData %>% pull(Total)
  IntData <- IntData %>% select(-Total)
  NotZero <- IntData[, apply(IntData, 2, function(col) any(col != 0)), drop = FALSE]
  TheseFluorophores <- NotZero %>% colnames(.)

  if (length(TheseFluorophores) >= 1){
  #x <- TheseFluorophores[1]
  FluorMatrices <- map(.x=TheseFluorophores, .f=FluorophoreVariance,
   LocalNumber=LocalNumber, ScaledData=ScaledData, distribution=distribution)

  if (length(FluorMatrices) > 1){
    SummedContents <- reduce(FluorMatrices, `+`)
  } else {SummedContents <- FluorMatrices[[1]]}

  return(SummedContents)
  } else {message("Institute blank data protocol here")}
}
