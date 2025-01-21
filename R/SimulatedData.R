#' Produces a simulated .fcs file according to your specificiations. 
#' 
#' @param populations A data.frame containing
#' @param abundance A data.frame containing
#' @param totalevents The number of desired total events
#' @param targets A data.frame containing Fluorophores and MFI
#' @param distribution A data.frame denoting 1, 2, 3 for respective pops
#' @param NumberDetectors Aurora number of detectors, used for reference signatures
#' 
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' 
#' @return Currently a matrix, flowframe implementation in progress
#' 
#' @noRd
SimulatedData <- function(populations, abundance, totalevents, 
  targets, distribution, NumberDetectors){

TheFluorophores <- targets %>% pull(Fluorophores)
ReferenceData <- Luciernaga:::InstrumentReferences(NumberDetectors=NumberDetectors)
Data <- ReferenceData %>% select(-Instrument) %>% filter(Fluorophore %in% TheFluorophores)
ScaledData <- map(.x=TheFluorophores, .f=FluorScaling, data=Data, targets=targets) %>% bind_rows()

IntAbund <- abundance %>% mutate(Total=Ratio*totalevents) %>% select(-Ratio)
ToAssemble <- left_join(IntAbund, populations, by="Pops")

Pops <- ToAssemble %>% pull(Pops)

ThePopulations <- map(.x=Pops, .f=DataSimulation,
 ToAssemble=ToAssemble, ScaledData=ScaledData, distribution=distribution)

#lapply(ThePopulations, dim)
#str(ThePopulations)
Dataset <- do.call(rbind, ThePopulations)
#ncol(Dataset)
#nrow(Dataset)
  
return(Dataset)
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