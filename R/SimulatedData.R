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