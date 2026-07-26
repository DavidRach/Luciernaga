#' Carries out the initial unmixing using just the single-color controls own reference signature. 
#' 
#' @param matrix The signature matrix for the respective single-color unmixing controls
#' @param matrix_columnname The column name used to match respective metadata in the GatingSet. 
#' Default is Fluorophore.
#' @param gs The GatingSet containing the raw single-color unmixing controls
#' @param subset The gate with the events you want to unmix. Default is root. 
#' @param outpath The location to store the unmixed .fcs files. 
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#' 
#' @importFrom flowWorkspace pData
#' @importFrom dplyr select pull mutate row_number
#' @importFrom tidyselect all_of
#' @importFrom purrr walk
#' 
#' @return Unmixed single-color unmixing controls for use in stain reduction calculations. 
#' 
#' @export
#' 
#' @examples A <- 2 + 2
#' 
SingleUnmix <- function(matrix, matrix_columnname="Fluorophore", gs, subset="root",
 outpath, inverse.transform=TRUE, sample.name="TUBENAME"){

    # Metadata Fluors
    pd <- pData(gs)
    Location <- pd |> select(all_of(matrix_columnname)) |> mutate(Row=row_number())
    TheseFluorophores <- Location |> pull(matrix_columnname)

    # Corrected Matrix Order
    RetainedMatrix <- matrix |> filter(Fluorophore %in% TheseFluorophores)
    new_order <- match(TheseFluorophores, RetainedMatrix$Fluorophore)
    RetainedMatrix <- RetainedMatrix[new_order,]

    FolderName <- file.path(outpath, "SingleUnmix")
    if(!dir.exists(FolderName)){dir.create(FolderName)}

    walk(.x=TheseFluorophores, .f=SingleUnmixIterator, gs=gs, matrix=RetainedMatrix, outpath=FolderName,
    subset=subset, inverse.transform=inverse.transform, sample.name=sample.name)

    FolderName2 <- file.path(outpath, "AllUnmix")
    if(!dir.exists(FolderName2)){dir.create(FolderName2)}  

    walk(.x=TheseFluorophores, .f=AllUnmixIterator, gs=gs, matrix=RetainedMatrix, outpath=FolderName2,
    subset=subset, inverse.transform=inverse.transform, sample.name=sample.name)

    FolderName3 <- file.path(outpath, "FMO")
    if(!dir.exists(FolderName3)){dir.create(FolderName3)}  

    walk(.x=TheseFluorophores, .f=FMOUnmixIterator, gs=gs, matrix=RetainedMatrix, outpath=FolderName3,
    subset=subset, inverse.transform=inverse.transform, sample.name=sample.name)
    
    message("Done!")
}


#' Internal for SingleUnmix, unmixes with all reference signatures
#' 
#' @param x The iterated in fluorophore name, matching metadata in gs and fluorophore in matrix
#' @param matrix The signature matrix for the respective single-color unmixing controls
#' @param gs The GatingSet containing the raw single-color unmixing controls
#' @param subset The gate with the events you want to unmix. Default is root. 
#' @param outpath The location to store the unmixed .fcs files. 
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#' 
#' @importFrom BiocGenerics subset
#' @importFrom dplyr filter select where 
#' @importFrom purrr walk
#' 
#' @return Unmixed single-color unmixing controls for use in stain reduction calculations. 
#' 
#' @noRd
AllUnmixIterator <- function(x, gs, matrix, outpath, subset, inverse.transform, sample.name){
    # x <- TheseFluorophores[1]
    InternalGS <- subset(gs, Fluorophore == x)
    Signature <- matrix |> dplyr::select(Fluorophore, Antigen, where(is.numeric))
    Panel <- matrix |> dplyr::select(Fluorophore, Antigen)
    colnames(Signature)[2] <- "Ligand"
    #colnames(Panel)[2] <- "Ligand"

    walk(.x=InternalGS, .f=Luciernaga_Unmix, controlData=Signature, sample.name=sample.name,
     addon="_AllUnmixed", subset=subset, removestrings="fcs", outpath=outpath, PanelPath=Panel,
    Verbose=FALSE, inverse.transform=inverse.transform)

    #devtools::load_all("/home/david/Documents/Luciernaga")
}

#' Internal for SingleUnmix, drops the given single-color signature, unmixes full panel without it. 
#' 
#' @param x The iterated in fluorophore name, matching metadata in gs and fluorophore in matrix
#' @param matrix The signature matrix for the respective single-color unmixing controls
#' @param gs The GatingSet containing the raw single-color unmixing controls
#' @param subset The gate with the events you want to unmix. Default is root. 
#' @param outpath The location to store the unmixed .fcs files. 
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#' 
#' @importFrom BiocGenerics subset
#' @importFrom dplyr filter select where 
#' @importFrom purrr walk
#' 
#' @return Unmixed single-color unmixing controls for use in stain reduction calculations. 
#' 
#' @noRd
FMOUnmixIterator <- function(x, gs, matrix, outpath, subset, inverse.transform, sample.name){
    # x <- TheseFluorophores[1]
    InternalGS <- subset(gs, Fluorophore != x)
    Signature <- matrix |> dplyr::filter(!Fluorophore %in% x) |> dplyr::select(Fluorophore, Antigen, where(is.numeric))
    Panel <- matrix |> dplyr::filter(!Fluorophore %in% x) |> dplyr::select(Fluorophore, Antigen)
    colnames(Signature)[2] <- "Ligand"
    #colnames(Panel)[2] <- "Ligand"

    TheName <- paste0("No", x)
    TheNameUnmixed <- paste0("_", TheName, "Unmixed")

    StashHere <- file.path(outpath, TheName)
    if(!dir.exists(StashHere)){dir.create(StashHere)}  

    walk(.x=InternalGS, .f=Luciernaga_Unmix, controlData=Signature, sample.name=sample.name,
     addon=TheNameUnmixed, subset=subset, removestrings="fcs", outpath=StashHere, PanelPath=Panel,
    Verbose=FALSE, inverse.transform=inverse.transform)

    #devtools::load_all("/home/david/Documents/Luciernaga")
}




#' Internal for SingleUnmix, unmixes with just single reference signature
#' 
#' @param x The iterated in fluorophore name, matching metadata in gs and fluorophore in matrix
#' @param matrix The signature matrix for the respective single-color unmixing controls
#' @param gs The GatingSet containing the raw single-color unmixing controls
#' @param subset The gate with the events you want to unmix. Default is root. 
#' @param outpath The location to store the unmixed .fcs files. 
#' @param inverse.transform Default equals TRUE
#' @param sample.name Keyword(s) to grab from the cytoset for renaming
#' 
#' @importFrom BiocGenerics subset
#' @importFrom dplyr filter select where 
#' @importFrom purrr walk
#' 
#' @return Unmixed single-color unmixing controls for use in stain reduction calculations. 
#' 
#' @noRd
SingleUnmixIterator <- function(x, gs, matrix, outpath, subset, inverse.transform, sample.name){
    # x <- TheseFluorophores[1]
    InternalGS <- subset(gs, Fluorophore == x)
    Signature <- matrix |> dplyr::filter(Fluorophore %in% x) |> dplyr::select(Fluorophore, Antigen, where(is.numeric))
    Panel <- matrix |> dplyr::filter(Fluorophore %in% x) |> dplyr::select(Fluorophore, Antigen)
    colnames(Signature)[2] <- "Ligand"
    #colnames(Panel)[2] <- "Ligand"

    walk(.x=InternalGS, .f=Luciernaga_Unmix, controlData=Signature, sample.name=sample.name,
     addon="_SingleUnmixed", subset=subset, removestrings="fcs", outpath=outpath, PanelPath=Panel,
    Verbose=FALSE, inverse.transform=inverse.transform)

    #devtools::load_all("/home/david/Documents/Luciernaga")
}


#' Takes the folder outputs from , parses the stain index for all the combinations,
#' returning as a large long data.frame
#' 
#' @param folder_location File.path to the parent folder where all the subfolders containing
#' the variant .fcs files were stored
#' @param outpath Where you want to store the processed data
#' @param excludeThese Used to exclude columns from transformation, default is "FSC|SSC|Time"
#' @param channelRange Default for biexponential transformation is 4096
#' @param maxValue Default for biexponential transformation is 4194304
#' @param pos Default for biexponential transformation is 5.62
#' @param neg Default for biexponential transformation is 0
#' @param widthBasis Default for biexponential transformation is -1000
#' @param inverse.transform Default is TRUE
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' 
#' @return A data.frame containing staining index for all the FMO folders fcs files. 
#' 
#' @export
#' 
#' @examples A <- 2+2
#' 
#' 
StainBrightnessIndexCalculator <- function(folder_location, outpath=NULL,
     excludeThese="FSC|SSC|Time", channelRange=4096, maxValue=4194304,
     pos=5.62, neg=0, widthBasis=-1000,
     inverse.transform=TRUE){

     AllUnmix <- list.files(folder_location, full.names=TRUE, pattern="AllUnmix")
     FMO <- list.files(folder_location, full.names=TRUE, pattern="FMO")
     FMO_Folders <- list.files(FMO, full.names=TRUE)
     SingleUnmix <- list.files(folder_location, full.names=TRUE, pattern="SingleUnmix")

     TheData <- purrr::map(.x=FMO_Folders, .f=BrightnessIndexIterator,
          excludeThese=excludeThese, channelRange=channelRange, maxValue=maxValue,
          pos=pos, neg=neg, widthBasis=widthBasis, inverse.transform=inverse.transform)

     TheData <- TheData |> bind_rows()
     return(TheData)
}

#' Internal for StainBrightnessIndexCalculator, creates positive and negative gates
#' using openCyto, before iterating through the .fcs files for stain index. 
#' 
#' @param x The file path to the respective FMO folder
#' @param excludeThese Used to exclude columns from transformation, default is "FSC|SSC|Time"
#' @param channelRange Default for biexponential transformation is 4096
#' @param maxValue Default for biexponential transformation is 4194304
#' @param pos Default for biexponential transformation is 5.62
#' @param neg Default for biexponential transformation is 0
#' @param widthBasis Default for biexponential transformation is -1000
#' @param inverse.transform Default is TRUE
#' 
#' @importFrom flowWorkspace load_cytoset_from_fcs GatingSet flowjo_biexp_trans transformerList transform
#' @importFrom data.table fread
#' @importFrom stringr str_detect
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom flowCore filterList
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#' 
#' @return A data.frame with staining index for the files within the FMO folder
#' 
#' @noRd
#' 
BrightnessIndexIterator <- function(x, excludeThese, channelRange, maxValue,
     pos, neg, widthBasis, inverse.transform){

     files <- list.files(x, pattern=".fcs", full.names=TRUE)
     theCytoset <-load_cytoset_from_fcs(files,
      truncate_max_range = FALSE, transformation = FALSE)
     theGatingSet <- GatingSet(theCytoset)
     SFC_Parameters <- colnames(theGatingSet)
     FluorophoresOnly <- SFC_Parameters[!stringr::str_detect(SFC_Parameters, excludeThese)]
     Biexponential <-  flowjo_biexp_trans(channelRange=channelRange,
      maxValue=maxValue, pos=pos, neg=neg, widthBasis=widthBasis)
     MyBiexTransform <- transformerList(FluorophoresOnly, Biexponential)
     transform(theGatingSet, MyBiexTransform)

     FileLocation <- system.file("extdata", package = "Luciernaga")
     UnmixedGates <- data.table::fread(file.path(path = FileLocation,
                                   pattern = 'GatesUnmixed.csv'))
     
     Example <- UnmixedGates[6]
     Example[1,1] <- FluorophoresOnly[1]
     Example[1,2] <- "+/-"
     Example[1,3] <- "root"
     Example[1,4] <- FluorophoresOnly[1]

     Template <- Example

     for (Fluorophore in FluorophoresOnly[2:length(FluorophoresOnly)]){
          Template1 <- Template
          Template1[1,1] <-Fluorophore
          Template1[1,4] <- Fluorophore
          Example <- rbind(Example, Template1)
     }

     UnmixedGating <- gatingTemplate(Example)
     gt_gating(UnmixedGating, theGatingSet) #flowCore filterList

     Data <- purrr::map(.x=theGatingSet, .f=StainIndexLocal, inverse.transform=inverse.transform)
     Data <- Data |> bind_rows()
     return(Data)
}

#' Internal for BrightnessIndexIterator (also an internal). When handed the gated flowframe,
#' proceeds to extract the correct positive and negative gate values for the respective single-color
#' and process the staining index. 
#' 
#' @param x The individual GatingHierarchy being iterated in
#' @param inverse.transform Whether to reverse the transformation 
#' before calculating the staining index. 
#' 
#' @importFrom flowWorkspace pData gs_pop_get_data
#' @importFrom dplyr pull summarise where across select
#' @importFrom tidyselect all_of
#' @importFrom stats quantile
#' @importFrom stringr str_extract str_match
#' @importFrom flowCore exprs
#' 
StainIndexLocal <- function(x, inverse.transform){

     NameString <- pData(x) |> pull(name)
     TheAbsentFluorophore <- str_extract(NameString, "(?<=_No).*(?=Unmixed\\.fcs)")

     TheFluorophore <- str_match(NameString, "^\\S+ (.*?) \\(Beads\\)")[,2]
     
     PositiveGate <- paste0(TheFluorophore, stringAppend, "+")
     NegativeGate <- paste0(TheFluorophore, stringAppend, "-")

     Positive <- gs_pop_get_data(x, PositiveGate, inverse.transform = inverse.transform)
     Positive <- exprs(Positive[[1]])
     Positive <- data.frame(Positive, check.names=FALSE)
     Positive <- Positive[,-grep(excludeThese, names(Positive))]
     TheMedian <- Positive |>
     summarise(across(where(is.numeric), \(x) quantile(x, probs = 0.95, na.rm = TRUE)))
     KeptMarkers <- colnames(TheMedian)
     Fluorophore <- names(TheMedian)[which.max(TheMedian)]
     Positive <- Positive |> select(all_of(Fluorophore))

     Negative <- gs_pop_get_data(x, NegativeGate, inverse.transform = inverse.transform)
     Negative <- exprs(Negative[[1]])
     Negative <- data.frame(Negative, check.names=FALSE)
     Negative <- Negative[,-grep(excludeThese, names(Negative))]
     Negative <- Negative |> select(all_of(Fluorophore))
  
     PositiveVals <- Positive |> pull(1)
     MFI_Pos <- quantile(PositiveVals, probs=0.5, na.rm =TRUE)
     MFI_Pos <- MFI_Pos[[1]]
     NegativeVals <- Negative |> pull(1)
     MFI_Neg <- quantile(NegativeVals, probs=0.5, na.rm =TRUE)
     MFI_Neg <- MFI_Neg[[1]]
     #hist(NegativeVals)
     MinNeg <- quantile(NegativeVals, probs=0.03, na.rm =TRUE)
     MinNeg <- MinNeg[[1]]
     MaxNeg <- quantile(NegativeVals, probs=0.97, na.rm =TRUE)
     MaxNeg <- MaxNeg[[1]]
     RSD <- (MaxNeg - MinNeg)/3.29
     StainingIndex <- (MFI_Pos - MFI_Neg)/(2*RSD)

     TheData <- data.frame(FMO=TheAbsentFluorophore,
                         Fluorophore=TheFluorophore,
                         StainIndex=StainingIndex)
     
     return(TheData)
}
