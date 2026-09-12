#' Takes AdditionalPageHandler outputs and consolidates them
#' 
#' @param x The list of data.frames from AdditionalPageHandler
#' @param Metadata The Metadata data.frame from QC_Chorus
#' @param FirstPage The First page data.frame from QC_Chorus
#' @param returnPreference Whether to return Imaging or HighSpeed data
#' 
#' @importFrom dplyr bind_rows
#' @importFrom dplyr rename
#' @importFrom tidyr pivot_wider
#' 
#' @noRd
Consolidator <- function(x, Metadata, FirstPage, returnPreference){
    #Metadata
    MainColumns <- ncol(FirstPage)
    MainColumnsNames <-colnames(FirstPage)
    TheList <- x

    FlatList <- unlist(lapply(TheList, function(x) {
        if (is.data.frame(x)){list(x)
            } else if (is.list(x)) {x
            } else {NULL}}), recursive = FALSE)

    Ncols <- sapply(FlatList, ncol)
    RunID <- cumsum(c(TRUE, diff(Ncols) != 0))
    Grouped <- split(FlatList, RunID)
    CombinedList <- lapply(Grouped, function(grp) do.call(rbind, grp))

    FirstData <- CombinedList[[1]]
    colnames(FirstData) <- MainColumnsNames
    MainData <- bind_rows(FirstPage, FirstData)

    SecondData <- CombinedList[[2]]
    colnames(SecondData) <- as.character(unlist(SecondData[1, ]))
    SecondData <- SecondData[-1, , drop = FALSE]

    if (length(CombinedList) == 4){
        ThirdData <- CombinedList[[3]]
        colnames(ThirdData) <- as.character(unlist(ThirdData[1, ]))
        ThirdData <- ThirdData[-1, , drop = FALSE]

        FourthData <- CombinedList[[2]]
        colnames(FourthData) <- as.character(unlist(FourthData[1, ]))
        FourthData <- FourthData[-1, , drop = FALSE]
    }
    
    if (length(CombinedList) == 4 && returnPreference != "Imaging"){
    MainData <- ThirdData

    FourthDataWide <- FourthData |>
        rename(LaserDelay = `Laser Delay`, 
        LaserPowerWithinSpec = `Laser Power within Spec`) |>
        pivot_wider(names_from = Laser,
             values_from = c(LaserDelay, LaserPowerWithinSpec),
             names_glue = "{Laser}_{.value}")

    MetaFull <- cbind(Metadata, FourthDataWide)
    MainData <- MainDataWiden(MainData)
    MetaExpanded <- MetaFull
    #MetaExpanded <- MetaFull[rep(1, nrow(MainData)), , drop = FALSE]
    FinalData <- cbind(MetaExpanded, MainData)
    } else if (returnPreference== "Imaging"){

    SecondDataWide <- SecondData |>
        rename(LaserDelay = `Laser Delay`, 
        LaserPowerWithinSpec = `Laser Power within Spec`) |>
        pivot_wider(names_from = Laser,
             values_from = c(LaserDelay, LaserPowerWithinSpec),
             names_glue = "{Laser}_{.value}")
    
    MainData <- MainDataWiden(MainData)
    MetaExpanded <- MetaFull

    MetaFull <- cbind(Metadata, SecondDataWide)
    #MetaExpanded <- MetaFull[rep(1, nrow(MainData)), , drop = FALSE]
    FinalData <- cbind(MetaExpanded, MainData)
    }



   return(FinalData)
}