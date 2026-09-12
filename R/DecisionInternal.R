
#' Internal for Wetlab_Decisions
#'
#' @param x An iterated specimen
#' @param data The passed data.frame
#' @param FinalConcentration_MillionperML Desired outcome
#' @param MillionCellsPerTube Desired Outcome
#' @param TheConditions Passed csv data
#' @param ReturnLeftover Passed Argument
#'
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom purrr map2
#'
#' @return An internal value
#'
#' @noRd
DecisionInternal <- function(x, data, FinalConcentration_MillionperML,
                             MillionCellsPerTube, TheConditions, ReturnLeftover){

  Conditions <- nrow(TheConditions)
  ConditionsList <- TheConditions %>% pull(Conditions)

  Interest <- data %>% dplyr::filter(name %in% x)

  name <- Interest %>% pull(name)
  Date <- Interest %>% pull(Date)
  TotalCells <- Interest %>% pull(TotalCells) %>% as.double(.)
  RestConcentration <- Interest %>% pull(DesiredConcentration) %>% as.double(.)
  FinalConcentration <- FinalConcentration_MillionperML*1000000
  CellsPerTube <- MillionCellsPerTube*1000000

  if(RestConcentration < FinalConcentration){
    warning("Currently Wetlab_Decisions is only coded going from high to low concentration")
  }

  # Checking Available Cells To Work With
  NeededCells <- CellsPerTube*Conditions
  Result <- TotalCells/NeededCells

  if (Result >= 1){
    TheTubes <- map(.x=ConditionsList, .f=ConditionWrap, name=name, Date=Date,
                    RestConcentration=RestConcentration,  FinalConcentration=FinalConcentration,
                    TheConditions=TheConditions, CellsPerTube=CellsPerTube,
                    TotalCells=TotalCells) %>% bind_rows()
    Insufficient <- FALSE

  } else if (Result < 1){
    #TotalCells
    Tentative <- TheConditions %>%
      mutate(MinCells = LowerLimit_MillionCells*CellRatio_Retained*1000000) %>%
      mutate(Theoretical = Partition*TotalCells) %>%
      mutate(Enough = case_when(Theoretical > MinCells ~ TRUE,  TRUE ~ FALSE))

    if(any(Tentative$Enough == FALSE)){
      NewTentative <- Tentative %>% arrange(desc(RankImportance)) %>% slice(-1)
      Residual <- sum(NewTentative$Partition)
      NewTentative <- NewTentative %>% mutate(Partition = Partition/Residual)
      NewTentative <- NewTentative %>% mutate(Theoretical = Partition*TotalCells) %>%
        mutate(Enough = case_when(Theoretical > MinCells ~ TRUE,  TRUE ~ FALSE))

      if(any(NewTentative$Enough == FALSE)){
        NewestTentative <- NewTentative %>% arrange(desc(RankImportance)) %>% slice(-1)
        Residual <- sum(NewestTentative$Partition)
        NewestTentative <- NewestTentative %>% mutate(Partition = Partition/Residual)
        NewestTentative <- NewestTentative %>% mutate(Theoretical = Partition*TotalCells) %>%
          mutate(Enough = case_when(Theoretical > MinCells ~ TRUE,  TRUE ~ FALSE))

        if(any(NewestTentative$Enough == FALSE)){
          Insufficient <- TRUE
          message("Insufficient cells for ", x)
        } else {#message("Proceed to sliding scale")
          SlidingScale <- NewestTentative
          Insufficient <- FALSE
        }
      } else{#message("Proceed to sliding Scale")
        SlidingScale <- NewTentative
        Insufficient <- FALSE
      }
    } else {#message("Proceed to sliding scale")
      SlidingScale <- Tentative
      Insufficient <- FALSE
    }

    if (Insufficient == FALSE){
      Sliding <- SlidingScale %>% select(-c("LowerLimit_MillionCells", "CellRatio_Retained",
                                            "MinCells", "Enough"))

      SliderConditions <- Sliding %>% pull(Conditions)
      SliderCellsPerTube <- Sliding %>% pull(Theoretical)
      SliderCellsPerTube[SliderCellsPerTube > CellsPerTube] <- CellsPerTube
      NeededCells <- sum(SliderCellsPerTube)

      TheTubes <- map2(.x=SliderConditions, .y=SliderCellsPerTube, .f=SliderConditionWrap,
                       name=name, Date=Date, RestConcentration=RestConcentration,
                       FinalConcentration=FinalConcentration,TheConditions=TheConditions,
                       TotalCells=TotalCells) %>% bind_rows()
    } else if (Insufficient == TRUE){
      message("Insufficient cells for ", x)
      Condition <- "Insufficient"
      TheTubes <- TheInsufficient(x=Condition, name=name, Date=Date,
                                  TotalCells=TotalCells, RestConcentration=RestConcentration)
    }
  }

  if (ReturnLeftover == TRUE && Insufficient == FALSE){
    LeftoverCells <- TotalCells-NeededCells
    Frame <- TheTubes %>% slice(1)
    Frame$Condition <- "Leftover"

    TotalCells <- format(LeftoverCells, scientific=TRUE, digits=2)
    Frame$TotalCells <- TotalCells

    CellsInTube <- Frame %>% pull(CellsPerTube) %>% as.double()
    TotalTubes <- round(LeftoverCells/CellsInTube, 1)
    TotalTubes <- data.frame(TotalTubes, check.names=FALSE)
    Frame <- cbind(Frame, TotalTubes)

    TheTubes <- TheTubes %>% mutate(TotalTubes=1)
    TheAssembly <- bind_rows(TheTubes, Frame)
  } else {TheAssembly <- TheTubes}

  return(TheAssembly)

}