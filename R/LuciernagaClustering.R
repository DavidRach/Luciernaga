#' Internal for UnstainedSignatures And SingleStainSignatures
#'
#' @param MyData The data.frame of raw .fcs exprs values
#' @param These A list of columns corresponding to local maxima peaks
#' @param DetectorName The brightest peak everything is normalized to.
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#'
#' @keywords internal
LuciernagaClustering <- function(MyData, These, DetectorName){
  MyData$Cluster <- paste(DetectorName, "10-", sep = "_")

  if (length(These) > 7){stop(
    "Only currently set up to handle up to 8 fluorescence peaks per fluorophore")
  } else if (length(These) == 7){second <- These[[1]]
                                 third <- These[[2]]
                                 fourth <- These[[3]]
                                 fifth <- These[[4]]
                                 sixth <- These[[5]]
                                 seventh <- These[[6]]
                                 eight <- These[[7]]
  } else if (length(These) == 6){second <- These[[1]]
                                 third <- These[[2]]
                                 fourth <- These[[3]]
                                 fifth <- These[[4]]
                                 sixth <- These[[5]]
                                 seventh <- These[[6]]
  } else if (length(These) == 5){second <- These[[1]]
                                 third <- These[[2]]
                                 fourth <- These[[3]]
                                 fifth <- These[[4]]
                                 sixth <- These[[5]]
  } else if (length(These) == 4){second <- These[[1]]
                                third <- These[[2]]
                                fourth <- These[[3]]
                                fifth <- These[[4]]
  } else if (length(These) == 3){second <- These[[1]]
                                 third <- These[[2]]
                                 fourth <- These[[3]]
  } else if (length(These) == 2){second <- These[[1]]
                                 third <- These[[2]]
  } else if (length(These) == 1){second <- These[[1]]
  } else if (length(These) == 0){
    if (Verbose == TRUE) {message("No second peak")}}

  if (length(These) >= 1){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[second]], 0.0) ~ paste0(MyData$Cluster, second, "_00-"),
    near(MyData[[second]], 0.1) ~ paste0(MyData$Cluster, second, "_01-"),
    near(MyData[[second]], 0.2) ~ paste0(MyData$Cluster, second, "_02-"),
    near(MyData[[second]], 0.3) ~ paste0(MyData$Cluster, second, "_03-"),
    near(MyData[[second]], 0.4) ~ paste0(MyData$Cluster, second, "_04-"),
    near(MyData[[second]], 0.5) ~ paste0(MyData$Cluster, second, "_05-"),
    near(MyData[[second]], 0.6) ~ paste0(MyData$Cluster, second, "_06-"),
    near(MyData[[second]], 0.7) ~ paste0(MyData$Cluster, second, "_07-"),
    near(MyData[[second]], 0.8) ~ paste0(MyData$Cluster, second, "_08-"),
    near(MyData[[second]], 0.9) ~ paste0(MyData$Cluster, second, "_09-"),
    near(MyData[[second]], 1.0) ~ paste0(MyData$Cluster, second, "_10-")))
  }

  if(length(These) >= 2){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[third]], 0.0) ~ paste0(MyData$Cluster, third, "_00"),
    near(MyData[[third]], 0.1) ~ paste0(MyData$Cluster, third, "_01"),
    near(MyData[[third]], 0.2) ~ paste0(MyData$Cluster, third, "_02"),
    near(MyData[[third]], 0.3) ~ paste0(MyData$Cluster, third, "_03"),
    near(MyData[[third]], 0.4) ~ paste0(MyData$Cluster, third, "_04"),
    near(MyData[[third]], 0.5) ~ paste0(MyData$Cluster, third, "_05"),
    near(MyData[[third]], 0.6) ~ paste0(MyData$Cluster, third, "_06"),
    near(MyData[[third]], 0.7) ~ paste0(MyData$Cluster, third, "_07"),
    near(MyData[[third]], 0.8) ~ paste0(MyData$Cluster, third, "_08"),
    near(MyData[[third]], 0.9) ~ paste0(MyData$Cluster, third, "_09"),
    near(MyData[[third]], 1.0) ~ paste0(MyData$Cluster, third, "_10")))
  }

  if(length(These) >= 3){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[fourth]], 0.0) ~ paste0(MyData$Cluster, fourth, "_00"),
    near(MyData[[fourth]], 0.1) ~ paste0(MyData$Cluster, fourth, "_01"),
    near(MyData[[fourth]], 0.2) ~ paste0(MyData$Cluster, fourth, "_02"),
    near(MyData[[fourth]], 0.3) ~ paste0(MyData$Cluster, fourth, "_03"),
    near(MyData[[fourth]], 0.4) ~ paste0(MyData$Cluster, fourth, "_04"),
    near(MyData[[fourth]], 0.5) ~ paste0(MyData$Cluster, fourth, "_05"),
    near(MyData[[fourth]], 0.6) ~ paste0(MyData$Cluster, fourth, "_06"),
    near(MyData[[fourth]], 0.7) ~ paste0(MyData$Cluster, fourth, "_07"),
    near(MyData[[fourth]], 0.8) ~ paste0(MyData$Cluster, fourth, "_08"),
    near(MyData[[fourth]], 0.9) ~ paste0(MyData$Cluster, fourth, "_09"),
    near(MyData[[fourth]], 1.0) ~ paste0(MyData$Cluster, fourth, "_10")))
  }

  if(length(These) >= 4){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[fifth]], 0.0) ~ paste0(MyData$Cluster, fifth, "_00"),
    near(MyData[[fifth]], 0.1) ~ paste0(MyData$Cluster, fifth, "_01"),
    near(MyData[[fifth]], 0.2) ~ paste0(MyData$Cluster, fifth, "_02"),
    near(MyData[[fifth]], 0.3) ~ paste0(MyData$Cluster, fifth, "_03"),
    near(MyData[[fifth]], 0.4) ~ paste0(MyData$Cluster, fifth, "_04"),
    near(MyData[[fifth]], 0.5) ~ paste0(MyData$Cluster, fifth, "_05"),
    near(MyData[[fifth]], 0.6) ~ paste0(MyData$Cluster, fifth, "_06"),
    near(MyData[[fifth]], 0.7) ~ paste0(MyData$Cluster, fifth, "_07"),
    near(MyData[[fifth]], 0.8) ~ paste0(MyData$Cluster, fifth, "_08"),
    near(MyData[[fifth]], 0.9) ~ paste0(MyData$Cluster, fifth, "_09"),
    near(MyData[[fifth]], 1.0) ~ paste0(MyData$Cluster, fifth, "_10")))
  }

  if(length(These) >= 5){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[sixth]], 0.0) ~ paste0(MyData$Cluster, sixth, "_00"),
    near(MyData[[sixth]], 0.1) ~ paste0(MyData$Cluster, sixth, "_01"),
    near(MyData[[sixth]], 0.2) ~ paste0(MyData$Cluster, sixth, "_02"),
    near(MyData[[sixth]], 0.3) ~ paste0(MyData$Cluster, sixth, "_03"),
    near(MyData[[sixth]], 0.4) ~ paste0(MyData$Cluster, sixth, "_04"),
    near(MyData[[sixth]], 0.5) ~ paste0(MyData$Cluster, sixth, "_05"),
    near(MyData[[sixth]], 0.6) ~ paste0(MyData$Cluster, sixth, "_06"),
    near(MyData[[sixth]], 0.7) ~ paste0(MyData$Cluster, sixth, "_07"),
    near(MyData[[sixth]], 0.8) ~ paste0(MyData$Cluster, sixth, "_08"),
    near(MyData[[sixth]], 0.9) ~ paste0(MyData$Cluster, sixth, "_09"),
    near(MyData[[sixth]], 1.0) ~ paste0(MyData$Cluster, sixth, "_10")))
  }

  if(length(These) >= 6){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[seventh]], 0.0) ~ paste0(MyData$Cluster, seventh, "_00"),
    near(MyData[[seventh]], 0.1) ~ paste0(MyData$Cluster, seventh, "_01"),
    near(MyData[[seventh]], 0.2) ~ paste0(MyData$Cluster, seventh, "_02"),
    near(MyData[[seventh]], 0.3) ~ paste0(MyData$Cluster, seventh, "_03"),
    near(MyData[[seventh]], 0.4) ~ paste0(MyData$Cluster, seventh, "_04"),
    near(MyData[[seventh]], 0.5) ~ paste0(MyData$Cluster, seventh, "_05"),
    near(MyData[[seventh]], 0.6) ~ paste0(MyData$Cluster, seventh, "_06"),
    near(MyData[[seventh]], 0.7) ~ paste0(MyData$Cluster, seventh, "_07"),
    near(MyData[[seventh]], 0.8) ~ paste0(MyData$Cluster, seventh, "_08"),
    near(MyData[[seventh]], 0.9) ~ paste0(MyData$Cluster, seventh, "_09"),
    near(MyData[[seventh]], 1.0) ~ paste0(MyData$Cluster, seventh, "_10")))
  }

  if(length(These) >= 7){MyData <- MyData %>% mutate(Cluster = case_when(
    near(MyData[[eight]], 0.0) ~ paste0(MyData$Cluster, eight, "_00"),
    near(MyData[[eight]], 0.1) ~ paste0(MyData$Cluster, eight, "_01"),
    near(MyData[[eight]], 0.2) ~ paste0(MyData$Cluster, eight, "_02"),
    near(MyData[[eight]], 0.3) ~ paste0(MyData$Cluster, eight, "_03"),
    near(MyData[[eight]], 0.4) ~ paste0(MyData$Cluster, eight, "_04"),
    near(MyData[[eight]], 0.5) ~ paste0(MyData$Cluster, eight, "_05"),
    near(MyData[[eight]], 0.6) ~ paste0(MyData$Cluster, eight, "_06"),
    near(MyData[[eight]], 0.7) ~ paste0(MyData$Cluster, eight, "_07"),
    near(MyData[[eight]], 0.8) ~ paste0(MyData$Cluster, eight, "_08"),
    near(MyData[[eight]], 0.9) ~ paste0(MyData$Cluster, eight, "_09"),
    near(MyData[[eight]], 1.0) ~ paste0(MyData$Cluster, eight, "_10")))
  }

  return(MyData)
}
