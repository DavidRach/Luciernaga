#' Creates a Description List for .fcs file from scratch.
#'
#' @param x The data.frame of exprs values (including time and scatters)
#'
#' @importFrom purrr map flatten
#' @importFrom stringr str_detect
#' @importFrom dplyr filter pull mutate relocate select case_when
#'
#' @return A TBD product
#' @noRd
DescriptionGenesis <- function(x){
 Parameters <- ParameterPrep(x)

 Description <- list(
   FCSversion = '3',
   '$BEGINANALYSIS' = '0',
   '$BEGINDATA' = '45270',
   '$BEGINSTEXT' = '0',
   '$BTIM' = '04:00:00.00',
   '$BYTEORD' = '4,3,2,1',
   '$CYT' = 'Aurora',
   '$CYTOLIB_VERSION' = '2.16.0',
   '$CYTSN' = 'U010524',
   '$DATATYPE' = 'F',
   '$DATE' = '05-Jan-2024',
   '$ENDANALYSIS' = '0',
   '$ENDDATA' = '3005269',
   '$ENDSTEXT' = '0',
   '$ETIM' = '04:01:00.00',
   '$FIL' = 'Example.fcs',
   '$FLOWRATE' = 'Medium',
   '$INST' = '0',
   '$MODE' = 'L',
   '$NEXTDATA' = '0',
   '$OP' = 'FulanoDeTal'
 )
 #View(Description)

 ParameterAdjust <- Parameters
 ParameterAdjust <- RowNameReorder(ParameterAdjust)
 NewParameters <- rownames(ParameterAdjust)
 TheParamList <- map(.x=NewParameters, .f=ParameterDictate, data=ParameterAdjust)
 TheParamList <- flatten(TheParamList)
 #View(TheParamList)

 TheNames <- Parameters[['name']]
 TheNames <- TheNames[!grepl("Time|FSC|SSC", TheNames, ignore.case = TRUE)]
 Diagonally <- diag(1, nrow = length(TheNames), ncol = length(TheNames))
 colnames(Diagonally) <- TheNames
 Diagonally[2, 1] <- 0.000001

 # SpillOver Matrix Creation

 Description_2 <- list(
   '$PAR'='74',
   '$PROJ'='2025_Example',
   '$SPILLOVER'= Diagonally,
   '$TIMESTEP'='0.0001',
   '$TOT'='10000',
   '$VOL'='30.28',
   'APPLY COMPENSATION'='FALSE',
   'CHARSET'='utf-8',
   'CREATOR'='Luciernaga 0.99.1',
   'FILENAME'='C:\\Users\\FulanoDeTal\\Desktop\\ExampleData.fcs',
   'FSC ASF'='1.04',
   'GROUPNAME'='ExampleGroup',
   'GUID'='ExampleData.fcs'
 )

 #View(Description_2)

 Lasers <- list()
 if(any(str_detect(Parameters$name, "UV"))){Lasers <- c(Lasers, "UV")
 }
 if(any(str_detect(Parameters$name, "^V"))){Lasers <- c(Lasers, "V")
 }
 if(any(str_detect(Parameters$name, "^B"))){Lasers <- c(Lasers, "B")
 }
 if(any(str_detect(Parameters$name, "^YG"))){Lasers <- c(Lasers, "YG")
 }
 if(any(str_detect(Parameters$name, "^R"))){Lasers <- c(Lasers, "R")
 }
 Lasers <- unlist(Lasers)

 LaserData <- data.frame(LaserOrder=c("YG", "V", "B", "R", "UV"),
                         ASF=c('1.11', '1.1', '1.03', '1.13', '1.19'),
                         DELAY=c('-40.875', '-20.65', '0', '19.2', '38.975'),
                         NAME=c("YellowGreen", "Violet",
                                     "Blue", "Red", "UV"))
 LaserData_subset <- LaserData |> filter(LaserOrder %in% Lasers) |>
   mutate(LaserNumber=row_number()) |> relocate(LaserNumber, .before=ASF)
 Lasers <- LaserData_subset |> pull(LaserOrder)
 AllLasers <- map(.x=Lasers, .f=LaserList, data=LaserData_subset)
 AllLasers <- flatten(AllLasers)
 #View(AllLasers)

 DisplaySetup <- ParameterAdjust |> select(name)
 DisplaySetup <- DisplaySetup |> mutate(Display="LOG")
 DisplaySetup <- DisplaySetup |> mutate(
   Display = case_when(name == "Time" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "FSC-A" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "FSC-H" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "FSC-W" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-A" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-H" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-W" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-B-A" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-B-H" ~ 'LIN', TRUE ~ Display),
   Display = case_when(name == "SSC-B-W" ~ 'LIN', TRUE ~ Display))

 TheDisplayNames <- rownames(DisplaySetup)

 DisplayList <- map(.x=TheDisplayNames, .f=DisplayInternal, data=DisplaySetup)
 Display <- flatten(DisplayList)

 #View(Display)

 Description_3 <- list(
   'THRESHOLD'='(FSC,150000)',
   'TUBENAME'='ExampleData',
   'USERSETTINGNAME'='2025_Simulation',
   'WINDOW EXTENSION'='3',
   'ORIGINALGUID'='ExampleData.fcs'
 )
 #View(Description_3)

 ###################
 # List, Assemble! #
 ###################

 TheMegaList <- c(Description, TheParamList, Description_2,
                  AllLasers, Display, Description_3)

 return(TheMegaList)
}