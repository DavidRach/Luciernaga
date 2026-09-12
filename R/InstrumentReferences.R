
#' Internal for QC_ProspectiveAdditions
#'
#' @param NumberDetectors A number corresponding to the number of detectors for
#' the spectral instrument
#'
#' @importFrom utils read.csv
#'
#' @return The stored reference data for that instrument contained within
#' Luciernaga
#' @noRd
InstrumentReferences <- function(NumberDetectors){
  if (NumberDetectors == 184){instrument <- "Sony_7L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_7L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "182_DUV"){instrument <- "Sony_6L_DUV"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_6LDUV_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 147){instrument <- "Sony_5L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_5L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 112){instrument <- "Sony 4L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_4L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 86){instrument <- "Sony 3L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "SonyID700_3L_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "52_7L"){instrument <- "ThermoFisher 7L_532-594"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot7L_532_594_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 55){instrument <- "ThermoFisher 7L_488-561"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot7L_488_561_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 51){instrument <- "ThermoFisher 6L_785"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot6L_785_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "52_6L"){instrument <- "ThermoFisher 6L_445"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "ThermoFisherBigfoot6L_445_ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 78){instrument <- "BD S8"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "BDFACSDiscoverS8ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "48_S"){instrument <- "BD S6"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "BDFACSSymphonyS6ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == "48_A"){instrument <- "BD A5"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "BDFACSSymphonyA5ReferenceLibrary.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 64){instrument <- "Cytek 5L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary5L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 54){instrument <- "Cytek 4L UV"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary4LUV.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 48){instrument <- "Cytek 4L YG"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary4LYG.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 38){instrument <- "Cytek 3L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary3L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 30){instrument <- "Cytek 2L VB"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary2LVB.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 22){instrument <- "Cytek 2L BR"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary2LBR.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else if (NumberDetectors == 14){instrument <- "Cytek 1L"
  FileLocation <- system.file("extdata", package = "Luciernaga")
  TheFile <- file.path(FileLocation, "CytekReferenceLibrary1L.csv")
  ReferenceData <- read.csv(TheFile, check.names = FALSE)
  } else {message("No References Found")}
}