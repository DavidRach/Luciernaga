  tmp <- withr::local_tempdir(pattern = "Luciernaga")
  withr::local_dir(tmp)

  #StorageLocation <- tmp

  # Load Libraries

  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(Luciernaga)
  library(stringr)
  library(dplyr)
  library(purrr)

  # Miscellaneous
  removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")

  # Find FCS Files
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)

  # Unstained Setup
  UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
  UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
  MyUnstainedCytoSet <- load_cytoset_from_fcs(UnstainedCells[c(1,3,5)],
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MyUnstainedGatingSet <- GatingSet(MyUnstainedCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyUnstainedGatingSet)

# SingleColor Setup
  CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
  CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
  MySCsCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MySCsGatingSet <- GatingSet(MySCsCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MySCsGatingSet)
