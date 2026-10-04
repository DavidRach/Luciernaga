#' Personal convenience function for initial repository setup for our projects. 
#' 
#' @noRd
FolderSetup <- function(){
  WorkingDirectory <- getwd()
  TheRepoName <- basename(WorkingDirectory)
  
  Current <- list.files(WorkingDirectory, include.dirs=TRUE)

  if (!"README.md" %in% Current){
    writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, "README.md"))
  }

  if (!"LICENSE.md" %in% Current){
    download.file("https://raw.githubusercontent.com/DavidRach/Luciernaga/main/LICENSE.md", "LICENSE.md")
  }

  if (!"data" %in% Current){
    dir.create("data")
    StoreHere <- file.path("data", "README.md")
    writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, StoreHere))
  } else{
    Objects <- list.files("data", pattern="README", full.names=TRUE)
    if (length(Objects) == 0){
      StoreHere <- file.path("data", "README.md")
      writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, StoreHere))
    }
  }

  if (!"images" %in% Current){
    dir.create("images")
    StoreHere <- file.path("images", "README.md")
    writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, StoreHere))
  } else{
    Objects <- list.files("images", pattern="README", full.names=TRUE)
    if (length(Objects) == 0){
      StoreHere <- file.path("images", "README.md")
      writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, StoreHere))
    }
  }

  if (!"outputs" %in% Current){
    dir.create("outputs")
    StoreHere <- file.path("outputs", "README.md")
    writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, StoreHere))
  } else{
    Objects <- list.files("outputs", pattern="README", full.names=TRUE)
    if (length(Objects) == 0){
      StoreHere <- file.path("outputs", "README.md")
      writeLines(paste0("# ", TheRepoName), file.path(WorkingDirectory, StoreHere))
    }
  }

  message("Done")

}
