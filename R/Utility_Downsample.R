#' Downsample a gs based on subset
#'
#' @param gs A gating set object
#' @param sample.name Keyword specifying sample name
#' @param removestrings Value to be removed from sample name
#' @param subsets The gating hierarchy subset you want to include
#' @param subsample Total number of events to subsample from each specimen
#' @param newName File Name for the Concatenate File
#' @param outpath Location to store the concatenated file
#'
#' @return A concatenated .fcs file to new location
#' @export
#'
#' @examples NULL


Utility_Downsample <- function(x, sample.name, removestrings,
                                  subsets, subsample){
  name <- keyword(x, sample.name)
  alternatename <- NameCleanUp(name, removestrings)

  #Retrieving the exprs data for my subset population of interest
  ff <- gs_pop_get_data(x, subsets)
  #newff <- realize_view(ff)

  df <- exprs(ff[[1]])
  DF <- as.data.frame(df, check.names = FALSE)

  # If down-sampling is specified
  if(!is.null(subsample)){DF <- slice_sample(DF, n = subsample,
                                             replace = FALSE)
  } else{DF <- DF}

  DF <- DF %>% mutate(specimen = alternatename)

  return(DF)
}
