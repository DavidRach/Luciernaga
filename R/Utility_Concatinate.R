#' Concatinate a gs based on subset and subsample
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

Utility_Concatinate <- function(gs, sample.name, removestrings, subsets, subsample, newName, outpath){
gs <- gs
sample.name <- sample.name
removestrings <- removestrings
subsample <- subsample

ConcatenatedFile <- map(gs, .f=.Internal_Downsample, sample.name=sample.name, removestrings=removestrings,
                        subsets=subsets, subsample=subsample) %>% bind_rows

#ConcatenatedMatrix <- as.matrix(ConcatenatedFile)

ConcatenatedExtra <- ConcatenatedFile %>% select(specimen)
ConcatenatedMain <- ConcatenatedFile %>% select(-specimen)
ExtraMatrix <- as.matrix(ConcatenatedExtra)
MainMatrix <- as.matrix(ConcatenatedFile)

ff <- gs_pop_get_data(gs, subsets)
FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
original_p <- parameters(FlowFrameTest)
original_d <- keyword(FlowFrameTest)

new_fcs <- new("flowFrame", exprs=MainMatrix, parameters=original_p, description=original_d)
#new_fcs <- new("flowFrame", exprs=ConcatenatedMatrix, parameters=original_p, description=original_d)

fr <- new_fcs
cols <- ExtraMatrix

#new_pd <- flowCore:::cols_to_pd(fr=fr, cols=cols) #Using Internal Function :( Bioconductor?

#pd <- pData(parameters(fr))
#pd <- rbind(pd, new_pd)
#fr@exprs <- cbind(exprs(fr), cols)
#pData(parameters(fr)) <- pd

#new_pid <- rownames(new_pd)
#new_kw <- fr@description

#for (i in new_pid){
#  new_kw[paste0(i,"B")] <- new_kw["$P1B"] #Unclear Purpose
#  new_kw[paste0(i,"E")] <- "0,0"
#  new_kw[paste0(i,"N")] <- new_pd[[i,1]]
#  #new_kw[paste0(i,"V")] <- new_kw["$P1V"] # Extra Unclear Purpose
#  new_kw[paste0(i,"R")] <- new_pd[[i,5]]
#  new_kw[paste0(i,"DISPLAY")] <- "LIN"
#  new_kw[paste0(i,"TYPE")] <- "Dimensionality"
#}

#new_kw
#UpdatedParameters <- parameters(fr)
#UpdatedExprs <- exprs(fr)

#new_fcs <- new("flowFrame", exprs=UpdatedExprs, parameters=UpdatedParameters, description=new_kw)

TheFileName <- paste0(newName, ".fcs")

fileSpot <- file.path(outpath, TheFileName)

write.FCS(new_fcs, filename = fileSpot, delimiter="#")
}

.Internal_Downsample <- function(x, sample.name, removestrings,
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
