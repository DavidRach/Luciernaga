#' Add Dimensionality Visualized parameters to raw .fcs files
#'
#' @param ff A realized_view object from flowWorkspace
#' @param DF The maybe downsampled exprs
#' @param columnframe The dimensionality visualized data.frame object to be added
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_cols
#' @importFrom flowCore parameters
#' @importFrom flowCore keyword
#' @importFrom flowCore exprs
#' @importFrom Biobase pData
#'
#' @return A new flow_frame object.
#' @export
#'
#' @examples NULL

Utility_ColAppend <- function(ff, DF, columnframe){
  TheColNames <- colnames(columnframe)

  #x <- TheColNames[1]

  ShiftedColumns <- map(.x=TheColNames, .f=InternalShift) %>% bind_cols
  Shifted <- as.matrix(ShiftedColumns)
  FCSSubset <- as.matrix(DF)

  #Updating in case of Downsample
  FlowFrameTest <- ff[[1, returnType = "flowFrame"]]
  original_p <- parameters(FlowFrameTest)
  original_d <- keyword(FlowFrameTest)
  new_fcs <- new("flowFrame", exprs=FCSSubset, parameters=original_p, description=original_d)

  #Adding new columns (modified flowCore utilities)
  fr <- new_fcs
  cols <- Shifted
  new_pd <- flowCore:::cols_to_pd(fr=fr, cols=cols) #Using Internal Function :( Bioconductor?

  pd <- pData(parameters(fr))
  pd <- rbind(pd, new_pd)
  fr@exprs <- cbind(exprs(fr), cols) ### Another Bioconductor :( for @
  pData(parameters(fr)) <- pd

  new_pid <- rownames(new_pd)
  new_kw <- fr@description ### Another Bioconductor :( for @

  for (i in new_pid){
    new_kw[paste0(i,"B")] <- new_kw["$P1B"] #Unclear Purpose
    new_kw[paste0(i,"E")] <- "0,0"
    new_kw[paste0(i,"N")] <- new_pd[[i,1]]
    #new_kw[paste0(i,"V")] <- new_kw["$P1V"] # Extra Unclear Purpose
    new_kw[paste0(i,"R")] <- new_pd[[i,5]]
    new_kw[paste0(i,"DISPLAY")] <- "LIN"
    new_kw[paste0(i,"TYPE")] <- "Dimensionality"
    new_kw[paste0("flowCore_", i, "Rmin")] <- new_pd[[i,4]]
    new_kw[paste0("flowCore_", i, "Rmax")] <- new_pd[[i,5]]
  }

  new_kw
  UpdatedParameters <- parameters(fr)
  UpdatedExprs <- exprs(fr)

  new_fcs <- new("flowFrame", exprs=UpdatedExprs, parameters=UpdatedParameters, description=new_kw)

  return(new_fcs)
}

#' Internal for Column Append
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @noRd

InternalShift <- function(x){
  TheColumn <- columnframe %>% select(all_of(x))
  ShiftedColumn <- TheColumn + abs(min(TheColumn))+1
  return(ShiftedColumn)
}
