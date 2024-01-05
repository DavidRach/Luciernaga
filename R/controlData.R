#' Prepare Luciernaga .fcs outputs as Control Inputs for Unmixing
#'
#' @param cs A flowSet (ex. cs or cs[[1]])
#' @param removestrings Portions of the single color names you want removed (Ex. c("DR_ILT_2023_", ".FCS")
#'
#' @return NULL
#' @export NULL
#'
#' @examples NULL
controlData <- function(cs, removestrings) {
  Control_Spectrums<-fsApply(cs,each_col,median)
  Control_Spectrums<-as.data.frame(Control_Spectrums)
  Control_Spectrums<-Control_Spectrums[,grep("-A", names(Control_Spectrums))]
  Control_Spectrums<-Control_Spectrums[,grep("SC|SS|FS", names(Control_Spectrums), invert=TRUE)]
  row.names(Control_Spectrums) <- NameCleanUp(row.names(Control_Spectrums), removestrings)
  row.names(Control_Spectrums) <- sub("(.*)_(.*)", "\\2_\\1", row.names(Control_Spectrums)) #Swaps order around the _
  Control_Spectrums <- Control_Spectrums[order(rownames(Control_Spectrums)), ]
  return(Control_Spectrums)
}
