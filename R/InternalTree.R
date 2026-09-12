
#' Internal for Luciernaga_Tree
#'
#' @param x The iterated Fluorophore to be filtered
#' @param TheData The Data for all fluorophores
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom dplyr row_number
#'
#' @return An internal value
#'
#' @noRd
InternalTree <- function(x, TheData){
  OriginalX <- x

  Internal <- TheData %>% dplyr::filter(str_detect(sample, fixed(x, ignore_case = TRUE)))

  if (x %in% c("PE", "APC")){

    Internal <- Internal %>% dplyr::filter(!str_detect(sample, "PE-|APC-|Per"))
  } #ExceptionHandling


  Total <- sum(Internal$Count, na.rm = TRUE)
  Internal <- Internal %>% mutate(Ratio = Count / Total) %>% relocate(Ratio, .after=Count)

  if(nrow(Internal) == 0){message("No Fluorophore for ", x, " was found")}

  if(nrow(Internal)>1){
    Internal <- Internal %>% arrange(desc(Detector1Raw)) #First Arrange
    MaxVal <- Internal %>% filter(row_number() == 1) %>% pull(Detector1Raw)
    Internal1 <- Internal %>% filter(Detector1Raw > MaxVal*0.8)

    if (nrow(Internal1)>1){Abundance <- Internal1 %>%
      filter(row_number() == 1) %>%  pull(Ratio)

    if(Abundance < 0.5){Internal2 <- Internal1 %>% arrange(desc(Ratio))
    #Second Arrange
    Proportion <- Internal2 %>% select(Ratio) %>% sum(.)
    Top <- Internal2 %>% filter(row_number() == 1) %>% pull(Ratio)

    if((Top/Proportion) < 0.5){MainComplexity <- Internal2 %>%
      filter(row_number() == 1) %>% pull(Brightness)
    Complexity <- Internal2 %>% arrange(desc(Brightness)) %>% pull(Brightness)
    ComplexityLowerBound <- Complexity[1]*0.6

    if(MainComplexity > ComplexityLowerBound){Internal3 <- Internal2 %>%
      arrange(Brightness)
    SubsetData <- Internal3 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Contested") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)

    } else {SubsetData <- Internal2 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Fifth Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
    } else {SubsetData <- Internal2 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Fourth Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
    } else {SubsetData <- Internal1 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Third Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
    } else {SubsetData <- Internal1 %>% filter(row_number() == 1)
    SubsetData <- SubsetData %>% mutate(Decision = "Second Level") %>% relocate(
      Decision, .after = Cluster)
    return(SubsetData)}
  } else {SubsetData <- Internal %>% filter(row_number() == 1)
  SubsetData <- SubsetData %>% mutate(Decision = "First Level") %>% relocate(
    Decision, .after = Cluster)
  return(SubsetData)}
}