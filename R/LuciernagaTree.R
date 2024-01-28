
#' Select candidate Luciernaga output .fcs files for future use in unmixing.
#'
#' @param x A vector of the fluorophores found in data.
#' @param data  The data.frame of Luciernaga outputs used to candidate .fcs files
#'  for unmixing.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr row_number
#' @importFrom dplyr relocate
#'
#' @return A data.frame listing the candidate .fcs files for future unmixing use.
#' @export
#'
#' @examples NULL

LuciernagaTree <- function(x, data){
  Internal <- data %>% filter(sample %in% x)

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
    ComplexityLowerBound <- Complexity[1]*0.9

    if(MainComplexity < ComplexityLowerBound){Internal3 <- Internal2 %>%
      arrange(desc(Brightness))
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
