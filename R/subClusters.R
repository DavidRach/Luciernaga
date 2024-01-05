subClusters <- function(data){
  data$experiment <- factor(data$experiment)
  numbers <- data %>% select(2:ncol(data))

  if(nrow(numbers) > 1){
    Names <- data$experiment
    NumericsT <- t(numbers)
    rownames(NumericsT) <- NULL
    colnames(NumericsT) <- Names
    NumericsT <- data.matrix(NumericsT)
    CosineMatrix <- lsa::cosine(NumericsT)
    CosineMatrix <- round(CosineMatrix, 2)
    data <- data %>% mutate(subCosine = CosineMatrix[,1]) %>% relocate(subCosine, .after = experiment)
  } else {data <- data %>% mutate(subCosine = 1)}

  return(data)
}
