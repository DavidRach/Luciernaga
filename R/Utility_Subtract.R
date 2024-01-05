Utility_Subtract <- function(i, Data){
  Z <- Data[i,]
  Alpha <- Data %>% dplyr::filter(Experiment == Z$Experiment) %>% dplyr::filter(Type == Z$Type) %>% dplyr::filter(Negative == Z$Negative) %>% dplyr::filter(str_detect(Ligand, "Unstained"))
  Subtracted <- Z[,6:69] - Alpha[,6:69]
  Updated <- cbind(Z[,1:5], Subtracted)
  Updated %>% mutate(Processed = "Yes") %>% relocate(Processed, .after = Negative)
}
