Utility_Luciernaga <- function(LuciernagaData, LuciernagaVariables, LuciernagaClusters, outfolder, filename){
  Items <- data.frame(table(LuciernagaData[[LuciernagaVariables]])) %>% pull(Var1) %>% as.character(.)

  FirstDetectorColumn <- which(grepl("\\d", colnames(LuciernagaData)))[1]
  LastDetectorColumn <- tail(which(grepl("\\d", colnames(LuciernagaData))), 1)

  theplots <- list()

  for (i in Items){LuciernagaSubset <- LuciernagaData %>% filter(.[[LuciernagaVariables]] %in% c(i))
  colnames(LuciernagaSubset) <- gsub("-A", "", colnames(LuciernagaSubset))
  zero_sum_rows <- LuciernagaSubset %>% filter(rowSums(select(., all_of(FirstDetectorColumn:LastDetectorColumn)), na.rm = TRUE) == 0)
  zero_sums <- nrow(zero_sum_rows)
  print(paste0(i, sep = "  ", zero_sums))
  DemonSpawn <- LuciernagaSubset %>% filter(rowSums(select(., all_of(FirstDetectorColumn:LastDetectorColumn)), na.rm = TRUE) == 0)
  LuciernagaSubset <- LuciernagaSubset %>% filter(rowSums(select(., all_of(FirstDetectorColumn:LastDetectorColumn)), na.rm = TRUE) != 0) #Temporary Debugging WorkAround.

  LinePlotData <- LuciernagaSubset %>% select({{LuciernagaClusters}}, {{FirstDetectorColumn}}:{{LastDetectorColumn}})
  LineCols <- ncol(LinePlotData)
  detector_order <- colnames(LuciernagaSubset)[FirstDetectorColumn:LastDetectorColumn]

  Melted <- gather(LinePlotData, key = "Detector", value = "value", all_of(2:LineCols)) #Gather is my New Best Friend

  Melted$Detector <- factor(Melted$Detector, levels = detector_order)
  Melted$Cluster <- factor(Melted$Cluster)

  #Change this in case raw values provided instead of normalized?
  Low <- 0
  High <- 1.1

  Melted1 <- data.frame(Melted)

  TheNormPlot <- ggplot(Melted1, aes(x = Detector, y = value, group = Cluster, color = Cluster)) + geom_line() + ylim(min = Low, max = High) +
    labs(title = "Fluorophores", x = "Detectors", y = "Normalized Values") + theme_bw() + scale_color_hue(direction = 1) + theme_linedraw() +
    theme(plot.title = element_text(size = 16L, face = "plain", hjust = 0.5), axis.title.y = element_text(size = 11L, face = "plain"),
          axis.title.x = element_text(size = 11L, face = "plain"), panel.grid.major = element_line(colour = "gray95", linetype = "twodash"),
          panel.grid.minor = element_line(colour = "gray95", linetype = "longdash"), panel.background = element_rect(fill = NA),
          plot.background = element_rect(colour = NA), legend.background = element_rect(fill = NA),
          axis.text.x = element_text(size = 5, angle = 45, hjust = 1))

  Names <- LinePlotData$Cluster
  Numbers <- LinePlotData %>% select_if(is.numeric)
  NumericsT <- t(Numbers)
  rownames(NumericsT) <- NULL
  colnames(NumericsT) <- Names
  NumericsT <- data.matrix(NumericsT)

  if(!ncol(NumericsT) == 1){
    CosineMatrix <- lsa::cosine(NumericsT)
    CosineMatrix <- round(CosineMatrix, 2)

    # Reorder By Similarity
    reorder_cormat <- function(cormat){
      dd <- as.dist((1-cormat)/2)
      hc <- hclust(dd)
      cormat <-cormat[hc$order, hc$order]
    }

    # Reorder the correlation matrix
    cormat <- reorder_cormat(CosineMatrix)
    melted_cormat <- reshape2::melt(cormat)

    #Generate a Red to Blue Heatmap
    CosineHeatMap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
      scale_fill_gradient2(low = "lightblue", high = "orange", mid = "white", midpoint = 0.7, limit = c(0.4,1), space = "Lab", name="Cosine\nSimilarity") +
      theme_bw() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + coord_fixed(ratio = 1.3) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),
            panel.background = element_blank(), axis.ticks = element_blank(), legend.position = c(1.2, 0.5),
            legend.direction = "vertical", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6), axis.text.y = element_text(size = 6),legend.key.size = unit(0.4, "cm")) #+ labs(title = "Cosine Similarity: Fluorophores")

    #print(CosineHeatMap)

    CosineOrder <- data.frame(table(melted_cormat$Var1)) %>% pull(Var1) %>% as.character(.)
    } else{CosineHeatMap <- NULL}

    #CosineOrder

    BarChartData <- LuciernagaSubset %>% mutate(Ratio = round(Ratio, 2))

    if(exists("CosineOrder")){BarChartData$Cluster <- factor(BarChartData$Cluster, levels = unique(BarChartData$Cluster)[order(match(unique(BarChartData$Cluster), CosineOrder))])}

  StackedBarChart <- ggplot(BarChartData, aes(x= sample, y = Ratio, fill = Cluster)) + geom_col() + theme_bw() + scale_fill_viridis(discrete = TRUE, option = "inferno", direction = -1) + labs(title = i) + theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10), axis.title.x = element_blank(), legend.key.size = unit(0.4, "cm")) + coord_fixed(ratio = 2)

HeatMapChart <- ggplot(BarChartData, aes(x= sample, y = Cluster, fill = Ratio)) + geom_tile() + geom_text(aes(label = Ratio)) + theme_bw() + scale_fill_gradient(name = "Ratio", low = "#FFFFFF", high = "#FF0000", limits = c(0, NA)) + theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), legend.key.size = unit(0.4, "cm"))  + coord_fixed(ratio = 1.1)

   theNegplots <- list(TheNormPlot, HeatMapChart, CosineHeatMap)
   #theNegplots <- list(TheNormPlot, StackedBarChart, HeatMapChart, CosineHeatMap)

   theplots[[length(theplots) +1]] <- theNegplots
   }

   theflattenedplots <- purrr::flatten(theplots)
   #theflattenedplots <- Filter(Negate(is.null), theflattenedplots)
   #theflattestplots <- purrr::flatten(theflattenedplots)

Utility_SCPatchwork(x = theflattenedplots, filename = filename, outfolder = outfolder)
}
