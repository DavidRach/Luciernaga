#' Plot Individual Fluorophore Gates, take the MFI, return the values.
#'
#' @param x A GatingSet type object (ex. gs or gs[[1]])
#' @param sample.name Keyword variable under which samples indentity is stored (ex. "GUID")
#' @param experiment Provide directly the experiment label (ex. "Jan2024")
#' @param experiment.date Keyword variable under which experiment is stored in absence of Experiment
#' @param rootlevel Gating Hierarchy node under which new gates will be constructed
#' @param samplinglevel Gating Hierarchy level from which MFI will be taken
#' @param bins When plotting, how many bins resolution
#' @param stats Whether to use "mean" or "median"
#' @param outpath Location which to save output
#' @param source  Whether to export for plotting
#' @param sourcelocation Location where source plotting file is stored.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
Utility_SingleColors <- function(x, sample.name, experiment = NULL, experiment.date, rootlevel, samplinglevel, bins, stats = NULL, outpath, source, sourcelocation){
  x <- x
  name <- keyword(x, sample.name)

  if(str_detect(name, "(Cells)")){Type <- "Cells"} else if(str_detect(name, "(Beads)")){Type <- "Beads"} else(Type <- NULL)

  name <- gsub(".fcs", "", gsub(" (Cells)", "", fixed = TRUE, gsub(" (Beads)", "", fixed = TRUE,name)))
  alternate.name <- name
  alternate.name <- gsub(".fcs", "", gsub(" ", "", gsub("(", "", fixed = TRUE, gsub(")", "", fixed = TRUE, gsub("_", "", fixed = TRUE, gsub("-", "", fixed = TRUE, alternate.name))))))

  if(!is.null(experiment)){Experiment <- experiment
  } else {experiment <- keyword(x, experiment.date)
  experiment <- gsub("-", "", fixed = TRUE, experiment)
  Experiment <- experiment}

  AggregateName <- alternate.name
  StorageLocation <- paste(outpath, AggregateName, sep = "/", collapse = NULL)

  #Unstained
  if (str_detect(name, "Unstained")){
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "FSC-A", y = "SSC-A"), subset = "nonDebris") + geom_hex(bins=bins) + geom_gate("lymph") + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    cells <- "lymph"
    if(source == TRUE){source(sourcelocation, local = TRUE)}

  } else if (str_detect(name, "BUV395")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BUV395")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV395", gating_method = "flowClust.1d", gating_args = "K=2", dims = "UV2-A")}
    cells <- "BUV395"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV2-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BUV496")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BUV496")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV496", gating_method = "flowClust", gating_args = "K=2", dims = "UV7-A")}
    cells <- "BUV496"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV7-A", y = "FSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BUV563")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BUV563")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV563", gating_method = "flowClust.1d", gating_args = "K=2", dims = "UV9-A")}
    cells <- "BUV563"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV9-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if(str_detect(name, "BUV615")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "/BUV615")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV615", gating_method = "flowClust.1d", gating_args = "K=2", dims = "UV10-A")}
    cells <- "BUV615"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV10-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BUV661")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BUV661")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV661", gating_method = "gate_quantile", gating_args = "probs = 0.94", dims = "UV11-A")}
    cells <- "BUV661"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV11-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BUV737")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BUV737")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV737", gating_method = "gate_quantile", gating_args = "probs = 0.94", dims = "UV14-A")}
    cells <- "BUV737"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV14-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BUV805")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BUV805")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BUV805", gating_method = "mindensity", dims = "UV16-A")}
    cells <- "BUV805"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "UV16-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV421")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV421")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV421", gating_method = "mindensity", dims = "V1-A")}
    cells <- "BV421"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V1-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "Pacific Blue")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "PacificBlue")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "PacificBlue", gating_method = "mindensity", dims = "V3-A")}
    cells <- "PacificBlue"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V3-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV480")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV480")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV480", gating_method = "flowClust.1d", gating_args = "K=2", dims = "V5-A")}
    cells <- "BV480"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V5-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV510")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV510")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV510", gating_method = "flowClust.1d", gating_args = "K=2", dims = "V7-A")}
    cells <- "BV510"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V7-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}

  } else if (str_detect(name, "BV605")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV605")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV605", gating_method = "flowClust.1d", gating_args = "K=2", dims = "V10-A")}
    cells <- "BV605"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V10-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV650")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV650")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV650", gating_method = "flowClust.1d", gating_args = c("K=2"), dims = "V11-A")}
    cells <- "BV650"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V11-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}

  } else if (str_detect(name, "BV711")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV711")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV711", gating_method = "flowClust.1d", gating_args = "K=2", dims = "V13-A")}
    cells <- "BV711"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V13-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV750")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV750")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV750", gating_method = "mindensity", dims = "V14-A")}
    cells <- "BV750"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V14-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV785")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV785")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV785", gating_method = "flowClust.1d", gating_args = c("K=2"), dims = "V15-A")}
    cells <- "BV785"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V15-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "BV786")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "BV786")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "BV786", gating_method = "flowClust.1d", gating_args = c("K=2"), dims = "V15-A")}
    cells <- "BV786"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "V15-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "FITC")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "FITC")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "FITC", gating_method = "mindensity", dims = "B2-A")}
    cells <- "FITC"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "B2-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}

  } else if (str_detect(name, "Alexa Fluor 488")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "AlexaFluor488")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "AlexaFluor488" , gating_method = "mindensity",  gating_args = "gate_range=c(0.3e5,1e6)", dims = "B2-A")}
    cells <- "AlexaFluor488"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "B2-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "Spark Blue 550")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "SparkBlue550")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "SparkBlue550", gating_method = "mindensity", dims = "B3-A")}
    cells <- "SparkBlue550"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "B3-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "PerCP-Cy5.5")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "PerCPCy55")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "PerCPCy55", gating_method = "gate_quantile", gating_args = "probs = 0.8", dims = "B9-A")}
    cells <- "PerCPCy55"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "B9-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "PE$")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "PE$")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "PE", gating_method = "flowClust.1d", gating_args = "K=2", dims = "YG1-A")}
    cells <- "PE"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "YG1-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "PE-Dazzle594")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "PEDazzle594")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "PEDazzle594", gating_method = "flowClust.1d", gating_args = "K=2", dims = "YG3-A")}
    cells <- "PEDazzle594"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "YG3-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "PE-Cy5")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "PECy5")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "PECy5", gating_method = "gate_quantile", gating_args = "probs = 0.9", dims = "YG5-A")}
    cells <- "PECy5"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "YG5-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "PE-Vio770")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "PEVio770")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "PEVio770", gating_method = "flowClust.1d", gating_args = "K=2", dims = "YG9-A")}
    cells <- "PEVio770"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "YG9-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}



  } else if (str_detect(name, "APC$")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "APC$")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "APC", gating_method = "flowClust.1d", gating_args = "K=2", dims = "R1-A")}
    cells <- "APC"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "R1-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "Alexa Fluor 647")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "AlexaFluor647")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "AlexaFluor647", gating_method = "mindensity", gating_args = "gate_range=c(0,3e5)", dims = "R2-A")}
    cells <- "AlexaFluor647"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "R2-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "APC-R700")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "APCR700")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "APCR700", gating_method = "flowClust.1d", gating_args = "K=2", dims = "R4-A")}
    cells <- "APCR700"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "R4-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "Zombie NIR")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "ZombieNIR")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "ZombieNIR", gating_method = "mindensity", dims = "R6-A")}
    cells <- "ZombieNIR"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "R6-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "APC-Fire 750")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "APCFire750")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "APCFire750", gating_method = "flowClust.1d", gating_args = "K=2", dims = "R7-A")}
    cells <- "APCFire750"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "R7-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}


  } else if (str_detect(name, "APC-Fire 810")){
    testthis <- gs_get_leaf_nodes(x)
    if(!str_detect(testthis, "APCFire810")){gs_add_gating_method(x, parent = rootlevel, pop = "+", alias = "APCFire810", gating_method = "flowClust.1d", gating_args = "K=2", dims = "R8-A")}
    cells <- "APCFire810"
    SingleColor <- as.ggplot(ggcyto(x, aes(x = "R8-A", y = "SSC-A"), subset = rootlevel) + geom_hex(bins=bins) + geom_gate(cells) + scale_x_flowjo_biexp() + theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), axis.title = element_text(size = 10, face = "bold"), legend.position = "none"))
    if(source == TRUE){source(sourcelocation, local = TRUE)}
  } else{newname <- paste0("WHERE DID YOU GO", name)
  print(newname)}

  if(samplinglevel == "cells"){data <- gs_pop_get_data(x, cells)
  } else {data <- gs_pop_get_data(x, samplinglevel)}

  #nrow(data)
  df <- flowCore::exprs(data[[1]])
  DF <- as.data.frame(df, check.names = FALSE)
  CleanedDF <- DF[,-grep("Time|FS|SC|SS|Original", names(DF))]

  if(stats == "mean"){Samples <- CleanedDF %>% summarize_all(mean)
  } else if (stats == "median"){Samples <- CleanedDF %>% summarize_all(median)
  } else(print("NA"))

  Specimen <- name
  Output <- cbind(Specimen, Experiment, Type, Samples)
  Output
}
