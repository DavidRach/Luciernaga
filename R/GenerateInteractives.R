#' Makes the website
#' 
#' @param x The file path to the parent folder
#' @param outpath The desired storage location
#' @param filename The desired file name
#' 
#' @importFrom stringr str_detect str_extract
#' @importFrom dplyr filter select group_by slice desc arrange mutate pull
#' @importFrom purrr flatten map
#' @importFrom htmltools tagList h1 save_html
#' 
#' @return A htmlwebpage to desired location
#' 
#' @noRd
#' 
#' @examples
#' A <- 2 + 2 
GenerateInteractives <- function(x, outpath = NULL, filename = NULL) {
  TheHoller <- file.path(x, "Luciernaga")
  TheFiles <- list.files(TheHoller, full.names = TRUE)
  TheFiles <- TheFiles[!str_detect(TheFiles, "GatingReport.pdf")]

  TheUnstaineds <- TheFiles[grep("Unstained", TheFiles)]
  RemoveThese <- c("Unstained (Cells)_", ".fcs", "PMA")
  TheseTypes <- NameCleanUp(basename(TheUnstaineds),
                             removestrings = RemoveThese)
  MainAFs <- sub("10.*", "", TheseTypes) |> unique()

  if ("V" %in% MainAFs) {
    MainAFs[MainAFs == "V"] <- "V10"
  }
  if ("UV" %in% MainAFs) {
    MainAFs[MainAFs == "UV"] <- "UV10"
  }
  if ("B" %in% MainAFs) {
    MainAFs[MainAFs == "B"] <- "B10"
  }
  if ("YG" %in% MainAFs) {
    MainAFs[MainAFs == "YG"] <- "YG10"
  }

  TheStained <- TheFiles[!str_detect(TheFiles, "Unstained")]
  WatchForThese <- paste0("_", MainAFs, collapse = "|")

  TheCleanStained <- TheStained[!str_detect(TheStained, WatchForThese)]
  TheseDudes <- sub(" \\(Cells\\).*", "", basename(TheCleanStained)) |>
    unique()
  TheseDudettes <- sub("^[^ ]+ ", "", TheseDudes)
  TheseDudettes <- TheseDudettes[!str_detect(TheseDudettes, "Unstim")] |>
    unique()

  TheRefs <- Luciernaga:::InstrumentReferences(NumberDetectors = 64)

  Ranking <- TheRefs |>
    filter(Fluorophore %in% TheseDudettes) |>
    group_by(Fluorophore) |>
    arrange(desc(AdjustedY)) |>
    slice(1) |>
    select(-Instrument, -AdjustedY)

  Order <- c("UV", "V", "B", "YG", "R")

  Sequence <- Ranking |>
    mutate(
      prefix = str_extract(Detector, "^[A-Z]+"),
      num = as.numeric(str_extract(Detector, "\\d+")),
      group_order = match(prefix, Order)
    ) |>
    arrange(group_order, num) |>
    pull(Fluorophore)

  ThePlots <- map(.x = Sequence, .f = LuciernagaCheck,
                   TheCleanStained = TheCleanStained, returnType = "plotly")

  ThePlots <- flatten(ThePlots)

  if (is.null(outpath)) {
    outpath <- getwd()
  }
  if (is.null(filename)) {
    filename <- basename(x)
  }

  html_content <- tagList(
    h1("Experiment ", filename),
    grid_layout(ThePlots, ncol = 2)
  )

  Assembled <- paste0(filename, ".html")
  StorageLocation <- file.path(outpath, Assembled)

  # Save to an HTML file
  save_html(html_content, file = StorageLocation)
}