#' Internal for CartoonSignatures, returns plot of all similar fluorophores
#'
#' @param TheseFluorophores The similar fluorophores identified by cosine
#' @param TheFluorophore The one we were originally interested in
#' @param data The reference data of fluorophore signatures
#' @param legend Default TRUE, alternately removes plot legend
#' @param plotname Default NULL, alternately specifies the plot title
#' @param plotlinecolor Expects NULL, otherwise if single line provide
#'  desired color
#' @param unstained Default NULL, alternatively adds corresponding
#'  unstained signature
#' @param block_size Numeric increment for LEGO blocks. Default NULL
#'  automatically chooses 0.1 for normalized data (<= 2) or 1000 for raw data (> 2).
#' @param ylim Numeric scalar (upper limit, e.g. 100000) or vector of length 2 (e.g. c(0, 100000)).
#'  Default NULL scales y-axis automatically to fit data.
#'
#' @importFrom dplyr filter rename select mutate relocate
#'  bind_cols bind_rows rowwise ungroup case_when
#' @importFrom tidyselect where everything
#' @importFrom tidyr pivot_longer unnest
#' @importFrom ggplot2 ggplot aes geom_rect geom_blank theme_bw labs
#'  geom_hline theme element_text element_blank scale_x_continuous
#'  scale_fill_identity coord_cartesian
#'
#' @return An internal value
#'
#' @noRd
CartoonInternal <- function(TheseFluorophores, TheFluorophore, data,
                            legend = TRUE, plotname = FALSE, plotlinecolor,
                            unstained = NULL, block_size = NULL, ylim = NULL){

      TheseFluorophores <- as.character(TheseFluorophores)

      if(!is.null(unstained)){
        These <- c(TheFluorophore, TheseFluorophores, "Unstained")

        if (any(unstained |> select(where(is.numeric)) > 1)){
          Metadata <- unstained |> select(!where(is.numeric))
          Numerics <- unstained |> select(where(is.numeric))
          n <- Numerics
          n[n < 0] <- 0
          A <- do.call(pmax, n)
          Normalized <- n/A
          unstained <- bind_cols(Metadata, Normalized)
        }

        UnstainedData <- unstained |>
          select(where(is.numeric)) |>
          pivot_longer(
            cols=everything(),
            names_to="Detector",
            values_to="value"
          ) |>
          mutate(
            Instrument="Existing",
            Fluorophore="Unstained"
          ) |>
          relocate(
            Instrument,
            Fluorophore,
            .before="Detector"
          )

      } else {
        These <- c(TheFluorophore, TheseFluorophores)
      }

      TheData <- data |>
        filter(Fluorophore %in% These) |>
        rename(value=AdjustedY)

      if(!is.null(unstained)){
        TheData <- bind_rows(TheData, UnstainedData)
      }

      TheData$Detector <- gsub("-A", "", TheData$Detector)

      Iterations <- TheData |>
        filter(Fluorophore %in% These[[1]]) |>
        nrow()

      if (is.character(TheData$Detector)) {
        MyVector <- TheData |>
          filter(Fluorophore %in% These[[1]]) |>
          pull(Detector)
      }

      if (is.numeric(TheData$Detector)) {
        MyVector <- seq_len(Iterations)
      }

      # --- Determine scale type & block step size ---
      max_val <- max(TheData$value, na.rm = TRUE)

      if (max_val > 2){
        YAxisLabel <- "Raw MFI"
        step_inc <- if (is.null(block_size)) 1000 else block_size
      } else {
        YAxisLabel <- "Normalized Value"
        step_inc <- if (is.null(block_size)) 0.1 else block_size
      }

      TheData$Detector <- factor(
        TheData$Detector,
        levels=MyVector
      )

      TheData$Fluorophore <- factor(
        TheData$Fluorophore,
        levels=These
      )

      # --- Assign block color based on laser designation ---
      TheData <- TheData |>
        mutate(
          BlockColor=case_when(
            grepl("UV", as.character(Detector), fixed=TRUE) ~ "purple",
            grepl("V", as.character(Detector), fixed=TRUE) ~ "violet",
            grepl("B", as.character(Detector), fixed=TRUE) ~ "blue",
            grepl("YG", as.character(Detector), fixed=TRUE) ~ "green",
            grepl("R", as.character(Detector), fixed=TRUE) ~ "red",
            TRUE ~ "gray"
          )
        )

      # --- Calculate horizontal position & scale values ---
      TheData <- TheData |>
        mutate(
          # Apply ceiling approach: values strictly between 0 and step_inc are bumped up to step_inc
          value = if (max_val <= 2) {
            v <- round(value, 1)
            if_else(v > 0 & v < step_inc, step_inc, v)
          } else {
            v <- round(value)
            if_else(v > 0 & v < step_inc, step_inc, v)
          },
          DetectorNum=as.integer(Detector),
          FluorNum=as.integer(Fluorophore),
          nFluor=n_distinct(Fluorophore),
          BlockWidth=0.8 / nFluor,
          XCenter=DetectorNum + (FluorNum - (nFluor + 1) / 2) * BlockWidth,
          xmin=XCenter - BlockWidth / 2,
          xmax=XCenter + BlockWidth / 2
        )

      # --- Expand positive values into dynamic LEGO blocks ---
      TheData <- TheData |>
        rowwise() |>
        mutate(
          Block=list(
            if (!is.na(value) && value >= step_inc) {
              seq(step_inc, value, by=step_inc)
            } else {
              NA_real_
            }
          )
        ) |>
        unnest(Block) |>
        ungroup() |>
        mutate(
          ymin=if_else(is.na(Block), 0, Block - step_inc),
          ymax=if_else(is.na(Block), 0, Block)
        )

      if (is.null(plotname)){
        TheTitle <- paste0(TheFluorophore)
      } else {
        TheTitle <- plotname
      }

      BaseTheme <- list(
        theme_bw(),
        labs(
          title=TheTitle,
          x=NULL,
          y=YAxisLabel
        ),
        # Fix y-axis limits if specified (handles single upper bound or two-element range)
        if (!is.null(ylim)) {
          if (length(ylim) == 1) {
            coord_cartesian(ylim = c(0, ylim))
          } else {
            coord_cartesian(ylim = ylim)
          }
        },
        # Draw threshold reference line conditionally if using 0-1 scale
        if (max_val <= 2) {
          geom_hline(yintercept=1, linetype="dashed", color="red")
        },
        theme(
          plot.title=element_text(size=8),
          legend.position="none",
          axis.text.x=element_text(size=6, angle=45),
          panel.grid=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=8)
        ),
        scale_x_continuous(
          breaks=seq(1, Iterations, by=5),
          labels=MyVector[seq(1, Iterations, by=5)]
        )
      )

      if (!is.null(plotlinecolor)){
        ThePlot <- ggplot(TheData) +
          geom_blank(aes(x=XCenter, y=0)) +
          geom_rect(
            data=TheData |> filter(!is.na(Block)),
            aes(
              xmin=xmin,
              xmax=xmax,
              ymin=ymin,
              ymax=ymax
            ),
            fill=plotlinecolor,
            color="white",
            linewidth=0.4
          ) +
          BaseTheme

      } else {
        ThePlot <- ggplot(TheData) +
          geom_blank(aes(x=XCenter, y=0)) +
          geom_rect(
            data=TheData |> filter(!is.na(Block)),
            aes(
              xmin=xmin,
              xmax=xmax,
              ymin=ymin,
              ymax=ymax,
              fill=BlockColor
            ),
            color="white",
            linewidth=0.4
          ) +
          scale_fill_identity() +
          BaseTheme
      }

      return(ThePlot)
}