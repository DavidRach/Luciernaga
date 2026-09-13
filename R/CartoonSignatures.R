#' Visualizes fluorophore signatures with smooth density-colored ribbons
#'
#' @param x Default NULL, the identities to filter for in columnname.
#'   When NULL, all unique values of columnname are used.
#' @param columnname Default "Sample", the column identifying each
#'   fluorophore/sample row.
#' @param data The data.frame containing the signature data.
#' @param characterColumns Default NULL, character columns to combine
#'   into a single Fluorophore identity when more than one is present.
#' @param Expand Default NULL, number of simulated replicate rows to
#'   generate per original row when adding variance.
#' @param variance Default "medium", either a preset string ("none",
#'   "low", "medium", "high") or a numeric vector of length 1 or 2
#'   controlling simulated row-to-row variance.
#' @param Normalize Default TRUE, whether to normalize detector values
#'   to a 0-1 scale when values exceed 1.
#' @param legend Default TRUE, whether to display a plot legend.
#' @param plotname Default NULL, title for the returned plot.
#' @param plotlinecolor Default NULL, color override for the plotted
#'   line(s).
#' @param ylim Default NULL, numeric vector of length 1 or 2 setting
#'   the y-axis limits.
#' @param returnType Default "Cartoon", one of "Cartoon", "Amalgamated",
#'   "Density", or "Ribbon", selecting the plot style returned.
#' @param ribbon_ci TBD
#' @param theme Default "viridis", the color palette/theme used for
#'   "Density" and "Ribbon" return types.
#' @param bins Default 100, number of bins used for "Density" return
#'   type.
#' @param logScale Default TRUE, whether to plot the y-axis on a log
#'   scale.
#'
#' @importFrom dplyr pull filter select all_of rename mutate slice
#'   row_number group_by summarise across bind_rows
#' @importFrom tidyselect where everything
#' @importFrom tidyr unite pivot_longer
#' @importFrom rlang .data
#' @importFrom stats approx median quantile runif
#' @importFrom viridis viridis
#' @importFrom ggplot2 ggplot geom_ribbon aes geom_line
#'   scale_x_continuous theme_minimal labs theme element_text
#'   scale_y_log10 coord_cartesian
#'
#' @return A ggplot2 object of the selected returnType
#'
#' @export
CartoonSignatures <- function(x = NULL, columnname = "Sample", data,
                              characterColumns = NULL, Expand = NULL,
                              variance = "medium", Normalize = TRUE,
                              legend = TRUE, plotname = NULL,
                              plotlinecolor = NULL, ylim = NULL,
                              returnType = "Cartoon", ribbon_ci = 0.95,
                              theme = "viridis", bins = 100, logScale = TRUE) {

  if (is.null(Expand) &&
      returnType %in% c("Density", "Ribbon", "Amalgamated")) {
    Expand <- 100
  }

  if (!is.null(ylim)) {
    if (!is.numeric(ylim)) {
      stop("`ylim` must be a numeric vector.")
    }
    if (length(ylim) == 1) {
      ylim <- c(0, ylim)
    } else if (length(ylim) != 2) {
      stop("`ylim` must be a numeric vector of length 2, e.g., c(0, 100).")
    }
  }

  if (is.null(x)) {
    x <- data |> dplyr::pull(columnname)
  }

  StartingData <- data |> dplyr::filter(.data[[columnname]] %in% x)
  CharacterLength <- StartingData |>
    dplyr::select(!tidyselect::where(is.numeric)) |> length()

  if (CharacterLength == 0) {
    stop("Please add a non-numeric column, and provide its columnname")
  }

  if (CharacterLength > 1 && !is.null(characterColumns)) {
    Identity <- StartingData |>
      dplyr::select(dplyr::all_of(characterColumns)) |>
      tidyr::unite("combined", tidyselect::everything(), sep = "_") |>
      dplyr::pull()
    Identity <- data.frame(Fluorophore = Identity) |>
      dplyr::rename("Fluorophore" = 1) |>
      dplyr::mutate(Fluorophore = paste0("ID_", Fluorophore))
  } else {
    Identity <- StartingData |>
      dplyr::select(!tidyselect::where(is.numeric)) |>
      dplyr::rename("Fluorophore" = 1) |>
      dplyr::mutate(Fluorophore = paste0("ID_", Fluorophore))
  }

  DetectorCols <- StartingData |> dplyr::select(tidyselect::where(is.numeric))

  if (!is.null(Expand) && is.numeric(Expand) && Expand > 0) {
    if (is.character(variance)) {
      variance <- match.arg(tolower(variance),
        choices = c("none", "low", "medium", "high"))
      var_range <- switch(variance,
        "none"   = c(0.999, 1.001),
        "low"    = c(0.97, 1.03),
        "medium" = c(0.90, 1.10),
        "high"   = c(0.75, 1.25)
      )
    } else if (is.numeric(variance)) {
      if (length(variance) == 1) {
        var_range <- if (variance == 0) {
          c(0.999, 1.001)
        } else {
          c(1 - abs(variance), 1 + abs(variance))
        }
      } else if (length(variance) == 2) {
        var_range <- variance
      } else {
        stop("`variance` numeric vector must be of length 1 or 2.")
      }
    } else {
      stop("`variance` must be a preset string ('none', 'low', 'medium',
        'high') or numeric.")
    }

    idx <- rep(seq_len(nrow(DetectorCols)), each = Expand)
    Identity <- Identity |> dplyr::slice(idx) |>
      dplyr::mutate(Fluorophore = paste0(Fluorophore, "_",
        dplyr::row_number()))
    DetectorCols <- DetectorCols |> dplyr::slice(idx)

    row_scales <- runif(nrow(DetectorCols), min = var_range[1],
      max = var_range[2])
    channel_dev <- (var_range[2] - 1) * 0.3
    channel_noise <- matrix(
      runif(prod(dim(DetectorCols)), min = 1 - channel_dev,
        max = 1 + channel_dev),
      nrow = nrow(DetectorCols), ncol = ncol(DetectorCols)
    )
    DetectorCols <- DetectorCols * row_scales * channel_noise
  }

  if (Normalize) {
    if (any(DetectorCols > 1)) {
      message("Normalizing Data for Signature Comparison")
      n <- DetectorCols
      row_maxes <- apply(n, 1, max)
      DetectorCols <- n / row_maxes
    }
  }

  WhoseThis <- cbind(Identity, DetectorCols)
  TheseFluorophores <- WhoseThis |> dplyr::pull(Fluorophore)
  detector_names <- colnames(DetectorCols)

  WhoseThis1 <- WhoseThis |>
    tidyr::pivot_longer(cols = tidyselect::where(is.numeric),
      names_to = "Detector", values_to = "AdjustedY")

  amalgamate_input <- WhoseThis |> dplyr::mutate(Count = 1) |>
    dplyr::relocate(Count, .after = Fluorophore)

  if (returnType == "Cartoon") {
    ThePlot <- CartoonInternal(TheseFluorophores = TheseFluorophores,
      TheFluorophore = NULL, data = WhoseThis1, legend = legend,
      plotname = plotname, plotlinecolor = plotlinecolor, ylim = ylim)
  } else if (returnType == "Amalgamated") {
    line_col <- if (!is.null(plotlinecolor)) plotlinecolor else "red"
    ThePlot <- QC_Amalgamate(data = amalgamate_input,
      samplecolumn = "Fluorophore", countcolumn = "Count",
      normalize = FALSE, returnType = "plot", titlename = plotname,
      linecolor = line_col, legend = legend)
  } else if (returnType == "Density") {
    ThePlot <- CartoonDensity2D(data = amalgamate_input, theme = theme,
      bins = bins, plotname = plotname, logScale = logScale)
  } else if (returnType == "Ribbon") {

    plot_df <- WhoseThis1 |>
      dplyr::mutate(DetectorIdx = as.numeric(factor(Detector,
        levels = detector_names)))

    if (isTRUE(logScale)) {
      plot_df <- plot_df |> dplyr::filter(AdjustedY > 0)
    }

    num_detectors <- length(detector_names)

    # 1. Define quantiles to construct stacked, continuous ribbon shells
    n_shells <- 10
    prob_steps <- seq(0.01, 0.49, length.out = n_shells)

    # Calculate quantiles per detector
    quantiles_df <- plot_df |>
      dplyr::group_by(DetectorIdx) |>
      dplyr::summarise(
        mid = median(AdjustedY, na.rm = TRUE),
        dplyr::across(
          .cols = tidyselect::everything(),
          .fns = list(),
          .names = "temp"
        ),
        .groups = "drop"
      )

    # Build a interpolated data frame for each quantile pair
    interp_x <- seq(1, num_detectors, length.out = (num_detectors - 1) * 10 + 1)

    ribbon_layers <- list()
    for (s in seq_len(n_shells)) {
      p_low <- prob_steps[s]
      p_high <- 1 - p_low

      q_per_det <- plot_df |>
        dplyr::group_by(DetectorIdx) |>
        dplyr::summarise(
          y_low = quantile(AdjustedY, probs = p_low, na.rm = TRUE),
          y_high = quantile(AdjustedY, probs = p_high, na.rm = TRUE),
          .groups = "drop"
        )

      y_low_interp <- stats::approx(q_per_det$DetectorIdx, q_per_det$y_low,
        xout = interp_x)$y
      y_high_interp <- stats::approx(q_per_det$DetectorIdx, q_per_det$y_high,
        xout = interp_x)$y

      ribbon_layers[[s]] <- data.frame(
        x = interp_x,
        ymin = y_low_interp,
        ymax = y_high_interp,
        level = s
      )
    }

    full_ribbon_df <- dplyr::bind_rows(ribbon_layers)

    # Interpolated center line
    median_per_det <- plot_df |>
      dplyr::group_by(DetectorIdx) |>
      dplyr::summarise(mid = median(AdjustedY, na.rm = TRUE),
        .groups = "drop")

    median_interp <- data.frame(
      x = interp_x,
      y = stats::approx(median_per_det$DetectorIdx, median_per_det$mid,
        xout = interp_x)$y
    )

    ThePlot <- ggplot2::ggplot()

    # Layer outer to inner ribbons for smooth gradient shading
    palette_colors <- viridis::viridis(n_shells, option = theme)
    for (s in seq_len(n_shells)) {
      layer_data <- full_ribbon_df |> dplyr::filter(level == s)
      ThePlot <- ThePlot +
        ggplot2::geom_ribbon(
          data = layer_data,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          fill = palette_colors[s],
          alpha = 0.35
        )
    }

    ThePlot <- ThePlot +
      ggplot2::geom_line(
        data = median_interp,
        ggplot2::aes(x = x, y = y),
        color = if (!is.null(plotlinecolor)) plotlinecolor else "white",
        linewidth = 0.9
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(detector_names),
        labels = detector_names,
        expand = c(0, 0)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = if (!is.null(plotname)) plotname else
          "Fluorophore Density Ribbon",
        x = "Detector",
        y = "Intensity"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
        vjust = 0.5, hjust = 1))

    if (!legend) {
      ThePlot <- ThePlot + ggplot2::theme(legend.position = "none")
    }
    if (isTRUE(logScale)) {
      ThePlot <- ThePlot + ggplot2::scale_y_log10()
    }
    if (!is.null(ylim)) {
      ThePlot <- ThePlot + ggplot2::coord_cartesian(ylim = ylim)
    }
  } else {
    stop("Invalid returnType. Options are: 'Cartoon', 'Amalgamated',
      'Density', or 'Ribbon'.")
  }

  return(ThePlot)
}