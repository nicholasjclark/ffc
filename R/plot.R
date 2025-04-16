#' Plot time series of fts basis coefficients
#'
#' Produces a time series plot of basis function coefficients from `fts_ts` objects
#'
#' @importFrom ggplot2 autoplot ggplot aes facet_wrap geom_line ylab
#' @importFrom ggplot2 scale_colour_viridis_d labeller label_wrap_gen
#' @importFrom rlang sym expr eval_tidy
#' @param object An object of class `fts_ts` containing time-varying
#' basis function coefficients extracted from an `ffc_gam` object
#' @param ... Ignored
#' @author Nicholas J Clark
#' @export
autoplot.fts_ts <- function(object, ...) {
  n_basis <- length(unique(object$.basis))

  # Construct mapping aesthetics
  aes_spec <- list(
    x = sym(attr(object, "time_var")),
    y = sym(".estimate"),
    colour = sym(".basis")
  )

  if (!attr(object, "summarized")) {
    aes_spec["group"] <- list(sym(".realisation"))
  }

  # Improve printing of basis names
  object$.basis <- gsub(
    "fts_bs_",
    "basis ",
    object$.basis
  )
  object$.basis <- gsub(
    "fts_",
    "basis ",
    object$.basis
  )

  # Construct the plot
  p <- ggplot(
    object,
    eval_tidy(expr(aes(!!!aes_spec)))
  ) +
    ylab("Coefficient estimate")

  # Add lines
  if (!requireNamespace("ggborderline", quietly = TRUE)) {
    rlang::inform(
      message = paste0(
        'Package "ggborderline" can enable more readable time series plots\n',
        "Please consider installing it"
      ),
      .frequency = "once",
      .frequency_id = "ggborderline_autoplot.fts_ts"
    )

    p <- p +
      geom_line(show.legend = FALSE) +
      scale_colour_viridis_d(
        option = "C",
        end = 0.85
      )
  } else {
    p <- p +
      do.call(
        what = `::`,
        args = list("ggborderline", "geom_borderline")
      )(bordercolour = "white",
        linewidth = 0.75,
        show.legend = FALSE) +
      scale_colour_viridis_d(
        option = "C",
        end = 0.85
      )
  }

  # Add facets and labels
  if (n_basis == 1L) {
    p <- p +
      ggplot2::labs(title = unique(object$.basis))
  }

  if (n_basis > 1L) {
    p <- p +
      facet_wrap(~.basis,
        scales = "free_y",
        ncol = min(4, n_basis),
        labeller = labeller(.basis = label_wrap_gen(10))
      )
  }

  return(p)
}
