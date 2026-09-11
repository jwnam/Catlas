# UMAP rendering is kept separate from app startup so it can be checked without
# loading the Seurat object. The caller provides the already filtered cells.
build_umap_plot <- function(df, group_var, color_mode, compare, threshold) {
  stopifnot(is.data.frame(df), length(threshold) == 1L, is.finite(threshold),
            color_mode %in% c("expr", "group"))

  group_labels <- as.character(df[[group_var]])
  group_labels[is.na(group_labels) | !nzchar(group_labels)] <- "Unknown"
  condition_labels <- as.character(df$Condition)
  condition_labels[is.na(condition_labels) | !nzchar(condition_labels)] <- "Unknown"
  df$.g <- group_labels
  df$.condition <- condition_labels
  df$.positive <- is.finite(df$.expr) & df$.expr > threshold
  groups <- sort(unique(group_labels))
  group_colors <- if (identical(group_var, "Condition")) cond_color(groups) else pal_for(groups)

  conditions <- sort(unique(condition_labels))
  split_panels <- isTRUE(compare) && length(conditions) > 1L
  panels <- if (split_panels) {
    lapply(conditions, function(cc) df[df$.condition == cc, , drop = FALSE])
  } else list(df)
  panel_titles <- if (split_panels) conditions else if (length(conditions) == 1L) conditions else "All conditions"

  # Every panel uses identical bounds, including when a condition has no
  # positive cells. Padding also handles an embedding with one unique value.
  padded_range <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(-1, 1))
    bounds <- range(x)
    padding <- max(diff(bounds) * 0.04, 0.1)
    bounds + c(-padding, padding)
  }
  x_range <- padded_range(df$UMAP_1)
  y_range <- padded_range(df$UMAP_2)
  n_panels <- length(panels)
  panel_gap <- if (n_panels > 1L) min(0.07, 0.25 / n_panels) else 0
  panel_width <- (1 - panel_gap * (n_panels - 1L)) / n_panels
  has_positive <- any(df$.positive)
  expression_values <- df$.expr[is.finite(df$.expr)]
  color_min <- min(c(0, expression_values))
  color_max <- max(c(0, expression_values))
  if (color_max <= color_min) color_max <- color_min + 1

  p <- plotly::plot_ly()
  plot_layout <- list(
    showlegend = identical(color_mode, "group"),
    margin = list(l = 58, r = if (color_mode == "expr") 100 else 170, t = 66, b = 52),
    hovermode = "closest",
    legend = list(x = 1.02, y = 1, xanchor = "left", yanchor = "top",
                  itemsizing = "trace", groupclick = "togglegroup", tracegroupgap = 0)
  )
  annotations <- list()

  escape_hover <- function(x) as.character(htmltools::htmlEscape(as.character(x)))
  hover_text <- function(d) {
    paste0(escape_hover(d$Major_celltype), " / ", escape_hover(d$Subtype),
           "<br>", escape_hover(d$Sample), " (", escape_hover(d$Condition), ")",
           "<br>expr: ", format(round(d$.expr, 3), trim = TRUE),
           "<br>", ifelse(d$.positive, "Above threshold", "At/below threshold"),
           " (", format(threshold, trim = TRUE), ")")
  }

  # Numeric marker.color vectors refer to ONE layout coloraxis. Per-trace
  # colorbars are disabled, so side-by-side panels cannot overlap their scales.
  if (color_mode == "expr" && has_positive) {
    plot_layout$coloraxis <- list(
      colorscale = "Viridis", cauto = FALSE, cmin = color_min, cmax = color_max,
      showscale = TRUE,
      colorbar = list(title = list(text = "log-norm", side = "top"),
                      x = 1.02, xanchor = "left", y = 0.5, len = 0.9, thickness = 16)
    )
  }

  for (i in seq_along(panels)) {
    d <- panels[[i]]
    suffix <- if (i == 1L) "" else as.character(i)
    x_axis <- paste0("x", suffix)
    y_axis <- paste0("y", suffix)
    domain_start <- (i - 1L) * (panel_width + panel_gap)
    domain <- c(domain_start, domain_start + panel_width)
    plot_layout[[paste0("xaxis", suffix)]] <- list(
      title = "UMAP_1", domain = domain, anchor = y_axis,
      range = x_range, zeroline = FALSE, matches = if (i > 1L) "x" else NULL
    )
    plot_layout[[paste0("yaxis", suffix)]] <- list(
      title = if (i == 1L) "UMAP_2" else "", domain = c(0, 1), anchor = x_axis,
      range = y_range, zeroline = FALSE, showticklabels = i == 1L,
      matches = if (i > 1L) "y" else NULL
    )
    annotations[[length(annotations) + 1L]] <- list(
      text = escape_hover(panel_titles[[i]]), x = mean(domain), y = 1.055,
      xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom",
      showarrow = FALSE, font = list(size = 14)
    )

    add_cells <- function(p, cells, marker, name, legendgroup = NULL) {
      if (!nrow(cells)) return(p)
      plotly::add_trace(
        p, x = cells$UMAP_1, y = cells$UMAP_2, type = "scattergl", mode = "markers",
        marker = marker, text = hover_text(cells), hoverinfo = "text",
        name = name, legendgroup = legendgroup, showlegend = FALSE,
        xaxis = x_axis, yaxis = y_axis, inherit = FALSE
      )
    }

    if (color_mode == "expr") {
      p <- add_cells(p, d[!d$.positive, , drop = FALSE],
                     list(size = 3, color = "#C9CDD3", showscale = FALSE),
                     "At/below threshold")
      positive <- d[d$.positive, , drop = FALSE]
      p <- add_cells(p, positive,
                     # Keep a JSON array even when only one cell is positive:
                     # a scalar numeric color is not a continuous color vector.
                     list(size = 3, color = I(unname(positive$.expr)), coloraxis = "coloraxis",
                          showscale = FALSE), "Expression")
    } else {
      # Add all gray cells first, then colored cells, so the latter remain
      # visible. The legend group includes both traces in every condition.
      for (positive in c(FALSE, TRUE)) {
        for (g in groups) {
          group_id <- paste0("group_", match(g, groups))
          selected <- d[d$.g == g & d$.positive == positive, , drop = FALSE]
          marker_color <- if (positive) unname(group_colors[[g]]) else "#C9CDD3"
          p <- add_cells(p, selected,
                         list(size = 3, color = marker_color, showscale = FALSE),
                         g, group_id)
        }
      }
    }

    if (!any(d$.positive)) {
      annotations[[length(annotations) + 1L]] <- list(
        text = if (nrow(d)) "No cells above threshold" else "No cells match the current filters",
        x = mean(domain), y = 0.98, xref = "paper", yref = "paper",
        xanchor = "center", yanchor = "top", showarrow = FALSE,
        bgcolor = "rgba(255,255,255,0.85)", font = list(size = 11, color = "#666666")
      )
    }
  }

  if (color_mode == "group") {
    # A null-position trace supplies only the legend symbol. Its size is 6,
    # twice the actual cell marker size of 3, without enlarging UMAP points.
    # Matching legendgroup values keep toggles linked across both panels.
    for (g in groups) {
      p <- plotly::add_trace(
        # Lists retain JSON null positions while preventing the R Plotly
        # builder from dropping this trace's name and legend group as NA rows.
        p, x = list(NA_real_), y = list(NA_real_), type = "scatter", mode = "markers",
        marker = list(size = 6, color = unname(group_colors[[g]])),
        name = g, legendgroup = paste0("group_", match(g, groups)),
        showlegend = TRUE, hoverinfo = "skip", xaxis = "x", yaxis = "y",
        inherit = FALSE
      )
    }
  }
  if (!nrow(df)) {
    p <- plotly::add_trace(p, x = list(NA_real_), y = list(NA_real_),
                          type = "scatter", mode = "markers",
                          marker = list(size = 3), showlegend = FALSE,
                          hoverinfo = "skip", inherit = FALSE)
  }
  plot_layout$annotations <- annotations
  p <- do.call(plotly::layout, c(list(p = p), plot_layout))
  plotly::config(p, displaylogo = FALSE, responsive = TRUE)
}
