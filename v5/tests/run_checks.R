#!/usr/bin/env Rscript
# Synthetic regression checks; never reads the atlas RDS or starts a server.
# Run with: /home/minwook/miniconda3/envs/crc_shiny/bin/Rscript tests/run_checks.R
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "tests/run_checks.R"
app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(DT)
  library(plotly)
})

checks <- 0L
check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  checks <<- checks + 1L
}
same <- function(a, b) isTRUE(all.equal(unname(a), unname(b), check.attributes = FALSE))
numbers <- function(x) as.numeric(unlist(x, use.names = FALSE))
strings <- function(x) as.character(unlist(x, use.names = FALSE))

# Read only definitions after startup; evaluating these does not load data.
app_code <- parse(file.path(app_dir, "app.R"))
assignment_name <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]))
    as.character(expr[[2]]) else ""
}
assignment_names <- vapply(app_code, assignment_name, character(1L))
plot_env <- new.env(parent = globalenv())
for (i in which(assignment_names %in% c("PAL", "pal_for", "COND_PAL", "cond_color")))
  eval(app_code[[i]], plot_env)
sys.source(file.path(app_dir, "R", "umap.R"), envir = plot_env)
sys.source(file.path(app_dir, "R", "violin.R"), envir = plot_env)
sys.source(file.path(app_dir, "R", "llm.R"), envir = plot_env)

df <- data.frame(
  cell_id = paste0("cell", seq_len(12)),
  UMAP_1 = seq_len(12), UMAP_2 = rep(c(-1, 0, 1), 4),
  Condition = rep(c("Normal", "Tumor"), each = 6),
  Major_celltype = rep(rep(c("B cells", "T cells"), each = 3), 2),
  Subtype = rep(c("subA", "subB"), 6),
  Sample = rep(c("N1", "T1"), each = 6),
  Patient = "P1", nCount_RNA = 100, nFeature_RNA = 20, percent.mt = 1,
  .expr = rep(c(0, 0.5, 1, 1.5, 2, 3), 2),
  stringsAsFactors = FALSE
)

# Check the serialized Plotly figure consumed by the browser, not snapshots.
figure <- function(widget) {
  jsonlite::fromJSON(
    plotly::plotly_json(plotly::plotly_build(widget), jsonedit = FALSE, pretty = FALSE),
    simplifyVector = FALSE
  )
}
umap <- function(data = df, threshold = 1, compare = TRUE, color_mode = "expr") {
  figure(plot_env$build_umap_plot(data, "Major_celltype", color_mode, compare, threshold))
}
violin <- function(data = df, threshold = 1, compare = TRUE, group_var = "Major_celltype") {
  figure(plot_env$build_violin_plot(data, group_var, compare, threshold))
}
point_count <- function(trace) sum(is.finite(numbers(trace$x)) & is.finite(numbers(trace$y)))
is_scatter <- function(trace) trace$type %in% c("scatter", "scattergl")
annotation_texts <- function(fig) vapply(fig$layout$annotations, function(a) a$text, character(1L))
color_traces <- function(fig) Filter(function(t) !is.null(t$marker$coloraxis), fig$data)
legend_traces <- function(fig) Filter(function(t) isTRUE(t$showlegend), fig$data)

expression <- umap()
condition_titles <- Filter(function(a) a$text %in% c("Normal", "Tumor"), expression$layout$annotations)
check(length(condition_titles) == 2L, "Expression UMAP must label both conditions.")
check(length(unique(vapply(condition_titles, function(a) a$x, numeric(1L)))) == 2L,
      "Condition labels must occupy separate panel positions.")
bars <- expression$layout[grep("^coloraxis[0-9]*$", names(expression$layout))]
check(length(bars) == 1L && isTRUE(bars[[1]]$showscale),
      "Expression panels must share exactly one colorbar.")
limits <- color_traces(expression)
check(length(limits) >= 2L && all(vapply(limits, function(t) identical(t$marker$coloraxis, "coloraxis"), logical(1L))) &&
        same(expression$layout$coloraxis$cmin, 0) && same(expression$layout$coloraxis$cmax, max(df$.expr)) &&
        !any(vapply(expression$data, function(t) isTRUE(t$marker$showscale), logical(1L))),
      "Expression colors must use common bounds across conditions.")
check(all(vapply(Filter(function(t) is_scatter(t) && point_count(t) > 0, expression$data),
                 function(t) same(numbers(t$marker$size), 3), logical(1L))),
      "UMAP data markers must retain their original size.")

# Gray traces use a literal color; expression traces carry numeric colors.
gray_count <- function(fig) {
  sum(vapply(Filter(function(t) {
    is_scatter(t) && is.character(t$marker$color) &&
      length(t$marker$color) == 1L && !isTRUE(t$showlegend)
  }, fig$data), point_count, numeric(1L)))
}
check(gray_count(umap(threshold = 0)) == 2, "Zero expression must be gray at threshold zero.")
check(gray_count(expression) == 6, "Cells exactly at the threshold must also be gray.")
check(sum(vapply(expression$data, point_count, numeric(1L))) == nrow(df),
      "Thresholding must preserve all UMAP cells.")

groups <- umap(color_mode = "group")
legends <- legend_traces(groups)
check(length(legends) == 2L && length(unique(vapply(legends, function(t) t$name, character(1L)))) == 2L,
      "Group UMAP must show one legend entry per cell type across panels.")
check(all(vapply(legends, function(t) same(numbers(t$marker$size), 6), logical(1L))),
      "Cell-type legend bullets must be twice the size of original markers.")
check(all(vapply(Filter(function(t) is_scatter(t) && point_count(t) > 0, groups$data),
                 function(t) same(numbers(t$marker$size), 3), logical(1L))),
      "Enlarging legends must not enlarge UMAP data markers.")
check(all(vapply(legends, function(leg) {
  linked <- Filter(function(t) identical(t$legendgroup, leg$legendgroup) && point_count(t) > 0, groups$data)
  panels <- unique(vapply(linked, function(t) if (is.null(t$xaxis)) "x" else t$xaxis, character(1L)))
  length(panels) == 2L
}, logical(1L))) && identical(groups$layout$legend$groupclick, "togglegroup"),
"Each legend entry must toggle that cell type in both panels.")

for (edge in list(transform(df, .expr = 0), df, df[df$Condition == "Normal", ])) {
  fig <- umap(edge, threshold = 10)
  check(gray_count(fig) == nrow(edge), "All cells must remain visible when no expression exceeds the threshold.")
  check(length(Filter(function(t) isTRUE(t$marker$showscale), fig$data)) <= 1L,
        "Empty positive-expression panels must not duplicate colorbars.")
}

# Plotly needs marker.color to remain a numeric JSON array even for one point.
# A singleton scalar can silently remove the expression colorbar in the browser.
one_positive <- transform(df, .expr = 0)
one_positive$.expr[6] <- 2
one_per_panel <- one_positive
one_per_panel$.expr[12] <- 3
for (edge in list(one_positive, one_per_panel, df[6, , drop = FALSE])) {
  fig <- umap(edge, threshold = 1)
  colored <- color_traces(fig)
  check(length(colored) == sum(edge$.expr > 1) && all(vapply(colored, function(t) {
    is.list(t$marker$color) && length(t$marker$color) == 1L &&
      is.numeric(t$marker$color[[1]]) && point_count(t) == 1
  }, logical(1L))), "Singleton expression colors must serialize as numeric arrays, not scalars.")
  check(isTRUE(fig$layout$coloraxis$showscale) &&
          length(grep("^coloraxis[0-9]*$", names(fig$layout))) == 1L,
        "Singleton expression traces must retain one shared colorbar.")
  check(gray_count(fig) == sum(edge$.expr <= 1) &&
          sum(vapply(fig$data, point_count, numeric(1L))) == nrow(edge),
        "Singleton positive panels must preserve every positive and gray cell.")
}

split <- violin()
halves <- Filter(function(t) identical(t$type, "violin"), split$data)
check(length(halves) == 2L && setequal(vapply(halves, function(t) t$side, character(1L)), c("negative", "positive")),
      "Two-condition comparison must retain split violin halves.")
negative <- Filter(function(t) identical(t$side, "negative"), halves)[[1]]
positive <- Filter(function(t) identical(t$side, "positive"), halves)[[1]]
centers_left <- sort(unique(numbers(negative$x)))
centers_right <- sort(unique(numbers(positive$x)))
ticks <- numbers(split$layout$xaxis$tickvals)
check(all(centers_right > centers_left) && all(centers_right - centers_left < 0.3),
      "Normal and Tumor violin halves need a small visible horizontal gap.")
check(same((centers_right + centers_left) / 2, ticks),
      "Group labels must remain centered between the violin halves.")
check(setequal(strings(split$layout$xaxis$ticktext), unique(df$Major_celltype)),
      "Numeric violin positions must retain readable cell-type labels.")
check(sum(vapply(halves, function(t) length(t$y), integer(1L))) == nrow(df),
      "The threshold must not filter the violin distribution.")

# v5: median + rank-sum p-value annotations, one per group, above the split violins.
`%||%` <- function(a, b) if (is.null(a)) b else a
stat_ann <- Filter(function(a) grepl("^med ", a$text %||% "") && grepl("p=|p<", a$text),
                   split$layout$annotations)
check(length(stat_ann) == length(unique(df$Major_celltype)),
      "Each split-violin group needs a median + rank-sum p-value annotation.")
check(all(vapply(stat_ann, function(a) identical(a$xref, "x") && identical(a$yref, "paper"),
                 logical(1L))),
      "Violin stat annotations must sit above each group on the shared x axis.")

extra <- df[df$Condition == "Normal", ]
extra$Condition <- "Adjacent"
extra$cell_id <- paste0("extra", seq_len(nrow(extra)))
df_three <- rbind(df, extra)
for (case in list(list(data = df, compare = TRUE), list(data = df, compare = FALSE),
                 list(data = df, group_var = "Condition"), list(data = df_three),
                 list(data = transform(df, .expr = 0)), list(data = df[df$Condition == "Normal", ]))) {
  threshold <- 4
  fig <- do.call(violin, c(case, list(threshold = threshold)))
  lines <- Filter(function(s) identical(s$type, "line") && same(s$y0, threshold) && same(s$y1, threshold),
                  fig$layout$shapes)
  panels <- unique(vapply(Filter(function(t) identical(t$type, "violin"), fig$data),
                          function(t) if (is.null(t$xaxis)) "x" else t$xaxis, character(1L)))
  check(length(lines) == length(panels), "Every violin panel must include the threshold line.")
  check(all(vapply(lines, function(s) grepl("domain$", s$xref) && identical(s$line$dash, "dash"), logical(1L))),
        "Threshold lines must span the panel width and use a dashed style.")
  yaxes <- fig$layout[grep("^yaxis[0-9]*$", names(fig$layout))]
  check(all(vapply(yaxes, function(a) {
    range <- numbers(a$range)
    length(range) == 2L && range[1] < threshold && range[2] > threshold
  }, logical(1L))), "Threshold must remain within the visible violin y range, even for an all-zero gene.")
}

# Exercise the actual reactive server with a tiny in-memory matrix. Startup
# expressions (including readRDS) and the final shinyApp() call are not evaluated.
app_env <- new.env(parent = globalenv())
app_env$app_dir <- app_dir
app_env$plot_df_base <- df[, setdiff(names(df), ".expr"), drop = FALSE]
app_env$expr_mat <- rbind(SNHG16 = df$.expr, ZERO = rep(0, nrow(df)))
colnames(app_env$expr_mat) <- df$cell_id
start <- match("available_genes", assignment_names)
stopifnot(!is.na(start))
for (expr in app_code[seq.int(start, length(app_code))]) {
  if (is.call(expr) && identical(expr[[1]], as.name("shinyApp"))) next
  eval(expr, app_env)
}
shiny::testServer(app_env$server, {
  session$setInputs(gene = "SNHG16", f_condition = c("Normal", "Tumor"),
                    f_celltype = c("B cells", "T cells"), f_subtype = character(),
                    f_sample = character(), split_by = "Major_celltype", minexpr = 0,
                    compare = TRUE, umap_color = "expr")
  cells_zero <- cells()
  summary_zero <- summary_df()
  session$setInputs(minexpr = 1)
  cells_one <- cells()
  summary_one <- summary_df()
  check(nrow(cells_one) == 12L && identical(cells_one$cell_id, cells_zero$cell_id),
        "Moving the threshold must preserve filtered cells and their denominator.")
  check(sum(cells_zero$.positive) == 10L && sum(cells_one$.positive) == 6L &&
          all(cells_one$.positive == (cells_one$.expr > 1)),
        "The slider must invalidate cells() and recompute a strict positive mask.")
  check(sum(summary_zero$n_cells) == 12L && sum(summary_one$n_cells) == 12L &&
          sum(summary_zero$n_positive) == 10L && sum(summary_one$n_positive) == 6L,
        "Summary positivity must update without changing the cell denominator.")
  check(all(summary_one$threshold == 1) &&
          same(summary_one$pct_pos, round(100 * summary_one$n_positive / summary_one$n_cells, 1)),
        "Summary percentages and exported thresholds must agree.")
  # v5: per-group Tumor/Normal breakdown + rank-sum p-value.
  new_cols <- c("n_Normal", "mean_Normal", "median_Normal",
                "n_Tumor", "mean_Tumor", "median_Tumor", "p_value")
  check(all(new_cols %in% names(summary_zero)),
        "Summary must report Tumor/Normal mean, median and a p-value per group.")
  check(all(summary_zero$n_Normal + summary_zero$n_Tumor == summary_zero$n_cells),
        "Per-condition counts must add up to each group's cell count.")
  check(is.numeric(summary_zero$p_value),
        "The rank-sum p-value column must be numeric (NA when a side is too small).")
})

# v5: natural-language parser validates against the real vocabulary.
nl <- plot_env$parse_nl_response(
  'ok {"gene":"snhg16","split_by":"Major_celltype","umap_color":"group","compare":true,"conditions":["Tumor","Normal"],"celltypes":["T cells","Made up"],"minexpr":0.5}',
  genes = c("SNHG16", "GAS5"), celltypes = c("T cells", "B cells"),
  conditions = c("Normal", "Tumor"), split_axes = c("Major_celltype", "Subtype", "Sample", "Patient", "Condition"))
check(identical(nl$settings$gene, "SNHG16"), "NL parser must resolve gene case-insensitively.")
check(identical(nl$settings$celltypes, "T cells"), "NL parser must drop invalid cell types.")
check(setequal(nl$settings$conditions, c("Normal", "Tumor")) &&
        isTRUE(nl$settings$compare) && identical(nl$settings$umap_color, "group") &&
        same(nl$settings$minexpr, 0.5), "NL parser must apply valid settings.")
check(isFALSE(plot_env$parse_nl_response("no json here", character(), character(),
                                         character(), character())$ok),
      "NL parser must fail cleanly when no JSON is present.")

cat(sprintf("PASS: %d synthetic regression checks; no RDS loaded and no service started.\n", checks))
