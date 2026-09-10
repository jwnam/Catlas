# CRC scRNA-seq Atlas Viewer (v4) — interactive lncRNA explorer
#
# Successor to v3. Same read-only RDS, same startup contract
# (start_catlas_v2.sh with CATLAS_APP_DIR pointing at v4), same data access.
#
# Changes vs v3
#   1. UMAP expression compare: single shared colorbar; Normal (left) / Tumor
#      (right) always drawn; larger square legend swatches in Group mode;
#      sub-type colors share the parent cell type's hue family.
#   2. "Split by" moved directly above the violin (it drives the violin/summary).
#      Violin always shows Normal (blue) vs Tumor (red) dodged with a small gap,
#      median + Wilcoxon rank-sum p-value annotated on top of each group.
#   3. Per-group summary reports Tumor/Normal mean, median and the p-value.
#   4. Min-expression filter now actually filters the cells used everywhere.
#   5. Natural-language search powered by a user-supplied LLM API key
#      (Anthropic Claude / OpenAI ChatGPT / Google Gemini) — configured in-app,
#      keys kept in the session only (never stored or logged).
#
# New hard dependency vs v2: plotly.  Optional (only for NL search): httr, jsonlite.

# --------------------------------------------------------------------------- #
#  Startup checks                                                              #
# --------------------------------------------------------------------------- #
app_dir_value <- Sys.getenv("CATLAS_APP_DIR", unset = "")
if (!nzchar(app_dir_value)) app_dir_value <- getwd()
app_dir <- normalizePath(app_dir_value, winslash = "/", mustWork = TRUE)

project_dir_value <- Sys.getenv("CATLAS_PROJECT_DIR", unset = file.path(app_dir, ".."))
project_dir <- normalizePath(project_dir_value, winslash = "/", mustWork = TRUE)

data_path_value <- Sys.getenv(
  "CATLAS_DATA_PATH",
  unset = file.path(project_dir, "crc_shiny_app_seurat.rds")
)
if (!file.exists(data_path_value)) stop("CRC Atlas input RDS was not found: ", data_path_value)
data_path <- normalizePath(data_path_value, winslash = "/", mustWork = TRUE)

required_packages <- c("shiny", "Seurat", "Matrix", "ggplot2", "dplyr", "DT", "plotly")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages) > 0L)
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "),
       ". Install with install.packages() in the crc_shiny environment.")

suppressPackageStartupMessages({
  library(shiny); library(Seurat); library(Matrix); library(ggplot2)
  library(dplyr); library(DT); library(plotly)
})

# --------------------------------------------------------------------------- #
#  Load data (same access pattern as v2/v3)                                   #
# --------------------------------------------------------------------------- #
message("CRC Atlas v4 loading: ", data_path)
obj <- readRDS(data_path)
DefaultAssay(obj) <- "RNA"
expr_mat <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                     error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))

umap_df <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- rownames(umap_df)

meta_df <- obj@meta.data
meta_df$cell_id        <- rownames(meta_df)
meta_df$Condition      <- as.character(meta_df$Condition)      # Tumor / Normal
meta_df$Sample         <- as.character(meta_df$Library)
meta_df$Patient        <- as.character(meta_df$Patient)
meta_df$Major_celltype <- as.character(meta_df$Cell_type_v2)   # incl. Malignant / Epithelial
meta_df$Subtype        <- as.character(meta_df$Cell_subtype_v2)

plot_df_base <- dplyr::left_join(umap_df, meta_df, by = "cell_id")
rm(obj); invisible(gc())

available_genes <- rownames(expr_mat)
gene_lower <- tolower(available_genes)

# Axes: UMAP can colour by any; the violin/summary "Split by" excludes Condition
# (Tumor vs Normal is always the within-group split of the violin itself).
UMAP_AXES   <- c("Cell type" = "Major_celltype", "Sub-type" = "Subtype",
                 "Sample" = "Sample", "Condition" = "Condition", "Patient" = "Patient")
SPLIT_AXES  <- c("Cell type" = "Major_celltype", "Sub-type" = "Subtype",
                 "Sample" = "Sample", "Patient" = "Patient")

conditions_all <- sort(unique(plot_df_base$Condition))
celltypes_all  <- sort(unique(plot_df_base$Major_celltype))
samples_all    <- sort(unique(plot_df_base$Sample))

# Show Normal on the left, Tumor on the right; others after.
cond_order <- function(x) {
  x <- unique(as.character(x)); c(intersect(c("Normal", "Tumor"), x),
                                  setdiff(sort(x), c("Normal", "Tumor")))
}

# --------------------------------------------------------------------------- #
#  Colour system: cell types get distinct hues; sub-types inherit the parent   #
#  cell type's hue family (light -> dark shades).                              #
# --------------------------------------------------------------------------- #
PAL <- c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860",
         "#DA8BC3","#8C8C8C","#CCB974","#64B5CD","#E377C2","#7F7F7F",
         "#BCBD22","#17BECF","#AEC7E8","#FFBB78","#98DF8A","#FF9896")
pal_for <- function(levels) setNames(rep(PAL, length.out = length(levels)), levels)

mix_col <- function(col, with, p)
  grDevices::rgb(t((1 - p) * grDevices::col2rgb(col) + p * grDevices::col2rgb(with)) / 255)
shade_palette <- function(base, n) {
  if (n <= 1) return(base)
  grDevices::colorRampPalette(
    c(mix_col(base, "white", 0.55), base, mix_col(base, "black", 0.35)))(n)
}

CELLTYPE_COLORS <- setNames(rep(PAL, length.out = length(celltypes_all)), celltypes_all)
SUBTYPE_COLORS  <- unlist(lapply(celltypes_all, function(ct) {
  subs <- sort(unique(plot_df_base$Subtype[plot_df_base$Major_celltype == ct]))
  subs <- subs[!is.na(subs)]
  if (!length(subs)) return(NULL)
  setNames(shade_palette(CELLTYPE_COLORS[[ct]], length(subs)), subs)
}))
COND_PAL <- c(Normal = "#4C72B0", Tumor = "#C44E52")
cond_color <- function(levels) {
  out <- COND_PAL[levels]; miss <- is.na(out)
  if (any(miss)) out[miss] <- PAL[seq_len(sum(miss))]
  setNames(unname(out), levels)
}

# Colour lookup for a grouping axis, always falling back to PAL for unknowns.
colors_for <- function(gv, levels) {
  base <- switch(gv,
    Major_celltype = CELLTYPE_COLORS[levels],
    Subtype        = SUBTYPE_COLORS[levels],
    Condition      = cond_color(levels),
    pal_for(levels))
  base <- unname(base); miss <- is.na(base)
  if (any(miss)) base[miss] <- pal_for(levels)[miss]
  setNames(base, levels)
}

fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("p=NA")
  if (p < 2.2e-16) return("p<2.2e-16")
  paste0("p=", signif(p, 2))
}

# --------------------------------------------------------------------------- #
#  LLM providers for natural-language search (optional)                        #
# --------------------------------------------------------------------------- #
LLM_PROVIDERS <- c("Anthropic (Claude)" = "anthropic",
                   "OpenAI (ChatGPT)"   = "openai",
                   "Google (Gemini)"    = "gemini")
LLM_MODELS <- list(
  anthropic = c("claude-opus-4-8", "claude-sonnet-4-6", "claude-haiku-4-5-20251001"),
  openai    = c("gpt-4o", "gpt-4o-mini", "o4-mini"),
  gemini    = c("gemini-2.0-flash", "gemini-1.5-pro", "gemini-1.5-flash"))

# Each returns the model's raw text answer, or throws with a readable message.
llm_call <- function(provider, model, key, prompt) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE))
    stop("Natural-language search needs the 'httr' and 'jsonlite' packages. ",
         "Install them in the crc_shiny environment.")
  if (!nzchar(key)) stop("No API key set. Open Settings and paste your key.")
  if (provider == "anthropic") {
    r <- httr::POST("https://api.anthropic.com/v1/messages",
      httr::add_headers(`x-api-key` = key, `anthropic-version` = "2023-06-01",
                        `content-type` = "application/json"),
      body = jsonlite::toJSON(list(model = model, max_tokens = 512,
        messages = list(list(role = "user", content = prompt))), auto_unbox = TRUE),
      encode = "raw")
    ct <- httr::content(r, as = "parsed", type = "application/json")
    if (!is.null(ct$error)) stop(ct$error$message %||% "Anthropic API error")
    ct$content[[1]]$text
  } else if (provider == "openai") {
    r <- httr::POST("https://api.openai.com/v1/chat/completions",
      httr::add_headers(Authorization = paste("Bearer", key),
                        `content-type` = "application/json"),
      body = jsonlite::toJSON(list(model = model,
        messages = list(list(role = "user", content = prompt))), auto_unbox = TRUE),
      encode = "raw")
    ct <- httr::content(r, as = "parsed", type = "application/json")
    if (!is.null(ct$error)) stop(ct$error$message %||% "OpenAI API error")
    ct$choices[[1]]$message$content
  } else if (provider == "gemini") {
    url <- paste0("https://generativelanguage.googleapis.com/v1beta/models/",
                  model, ":generateContent?key=", key)
    r <- httr::POST(url, httr::add_headers(`content-type` = "application/json"),
      body = jsonlite::toJSON(list(contents = list(list(
        parts = list(list(text = prompt))))), auto_unbox = TRUE), encode = "raw")
    ct <- httr::content(r, as = "parsed", type = "application/json")
    if (!is.null(ct$error)) stop(ct$error$message %||% "Gemini API error")
    ct$candidates[[1]]$content$parts[[1]]$text
  } else stop("Unknown provider: ", provider)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# --------------------------------------------------------------------------- #
#  UI                                                                          #
# --------------------------------------------------------------------------- #
ui <- fluidPage(
  tags$head(tags$style(HTML(
    ".ctrl-bar{background:#f5f7fa;border:1px solid #e3e8ee;border-radius:8px;padding:8px 12px;margin-bottom:10px;}
     .gene-note{font-size:12px;color:#666;margin:-6px 0 8px;}
     .nl-box{background:#eef4fb;border:1px solid #cfe0f3;border-radius:8px;padding:8px;margin-bottom:8px;}
     .split-bar{background:#fbf6ef;border:1px solid #eadfce;border-radius:8px;padding:6px 10px;margin:6px 0;}"))),
  titlePanel("CATLAS · CRC lncRNA Atlas (v4)"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "nl-box",
        tags$b("🔎 Natural-language search"),
        textInput("nl_query", NULL, placeholder = "e.g. SNHG16 in T cells, tumor vs normal"),
        fluidRow(
          column(7, actionButton("nl_go", "Search", class = "btn-sm btn-primary", width = "100%")),
          column(5, actionLink("nl_settings", "⚙ API settings"))),
        uiOutput("nl_status")),
      selectizeInput("gene", "Gene (lncRNA / any feature)", choices = NULL, selected = NULL,
                     options = list(placeholder = "e.g. SNHG16, NORAD, LINC...", maxOptions = 100)),
      uiOutput("gene_note"),
      tags$hr(),
      tags$b("Filters"),
      checkboxGroupInput("f_condition", "Condition", choices = conditions_all,
                         selected = conditions_all, inline = TRUE),
      selectizeInput("f_celltype", "Cell type", choices = celltypes_all, selected = celltypes_all,
                     multiple = TRUE, options = list(plugins = list("remove_button"))),
      selectizeInput("f_subtype", "Sub-type (optional)", choices = NULL, selected = NULL,
                     multiple = TRUE, options = list(plugins = list("remove_button"),
                                                     placeholder = "all sub-types")),
      selectizeInput("f_sample", "Sample (optional)", choices = samples_all, selected = NULL,
                     multiple = TRUE, options = list(plugins = list("remove_button"),
                                                     placeholder = "all samples")),
      sliderInput("minexpr", "Min expression filter (log-norm)", min = 0, max = 3,
                  value = 0, step = 0.1),
      helpText("Cells with expression below this value are excluded from all views."),
      tags$hr(),
      downloadButton("dl_csv", "Summary (CSV)", class = "btn-sm btn-outline-primary")
    ),

    mainPanel(
      width = 9,
      div(class = "ctrl-bar",
        fluidRow(
          column(5, radioButtons("umap_color", "UMAP color",
                                 choices = c("Expression" = "expr", "Group" = "group"),
                                 selected = "expr", inline = TRUE)),
          column(4, conditionalPanel("input.umap_color == 'group'",
                     selectInput("umap_group", "Color by", choices = UMAP_AXES,
                                 selected = "Major_celltype"))),
          column(3, checkboxInput("compare", "UMAP: Normal | Tumor side-by-side", TRUE)))),
      tabsetPanel(
        tabPanel("Gene expression",
          br(),
          h5(textOutput("umap_title")),
          plotlyOutput("umap", height = "440px"),
          div(class = "split-bar",
            fluidRow(
              column(8, radioButtons("split_by", "Split violin by",
                                     choices = SPLIT_AXES, selected = "Major_celltype", inline = TRUE)),
              column(4, div(style = "padding-top:6px;color:#666;font-size:12px;",
                            "Blue = Normal · Red = Tumor · p = Wilcoxon rank-sum")))),
          h5(textOutput("violin_title")),
          plotlyOutput("violin", height = "420px"),
          br(),
          h5("Per-group summary (Tumor vs Normal)"),
          DTOutput("summary_tbl")),
        tabPanel("Composition", br(), plotlyOutput("composition", height = "600px")),
        tabPanel("Metadata", br(), DTOutput("metadata_tbl")))
    )
  )
)

# --------------------------------------------------------------------------- #
#  Server                                                                      #
# --------------------------------------------------------------------------- #
server <- function(input, output, session) {

  updateSelectizeInput(session, "gene", choices = available_genes,
                       selected = if ("SNHG16" %in% available_genes) "SNHG16" else available_genes[1],
                       server = TRUE)

  observeEvent(input$f_celltype, {
    subs <- plot_df_base %>% filter(Major_celltype %in% input$f_celltype) %>%
      pull(Subtype) %>% unique() %>% sort()
    updateSelectizeInput(session, "f_subtype", choices = subs,
                         selected = intersect(input$f_subtype, subs))
  }, ignoreNULL = FALSE)

  # ---- cells passing all filters (incl. min-expression) --------------------
  cells <- reactive({
    req(input$gene)
    validate(need(input$gene %in% available_genes, "Gene not found."))
    df <- plot_df_base
    if (length(input$f_condition)) df <- df[df$Condition %in% input$f_condition, , drop = FALSE]
    if (length(input$f_celltype))  df <- df[df$Major_celltype %in% input$f_celltype, , drop = FALSE]
    if (length(input$f_subtype))   df <- df[df$Subtype %in% input$f_subtype, , drop = FALSE]
    if (length(input$f_sample))    df <- df[df$Sample %in% input$f_sample, , drop = FALSE]
    df$.expr <- as.numeric(expr_mat[input$gene, df$cell_id])
    if (input$minexpr > 0) df <- df[df$.expr >= input$minexpr, , drop = FALSE]  # (4) real filter
    validate(need(nrow(df) > 0, "No cells match the current filters."))
    df
  })

  split_label <- reactive(names(SPLIT_AXES)[SPLIT_AXES == input$split_by])

  output$gene_note <- renderUI(div(class = "gene-note", {
    df <- cells()
    sprintf("%s — %s cells shown · mean %.2f · %.1f%% > 0 (log-norm)",
            input$gene, format(nrow(df), big.mark = ","),
            mean(df$.expr), 100 * mean(df$.expr > 0))
  }))

  output$umap_title <- renderText(
    if (input$umap_color == "expr") paste0("UMAP — ", input$gene, " expression")
    else paste0("UMAP — colored by ", names(UMAP_AXES)[UMAP_AXES == input$umap_group]))
  output$violin_title <- renderText(
    paste0(input$gene, " expression by ", split_label(), " (Normal vs Tumor)"))

  # ---- UMAP ----------------------------------------------------------------
  output$umap <- renderPlotly({
    df <- cells()
    if (input$umap_color == "expr") {
      cmax <- max(df$.expr, 1e-6)
      expr_panel <- function(d, title, showbar)
        plot_ly(d, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                marker = list(size = 3, color = ~.expr, colorscale = "Viridis",
                              cmin = 0, cmax = cmax, showscale = showbar,
                              colorbar = list(title = "log-norm", len = 0.9)),
                text = ~paste0(Major_celltype, " / ", Subtype, "<br>", Sample,
                               " (", Condition, ")<br>expr: ", round(.expr, 2)),
                hoverinfo = "text") %>%
          layout(title = list(text = title, font = list(size = 12)),
                 xaxis = list(title = "UMAP_1", zeroline = FALSE),
                 yaxis = list(title = "UMAP_2", zeroline = FALSE))
      if (isTRUE(input$compare) && length(unique(df$Condition)) > 1) {
        conds <- cond_order(df$Condition)                     # (1) Normal left, Tumor right
        n <- length(conds)
        plts <- lapply(seq_along(conds), function(i)           # (1) single shared colorbar
          expr_panel(df[df$Condition == conds[i], ], conds[i], showbar = (i == n)))
        subplot(plts, nrows = 1, shareY = TRUE, titleX = TRUE, margin = 0.04)
      } else expr_panel(df, NULL, TRUE)
    } else {
      gv <- input$umap_group
      df$.g <- factor(df[[gv]]); lv <- levels(df$.g)
      p <- plot_ly(df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                   color = ~.g, colors = colors_for(gv, lv),
                   marker = list(size = 3, symbol = "square"),          # (1) square swatches
                   text = ~paste0(.g, "<br>expr: ", round(.expr, 2)), hoverinfo = "text")
      p %>% layout(xaxis = list(title = "UMAP_1", zeroline = FALSE),
                   yaxis = list(title = "UMAP_2", zeroline = FALSE),
                   legend = list(itemsizing = "constant", font = list(size = 13),  # (1) bigger legend
                                 title = list(text = names(UMAP_AXES)[UMAP_AXES == gv])))
    }
  })

  # ---- Violin: Normal vs Tumor dodged, median + rank-sum p on top ----------
  output$violin <- renderPlotly({
    df <- cells(); gv <- input$split_by
    df$.g <- as.character(df[[gv]])
    lv <- sort(unique(df$.g))
    conds <- cond_order(df$Condition)
    ymax <- max(df$.expr, 1e-6)

    if (length(conds) == 2) {
      off <- 0.2                                   # (2) offset -> gap between N and T
      cc <- cond_color(conds)
      dA <- df[df$Condition == conds[1], ]; dB <- df[df$Condition == conds[2], ]
      p <- plot_ly() %>%
        add_trace(type = "violin", x = match(dA$.g, lv) - off, y = dA$.expr, name = conds[1],
                  width = 0.34, opacity = 0.65, line = list(color = cc[[conds[1]]]),
                  fillcolor = cc[[conds[1]]], points = FALSE, box = list(visible = TRUE),
                  meanline = list(visible = TRUE)) %>%
        add_trace(type = "violin", x = match(dB$.g, lv) + off, y = dB$.expr, name = conds[2],
                  width = 0.34, opacity = 0.65, line = list(color = cc[[conds[2]]]),
                  fillcolor = cc[[conds[2]]], points = FALSE, box = list(visible = TRUE),
                  meanline = list(visible = TRUE))
      ann <- lapply(seq_along(lv), function(i) {                     # (2) median + p per group
        g <- lv[i]
        ea <- df$.expr[df$.g == g & df$Condition == conds[1]]
        eb <- df$.expr[df$.g == g & df$Condition == conds[2]]
        p_val <- tryCatch(
          if (length(ea) > 1 && length(eb) > 1) wilcox.test(eb, ea)$p.value else NA_real_,
          error = function(e) NA_real_)
        list(x = i, y = ymax * 1.08, showarrow = FALSE, xref = "x", yref = "y",
             font = list(size = 9), align = "center",
             text = sprintf("med %s/%s %.2f/%.2f<br>%s", substr(conds[1],1,1),
                            substr(conds[2],1,1), median(ea), median(eb), fmt_p(p_val)))
      })
      p %>% layout(
        xaxis = list(tickmode = "array", tickvals = seq_along(lv), ticktext = lv,
                     range = c(0.5, length(lv) + 0.5), tickangle = -30),
        yaxis = list(title = "log-norm expression", range = c(-0.03 * ymax, ymax * 1.20)),
        annotations = ann, legend = list(orientation = "h", x = 0, y = 1.02))
    } else {
      cond1 <- if (length(conds)) conds[1] else "cells"
      col <- if (length(conds)) cond_color(conds)[[conds[1]]] else "#8C8C8C"
      plot_ly(df, x = ~.g, y = ~.expr, type = "violin", name = cond1, points = FALSE,
              line = list(color = col), fillcolor = col, opacity = 0.65,
              box = list(visible = TRUE), meanline = list(visible = TRUE)) %>%
        layout(xaxis = list(title = "", tickangle = -30),
               yaxis = list(title = "log-norm expression"))
    }
  })

  # ---- Per-group summary: Tumor/Normal mean, median, rank-sum p ------------
  summary_df <- reactive({
    df <- cells(); gv <- input$split_by
    df$.g <- as.character(df[[gv]])
    lv <- sort(unique(df$.g))
    rows <- lapply(lv, function(g) {
      en <- df$.expr[df$.g == g & df$Condition == "Normal"]
      et <- df$.expr[df$.g == g & df$Condition == "Tumor"]
      p_val <- tryCatch(
        if (length(en) > 1 && length(et) > 1) wilcox.test(et, en)$p.value else NA_real_,
        error = function(e) NA_real_)
      data.frame(group = g,
                 n_Normal = length(en), mean_Normal = round(mean0(en), 3), median_Normal = round(med0(en), 3),
                 n_Tumor = length(et),  mean_Tumor  = round(mean0(et), 3), median_Tumor  = round(med0(et), 3),
                 p_value = signif(p_val, 3), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    names(out)[1] <- split_label()
    out[order(-(out$mean_Tumor + out$mean_Normal)), ]
  })
  mean0 <- function(x) if (length(x)) mean(x) else NA_real_
  med0  <- function(x) if (length(x)) median(x) else NA_real_

  output$summary_tbl <- renderDT({
    d <- summary_df()
    datatable(d, rownames = FALSE, options = list(pageLength = 15, dom = "tip")) %>%
      formatStyle("mean_Tumor",  background = styleColorBar(c(0, max(d$mean_Tumor, na.rm = TRUE)),  "#f2c5c5")) %>%
      formatStyle("mean_Normal", background = styleColorBar(c(0, max(d$mean_Normal, na.rm = TRUE)), "#c5d4f2"))
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0("CATLAS_", input$gene, "_by_", input$split_by, ".csv"),
    content  = function(file) write.csv(summary_df(), file, row.names = FALSE))

  # ---- Composition ----------------------------------------------------------
  output$composition <- renderPlotly({
    df <- cells()
    comp <- df %>% count(Condition, Major_celltype) %>%
      group_by(Condition) %>% mutate(freq = n / sum(n)) %>% ungroup()
    lv <- sort(unique(comp$Major_celltype))
    plot_ly(comp, x = ~Condition, y = ~freq, color = ~Major_celltype,
            colors = colors_for("Major_celltype", lv), type = "bar",
            text = ~paste0(Major_celltype, ": ", round(100 * freq, 1), "%"),
            hoverinfo = "text") %>%
      layout(barmode = "stack", yaxis = list(title = "cell fraction", tickformat = ".0%"),
             xaxis = list(title = ""))
  })

  # ---- Metadata -------------------------------------------------------------
  output$metadata_tbl <- renderDT(datatable(
    cells() %>% select(cell_id, Patient, Sample, Condition, Major_celltype, Subtype,
                       nCount_RNA, nFeature_RNA, percent.mt),
    options = list(pageLength = 20)))

  # ---- Natural-language search ---------------------------------------------
  llm <- reactiveValues(provider = "anthropic",
                        model = LLM_MODELS$anthropic[1], key = "")

  observeEvent(input$nl_settings, {
    showModal(modalDialog(title = "LLM API settings", easyClose = TRUE,
      selectInput("llm_provider", "Provider", choices = LLM_PROVIDERS, selected = llm$provider),
      selectizeInput("llm_model", "Model", choices = LLM_MODELS[[llm$provider]],
                     selected = llm$model, options = list(create = TRUE)),
      passwordInput("llm_key", "API key", value = llm$key),
      helpText("Keys are used only for your requests in this session and are never stored or logged."),
      footer = tagList(modalButton("Cancel"), actionButton("llm_save", "Save", class = "btn-primary"))))
  })
  observeEvent(input$llm_provider,
    updateSelectizeInput(session, "llm_model", choices = LLM_MODELS[[input$llm_provider]],
                         selected = LLM_MODELS[[input$llm_provider]][1]))
  observeEvent(input$llm_save, {
    llm$provider <- input$llm_provider; llm$model <- input$llm_model; llm$key <- input$llm_key
    removeModal()
    output$nl_status <- renderUI(div(class = "gene-note",
      sprintf("Provider set: %s / %s", names(LLM_PROVIDERS)[LLM_PROVIDERS == llm$provider], llm$model)))
  })

  observeEvent(input$nl_go, {
    q <- trimws(input$nl_query %||% "")
    if (!nzchar(q)) return(NULL)
    prompt <- paste0(
      "You translate a request about a colorectal-cancer single-cell atlas into UI settings. ",
      "Return ONLY a compact JSON object (no prose) with keys: ",
      "gene (a gene symbol string or null), ",
      "split_by (one of Major_celltype, Subtype, Sample, Patient), ",
      "umap_color ('expr' or 'group'), compare (true/false), ",
      "conditions (array subset of [\"Normal\",\"Tumor\"]), ",
      "celltypes (array subset of the provided list, or []), ",
      "minexpr (number 0-3). ",
      "Available cell types: ", paste(celltypes_all, collapse = ", "), ". ",
      "User request: \"", q, "\". JSON:")
    output$nl_status <- renderUI(div(class = "gene-note", "Querying LLM..."))
    res <- tryCatch(llm_call(llm$provider, llm$model, llm$key, prompt),
                    error = function(e) structure(conditionMessage(e), class = "llm_err"))
    if (inherits(res, "llm_err")) {
      output$nl_status <- renderUI(div(class = "gene-note", style = "color:#b00;",
                                       paste("Error:", as.character(res)))); return(NULL)
    }
    js <- regmatches(res, regexpr("\\{.*\\}", res, perl = TRUE))
    parsed <- tryCatch(jsonlite::fromJSON(js), error = function(e) NULL)
    if (is.null(parsed)) {
      output$nl_status <- renderUI(div(class = "gene-note", style = "color:#b00;",
                                       "Could not parse LLM response.")); return(NULL)
    }
    applied <- character(0)
    if (!is.null(parsed$gene) && length(parsed$gene) == 1 && !is.na(parsed$gene)) {
      hit <- available_genes[match(tolower(parsed$gene), gene_lower)]
      if (!is.na(hit)) { updateSelectizeInput(session, "gene", choices = available_genes,
                                              selected = hit, server = TRUE)
                         applied <- c(applied, paste0("gene=", hit)) }
    }
    if (!is.null(parsed$split_by) && parsed$split_by %in% SPLIT_AXES) {
      updateRadioButtons(session, "split_by", selected = parsed$split_by)
      applied <- c(applied, paste0("split_by=", parsed$split_by)) }
    if (!is.null(parsed$umap_color) && parsed$umap_color %in% c("expr", "group")) {
      updateRadioButtons(session, "umap_color", selected = parsed$umap_color)
      applied <- c(applied, paste0("umap=", parsed$umap_color)) }
    if (!is.null(parsed$compare) && is.logical(parsed$compare)) {
      updateCheckboxInput(session, "compare", value = parsed$compare) }
    if (!is.null(parsed$conditions) && length(parsed$conditions)) {
      cs <- intersect(as.character(parsed$conditions), conditions_all)
      if (length(cs)) { updateCheckboxGroupInput(session, "f_condition", selected = cs)
                        applied <- c(applied, paste0("conditions=", paste(cs, collapse = "/"))) } }
    if (!is.null(parsed$celltypes) && length(parsed$celltypes)) {
      ct <- intersect(as.character(parsed$celltypes), celltypes_all)
      if (length(ct)) { updateSelectizeInput(session, "f_celltype", selected = ct)
                        applied <- c(applied, paste0("celltypes=", paste(ct, collapse = "/"))) } }
    if (!is.null(parsed$minexpr) && is.numeric(parsed$minexpr)) {
      updateSliderInput(session, "minexpr", value = max(0, min(3, parsed$minexpr))) }
    output$nl_status <- renderUI(div(class = "gene-note",
      if (length(applied)) paste("Applied:", paste(applied, collapse = ", "))
      else "No actionable settings found."))
  })
}

shinyApp(ui, server)
