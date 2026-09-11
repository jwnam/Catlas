# CRC scRNA-seq Atlas Viewer (v5) — interactive lncRNA explorer
# Uses the original RDS read-only. Builds on v4 (condition panel labels/colorbars,
# legend marker size, split-violin spacing, visible expression thresholds) and adds:
#   - "Split by" moved directly above the violin (it drives violin + summary);
#   - median (per condition) + Wilcoxon rank-sum p-value annotated on the violin;
#   - per-group summary with Tumor/Normal n, mean, median and the rank-sum p-value;
#   - natural-language search via a user-supplied LLM API key (Anthropic / OpenAI
#     / Google), configured in-app, keys kept in the session only.
# New optional dependency (natural-language search only): httr, jsonlite.

# --------------------------------------------------------------------------- #
#  Startup checks (kept identical in spirit to v2/app.R)                       #
# --------------------------------------------------------------------------- #
app_dir_value <- Sys.getenv("CATLAS_APP_DIR", unset = "")
if (!nzchar(app_dir_value)) app_dir_value <- getwd()
app_dir <- normalizePath(app_dir_value, winslash = "/", mustWork = TRUE)

project_dir_value <- Sys.getenv("CATLAS_PROJECT_DIR", unset = file.path(app_dir, "..", ".."))
project_dir <- normalizePath(project_dir_value, winslash = "/", mustWork = TRUE)

data_path_value <- Sys.getenv(
  "CATLAS_DATA_PATH",
  unset = file.path(project_dir, "crc_shiny_app_seurat.rds")
)
if (!file.exists(data_path_value)) {
  stop("CRC Atlas input RDS was not found: ", data_path_value)
}
data_path <- normalizePath(data_path_value, winslash = "/", mustWork = TRUE)

required_packages <- c("shiny", "Seurat", "Matrix", "ggplot2",
                       "dplyr", "DT", "plotly")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "),
       ". Install with install.packages() in the crc_shiny environment.")
}

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(DT)
  library(plotly)
})

# --------------------------------------------------------------------------- #
#  Load data (same access pattern as v2)                                      #
# --------------------------------------------------------------------------- #
message("CRC Atlas v4 loading: ", data_path)
obj <- readRDS(data_path)
DefaultAssay(obj) <- "RNA"

expr_mat <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
)

umap_df <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- rownames(umap_df)

meta_df <- obj@meta.data
meta_df$cell_id        <- rownames(meta_df)
meta_df$Condition      <- as.character(meta_df$Condition)      # Tumor / Normal
meta_df$Sample         <- as.character(meta_df$Library)
meta_df$Patient        <- as.character(meta_df$Patient)
meta_df$Major_celltype <- as.character(meta_df$Cell_type_v2)
meta_df$Subtype        <- as.character(meta_df$Cell_subtype_v2)
if ("Class" %in% colnames(meta_df)) meta_df$Class <- as.character(meta_df$Class)

plot_df_base <- dplyr::left_join(umap_df, meta_df, by = "cell_id")
rm(obj); invisible(gc())

available_genes <- rownames(expr_mat)

# The grouping axes the UI can "split by" (label -> column in plot_df_base)
AXES <- c("Cell type" = "Major_celltype",
          "Sub-type"  = "Subtype",
          "Sample"    = "Sample",
          "Condition" = "Condition",
          "Patient"   = "Patient")

conditions_all <- sort(unique(plot_df_base$Condition))
celltypes_all  <- sort(unique(plot_df_base$Major_celltype))
samples_all    <- sort(unique(plot_df_base$Sample))

# Shared discrete palette so UMAP and violin agree on group colors
PAL <- c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860",
         "#DA8BC3","#8C8C8C","#CCB974","#64B5CD","#E377C2","#7F7F7F",
         "#BCBD22","#17BECF","#AEC7E8","#FFBB78","#98DF8A","#FF9896")
pal_for <- function(levels) setNames(rep(PAL, length.out = length(levels)), levels)

# Condition is Tumor / Normal. Give those stable colors; the N-category fallback
# below still applies should Condition ever carry more levels.
COND_PAL <- c(Tumor = "#C44E52", Normal = "#4C72B0")
cond_color <- function(levels) {
  out <- COND_PAL[levels]
  miss <- is.na(out)
  if (any(miss)) out[miss] <- PAL[seq_len(sum(miss))]
  setNames(unname(out), levels)
}

source(file.path(app_dir, "R", "umap.R"), local = TRUE)
source(file.path(app_dir, "R", "violin.R"), local = TRUE)
source(file.path(app_dir, "R", "llm.R"), local = TRUE)

# --------------------------------------------------------------------------- #
#  UI                                                                          #
# --------------------------------------------------------------------------- #
ui <- fluidPage(
  tags$head(tags$style(HTML(
    ".ctrl-bar{background:#f5f7fa;border:1px solid #e3e8ee;border-radius:8px;padding:8px 12px;margin-bottom:10px;}
     .gene-note{font-size:12px;color:#666;margin-top:-6px;margin-bottom:8px;}
     .nl-box{background:#eef4fb;border:1px solid #cfe0f3;border-radius:8px;padding:8px;margin-bottom:8px;}
     .split-bar{background:#fbf6ef;border:1px solid #eadfce;border-radius:8px;padding:6px 10px;margin:6px 0;}"
  ))),
  titlePanel("CATLAS · CRC lncRNA Atlas"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "nl-box",
        tags$b("\U0001F50E Natural-language search"),
        textInput("nl_query", NULL, placeholder = "e.g. SNHG16 in T cells, tumor vs normal"),
        fluidRow(
          column(7, actionButton("nl_go", "Search", class = "btn-sm btn-primary", width = "100%")),
          column(5, actionLink("nl_settings", "⚙ API settings"))),
        uiOutput("nl_status")),
      selectizeInput("gene", "Gene (lncRNA / any feature)",
                     choices = NULL, selected = NULL, multiple = FALSE,
                     options = list(placeholder = "e.g. SNHG16, NORAD, LINC...",
                                    maxOptions = 100)),
      uiOutput("gene_note"),
      tags$hr(),
      tags$b("Filters"),
      checkboxGroupInput("f_condition", "Condition",
                         choices = conditions_all, selected = conditions_all, inline = TRUE),
      selectizeInput("f_celltype", "Cell type", choices = celltypes_all,
                     selected = celltypes_all, multiple = TRUE,
                     options = list(plugins = list("remove_button"))),
      selectizeInput("f_subtype", "Sub-type (optional)", choices = NULL,
                     selected = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"),
                                    placeholder = "all sub-types")),
      selectizeInput("f_sample", "Sample (optional)", choices = samples_all,
                     selected = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"),
                                    placeholder = "all samples")),
      sliderInput("minexpr", "Positive-expression threshold", min = 0, max = 3,
                  value = 0, step = 0.1),
      helpText("Positive means expression > threshold. Cells at or below the threshold are gray on the UMAP; the dashed violin line marks the threshold. All selected cells remain in the distributions and summary denominator."),
      tags$hr(),
      downloadButton("dl_csv", "Summary (CSV)", class = "btn-sm btn-outline-primary")
    ),

    mainPanel(
      width = 9,
      div(class = "ctrl-bar",
        fluidRow(
          column(7, radioButtons("umap_color", "UMAP color",
                                 choices = c("Expression" = "expr", "Group" = "group"),
                                 selected = "expr", inline = TRUE)),
          column(5, checkboxInput("compare", "Compare conditions side-by-side", TRUE))
        )
      ),
      tabsetPanel(
        tabPanel("Gene expression",
          br(),
          h5(textOutput("umap_title")),
          plotlyOutput("umap", height = "440px"),
          div(class = "split-bar",
            fluidRow(
              column(8, radioButtons("split_by", "Split violin by",
                                     choices = AXES, selected = "Major_celltype", inline = TRUE)),
              column(4, div(style = "padding-top:6px;color:#666;font-size:12px;",
                            "Blue = Normal · Red = Tumor · p = Wilcoxon rank-sum")))),
          h5(textOutput("violin_title")),
          plotlyOutput("violin", height = "420px"),
          br(),
          h5("Per-group summary (Tumor vs Normal)"),
          DTOutput("summary_tbl")
        ),
        tabPanel("Composition", br(), plotlyOutput("composition", height = "600px")),
        tabPanel("Metadata", br(), DTOutput("metadata_tbl"))
      )
    )
  )
)

# --------------------------------------------------------------------------- #
#  Server                                                                      #
# --------------------------------------------------------------------------- #
server <- function(input, output, session) {

  updateSelectizeInput(session, "gene", choices = available_genes,
                       selected = if ("SNHG16" %in% available_genes) "SNHG16"
                                  else available_genes[1], server = TRUE)

  # Include the selected gene's full expression range instead of fixing the
  # slider at 3. Preserve the current threshold whenever the new range allows it.
  observeEvent(input$gene, {
    req(input$gene %in% available_genes)
    values <- as.numeric(expr_mat[input$gene, ])
    finite_values <- values[is.finite(values)]
    upper <- max(3, if (length(finite_values)) ceiling(max(finite_values) * 10) / 10 else 3)
    current <- isolate(input$minexpr)
    if (is.null(current) || !is.finite(current)) current <- 0
    updateSliderInput(session, "minexpr", max = upper, value = min(current, upper))
  }, ignoreNULL = TRUE)

  # sub-type choices follow the selected cell types
  observeEvent(input$f_celltype, {
    subs <- plot_df_base %>%
      filter(Major_celltype %in% input$f_celltype) %>%
      pull(Subtype) %>% unique() %>% sort()
    updateSelectizeInput(session, "f_subtype", choices = subs,
                         selected = intersect(input$f_subtype, subs))
  }, ignoreNULL = FALSE)

  # cells passing filters, with the current gene's expression attached
  cells <- reactive({
    req(input$gene)
    validate(need(input$gene %in% available_genes, "Gene not found."))
    df <- plot_df_base
    df <- df[df$Condition %in% input$f_condition, , drop = FALSE]
    df <- df[df$Major_celltype %in% input$f_celltype, , drop = FALSE]
    if (length(input$f_subtype))   df <- df[df$Subtype %in% input$f_subtype, , drop = FALSE]
    if (length(input$f_sample))    df <- df[df$Sample %in% input$f_sample, , drop = FALSE]
    validate(need(nrow(df) > 0, "No cells match the current filters."))
    df$.expr <- as.numeric(expr_mat[input$gene, df$cell_id])
    req(!is.null(input$minexpr), is.finite(input$minexpr))
    df$.positive <- df$.expr > input$minexpr
    df
  })

  axis_label <- reactive(names(AXES)[AXES == input$split_by])

  output$gene_note <- renderUI({
    df <- cells()
    div(class = "gene-note",
        sprintf("%s — %s / %s cells positive (%.1f%%; expression > %.1f; mean %.2f, log-norm)",
                input$gene, format(sum(df$.positive), big.mark = ","),
                format(nrow(df), big.mark = ","), 100 * mean(df$.positive),
                input$minexpr, mean(df$.expr)))
  })

  output$umap_title <- renderText({
    if (input$umap_color == "expr")
      paste0("UMAP — ", input$gene, " expression")
    else paste0("UMAP — colored by ", axis_label())
  })
  output$violin_title <- renderText(
    paste0(input$gene, " expression by ", axis_label()))

  # ---- UMAP / violin --------------------------------------------------------
  output$umap <- renderPlotly({
    build_umap_plot(cells(), input$split_by, input$umap_color,
                    isTRUE(input$compare), input$minexpr)
  })

  output$violin <- renderPlotly({
    build_violin_plot(cells(), input$split_by,
                      isTRUE(input$compare), input$minexpr)
  })

  # ---- Per-group summary (adds Tumor/Normal mean, median, rank-sum p) -------
  .mean0 <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  .med0  <- function(x) { x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
  .rsp   <- function(a, b) { a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) > 1L && length(b) > 1L)
      tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_)
    else NA_real_ }

  summary_df <- reactive({
    df <- cells(); thr <- input$minexpr; gv <- input$split_by
    df %>%
      group_by(.grp = .data[[gv]]) %>%
      summarise(n_cells = n(),
                n_positive = sum(.positive),
                pct_pos = round(100 * mean(.positive), 1),
                threshold = thr,
                mean    = round(mean(.expr), 3),
                median  = round(median(.expr), 3),
                n_Normal      = sum(Condition == "Normal"),
                mean_Normal   = round(.mean0(.expr[Condition == "Normal"]), 3),
                median_Normal = round(.med0(.expr[Condition == "Normal"]), 3),
                n_Tumor       = sum(Condition == "Tumor"),
                mean_Tumor    = round(.mean0(.expr[Condition == "Tumor"]), 3),
                median_Tumor  = round(.med0(.expr[Condition == "Tumor"]), 3),
                p_value = signif(.rsp(.expr[Condition == "Tumor"],
                                      .expr[Condition == "Normal"]), 3),
                .groups = "drop") %>%
      rename(!!axis_label() := .grp) %>%
      arrange(desc(mean))
  })

  output$summary_tbl <- renderDT({
    d <- summary_df()
    datatable(d, rownames = FALSE, options = list(pageLength = 15, dom = "tip",
                                                  scrollX = TRUE)) %>%
      formatStyle("mean_Tumor",  background =
                    styleColorBar(c(0, max(d$mean_Tumor, na.rm = TRUE)),  "#f2c5c5")) %>%
      formatStyle("mean_Normal", background =
                    styleColorBar(c(0, max(d$mean_Normal, na.rm = TRUE)), "#c5d4f2"))
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0("CATLAS_", input$gene, "_by_", input$split_by, ".csv"),
    content  = function(file) write.csv(summary_df(), file, row.names = FALSE)
  )

  # ---- Composition ----------------------------------------------------------
  output$composition <- renderPlotly({
    df <- cells()
    comp <- df %>% count(Condition, Major_celltype) %>%
      group_by(Condition) %>% mutate(freq = n / sum(n)) %>% ungroup()
    lv <- sort(unique(comp$Major_celltype))
    plot_ly(comp, x = ~Condition, y = ~freq, color = ~Major_celltype,
            colors = pal_for(lv), type = "bar",
            text = ~paste0(Major_celltype, ": ", round(100 * freq, 1), "%"),
            hoverinfo = "text") %>%
      layout(barmode = "stack", yaxis = list(title = "cell fraction", tickformat = ".0%"),
             xaxis = list(title = ""))
  })

  # ---- Metadata -------------------------------------------------------------
  output$metadata_tbl <- renderDT({
    datatable(
      cells() %>% select(cell_id, Patient, Sample, Condition,
                         Major_celltype, Subtype, nCount_RNA, nFeature_RNA, percent.mt),
      options = list(pageLength = 20))
  })

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
      footer = tagList(modalButton("Cancel"),
                       actionButton("llm_save", "Save", class = "btn-primary"))))
  })
  observeEvent(input$llm_provider,
    updateSelectizeInput(session, "llm_model", choices = LLM_MODELS[[input$llm_provider]],
                         selected = LLM_MODELS[[input$llm_provider]][1]))
  observeEvent(input$llm_save, {
    llm$provider <- input$llm_provider; llm$model <- input$llm_model; llm$key <- input$llm_key
    removeModal()
    output$nl_status <- renderUI(div(class = "gene-note",
      sprintf("Provider set: %s / %s",
              names(LLM_PROVIDERS)[LLM_PROVIDERS == llm$provider], llm$model)))
  })

  observeEvent(input$nl_go, {
    q <- trimws(input$nl_query %||% "")
    if (!nzchar(q)) return(NULL)
    output$nl_status <- renderUI(div(class = "gene-note", "Querying LLM..."))
    prompt <- build_nl_prompt(q, celltypes_all)
    res <- tryCatch(llm_call(llm$provider, llm$model, llm$key, prompt),
                    error = function(e) structure(conditionMessage(e), class = "llm_err"))
    if (inherits(res, "llm_err")) {
      output$nl_status <- renderUI(div(class = "gene-note", style = "color:#b00;",
                                       paste("Error:", as.character(res))))
      return(NULL)
    }
    pr <- parse_nl_response(res, available_genes, celltypes_all,
                            conditions_all, unname(AXES))
    s <- pr$settings
    if (!is.null(s$gene))
      updateSelectizeInput(session, "gene", choices = available_genes,
                           selected = s$gene, server = TRUE)
    if (!is.null(s$split_by))   updateRadioButtons(session, "split_by", selected = s$split_by)
    if (!is.null(s$umap_color)) updateRadioButtons(session, "umap_color", selected = s$umap_color)
    if (!is.null(s$compare))    updateCheckboxInput(session, "compare", value = s$compare)
    if (!is.null(s$conditions)) updateCheckboxGroupInput(session, "f_condition", selected = s$conditions)
    if (!is.null(s$celltypes))  updateSelectizeInput(session, "f_celltype", selected = s$celltypes)
    if (!is.null(s$minexpr))    updateSliderInput(session, "minexpr", value = s$minexpr)
    output$nl_status <- renderUI(div(class = "gene-note", pr$msg))
  })
}

shinyApp(ui, server)
