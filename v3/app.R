# CRC scRNA-seq Atlas Viewer (v3) — interactive lncRNA explorer
#
# Drop-in successor to v2/app.R. Reuses the same read-only RDS, the same startup
# environment contract (start_catlas_v2.sh with CATLAS_APP_DIR pointing at v3),
# and the same data-access pattern. The UI is rebuilt around a single
# gene-centric view: pick one gene, then explore its expression by
# cell type / sub-type / sample / condition with linked UMAP + violin plots.
#
# New dependency vs v2: `plotly` (for hover / zoom / linked interactivity).

# --------------------------------------------------------------------------- #
#  Startup checks (kept identical in spirit to v2/app.R)                       #
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
message("CRC Atlas v3 loading: ", data_path)
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

# --------------------------------------------------------------------------- #
#  UI                                                                          #
# --------------------------------------------------------------------------- #
ui <- fluidPage(
  tags$head(tags$style(HTML(
    ".ctrl-bar{background:#f5f7fa;border:1px solid #e3e8ee;border-radius:8px;padding:8px 12px;margin-bottom:10px;}
     .gene-note{font-size:12px;color:#666;margin-top:-6px;margin-bottom:8px;}"
  ))),
  titlePanel("CATLAS · CRC lncRNA Atlas (v3)"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
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
      tags$hr(),
      downloadButton("dl_csv", "Summary (CSV)", class = "btn-sm btn-outline-primary")
    ),

    mainPanel(
      width = 9,
      div(class = "ctrl-bar",
        fluidRow(
          column(5, radioButtons("split_by", "Split by",
                                 choices = AXES, selected = "Major_celltype", inline = TRUE)),
          column(4, radioButtons("umap_color", "UMAP color",
                                 choices = c("Expression" = "expr", "Group" = "group"),
                                 selected = "expr", inline = TRUE)),
          column(3, checkboxInput("compare", "Compare conditions side-by-side", TRUE))
        )
      ),
      tabsetPanel(
        tabPanel("Gene expression",
          br(),
          h5(textOutput("umap_title")),
          plotlyOutput("umap", height = "440px"),
          br(),
          h5(textOutput("violin_title")),
          plotlyOutput("violin", height = "380px"),
          br(),
          h5("Per-group summary"),
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
    if (length(input$f_condition)) df <- df[df$Condition %in% input$f_condition, , drop = FALSE]
    if (length(input$f_celltype))  df <- df[df$Major_celltype %in% input$f_celltype, , drop = FALSE]
    if (length(input$f_subtype))   df <- df[df$Subtype %in% input$f_subtype, , drop = FALSE]
    if (length(input$f_sample))    df <- df[df$Sample %in% input$f_sample, , drop = FALSE]
    validate(need(nrow(df) > 0, "No cells match the current filters."))
    df$.expr <- as.numeric(expr_mat[input$gene, df$cell_id])
    df
  })

  axis_label <- reactive(names(AXES)[AXES == input$split_by])

  output$gene_note <- renderUI({
    df <- cells()
    div(class = "gene-note",
        sprintf("%s — %.1f%% of %s filtered cells positive (mean %.2f, log-norm)",
                input$gene, 100 * mean(df$.expr > input$minexpr),
                format(nrow(df), big.mark = ","), mean(df$.expr)))
  })

  output$umap_title <- renderText({
    if (input$umap_color == "expr")
      paste0("UMAP — ", input$gene, " expression")
    else paste0("UMAP — colored by ", axis_label())
  })
  output$violin_title <- renderText(
    paste0(input$gene, " expression by ", axis_label()))

  # ---- UMAP (plotly / WebGL) ------------------------------------------------
  output$umap <- renderPlotly({
    df <- cells(); gv <- input$split_by
    df$.g <- factor(df[[gv]])

    one_panel <- function(d, title = NULL) {
      if (input$umap_color == "expr") {
        p <- plot_ly(d, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                marker = list(size = 3, color = ~.expr, colorscale = "Viridis",
                              showscale = TRUE, colorbar = list(title = "log-norm")),
                text = ~paste0(Major_celltype, " / ", Subtype,
                               "<br>", Sample, " (", Condition, ")",
                               "<br>expr: ", round(.expr, 2)),
                hoverinfo = "text")
      } else {
        p <- plot_ly(d, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                color = ~.g, colors = pal_for(levels(d$.g)),
                marker = list(size = 3),
                text = ~paste0(.g, "<br>expr: ", round(.expr, 2)), hoverinfo = "text")
      }
      p %>% layout(title = list(text = title, font = list(size = 12)),
                   xaxis = list(title = "UMAP_1", zeroline = FALSE),
                   yaxis = list(title = "UMAP_2", zeroline = FALSE))
    }

    if (isTRUE(input$compare) && length(unique(df$Condition)) > 1) {
      conds <- sort(unique(df$Condition))
      plts <- lapply(conds, function(cc) one_panel(df[df$Condition == cc, ], cc))
      subplot(plts, nrows = 1, shareY = TRUE, titleX = TRUE, margin = 0.03) %>%
        layout(showlegend = (input$umap_color == "group"))
    } else {
      one_panel(df) %>% layout(showlegend = (input$umap_color == "group"))
    }
  })

  # ---- Violin (plotly) ------------------------------------------------------
  output$violin <- renderPlotly({
    df <- cells(); gv <- input$split_by
    df$.g <- factor(df[[gv]])
    conds <- sort(unique(as.character(df$Condition)))

    # grouped violin across the split axis (used stand-alone and per facet)
    grouped_violin <- function(d) {
      lv <- levels(factor(d$.g))
      plot_ly(d, x = ~.g, y = ~.expr, type = "violin", color = ~.g,
              colors = pal_for(lv), points = FALSE,
              box = list(visible = TRUE), meanline = list(visible = TRUE),
              showlegend = FALSE) %>%
        layout(xaxis = list(title = ""), yaxis = list(title = "log-norm expression"))
    }

    if (isTRUE(input$compare) && gv != "Condition" && length(conds) == 2) {
      # exactly 2 conditions: split violin (back-to-back) per group.
      # (violinmode defaults to "overlay"; setting it in layout() errors on some
      #  plotly versions, so we rely on the default.)
      cc <- cond_color(conds)
      plot_ly() %>%
        add_trace(type = "violin", data = df[df$Condition == conds[1], ],
                  x = ~.g, y = ~.expr, name = conds[1], side = "negative", opacity = 0.6,
                  line = list(color = cc[[conds[1]]]), fillcolor = cc[[conds[1]]],
                  points = FALSE, meanline = list(visible = TRUE)) %>%
        add_trace(type = "violin", data = df[df$Condition == conds[2], ],
                  x = ~.g, y = ~.expr, name = conds[2], side = "positive", opacity = 0.6,
                  line = list(color = cc[[conds[2]]]), fillcolor = cc[[conds[2]]],
                  points = FALSE, meanline = list(visible = TRUE)) %>%
        layout(xaxis = list(title = ""),
               yaxis = list(title = "log-norm expression"))
    } else if (isTRUE(input$compare) && gv != "Condition" && length(conds) > 2) {
      # 3+ conditions: one grouped-violin panel per condition, side by side.
      plts <- lapply(conds, function(cc) {
        grouped_violin(df[df$Condition == cc, ]) %>%
          layout(annotations = list(text = cc, x = 0.5, y = 1.04, xref = "paper",
                                    yref = "paper", showarrow = FALSE,
                                    font = list(size = 12)))
      })
      subplot(plts, nrows = 1, shareY = TRUE, margin = 0.02)
    } else {
      grouped_violin(df)
    }
  })

  # ---- Per-group summary ----------------------------------------------------
  summary_df <- reactive({
    df <- cells(); thr <- input$minexpr; gv <- input$split_by
    df %>%
      group_by(.grp = .data[[gv]]) %>%
      summarise(n_cells = n(),
                pct_pos = round(100 * mean(.expr > thr), 1),
                mean    = round(mean(.expr), 3),
                median  = round(median(.expr), 3),
                .groups = "drop") %>%
      rename(!!axis_label() := .grp) %>%
      arrange(desc(mean))
  })

  output$summary_tbl <- renderDT({
    d <- summary_df()
    datatable(d, rownames = FALSE, options = list(pageLength = 15, dom = "tip")) %>%
      formatStyle("mean", background =
                    styleColorBar(c(0, max(d$mean, na.rm = TRUE)), "#c8e6c9"))
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
}

shinyApp(ui, server)
