# Natural-language search helpers. The prompt builder and response parser are
# pure functions (no network, no Seurat) so they can be checked with synthetic
# vocabularies. Only llm_call() touches the network.

LLM_PROVIDERS <- c("Anthropic (Claude)" = "anthropic",
                   "OpenAI (ChatGPT)"   = "openai",
                   "Google (Gemini)"    = "gemini")
LLM_MODELS <- list(
  anthropic = c("claude-opus-4-8", "claude-sonnet-4-6", "claude-haiku-4-5-20251001"),
  openai    = c("gpt-4o", "gpt-4o-mini", "o4-mini"),
  gemini    = c("gemini-2.0-flash", "gemini-1.5-pro", "gemini-1.5-flash"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Instruction that turns a free-text request into a strict settings JSON.
build_nl_prompt <- function(query, celltypes) {
  paste0(
    "You translate a request about a colorectal-cancer single-cell atlas into UI ",
    "settings. Return ONLY a compact JSON object (no prose) with keys: ",
    "gene (a gene symbol string or null), ",
    "split_by (one of Major_celltype, Subtype, Sample, Patient, Condition), ",
    "umap_color ('expr' or 'group'), compare (true/false), ",
    "conditions (array subset of [\"Normal\",\"Tumor\"]), ",
    "celltypes (array subset of the provided list, or []), ",
    "minexpr (number 0-3). ",
    "Available cell types: ", paste(celltypes, collapse = ", "), ". ",
    "User request: \"", query, "\". JSON:")
}

# Validate the model's answer against the real vocabulary; drop anything unknown
# so a bad LLM response can never inject invalid state into the app.
parse_nl_response <- function(text, genes, celltypes, conditions, split_axes) {
  gene_lower <- tolower(genes)
  fail <- function(msg) list(ok = FALSE, msg = msg, settings = list(), applied = character())
  if (is.null(text) || !nzchar(text)) return(fail("Empty LLM response."))
  js <- regmatches(text, regexpr("\\{.*\\}", text, perl = TRUE))
  if (!length(js)) return(fail("No JSON found in LLM response."))
  parsed <- tryCatch(jsonlite::fromJSON(js), error = function(e) NULL)
  if (is.null(parsed)) return(fail("Could not parse LLM JSON."))

  settings <- list(); applied <- character()
  if (!is.null(parsed$gene) && length(parsed$gene) == 1L && !is.na(parsed$gene)) {
    hit <- genes[match(tolower(as.character(parsed$gene)), gene_lower)]
    if (!is.na(hit)) { settings$gene <- hit; applied <- c(applied, paste0("gene=", hit)) }
  }
  if (!is.null(parsed$split_by) && parsed$split_by %in% split_axes) {
    settings$split_by <- parsed$split_by
    applied <- c(applied, paste0("split_by=", parsed$split_by))
  }
  if (!is.null(parsed$umap_color) && parsed$umap_color %in% c("expr", "group")) {
    settings$umap_color <- parsed$umap_color
    applied <- c(applied, paste0("umap=", parsed$umap_color))
  }
  if (!is.null(parsed$compare) && is.logical(parsed$compare)) {
    settings$compare <- isTRUE(parsed$compare)
    applied <- c(applied, paste0("compare=", settings$compare))
  }
  if (!is.null(parsed$conditions)) {
    cs <- intersect(as.character(parsed$conditions), conditions)
    if (length(cs)) { settings$conditions <- cs
      applied <- c(applied, paste0("conditions=", paste(cs, collapse = "/"))) }
  }
  if (!is.null(parsed$celltypes)) {
    ct <- intersect(as.character(parsed$celltypes), celltypes)
    if (length(ct)) { settings$celltypes <- ct
      applied <- c(applied, paste0("celltypes=", paste(ct, collapse = "/"))) }
  }
  if (!is.null(parsed$minexpr) && is.numeric(parsed$minexpr)) {
    settings$minexpr <- max(0, min(3, parsed$minexpr))
    applied <- c(applied, paste0("minexpr=", settings$minexpr))
  }
  list(ok = length(applied) > 0L,
       msg = if (length(applied)) paste("Applied:", paste(applied, collapse = ", "))
             else "No actionable settings found.",
       settings = settings, applied = applied)
}

# Network call. Returns the model's raw text answer, or throws a readable error.
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
