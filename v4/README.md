# CRC Atlas Shiny v4 — interactive lncRNA explorer

v4 builds on v3 with UI/plot refinements and a natural-language search.
Same read-only RDS, same startup contract (run with `CATLAS_APP_DIR=v4`),
same data access as v2/v3.

## Changes vs v3

1. **UMAP (expression, compare)** — a **single shared colorbar** (was two,
   overlapping); panels always drawn as **Normal (left) / Tumor (right)**.
   In **Group** mode the legend uses **larger square swatches** (constant
   sizing) so colours are easy to tell apart, and **sub-type colours inherit
   the parent cell type's hue family** (light→dark shades).
2. **"Split by" moved directly above the violin** (it drives the violin +
   summary). The violin always shows **Normal (blue) vs Tumor (red)** dodged
   with a small gap, and annotates **median (N/T) + Wilcoxon rank-sum p-value**
   on top of each group.
3. **Per-group summary** now reports, per group, **Tumor/Normal n, mean, median**
   and the **rank-sum p-value** (CSV download included).
4. **Min-expression filter** actually filters the cells used in every view
   (slider > 0 → cells with `expr < value` are dropped).
5. **Natural-language search** — type a request (e.g. *"SNHG16 in T cells,
   tumor vs normal"*); an LLM maps it to gene + filters + view settings.

## Natural-language search / LLM configuration

Bring your own API key. Click **⚙ API settings** to pick a provider + model and
paste a key. Supported providers:

| Provider | Endpoint | Example models |
| --- | --- | --- |
| Anthropic (Claude) | `api.anthropic.com/v1/messages` | claude-opus-4-8, claude-sonnet-4-6, claude-haiku-4-5 |
| OpenAI (ChatGPT) | `api.openai.com/v1/chat/completions` | gpt-4o, gpt-4o-mini, o4-mini |
| Google (Gemini) | `generativelanguage.googleapis.com` | gemini-2.0-flash, gemini-1.5-pro |

- Keys are held **in the session only** — never written to disk or logged.
- The model must return a small JSON object; the app **validates** every field
  against the real vocabulary (gene symbols, cell types, conditions) and ignores
  anything invalid, so a bad LLM answer can't break the app.
- Model names are editable (type a custom one in the Model box).

## Dependencies

- **Required:** v2's packages + `plotly` (same as v3).
- **Optional (NL search only):** `httr`, `jsonlite`. The app starts fine without
  them; NL search shows a message asking to install them if used.

```r
install.packages(c("plotly", "httr", "jsonlite"))
```

## Run

```bash
CATLAS_PROJECT_DIR="$PWD" CATLAS_APP_DIR="$PWD/v4" \
CATLAS_R_BIN="$(command -v R)" bash v2/start_catlas_v2.sh
```

## Validation status

Authored on a machine without the 2.53 GB Seurat RDS. Data-loading is reused
verbatim from v2/v3. The new logic was validated headlessly under R against a
30k-cell synthetic dataset that shares the plotting/stat code:

- syntax parse of `app.R`
- dodged Normal/Tumor violin with per-group median + rank-sum annotations
- single shared colorbar across the 2-panel expression UMAP
- sub-type shade palette (same hue family as parent cell type)
- per-group Tumor/Normal summary (mean/median/p) 
- natural-language JSON parsing + field validation (invalid values dropped)

Still needs one launch on the crc_shiny server to confirm end-to-end with the
real object and a live LLM key.
