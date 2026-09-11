# CATLAS v5

CATLAS v5 is an interactive R Shiny application for exploring colorectal cancer single-cell RNA sequencing data, with a focus on lncRNA expression. Users can search for any feature available in the supplied Seurat object.

v5 builds on v4 (shared expression colorbar, condition panel titles, split-violin spacing, legend marker size, and a visible positive-expression threshold) and adds a statistical annotation layer, a Tumor/Normal summary breakdown, and a natural-language search. It preserves the existing v2–v4 files and uses the original RDS as read-only input.

## Features

- Explore gene expression on interactive UMAP and violin plots.
- Compare Normal and Tumor conditions side by side, with **median and Wilcoxon rank-sum p-value** annotated per group on the violin.
- Filter cells by condition, cell type, subtype, and sample.
- Group results by cell type, subtype, sample, condition, or patient — the **Split by** control now sits directly above the violin it drives.
- View per-group expression summaries (**including Tumor/Normal mean, median, and rank-sum p-value**), cell composition, and metadata.
- **Natural-language search**: describe what you want (e.g. *"SNHG16 in T cells, tumor vs normal"*) and an LLM maps it to gene + filters + view settings, using your own API key.
- Download the current per-group summary as a CSV file.

## Changes from v4

| Area | v4 | v5 |
| --- | --- | --- |
| Split-by placement | In the top control bar above the UMAP. | Moved directly above the violin plot, since it drives the violin and summary. |
| Violin statistics | Split halves with a small gap and a threshold line. | Adds a per-group **median (Normal/Tumor) + Wilcoxon rank-sum p-value** label on top of each group. |
| Per-group summary | `n_cells`, `n_positive`, `pct_pos`, `threshold`, `mean`, `median`. | Adds `n_Normal`, `mean_Normal`, `median_Normal`, `n_Tumor`, `mean_Tumor`, `median_Tumor`, and `p_value`. |
| Search | Gene dropdown only. | Adds an optional natural-language search backed by a user-supplied LLM API key. |

## Natural-language search

Click **⚙ API settings** in the sidebar to choose a provider and model and paste an API key, then type a request and press **Search**.

| Provider | Endpoint | Example models |
| --- | --- | --- |
| Anthropic (Claude) | `api.anthropic.com/v1/messages` | claude-opus-4-8, claude-sonnet-4-6, claude-haiku-4-5 |
| OpenAI (ChatGPT) | `api.openai.com/v1/chat/completions` | gpt-4o, gpt-4o-mini, o4-mini |
| Google (Gemini) | `generativelanguage.googleapis.com` | gemini-2.0-flash, gemini-1.5-pro |

- Keys are held **in the session only** — never written to disk or logged.
- The model must return a small JSON object; `parse_nl_response()` **validates every field** against the real vocabulary (gene symbols, cell types, conditions, grouping axes) and silently drops anything invalid, so a bad LLM answer cannot break the app.
- Model names are editable (type a custom one in the Model box).

## Changes from v3

| Issue | Behavior in v3 | Improvement in v4 |
| --- | --- | --- |
| Overlapping expression colorbars | Each condition panel created a separate colorbar at the same default position, and expression color limits were determined separately. | Condition panels share one expression color scale and a single colorbar when cells above the threshold are present. |
| Missing Normal panel title | Individual `layout(title=...)` settings were combined into a single figure title by `subplot()`. | Normal and Tumor titles are placed separately above their respective panels. |
| No gap between split violin halves | Both condition halves used the same category center. | The halves are shifted slightly left and right to create a small gap, while group labels remain centered. |
| Small cell-type legend markers | Both UMAP points and legend markers used size 3. | UMAP points retain size 3, while legend markers use size 6. Each group has one legend entry, and clicking it toggles that group across both condition panels. |
| Threshold changes were not visible in the plots | The threshold affected positivity notes and summaries but was not passed to the UMAP or violin functions. | Changing the threshold immediately updates gray UMAP cells and the dashed threshold line on the violin plot. |

## Positive-expression threshold

A cell is positive when **`expression > threshold`**. Cells with **`expression <= threshold`**, including cells exactly at the threshold, appear gray on the UMAP.

This rule applies to both UMAP color modes:

- **Expression:** Cells above the threshold use the shared continuous expression scale.
- **Group:** Cells above the threshold use their group color.

The violin plot shows the current threshold as a horizontal dashed line. Changing the threshold does not remove cells from the selected population, so the violin distributions, total cell counts, and the cells used to calculate means and medians remain unchanged.

The summary table and CSV contain the selected grouping column and these fields:

| Field | Description |
| --- | --- |
| `n_cells` | Total number of selected cells in the group. |
| `n_positive` | Number of cells with expression strictly greater than the threshold. |
| `pct_pos` | `100 * n_positive / n_cells`, rounded to one decimal place. |
| `threshold` | The threshold used for the summary. |
| `mean` | Mean expression across all selected cells in the group, rounded to three decimal places. |
| `median` | Median expression across all selected cells in the group, rounded to three decimal places. |
| `n_Normal`, `n_Tumor` | Selected cell counts in the group split by condition (they sum to `n_cells`). |
| `mean_Normal`, `median_Normal` | Mean/median expression of the group's Normal cells. |
| `mean_Tumor`, `median_Tumor` | Mean/median expression of the group's Tumor cells. |
| `p_value` | Wilcoxon rank-sum p-value comparing Tumor vs Normal expression in the group (`NA` when either side has fewer than two cells). |

The slider maximum follows the selected gene's maximum expression, rounded upward to the next 0.1 increment, with a minimum upper limit of 3. When a different gene is selected, the current threshold is retained if it remains within the new range; otherwise, it is reduced to the new maximum.

Clearing all Condition or Cell type selections produces a message that no cells match the filters. Leaving the optional Sub-type or Sample selections empty applies no filter for that field.

## Files and data

| File | Purpose |
| --- | --- |
| [app.R](app.R) | Shiny interface, data loading, filters, summaries, natural-language search wiring, and downloads. |
| [R/umap.R](R/umap.R) | UMAP panels, shared expression color scale, group legends, and threshold coloring. |
| [R/violin.R](R/violin.R) | Violin plots, condition spacing, threshold overlays, and median + rank-sum p-value annotations. |
| [R/llm.R](R/llm.R) | Natural-language prompt builder, response parser/validator, and provider API calls. |
| [start_catlas_v5.sh](start_catlas_v5.sh) | Launcher that configures paths, temporary directories, host, and port before starting R. |
| [tests/run_checks.R](tests/run_checks.R) | Synthetic regression checks for plot construction, summary statistics, the NL parser, and Shiny reactive behavior. |

The original server layout uses:

- Application: `/home/minwook/Shiny_CRC_atlas/Catlas/v5/app.R`
- Read-only input data: `/home/minwook/Shiny_CRC_atlas/crc_shiny_app_seurat.rds`
- R executable: `/home/minwook/miniconda3/envs/crc_shiny/bin/R`

The application reads expression values from the RNA assay's `data` layer, with a fallback to the `data` slot for compatible Seurat versions.

## Requirements

Use a Linux environment with Bash and an R environment containing `shiny`, `Seurat`, `Matrix`, `ggplot2`, `dplyr`, `DT`, and `plotly`. The plotting code also uses `htmltools`, and the regression checks use `jsonlite`. Natural-language search is optional and additionally needs `httr` and `jsonlite`; the app starts without them and only reports the missing packages if you run a search.

```r
install.packages(c("httr", "jsonlite"))  # only needed for natural-language search
```

The launcher requires read access to the RDS file and write access to its temporary-directory parent. By default, temporary files are created under `v5/tmp_for_catlas`, with directory permissions set to `0700`.

## Preview on the original server

The launcher defaults to `127.0.0.1:4511`. Confirm the current service configuration and use a separate available port for a v5 preview, distinct from any running service.

```bash
CATLAS_PORT=4511 bash /home/minwook/Shiny_CRC_atlas/Catlas/v5/start_catlas_v5.sh
```

The preview is accessible at `http://127.0.0.1:4511` on the server. This loopback address refers to the machine running the application. Keep the preview port distinct from the port used by an existing service.

The v5 implementation does not change the production systemd entry point or Apache configuration. Publishing code to the `v5` branch does not switch the running web service to v5.

For another server layout, set `CATLAS_PROJECT_DIR`, `CATLAS_APP_DIR`, `CATLAS_DATA_PATH`, and `CATLAS_R_BIN` to the appropriate locations. `CATLAS_TMP_ROOT`, `CATLAS_HOST`, and `CATLAS_PORT` can also be overridden through environment variables.

## Validation

Run the synthetic regression checks with the configured R environment:

```bash
/home/minwook/miniconda3/envs/crc_shiny/bin/Rscript /home/minwook/Shiny_CRC_atlas/Catlas/v5/tests/run_checks.R
```

The v5 validation passed **63 synthetic regression checks** (the 54 inherited from v4 plus 9 new). These checks do not load the atlas RDS or start a web server. They cover:

- Condition panel titles and a shared expression color scale.
- UMAP point sizes, legend marker sizes, and linked legend groups.
- Gray-cell counts at different thresholds, including equality at the threshold.
- Split violin spacing, centered labels, and threshold lines.
- **Per-group median + rank-sum p-value annotations above the split violins.**
- Reactive threshold and summary updates in the actual Shiny server function.
- **The Tumor/Normal summary breakdown (per-condition counts summing to `n_cells`, numeric `p_value`).**
- **The natural-language parser: case-insensitive gene resolution, dropping invalid cell types, applying valid settings, and failing cleanly with no JSON.**
- Single-positive-cell cases, ensuring expression colors remain JSON arrays and retain the shared colorbar.

Because this development environment had no access to the 2.53 GB Seurat RDS, v5 has **not yet been launched against the real object or a live LLM key**; the data-loading path is unchanged from v4, and all new logic is covered by the synthetic checks above. A server launch is still recommended to confirm end-to-end behavior.
