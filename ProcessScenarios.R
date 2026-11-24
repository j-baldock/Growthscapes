# ProcessScenarios.R
#
# Combine the results of every experimental run (one .rds per entry in
# `exp_scenarios`) into a single .rds per key simulation output, tagged with the
# experimental-design metadata from `exp_design`. Downstream .qmd files read
# these combined files instead of the individual run results.
#
# Outputs (results/combined/<output>.rds), one row-block per scenario:
#   census, daily_patch, ibm_summary, spawn_log, surv_by_dsh, extinct,
#   fish_registry (parent fish only), stock_daily, habitat_df, harvest_log
# (harvest_log is only written if at least one run has harvest.)
#
# Two stages, so run time is linear in the number of scenarios and memory is
# bounded:
#   1. Per scenario: read the run, write one small intermediate .rds per output
#      to results/combined/tmp/<output>/<scenario>.rds, and free the run.
#   2. Per output: read that output's intermediates, bind them, save a single
#      combined .rds, and delete the intermediates. Peak memory is the size of
#      the largest single combined output.
# Intermediates are kept if stage 2 fails, and `resume <- TRUE` skips scenarios
# whose intermediates are already complete. Run after BatchRunScenarios.R.

# ── Setup ──────────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)

# Written by ScenarioDefs.qmd (re-render / qmd_to_r_script + source it to refresh)
exp_design    <- readRDS("source_files/exp_design.rds")
exp_scenarios <- readRDS("source_files/exp_scenarios.rds")

resume      <- FALSE   # TRUE: skip scenarios already fully processed in stage 1
results_dir <- "results"

if (!dir.exists(file.path(results_dir, "combined"))) dir.create(file.path(results_dir, "combined"))
out_dir     <- file.path(results_dir, "combined")
tmp_dir     <- file.path(out_dir, "tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

outputs <- c("census", "daily_patch", "ibm_summary", "spawn_log", "surv_by_dsh",
             "extinct", "fish_registry", "stock_daily", "habitat_df", "harvest_log")

# ── Metadata ───────────────────────────────────────────────────────────────────
# Design columns attached to every output. A_cold/A_warm/pcmax_* come from the
# design table (habitat_df also carries its own time-varying versions; the
# habitat_df copies are kept and the design values are dropped to avoid clashes).
meta_cols <- c("scenario", "habitat", "prop_cold", "harvest_rate",
               "pcmax_cold", "pcmax_warm", "A_cold", "A_warm", "harvestable_size",
               "rep", "sim_seed")

meta_tbl <- exp_design |>
  # harvestable_size is a scenario parameter, not a design column
  mutate(harvestable_size = vapply(exp_scenarios[scenario],
                                   \(p) p$harvestable_size, numeric(1))) |>
  select(all_of(meta_cols))

# Prepend metadata for one scenario to an output table
add_meta <- function(df, meta_row) {
  meta_row <- meta_row |> select(-any_of(setdiff(names(df), "scenario")))
  df |>
    mutate(scenario = meta_row$scenario) |>
    select(-any_of(setdiff(names(meta_row), "scenario"))) |>
    bind_cols(meta_row |> select(-scenario)) |>
    relocate(all_of(c("scenario", setdiff(names(meta_row), "scenario"))))
}

# Coerce a result element to a data frame (extinct is a scalar flag)
as_output_df <- function(res, name) {
  x <- res[[name]]
  if (name == "extinct") {
    return(tibble(extinct         = isTRUE(x),
                  extinction_date = as.Date(res$extinction_date %||% NA)))
  }
  if (is.null(x)) return(NULL)
  x <- as_tibble(x)

  # Only parent fish are used downstream (age and harvest-eligibility of spawners)
  if (name == "fish_registry") {
    parents <- unique(res$spawn_log$parent_pid)
    x <- filter(x, pid %in% parents)
    if (nrow(x) < length(parents)) {
      warning(length(parents) - nrow(x), " parent_pid(s) not found in fish_registry",
              call. = FALSE)
    }
  }
  x
}

# ── Stage 1: per-scenario intermediates ────────────────────────────────────────
for (out in outputs) dir.create(file.path(tmp_dir, out), showWarnings = FALSE)
done_dir <- file.path(tmp_dir, "_done")
dir.create(done_dir, showWarnings = FALSE)

processed     <- character()
missing_files <- character()

for (nm in names(exp_scenarios)) {
  path <- file.path(results_dir, paste0(nm, ".rds"))
  if (!file.exists(path)) {
    warning("No results file for '", nm, "'; skipping.", call. = FALSE)
    missing_files <- c(missing_files, nm)
    next
  }
  processed <- c(processed, nm)

  marker <- file.path(done_dir, nm)
  if (resume && file.exists(marker)) {
    message("Already processed (resume): ", nm)
    next
  }

  message("Processing: ", nm)
  res      <- readRDS(path)
  meta_row <- meta_tbl |> filter(scenario == nm)

  for (out in outputs) {
    df <- as_output_df(res, out)
    if (is.null(df)) {
      warning("'", out, "' missing in ", nm, call. = FALSE)
      next
    }
    # Extinct runs can have empty tables (even 0 columns); skip them so they
    # cannot cause type conflicts when binding. Loaders substitute zero-row tables.
    if (nrow(df) == 0 || ncol(df) == 0) {
      message("  '", out, "' is empty for ", nm, "; skipping.")
      next
    }
    saveRDS(add_meta(df, meta_row), file.path(tmp_dir, out, paste0(nm, ".rds")))
  }
  file.create(marker)

  rm(res, df)
  invisible(gc())
}

# ── Stage 2: combine per output ────────────────────────────────────────────────
combined_paths <- character()

for (out in outputs) {
  # Design order, not alphabetical
  files <- file.path(tmp_dir, out, paste0(processed, ".rds"))
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    message("No data for '", out, "'; not written.")
    next
  }

  message("Combining: ", out, " (", length(files), " scenario(s))")
  combined <- purrr::map(files, readRDS) |> bind_rows()
  out_path <- file.path(out_dir, paste0(out, ".rds"))
  saveRDS(combined, out_path)
  combined_paths[out] <- out_path

  unlink(file.path(tmp_dir, out), recursive = TRUE)
  rm(combined)
  invisible(gc())
}

# Remove the scratch folder once every output has been combined
# (only if no intermediate .rds files are left, e.g. after a failed combine)
if (length(list.files(tmp_dir, pattern = "\\.rds$", recursive = TRUE)) == 0) {
  unlink(tmp_dir, recursive = TRUE)
}

# ── Summary ────────────────────────────────────────────────────────────────────
message("\nProcessed ", length(processed), " of ", nrow(meta_tbl), " scenarios.")
if (length(missing_files) > 0) {
  message("Missing results: ", paste(missing_files, collapse = ", "))
}
message("Files written to ", out_dir, ":")
print(file.info(combined_paths)["size"] |>
        mutate(size_MB = round(size / 1e6, 1), size = NULL))
