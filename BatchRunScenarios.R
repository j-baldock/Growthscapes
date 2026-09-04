# BatchRunScenarios.R
#
# Run IBM simulations for one or all scenarios, saving results as .rds files
# in the results/ folder. Simulations should be run from this script rather
# than from within the book .qmd files so that results are not regenerated
# on every render.
#
# Prerequisites: ensure the following .R files are up to date by running
# qmd_to_r_script() in the console after any changes to their .qmd counterparts:
#   - quarto::qmd_to_r_script("Habitat.qmd")
#   - quarto::qmd_to_r_script("Functions.qmd")
#   - quarto::qmd_to_r_script("SimulationLoop.qmd")
#   - quarto::qmd_to_r_script("ScenarioDefs.qmd")

# ── Setup ──────────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(beepr)

source("Habitat.R")
source("Functions.R")
source("SimulationLoop.R")
source("ScenarioDefs.R")   # defines `scenarios` (list of named parameter sets)

load("data/wt.growth.array.RData")  # loads `wt.growth`

if (!dir.exists("results")) dir.create("results")

# ── run_scenario(): run a single named scenario ────────────────────────────────
# name      - must match a name in the `scenarios` list (e.g. "base")
# overwrite - if FALSE (default), skip silently when results already exist

run_scenario <- function(name, overwrite = FALSE) {
  out_path <- file.path("results", paste0(name, ".rds"))

  if (file.exists(out_path) && !overwrite) {
    message("Skipping '", name, "' — results file already exists. ",
            "Use overwrite = TRUE to re-run.")
    return(invisible(NULL))
  }

  if (!name %in% names(scenarios)) {
    stop("'", name, "' is not a defined scenario. ",
         "Available: ", paste(names(scenarios), collapse = ", "))
  }

  params     <- scenarios[[name]]
  habitat_df <- build_habitat(params)

  message("Running scenario: '", name, "' ...")
  t_start <- Sys.time()
  result <- run_simulation(habitat_df = habitat_df,
                           params     = params,
                           wt_growth  = wt.growth)
  elapsed_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  elapsed_fmt <- sprintf("%dh %02dm %02ds",
                         as.integer(elapsed_s %/% 3600),
                         as.integer((elapsed_s %% 3600) %/% 60),
                         as.integer(elapsed_s %% 60))

  saveRDS(result, out_path)
  message("Saved: ", out_path, " (run time: ", elapsed_fmt, ")")
  invisible(result)
}

# ── run_all_scenarios(): batch run over all defined scenarios ──────────────────
# Iterates over every entry in `scenarios`. Skips those that already have saved
# results unless overwrite = TRUE.

run_all_scenarios <- function(overwrite = FALSE) {
  nms <- names(scenarios)
  message("Batch run: ", length(nms), " scenario(s): ",
          paste(nms, collapse = ", "))
  for (name in nms) {
    run_scenario(name, overwrite = overwrite)
  }
  message("Batch run complete.")
}


# ── Individual scenario calls ──────────────────────────────────────────────────
# Uncomment a line below to run (or re-run) a single scenario:

# run_scenario("base")
run_scenario("temp_offset", overwrite = TRUE)
run_scenario("temp_offset_diffP", overwrite = TRUE)   # re-run and overwrite
run_scenario("temp_offset_fixedhab", overwrite = TRUE)   # re-run and overwrite

# fixed habitat: cold only (null)
run_scenario("TempOffset_ColdOnly_95percold", overwrite = TRUE)
run_scenario("TempOffset_ColdOnly_75percold", overwrite = TRUE)
run_scenario("TempOffset_ColdOnly_50percold", overwrite = TRUE)
run_scenario("TempOffset_ColdOnly_25percold", overwrite = TRUE)
run_scenario("TempOffset_ColdOnly_05percold", overwrite = TRUE)
# beep()

# fixed habitat: cold + warm (same pcmax)
run_scenario("TempOffset_ColdWarm_95percold", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_75percold", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_50percold", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_25percold", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_05percold", overwrite = TRUE)
beep()

# fixed habitat: cold + warm (same pcmax)
run_scenario("TempOffset_ColdWarm_95percold_highPwarm", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_75percold_highPwarm", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_50percold_highPwarm", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_25percold_highPwarm", overwrite = TRUE)
run_scenario("TempOffset_ColdWarm_05percold_highPwarm", overwrite = TRUE)

# ── Batch run ─────────────────────────────────────────────────────────────────
# Uncomment to run all scenarios (skips any with existing results):

st <- Sys.time()
# run_all_scenarios()
run_all_scenarios(overwrite = TRUE)   # re-run everything
Sys.time() - st
beep()
