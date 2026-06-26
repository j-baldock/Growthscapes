#' ---
#' title: Define Scenarios & Experiments
#' ---
#' 


library(tidyverse)
library(lubridate)
library(ggpubr)
library(egg)

#| cache: false
source("Habitat.R")


base_params <- list(

  # ── Timeline ──────────────────────────────────────────
  nyears_burnin     = 40,
  nyears_experiment = 60,
  ramp_years        = 60, # ramp years should almost always be set to the same as nyears_experiment
  mindate           = as.Date("2001-05-01"),

  # ── Temperature — warm patch ───────────────────────────
  base_temp_warm  = 11,
  amplitude_warm  = 14,
  peak_doy_warm   = 196,
  temp_min_warm   = 0.5,
  temp_noise_warm = 0,

  # ── Temperature — cold patch ───────────────────────────
  base_temp_cold  = 6.6,
  amplitude_cold  = 8.4,
  peak_doy_cold   = 196,
  temp_min_cold   = 0.3,
  temp_noise_cold = 0,

  # ── Patch burn-in values ───────────────────────────────
  pcmax_warm = 0.5,  pcmax_cold = 0.5,
  A_warm     = 1.0,  A_cold     = 1.0,
  K_warm     = 500,  K_cold     = 500,
  S_max_warm = 0.9994, S_max_cold = 0.9994,

  # ── Experiment targets (NA = no change) ───────────────
  pcmax_warm_target = NA,  pcmax_cold_target = NA,
  A_warm_target     = 2,   A_cold_target     = 0,
  K_warm_target     = NA,  K_cold_target     = NA,
  S_max_warm_target = NA,  S_max_cold_target = NA,

  # ── Behavioural / movement ─────────────────────────────
  move_stochastic   = "prob",
  food_densdepen    = "hyperbolic",
  sense_environment = "density_current",

  # ── Population ─────────────────────────────────────────
  n_fish_per_strat  = 50,
  start_wt          = 0.5,
  wt_sd             = 0.05,

  # ── Movement cost ──────────────────────────────────────
  movecost_c   = 0.1,
  movecost_b   = 0.3,
  sigma_bold   = 0.002,
  tau          = 0.001,

  # ── Survival ──────────────────────────────────────────
  T1_mort           = 30,
  T9_mort           = 25.8,
  K9_starv          = 0.55,
  K1_starv          = 0.45,
  MaxDensity4Growth = 50,

  # ── Competition ───────────────────────────────────────
  dominance_beta = 1,

  # ── Reproduction ──────────────────────────────────────
  egg_wt       = 0.07,
  repro_cost   = 0.2,
  sigma_fecund = 0.085,
  egg_surv     = 0.1,

  # ── Seeds ─────────────────────────────────────────────
  hab_seed = 123,
  sim_seed = 7843
)


scenarios <- list(
  null_cold                   = modifyList(base_params, list(A_warm = 0, 
                                                             A_warm_target = NA)),
  temp_mult                   = base_params,
  temp_offset                 = modifyList(base_params, list(base_temp_warm  = 15,
                                                             amplitude_warm  = 9.5,
                                                             temp_min_warm   = 0.5,
                                                             base_temp_cold  = 6.6,
                                                             amplitude_cold  = 8.4,
                                                             temp_min_cold   = 0)),
  temp_offset_diffP           = modifyList(base_params, list(base_temp_warm  = 15,
                                                             amplitude_warm  = 9.5,
                                                             temp_min_warm   = 0.5,
                                                             base_temp_cold  = 6.6,
                                                             amplitude_cold  = 8.4,
                                                             temp_min_cold   = 0,
                                                             pcmax_warm      = 0.6,
                                                             pcmax_cold      = 0.4))
)


plot_habitat(build_habitat(scenarios[["null_cold"]]), exp_start_date = scenarios[["null_cold"]]$mindate + years(scenarios[["null_cold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["temp_mult"]]), exp_start_date = scenarios[["temp_mult"]]$mindate + years(scenarios[["temp_mult"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["temp_offset"]]), exp_start_date = scenarios[["temp_offset"]]$mindate + years(scenarios[["temp_offset"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["temp_offset_diffP"]]), exp_start_date = scenarios[["temp_offset_diffP"]]$mindate + years(scenarios[["temp_offset_diffP"]]$nyears_burnin))


