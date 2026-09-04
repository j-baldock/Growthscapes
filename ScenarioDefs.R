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
  nyears_burnin     = 50,
  nyears_experiment = 50,
  ramp_years        = 0,
  mindate           = as.Date("2001-05-01"),

  # ── Temperature — warm patch ───────────────────────────
  base_temp_warm  = 15,
  amplitude_warm  = 9.5,
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
  K_warm     = 400,  K_cold     = 400, 
  S_max_warm = 0.9994, S_max_cold = 0.9994,

  # ── Experiment targets (NA = no change) ───────────────
  pcmax_warm_target = NA,  pcmax_cold_target = NA,
  A_warm_target     = NA,  A_cold_target     = NA,
  K_warm_target     = NA,  K_cold_target     = NA,
  S_max_warm_target = NA,  S_max_cold_target = NA,

  # ── Behavioural / movement ─────────────────────────────
  move_stochastic   = "prob",
  food_densdepen    = "hyperbolic",
  sense_environment = "density_all",

  # ── Population ─────────────────────────────────────────
  n_fish_per_strat  = 50,
  start_wt          = 0.5,
  wt_sd             = 0.05,

  # ── Movement cost ──────────────────────────────────────
  movecost_c   = 0.0, # previously 0.1. 0.0 is effective no cost
  movecost_b   = 0.3,
  sigma_bold   = 0.002,
  tau          = 0.001,

  # ── Survival ──────────────────────────────────────────
  s_min             = 0.96,
  s_w0              = 5,
  s_k               = 0.8,
  T1_mort           = 30,
  T9_mort           = 25.8,
  K9_starv          = 0.55,
  K1_starv          = 0.45,
  MaxDensity4Growth = 50,

  # Critical-period consumption-based survival (Elliott 1989)
  # Applies fncSurviveConsumption() to age-0 fry for their first crit_period_days.
  # Parameterised so that sustained pcmax_dd ≤ crit_pcmax_lo ≈ 1% 60-day survival;
  # sustained pcmax_dd ≥ crit_pcmax_hi ≈ 99% 60-day survival (negligible cost).
  crit_pcmax_lo    = 0.30,
  crit_pcmax_hi    = 0.40,
  crit_period_days = 60L,
  
  # ── Competition ───────────────────────────────────────
  dominance_beta             = 1,
  age_structured_competition = TRUE,

  # ── Reproduction ──────────────────────────────────────
  egg_wt       = 0.07,
  repro_cost   = 0.2,
  sigma_fecund = 0.085,
  egg_surv     = 0.15, 

  # ── Seeds ─────────────────────────────────────────────
  hab_seed = 123,
  sim_seed = 7843
)


scenarios <- list(
  # null_cold                   = modifyList(base_params, list(A_warm = 0, 
  #                                                            A_warm_target = NA)),
  # temp_mult                   = base_params,
  # temp_offset                 = modifyList(base_params, list(base_temp_warm  = 15,
  #                                                            amplitude_warm  = 9.5,
  #                                                            temp_min_warm   = 0.5,
  #                                                            base_temp_cold  = 6.6,
  #                                                            amplitude_cold  = 8.4,
  #                                                            temp_min_cold   = 0)),
  # temp_offset_diffP           = modifyList(base_params, list(base_temp_warm  = 15,
  #                                                            amplitude_warm  = 9.5,
  #                                                            temp_min_warm   = 0.5,
  #                                                            base_temp_cold  = 6.6,
  #                                                            amplitude_cold  = 8.4,
  #                                                            temp_min_cold   = 0,
  #                                                            pcmax_warm      = 0.6,
  #                                                            pcmax_cold      = 0.4)),
  
  # Fixed habitat scenarios: changing habitat availability is tested by comparing among scenarios
  ### Null (cold only)
  TempOffset_ColdOnly_95percold      = modifyList(base_params, list(A_cold    = 1.9,
                                                                    A_warm    = 0)),
  
  TempOffset_ColdOnly_75percold      = modifyList(base_params, list(A_cold    = 1.5,
                                                                    A_warm    = 0)),
  
  TempOffset_ColdOnly_50percold      = modifyList(base_params, list(A_cold    = 1.0,
                                                                    A_warm    = 0)),
  
  TempOffset_ColdOnly_25percold      = modifyList(base_params, list(A_cold    = 0.5,
                                                                    A_warm    = 0)),
  
  TempOffset_ColdOnly_05percold      = modifyList(base_params, list(A_cold    = 0.1,
                                                                    A_warm    = 0)),
  
  ### Cold + Warm, same Pcmax
  TempOffset_ColdWarm_95percold      = modifyList(base_params, list(A_cold    = 1.9,
                                                                    A_warm    = 0.1)),
  
  TempOffset_ColdWarm_75percold      = modifyList(base_params, list(A_cold    = 1.5,
                                                                    A_warm    = 0.5)),
  
  TempOffset_ColdWarm_50percold      = modifyList(base_params, list(A_cold    = 1.0,
                                                                    A_warm    = 1.0)),
  
  TempOffset_ColdWarm_25percold      = modifyList(base_params, list(A_cold    = 0.5,
                                                                    A_warm    = 1.5)),
  
  TempOffset_ColdWarm_05percold      = modifyList(base_params, list(A_cold    = 0.1,
                                                                    A_warm    = 1.9)),
  
  ### Cold + Warm, high Pcmax in warm
  TempOffset_ColdWarm_95percold_highPwarm      = modifyList(base_params, list(A_cold    = 1.9,
                                                                    A_warm    = 0.1,
                                                                    pcmax_warm = 0.6)),
  
  TempOffset_ColdWarm_75percold_highPwarm      = modifyList(base_params, list(A_cold    = 1.5,
                                                                    A_warm    = 0.5,
                                                                    pcmax_warm = 0.6)),
  
  TempOffset_ColdWarm_50percold_highPwarm      = modifyList(base_params, list(A_cold    = 1.0,
                                                                    A_warm    = 1.0,
                                                                    pcmax_warm = 0.6)),
  
  TempOffset_ColdWarm_25percold_highPwarm      = modifyList(base_params, list(A_cold    = 0.5,
                                                                    A_warm    = 1.5,
                                                                    pcmax_warm = 0.6)),
  
  TempOffset_ColdWarm_05percold_highPwarm      = modifyList(base_params, list(A_cold    = 0.1,
                                                                    A_warm    = 1.9,
                                                                    pcmax_warm = 0.6))
  
)


plot_habitat(build_habitat(scenarios[["TempOffset_ColdOnly_95percold"]]), exp_start_date = scenarios[["TempOffset_ColdOnly_95percold"]]$mindate + years(scenarios[["TempOffset_ColdOnly_95percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdOnly_75percold"]]), exp_start_date = scenarios[["TempOffset_ColdOnly_75percold"]]$mindate + years(scenarios[["TempOffset_ColdOnly_75percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdOnly_50percold"]]), exp_start_date = scenarios[["TempOffset_ColdOnly_50percold"]]$mindate + years(scenarios[["TempOffset_ColdOnly_50percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdOnly_25percold"]]), exp_start_date = scenarios[["TempOffset_ColdOnly_25percold"]]$mindate + years(scenarios[["TempOffset_ColdOnly_25percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdOnly_05percold"]]), exp_start_date = scenarios[["TempOffset_ColdOnly_05percold"]]$mindate + years(scenarios[["TempOffset_ColdOnly_05percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_95percold"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_95percold"]]$mindate + years(scenarios[["TempOffset_ColdWarm_95percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_75percold"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_75percold"]]$mindate + years(scenarios[["TempOffset_ColdWarm_75percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_50percold"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_50percold"]]$mindate + years(scenarios[["TempOffset_ColdWarm_50percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_25percold"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_25percold"]]$mindate + years(scenarios[["TempOffset_ColdWarm_25percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_05percold"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_05percold"]]$mindate + years(scenarios[["TempOffset_ColdWarm_05percold"]]$nyears_burnin))


