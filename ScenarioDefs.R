#' ---
#' title: Define Scenarios & Experiments
#' ---
#' 


library(tidyverse)
library(lubridate)
library(ggpubr)
library(egg)
library(patchwork)

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
  
  # ── Harvest ───────────────────────────────────────────
  harvest_rate     = 0,
  harvest_window   = c(91L, 91L),

  # ── Seeds ─────────────────────────────────────────────
  hab_seed = 123,
  sim_seed = 7843
)


scenarios <- list(

  ### Null (cold only)
  TempOffset_ColdOnly_95percold             = modifyList(base_params, list(A_cold = 1.9,  A_warm = 0)),
  TempOffset_ColdOnly_75percold             = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0)),
  TempOffset_ColdOnly_50percold             = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0)),
  TempOffset_ColdOnly_25percold             = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0)),
  TempOffset_ColdOnly_05percold             = modifyList(base_params, list(A_cold = 0.1,  A_warm = 0)),
  
  ### Cold + Warm, same Pcmax
  TempOffset_ColdWarm_95percold             = modifyList(base_params, list(A_cold = 1.9,  A_warm = 0.1)),
  TempOffset_ColdWarm_75percold             = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5)),
  TempOffset_ColdWarm_50percold             = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0)),
  TempOffset_ColdWarm_25percold             = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5)),
  TempOffset_ColdWarm_05percold             = modifyList(base_params, list(A_cold = 0.1,  A_warm = 1.9)),
  
  ### Cold + Warm, high Pcmax in warm
  TempOffset_ColdWarm_95percold_highPwarm   = modifyList(base_params, list(A_cold = 1.9,  A_warm = 0.1,  pcmax_warm = 0.6)),
  TempOffset_ColdWarm_75percold_highPwarm   = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  pcmax_warm = 0.6)),
  TempOffset_ColdWarm_50percold_highPwarm   = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  pcmax_warm = 0.6)),
  TempOffset_ColdWarm_25percold_highPwarm   = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  pcmax_warm = 0.6)),
  TempOffset_ColdWarm_05percold_highPwarm   = modifyList(base_params, list(A_cold = 0.1,  A_warm = 1.9,  pcmax_warm = 0.6)),
  
  
  
  # Harvested
  
  ## 75% cold (25% warm)
  #### Cold only
  Harvest_10_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.10,   harvest_window = c(152L, 273L))),
  Harvest_20_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.20,   harvest_window = c(152L, 273L))),
  Harvest_30_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.30,   harvest_window = c(152L, 273L))),
  Harvest_40_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.40,   harvest_window = c(152L, 273L))),
  Harvest_50_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.50,   harvest_window = c(152L, 273L))),
  Harvest_55_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.55,   harvest_window = c(152L, 273L))),
  Harvest_60_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.60,   harvest_window = c(152L, 273L))),
  Harvest_65_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.65,   harvest_window = c(152L, 273L))),
  Harvest_70_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.70,   harvest_window = c(152L, 273L))),
  Harvest_75_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.75,   harvest_window = c(152L, 273L))),
  Harvest_80_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.80,   harvest_window = c(152L, 273L))),
  Harvest_90_ColdOnly_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0,    harvest_rate = 0.90,   harvest_window = c(152L, 273L))),
  #### Cold + Warm
  Harvest_10_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.10,   harvest_window = c(152L, 273L))),
  Harvest_20_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.20,   harvest_window = c(152L, 273L))),
  Harvest_30_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.30,   harvest_window = c(152L, 273L))),
  Harvest_40_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.40,   harvest_window = c(152L, 273L))),
  Harvest_50_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.50,   harvest_window = c(152L, 273L))),
  Harvest_55_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.55,   harvest_window = c(152L, 273L))),
  Harvest_60_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.60,   harvest_window = c(152L, 273L))),
  Harvest_65_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.65,   harvest_window = c(152L, 273L))),
  Harvest_70_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.70,   harvest_window = c(152L, 273L))),
  Harvest_75_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.75,   harvest_window = c(152L, 273L))),
  Harvest_80_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.80,   harvest_window = c(152L, 273L))),
  Harvest_90_ColdWarm_75percold      = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.90,   harvest_window = c(152L, 273L))),
  #### Cold + Warm (high P in warm)
  Harvest_10_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.10,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_20_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.20,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_30_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.30,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_40_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.40,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_50_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.50,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_55_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.55,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_60_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.60,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_65_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.65,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_70_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.70,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_75_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.75,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_80_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.80,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_90_ColdWarmHigh_75percold  = modifyList(base_params, list(A_cold = 1.5,  A_warm = 0.5,  harvest_rate = 0.90,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  
  
  ## 50% cold (50% warm)
  #### Cold only
  Harvest_10_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.10,   harvest_window = c(152L, 273L))),
  Harvest_20_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.20,   harvest_window = c(152L, 273L))),
  Harvest_30_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.30,   harvest_window = c(152L, 273L))),
  Harvest_40_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.40,   harvest_window = c(152L, 273L))),
  Harvest_50_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.50,   harvest_window = c(152L, 273L))),
  Harvest_55_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.55,   harvest_window = c(152L, 273L))),
  Harvest_60_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.60,   harvest_window = c(152L, 273L))),
  Harvest_65_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.65,   harvest_window = c(152L, 273L))),
  Harvest_70_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.70,   harvest_window = c(152L, 273L))),
  Harvest_75_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.75,   harvest_window = c(152L, 273L))),
  Harvest_80_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.80,   harvest_window = c(152L, 273L))),
  Harvest_90_ColdOnly_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 0,    harvest_rate = 0.90,   harvest_window = c(152L, 273L))),
  #### Cold + Warm
  Harvest_10_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.10,   harvest_window = c(152L, 273L))),
  Harvest_20_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.20,   harvest_window = c(152L, 273L))),
  Harvest_30_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.30,   harvest_window = c(152L, 273L))),
  Harvest_40_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.40,   harvest_window = c(152L, 273L))),
  Harvest_50_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.50,   harvest_window = c(152L, 273L))),
  Harvest_55_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.55,   harvest_window = c(152L, 273L))),
  Harvest_60_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.60,   harvest_window = c(152L, 273L))),
  Harvest_65_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.65,   harvest_window = c(152L, 273L))),
  Harvest_70_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.70,   harvest_window = c(152L, 273L))),
  Harvest_75_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.75,   harvest_window = c(152L, 273L))),
  Harvest_80_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.80,   harvest_window = c(152L, 273L))),
  Harvest_90_ColdWarm_50percold      = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.90,   harvest_window = c(152L, 273L))),
  #### Cold + Warm (high P in warm)
  Harvest_10_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.10,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_20_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.20,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_30_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.30,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_40_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.40,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_50_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.50,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_55_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.55,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_60_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.60,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_65_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.65,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_70_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.70,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_75_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.75,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_80_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.80,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_90_ColdWarmHigh_50percold  = modifyList(base_params, list(A_cold = 1.0,  A_warm = 1.0,  harvest_rate = 0.90,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  
  
  ## 25% cold (75% warm)
  #### Cold only
  Harvest_10_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.10,   harvest_window = c(152L, 273L))),
  Harvest_20_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.20,   harvest_window = c(152L, 273L))),
  Harvest_30_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.30,   harvest_window = c(152L, 273L))),
  Harvest_40_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.40,   harvest_window = c(152L, 273L))),
  Harvest_50_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.50,   harvest_window = c(152L, 273L))),
  Harvest_55_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.55,   harvest_window = c(152L, 273L))),
  Harvest_60_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.60,   harvest_window = c(152L, 273L))),
  Harvest_65_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.65,   harvest_window = c(152L, 273L))),
  Harvest_70_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.70,   harvest_window = c(152L, 273L))),
  Harvest_75_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.75,   harvest_window = c(152L, 273L))),
  Harvest_80_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.80,   harvest_window = c(152L, 273L))),
  Harvest_90_ColdOnly_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 0,    harvest_rate = 0.90,   harvest_window = c(152L, 273L))),
  #### Cold + Warm
  Harvest_10_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.10,   harvest_window = c(152L, 273L))),
  Harvest_20_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.20,   harvest_window = c(152L, 273L))),
  Harvest_30_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.30,   harvest_window = c(152L, 273L))),
  Harvest_40_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.40,   harvest_window = c(152L, 273L))),
  Harvest_50_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.50,   harvest_window = c(152L, 273L))),
  Harvest_55_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.55,   harvest_window = c(152L, 273L))),
  Harvest_60_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.60,   harvest_window = c(152L, 273L))),
  Harvest_65_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.65,   harvest_window = c(152L, 273L))),
  Harvest_70_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.70,   harvest_window = c(152L, 273L))),
  Harvest_75_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.75,   harvest_window = c(152L, 273L))),
  Harvest_80_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.80,   harvest_window = c(152L, 273L))),
  Harvest_90_ColdWarm_25percold      = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.90,   harvest_window = c(152L, 273L))),
  #### Cold + Warm (high P in warm)
  Harvest_10_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.10,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_20_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.20,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_30_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.30,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_40_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.40,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_50_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.50,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_55_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.55,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_60_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.60,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_65_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.65,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_70_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.70,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_75_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.75,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_80_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.80,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6)),
  Harvest_90_ColdWarmHigh_25percold  = modifyList(base_params, list(A_cold = 0.5,  A_warm = 1.5,  harvest_rate = 0.90,   harvest_window = c(152L, 273L),  pcmax_warm = 0.6))
  
  
  
  
)


# ── Timeline ────────────────────────────────────────────────────────────────
nyears_burnin     <- base_params$nyears_burnin
nyears_experiment <- base_params$nyears_experiment
mindate           <- base_params$mindate
maxdate           <- mindate + years(nyears_burnin) + years(nyears_experiment)
dates             <- seq(mindate, maxdate, by = "day")
n                 <- length(dates)
doy               <- as.numeric(format(dates, "%j"))


sim_temp <- function(doy, base_temp, amplitude, peak_doy, temp_min, noise_sd, n) {
    phase <- peak_doy - 365 / 4   # shift so sin() peaks on peak_doy
    temps <- base_temp + amplitude * sin(2 * pi * (doy - phase) / 365)
    noise <- rnorm(n, 0, noise_sd)
    pmax(temp_min, pmin(25, temps + noise))
  }

simulated_temp      <- sim_temp(doy, base_params$base_temp_warm, base_params$amplitude_warm,
                                base_params$peak_doy_warm, base_params$temp_min_warm, base_params$temp_noise_warm, n)
simulated_temp_cold <- sim_temp(doy, base_params$base_temp_cold, base_params$amplitude_cold,
                                base_params$peak_doy_cold, base_params$temp_min_cold, base_params$temp_noise_cold, n)

habitat_df <- tibble(
  date     = dates,
  doy      = doy,
  dayofsim = seq_along(dates),
  temp_warm = round(simulated_temp,      digits = 1),
  temp_cold = round(simulated_temp_cold, digits = 1)
)
ref_year <- year(min(habitat_df$date)) + 1

#| cache: false
#| label: fig-habitat-conceptual
#| fig-cap: Habitat area composition as a function of the proportion of cold habitat.
#|   Total habitat area is fixed at 2 units.
#| fig-width: 9
#| fig-height: 3.5

season_shading <- data.frame(
  xmin   = as.Date(c(paste0(ref_year, "-01-01"),
                     paste0(ref_year, "-03-01"),
                     paste0(ref_year, "-06-01"),
                     paste0(ref_year, "-09-01"),
                     paste0(ref_year, "-12-01"))),
  xmax   = c(as.Date(paste0(ref_year,     "-03-01")) - 1,  # last day of Feb (leap-year safe)
             as.Date(paste0(ref_year,     "-06-01")) - 1,
             as.Date(paste0(ref_year,     "-09-01")) - 1,
             as.Date(paste0(ref_year,     "-12-01")) - 1,
             as.Date(paste0(ref_year + 1, "-01-01")) - 1),
  season = factor(c("Winter", "Spring", "Summer", "Autumn", "Winter"),
                  levels = c("Winter", "Spring", "Summer", "Autumn"))
)

p_temp_year <- habitat_df |>
  filter(year(date) == ref_year) |>
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_rect(data = season_shading,
            aes(xmin = xmin, xmax = xmax, fill = season),
            ymin = -Inf, ymax = Inf, alpha = 0.5, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(Winter = "#AED6F1", Spring = "#A9DFBF",
               Summer = "#F9E79F", Autumn = "#F0B27A"),
    name   = "Season"
  ) +
  geom_line(aes(x = date, y = temp_warm), color = "red", linewidth = 1) +
  geom_line(aes(x = date, y = temp_cold), color = "blue", linewidth = 1) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b",
    expand = expansion(0)
  ) +
  theme_bw() +
  theme(legend.position = "right", panel.grid = element_blank(), text = element_text(color = "black")) +
  xlab("Date") +
  ylab("Daily stream temperature (°C)") +
  ylim(0,25)


p_area_percent <- tibble(
  prop_cold = seq(0, 1, by = 0.1)
) |>
  mutate(
    cold_area = prop_cold * 2,
    warm_area = (1 - prop_cold) * 2
  ) |>
  pivot_longer(cols = c(cold_area, warm_area), names_to = "habitat", values_to = "area") |>
  mutate(
    habitat = factor(habitat, levels = c("warm_area", "cold_area"),
                     labels = c("Seasonally warm", "Perennially cold")),
    prop_cold_label = factor(prop_cold, labels = paste0(seq(0, 100, 10)))
  ) |>
  ggplot(aes(x = prop_cold_label, y = area, fill = habitat)) +
  geom_col(alpha = 0.8) +
  geom_col(
    data = ~ filter(.x, habitat == "Perennially cold"),
    aes(color = "Cold-only"),
    fill = NA, linewidth = 0.6
  ) +
  scale_x_discrete(limits = rev) +
  scale_fill_manual(values = c("Perennially cold" = "blue", "Seasonally warm" = "red")) +
  scale_color_manual(values = c("Cold-only" = "black"), name = "Habitat access") +
  labs(
    x = "Proportion cold habitat, % (decreasing)",
    y = "Habitat area (unitless)",
    fill = "Habitat type"
  ) +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(color = "black"))

p_temp_year + p_area_percent +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


plot_habitat(build_habitat(scenarios[["TempOffset_ColdOnly_75percold"]]), exp_start_date = scenarios[["TempOffset_ColdOnly_75percold"]]$mindate + years(scenarios[["TempOffset_ColdOnly_75percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_75percold"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_75percold"]]$mindate + years(scenarios[["TempOffset_ColdWarm_75percold"]]$nyears_burnin))


plot_habitat(build_habitat(scenarios[["TempOffset_ColdWarm_75percold_highPwarm"]]), exp_start_date = scenarios[["TempOffset_ColdWarm_75percold_highPwarm"]]$mindate + years(scenarios[["TempOffset_ColdWarm_75percold_highPwarm"]]$nyears_burnin))


