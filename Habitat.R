#' ---
#' title: Build Habitat
#' ---
#' 


library(tidyverse)
library(lubridate)
library(ggpubr)
library(egg)


build_habitat <- function(params = list()) {

  # ── Timeline ────────────────────────────────────────────────────────────────
  nyears_burnin     <- params$nyears_burnin     %||% 20
  nyears_experiment <- params$nyears_experiment %||% 20
  mindate           <- params$mindate           %||% as.Date("2001-05-01")

  maxdate <- mindate + years(nyears_burnin) + years(nyears_experiment)
  dates   <- seq(mindate, maxdate, by = "day")
  n       <- length(dates)
  doy     <- as.numeric(format(dates, "%j"))

  # ── Temperature ─────────────────────────────────────────────────────────────
  # Each patch is an independent annual sine wave with optional daily noise.
  # Defaults for the cold patch preserve the original cold_multiplier = 0.6 behaviour
  # (base_temp_cold = 11 * 0.6 = 6.6; amplitude_cold = 14 * 0.6 = 8.4).

  # Warm patch
  base_temp_warm  <- params$base_temp_warm  %||% 11     # mean annual temperature (°C)
  amplitude_warm  <- params$amplitude_warm  %||% 14     # seasonal amplitude (°C)
  peak_doy_warm   <- params$peak_doy_warm   %||% 213    # day of peak temperature (Aug 1)
  temp_min_warm   <- params$temp_min_warm   %||% 0.5    # minimum temperature floor (°C)
  temp_noise_warm <- params$temp_noise_warm %||% 0      # daily noise SD (0 = deterministic)

  # Cold patch
  base_temp_cold  <- params$base_temp_cold  %||% 6.6
  amplitude_cold  <- params$amplitude_cold  %||% 8.4
  peak_doy_cold   <- params$peak_doy_cold   %||% 213
  temp_min_cold   <- params$temp_min_cold   %||% 0.3
  temp_noise_cold <- params$temp_noise_cold %||% 0

  # Simulate one patch: sine wave + noise, floored at temp_min, capped at 25 °C
  sim_temp <- function(doy, base_temp, amplitude, peak_doy, temp_min, noise_sd, n) {
    phase <- peak_doy - 365 / 4   # shift so sin() peaks on peak_doy
    temps <- base_temp + amplitude * sin(2 * pi * (doy - phase) / 365)
    noise <- rnorm(n, 0, noise_sd)
    pmax(temp_min, pmin(25, temps + noise))
  }

  set.seed(params$hab_seed %||% 123)
  simulated_temp      <- sim_temp(doy, base_temp_warm, amplitude_warm,
                                  peak_doy_warm, temp_min_warm, temp_noise_warm, n)
  simulated_temp_cold <- sim_temp(doy, base_temp_cold, amplitude_cold,
                                  peak_doy_cold, temp_min_cold, temp_noise_cold, n)

  habitat_df <- tibble(
    date     = dates,
    doy      = doy,
    dayofsim = seq_along(dates),
    temp_warm = round(simulated_temp,      digits = 1),
    temp_cold = round(simulated_temp_cold, digits = 1)
  )

  # ── Patch parameters (burn-in values) ───────────────────────────────────────
  # Use [[ for exact matching — $ partial-matches e.g. params$A_cold to A_cold_target.
  # Stored as time series columns so they can vary during the experiment period.
  habitat_df <- habitat_df |>
    mutate(
      pcmax_warm = params[["pcmax_warm"]]  %||% 0.5,    # max proportion of C_max
      pcmax_cold = params[["pcmax_cold"]]  %||% 0.5,
      A_warm     = params[["A_warm"]]      %||% 1.0,    # patch area (arbitrary units)
      A_cold     = params[["A_cold"]]      %||% 1.0,
      K_warm     = params[["K_warm"]]      %||% 500,    # half-saturation density (fish/area)
      K_cold     = params[["K_cold"]]      %||% 500,
      S_max_warm = params[["S_max_warm"]]  %||% 0.9994, # max daily survival probability
      S_max_cold = params[["S_max_cold"]]  %||% 0.9994
    )

  # ── Experiment period ramps ──────────────────────────────────────────────────
  # After the burn-in period, one or more patch parameters can be ramped to new
  # target values. Set a target to NA to leave that parameter unchanged.
  exp_start_date <- mindate + years(nyears_burnin)
  d_exp_start    <- which(habitat_df$date == exp_start_date)
  ramp_years     <- params$ramp_years %||% nyears_experiment

  # Linear ramp from end-of-burn-in value to target_val over ramp_yrs years,
  # then holds target_val for the remainder. ramp_yrs = 0 gives a step change.
  apply_ramp <- function(col, d_start, target_val, ramp_yrs) {
    baseline  <- col[d_start - 1]
    n_ramp    <- max(round(ramp_yrs * 365.25), 1)
    exp_idx   <- d_start:length(col)
    ramp_frac <- pmin(seq_along(exp_idx) / n_ramp, 1)
    col[exp_idx] <- baseline + (target_val - baseline) * ramp_frac
    col
  }

  params_to_change <- list(
    pcmax_warm = params$pcmax_warm_target %||% NA_real_,
    pcmax_cold = params$pcmax_cold_target %||% NA_real_,
    A_warm     = params$A_warm_target     %||% NA_real_,
    A_cold     = params$A_cold_target     %||% NA_real_,
    K_warm     = params$K_warm_target     %||% NA_real_,
    K_cold     = params$K_cold_target     %||% NA_real_,
    S_max_warm = params$S_max_warm_target %||% NA_real_,
    S_max_cold = params$S_max_cold_target %||% NA_real_
  )

  for (param in names(params_to_change)) {
    target <- params_to_change[[param]]
    if (!is.na(target)) {
      habitat_df[[param]] <- apply_ramp(habitat_df[[param]], d_exp_start, target, ramp_years)
    }
  }

  habitat_df
}


plot_habitat <- function(habitat_df, exp_start_date = NULL) {
  
  # Second calendar year in the data — used for the single-year temperature view
  ref_year <- year(min(habitat_df$date)) + 1

  # Jun–Aug shading rectangles
  summer_shading <- data.frame(
    xmin = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-06-01")),
    xmax = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-08-31")),
    ymin = -Inf, ymax = Inf
  )

  # Helper: optionally add experiment-start vline
  vline <- if (!is.null(exp_start_date)) {
    geom_vline(xintercept = exp_start_date, linetype = "dashed", color = "grey40")
  } else {
    NULL
  }

  # ── Temperature: full time series ───────────────────────────────────────────
  p_temp_full <- habitat_df |> ggplot() +
    geom_rect(data = summer_shading,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
    vline +
    geom_line(aes(x = date, y = temp_warm), color = "red") +
    geom_line(aes(x = date, y = temp_cold), color = "blue") +
    theme_bw() +
    xlab("Date") +
    ylab("Simulated daily stream temperature (°C)") +
    labs(title = "Daily stream temperature with Jun–Aug shaded") +
    ylim(0,27)

  # ── Temperature: single-year view ───────────────────────────────────────────
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
    geom_rect(data = season_shading,
              aes(xmin = xmin, xmax = xmax, fill = season),
              ymin = -Inf, ymax = Inf, alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(
      values = c(Winter = "#AED6F1", Spring = "#A9DFBF",
                 Summer = "#F9E79F", Autumn = "#F0B27A"),
      name   = NULL
    ) +
    geom_line(aes(x = date, y = temp_warm), color = "red") +
    geom_line(aes(x = date, y = temp_cold), color = "blue") +
    theme_bw() +
    theme(legend.position = "top") +
    xlab("Date") +
    ylab("Simulated daily stream temperature (°C)") +
    labs(title = paste("Single-year view:", ref_year)) +
    ylim(0,27)

  # ── Patch parameters: 2×2 panel ─────────────────────────────────────────────
  p_pcmax <- habitat_df |> ggplot() +
    vline +
    geom_line(aes(x = date, y = pcmax_warm), color = "red") +
    geom_line(aes(x = date, y = pcmax_cold), color = "blue") +
    theme_bw() +
    xlab("Date") + ylab("Max. P_Cmax") +
    ylim(0,1)

  p_area <- habitat_df |> ggplot() +
    vline +
    geom_line(aes(x = date, y = A_warm), color = "red") +
    geom_line(aes(x = date, y = A_cold), color = "blue") +
    theme_bw() +
    xlab("Date") + ylab("Patch area (unitless)") +
    ylim(0,2)

  p_K <- habitat_df |> ggplot() +
    vline +
    geom_line(aes(x = date, y = K_warm), color = "red") +
    geom_line(aes(x = date, y = K_cold), color = "blue") +
    theme_bw() +
    xlab("Date") + ylab("Strength of density-dependence\n(half-saturation density)") +
    ylim(0,1000)

  p_smax <- habitat_df |> ggplot() +
    vline +
    geom_line(aes(x = date, y = S_max_warm), color = "red") +
    geom_line(aes(x = date, y = S_max_cold), color = "blue") +
    theme_bw() +
    xlab("Date") + ylab("Max. daily probability of survival") +
    ylim(0.998,1)

  print(p_temp_full)
  print(p_temp_year)
  egg::ggarrange(p_pcmax, p_area, p_K, p_smax, nrow = 2, ncol = 2)

  invisible(list(
    temp_full = p_temp_full,
    temp_year = p_temp_year,
    pcmax     = p_pcmax,
    area      = p_area,
    K         = p_K,
    S_max     = p_smax
  ))
}


hab_params <- list(
  # Timeline
  nyears_burnin     = 20,
  nyears_experiment = 20,
  mindate           = as.Date("2001-05-01"),

  # Temperature — warm patch
  base_temp_warm  = 15,
  amplitude_warm  = 9.5,
  peak_doy_warm   = 196,   # July 15
  temp_min_warm   = 0.5,
  temp_noise_warm = 0,

  # Temperature — cold patch
  base_temp_cold  = 6.6,
  amplitude_cold  = 8.4,
  peak_doy_cold   = 196,
  temp_min_cold   = 0.1,
  temp_noise_cold = 0,

  # Patch burn-in values
  pcmax_warm = 0.5,
  pcmax_cold = 0.5,
  A_warm     = 1.0,
  A_cold     = 1.0,
  K_warm     = 500,
  K_cold     = 500,
  S_max_warm = 0.9994,
  S_max_cold = 0.9994,

  # Experiment targets (NA = no change)
  pcmax_warm_target = NA,
  pcmax_cold_target = NA,
  A_warm_target     = 2,
  A_cold_target     = 0,
  K_warm_target     = NA,
  K_cold_target     = NA,
  S_max_warm_target = NA,
  S_max_cold_target = NA,
  ramp_years        = 20,

  hab_seed = 123
)

habitat_df     <- build_habitat(hab_params)
exp_start_date <- hab_params$mindate + years(hab_params$nyears_burnin)


plot_habitat(habitat_df, exp_start_date)


str(habitat_df)
head(habitat_df)


