#' ---
#' title: Functions to Plot Results
#' ---
#' 


plot_mechanistic <- function(ibm_long, habitat_df, sim_params,
                              spawn_log, fish_registry,
                              exp_start_date = NULL) {

  # ── Shared objects ──────────────────────────────────────────────────────────
  T1_mort   <- sim_params$T1_mort   %||% 30
  T9_mort   <- sim_params$T9_mort   %||% 25.8
  K9_starv  <- sim_params$K9_starv  %||% 0.55
  K1_starv  <- sim_params$K1_starv  %||% 0.45
  start_wt  <- sim_params$start_wt  %||% 0.5

  strategy_colors <- c(
    resident_warm = "red",
    resident_cold = "blue",
    optimal_mover = "darkgreen"
  )

  summer_shading <- data.frame(
    xmin = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-06-01")),
    xmax = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-08-31")),
    ymin = -Inf, ymax = Inf
  )

  vline <- if (!is.null(exp_start_date))
    geom_vline(xintercept = exp_start_date, linetype = "dashed", color = "grey40")
  else
    NULL

  wt_seq <- seq(0.1, 25, 0.1)


  # ── 1. BIOENERGETICS ────────────────────────────────────────────────────────
  bioe_output <- ibm_long |>
    filter(survived == 1) |>
    left_join(
      habitat_df |> select(dayofsim, pcmax_warm, pcmax_cold, K_warm, K_cold),
      by = "dayofsim"
    ) |>
    group_by(dayofsim, patch) |>
    mutate(density = n()) |>
    ungroup() |>
    mutate(
      func_temp         = fncTempDepend(temp),
      pcmax_baseline    = ifelse(patch == "warm", pcmax_warm, pcmax_cold),
      pcmax_adjusted    = pmin(pcmax_baseline, func_temp),
      cmax_allometric   = fncAllomCmax(weight),
      k                 = if_else(patch == "warm", K_warm, K_cold),
      pcmax_adjusted_dd = pcmax_adjusted * (k / (k + density - 1)),
      ration            = cmax_allometric * pcmax_adjusted_dd
    )

  pcmax_warm_ref <- habitat_df$pcmax_warm[1]
  pcmax_cold_ref <- habitat_df$pcmax_cold[1]

  # Allometric Cmax reference line
  constants <- fncReadConstants()
  xs <- seq(min(ibm_long$weight, na.rm = TRUE),
            max(ibm_long$weight, na.rm = TRUE),
            length.out = 100)
  ys <- constants$Consumption$CA * (xs ^ constants$Consumption$CB)

  p_ration_warm <- ggplot() +
    geom_point(data = bioe_output |> filter(patch == "warm"),
               aes(x = weight, y = ration, color = temp), alpha = 0.3) +
    geom_line(aes(x = xs, y = ys * pcmax_warm_ref), color = "black") +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "top") +
    xlab("Fish mass (g)") + ylab("Ration size (g/g/d)") +
    labs(title = "Bioenergetics — warm patch", color = "Temp (°C)")

  p_ration_cold <- ggplot() +
    geom_point(data = bioe_output |> filter(patch == "cold"),
               aes(x = weight, y = ration, color = temp), alpha = 0.3) +
    geom_line(aes(x = xs, y = ys * pcmax_cold_ref), color = "black") +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "top") +
    xlab("Fish mass (g)") + ylab("Ration size (g/g/d)") +
    labs(title = "Bioenergetics — cold patch", color = "Temp (°C)")

  p_pcmax_warm <- ggplot() +
    geom_point(data = bioe_output |> filter(patch == "warm"),
               aes(x = temp, y = pcmax_adjusted), alpha = 0.1) +
    geom_abline(slope = 0, intercept = pcmax_warm_ref,
                color = "red", linetype = "dashed") +
    theme_bw() +
    xlab("Temperature (°C)") + ylab("Temperature-adjusted P_Cmax") +
    labs(title = "P_Cmax — warm patch")

  p_pcmax_cold <- ggplot() +
    geom_point(data = bioe_output |> filter(patch == "cold"),
               aes(x = temp, y = pcmax_adjusted), alpha = 0.1) +
    geom_abline(slope = 0, intercept = pcmax_cold_ref,
                color = "red", linetype = "dashed") +
    theme_bw() +
    xlab("Temperature (°C)") + ylab("Temperature-adjusted P_Cmax") +
    labs(title = "P_Cmax — cold patch")


  # ── 2. TEMPERATURE EXPERIENCED ──────────────────────────────────────────────
  p_temp_full <- ibm_long |>
    filter(strategy == "optimal_mover") |>
    ggplot(aes(x = date, y = temp, color = strategy, group = pid)) +
    geom_line(alpha = 0.4) +
    scale_color_manual(values = strategy_colors) +
    theme_bw() +
    xlab("Date") + ylab("Experienced temperature (°C)") +
    labs(title = "Individual experienced temperature — optimal movers (full)")

  ref_year <- year(min(habitat_df$date)) + 1
  p_temp_year <- ibm_long |>
    filter(strategy == "optimal_mover", year(date) == ref_year) |>
    ggplot(aes(x = date, y = temp, color = strategy, group = pid)) +
    geom_line(alpha = 0.5) +
    scale_color_manual(values = strategy_colors) +
    theme_bw() +
    xlab("Date") + ylab("Experienced temperature (°C)") +
    labs(title = paste("Individual experienced temperature —", ref_year))


  # ── 3. GROWTH RATE ──────────────────────────────────────────────────────────
  p_ggd_full <- ibm_long |>
    ggplot(aes(x = date, y = ggd, color = strategy, group = pid)) +
    geom_line(alpha = 0.4) +
    scale_color_manual(values = strategy_colors) +
    theme_bw() +
    xlab("Date") + ylab("Specific growth rate (g/g/d)") +
    labs(title = "Individual growth rates")

  p_ggd_age <- ibm_long |>
    filter(year(date) == ref_year) |>
    mutate(age_f = as.factor(round(age))) |>
    ggplot(aes(x = date, y = ggd, color = age_f, group = age_f)) +
    geom_line(alpha = 0.6) +
    theme_bw() +
    xlab("Date") + ylab("Specific growth rate (g/g/d)") +
    labs(title = paste("Growth rate by age class —", ref_year), color = "Age")


  # ── 4. SIZE TRAJECTORIES ────────────────────────────────────────────────────
  p_wt_age <- ibm_long |>
    ggplot(aes(x = age, y = weight, color = strategy, group = pid)) +
    geom_line(alpha = 0.4) +
    scale_color_manual(values = strategy_colors) +
    scale_fill_manual(values = strategy_colors) +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Age (years)") + ylab("Individual weight (g)") +
    labs(title = "Weight trajectories over life span")

  p_wt_date <- ibm_long |>
    ggplot() +
    geom_rect(data = summer_shading,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
    vline +
    geom_line(aes(x = date, y = weight, color = strategy, group = pid), alpha = 0.4) +
    scale_color_manual(values = strategy_colors) +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Date") + ylab("Individual weight (g)") +
    labs(title = "Weight trajectories over time")


  # ── 5. MAXIMUM AGE ──────────────────────────────────────────────────────────
  max_age_df <- ibm_long |>
    filter(survived == 1) |>
    group_by(strategy, pid) |>
    summarise(max_age = max(age), .groups = "drop")

  p_max_age <- ggplot(max_age_df, aes(x = max_age, fill = strategy)) +
    geom_histogram(binwidth = 0.25, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = strategy_colors) +
    facet_wrap(~strategy, ncol = 1, scales = "free_y") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    xlab("Maximum age (years)") + ylab("Number of fish") +
    labs(title = "Maximum age by strategy")


  # ── 6. SURVIVAL ─────────────────────────────────────────────────────────────
  mycols <- RColorBrewer::brewer.pal(4, "Set2")
  surv_source_colors <- c(
    "Temperature (p_temp)" = mycols[1],
    "Starvation (p_starv)" = mycols[2],
    "Size (p_sg)"          = mycols[3]
  )

  all_surv_df <- ibm_long |>
    left_join(
      habitat_df |> select(dayofsim, S_max_warm, S_max_cold),
      by = "dayofsim"
    ) |>
    mutate(
      p_temp  = fncSurviveTemp(temp, T1 = T1_mort, T9 = T9_mort),
      p_starv = fncSurviveStarve(condition, K9 = K9_starv, K1 = K1_starv),
      p_sg    = fncSurviveSize(weight,
                               maxprob = ifelse(patch == "warm", S_max_warm, S_max_cold))[[1]]
    ) |>
    pivot_longer(c(p_temp, p_starv, p_sg),
                 names_to = "source", values_to = "probability") |>
    mutate(
      source = recode(source,
        p_temp  = "Temperature (p_temp)",
        p_starv = "Starvation (p_starv)",
        p_sg    = "Size (p_sg)"
      ),
      strategy = factor(strategy,
        levels = c("resident_warm", "resident_cold", "optimal_mover"))
    ) |>
    filter(strategy != "resident_warm")

  p_surv_movers <- all_surv_df |>
    filter(strategy == "optimal_mover") |>
    group_by(strategy, date, cohort, source) |>
    summarise(mean_p = mean(probability), .groups = "drop") |>
    ggplot(aes(x = date, y = mean_p, color = source)) +
    geom_rect(data = summer_shading,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = surv_source_colors) +
    facet_wrap(~cohort, scales = "free_y") +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Date") + ylab("Mean daily survival probability") +
    labs(title = "Mortality sources by cohort — optimal movers (group means)",
         color = "Source")

  p_condition <- ibm_long |>
    ggplot(aes(x = date, y = condition, color = strategy, fill = strategy)) +
    geom_line(alpha = 0.4) +
    geom_hline(yintercept = c(K9_starv, K1_starv), linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = strategy_colors) +
    scale_fill_manual(values = strategy_colors) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Date") + ylab("Relative condition (W / W_peak)") +
    labs(title = "Condition over time (dashed = starvation anchors)",
         color = "Strategy", fill = "Strategy") +
    facet_wrap(~cohort, scales = "free_x")


  # ── 7. MOVEMENT ─────────────────────────────────────────────────────────────
  p_occupancy <- ibm_long |>
    filter(strategy == "optimal_mover") |>
    group_by(cohort, date) |>
    summarise(prop_warm = mean(patch == "warm"), .groups = "drop") |>
    ggplot(aes(x = date, y = prop_warm)) +
    vline +
    geom_line() +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    facet_wrap(~cohort, scales = "free_x") +
    theme_bw() +
    xlab("Date") + ylab("Proportion of optimal movers in warm patch") +
    labs(title = "Patch occupancy by cohort — optimal movers")

  pid_order <- ibm_long |>
    group_by(pid) |>
    summarise(prop_warm = mean(patch == "warm"), .groups = "drop") |>
    arrange(prop_warm)

  p_raster <- ibm_long |>
    mutate(pid = factor(pid, levels = pid_order$pid)) |>
    ggplot(aes(x = date, y = pid, fill = patch)) +
    geom_tile() +
    scale_fill_manual(values = c(warm = "#d95f02", cold = "#1f78b4")) +
    facet_wrap(~cohort, scales = "free") +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank()) +
    xlab("Date") + ylab("Individual fish (ordered by % time in warm)") +
    labs(title = "Patch occupancy raster — optimal movers", fill = "Patch")


  # ── 8. REPRODUCTION ─────────────────────────────────────────────────────────
  spawn_days <- unique(spawn_log$dayofsim)

  spawner_traits <- spawn_log |>
    left_join(fish_registry |> select(pid, strategy), by = c("parent_pid" = "pid")) |>
    left_join(ibm_long |> select(pid, dayofsim, age),
              by = c("parent_pid" = "pid", "dayofsim")) |>
    transmute(pid = parent_pid, dayofsim, weight, condition, age, strategy,
              spawned = "Spawner")

  non_spawner_traits <- ibm_long |>
    filter(dayofsim %in% spawn_days, !is.na(weight), !is.na(condition)) |>
    anti_join(spawn_log, by = c("pid" = "parent_pid", "dayofsim")) |>
    transmute(pid, dayofsim, weight, condition, age, strategy, spawned = "Non-spawner")

  spawner_df <- bind_rows(spawner_traits, non_spawner_traits) |>
    filter(strategy != "resident_warm") |>
    mutate(strategy = factor(strategy,
      levels = c("optimal_mover", "resident_cold"),
      labels = c("Optimal mover", "Resident (cold)")))

  make_trait_plot <- function(data, x_var, x_label) {
    ggplot(data, aes(x = .data[[x_var]], fill = spawned, color = spawned)) +
      geom_density(aes(y = after_stat(scaled)), alpha = 0.35, linewidth = 0.6) +
      geom_rug(data = filter(data, spawned == "Spawner"),
               sides = "b", linewidth = 0.5, alpha = 0.8, color = "#d62728") +
      facet_wrap(~strategy, ncol = 1) +
      scale_fill_manual(values = c("Non-spawner" = "#7bafd4", "Spawner" = "#d62728")) +
      scale_color_manual(values = c("Non-spawner" = "#7bafd4", "Spawner" = "#d62728")) +
      labs(x = x_label, y = "Scaled density", fill = NULL, color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold"))
  }

  p_spawner_wt   <- make_trait_plot(spawner_df, "weight",    "Weight (g)")
  p_spawner_cond <- make_trait_plot(spawner_df, "condition", "Condition (K)")

  p_spawn_traits <- (p_spawner_wt + p_spawner_cond + patchwork::plot_layout(guides = "collect")) +
    patchwork::plot_annotation(
      title    = "Spawner vs. non-spawner traits on spawning days",
      subtitle = "Resident (warm) excluded"
    ) &
    theme(legend.position = "bottom")

  spawn_plot_df <- spawn_log |>
    left_join(fish_registry |> select(pid, strategy), by = c("parent_pid" = "pid")) |>
    left_join(ibm_long |> select(pid, dayofsim, date),
              by = c("parent_pid" = "pid", "dayofsim")) |>
    mutate(doy = yday(date), year = year(date)) |>
    filter(strategy != "resident_warm")

  p_spawn_timing <- ggplot(spawn_plot_df, aes(x = doy, fill = strategy)) +
    geom_histogram(bins = 52, alpha = 0.7) +
    scale_fill_manual(values = strategy_colors, name = NULL) +
    scale_x_continuous(
      breaks = c(1, 60, 121, 182, 244, 305),
      labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")
    ) +
    theme_bw() +
    xlab("Day of year") + ylab("Spawn events") +
    labs(title = "Spawn timing — day of year") +
    facet_wrap(~strategy, ncol = 1)


  # ── Print ────────────────────────────────────────────────────────────────────
  print(ggpubr::ggarrange(p_ration_cold, p_ration_warm,
                          nrow = 1, common.legend = TRUE, legend = "right"))
  print(ggpubr::ggarrange(p_pcmax_cold, p_pcmax_warm, nrow = 1))
  print(p_temp_full)
  print(p_temp_year)
  print(p_ggd_full)
  print(p_ggd_age)
  print(p_wt_age)
  print(p_wt_date)
  print(p_max_age)
  print(p_surv_movers)
  print(p_condition)
  print(p_occupancy)
  print(p_raster)
  print(p_spawn_traits)
  print(p_spawn_timing)

  invisible(list(
    ration_warm   = p_ration_warm,
    ration_cold   = p_ration_cold,
    pcmax_warm    = p_pcmax_warm,
    pcmax_cold    = p_pcmax_cold,
    temp_full     = p_temp_full,
    temp_year     = p_temp_year,
    ggd_full      = p_ggd_full,
    ggd_age       = p_ggd_age,
    wt_age        = p_wt_age,
    wt_date       = p_wt_date,
    max_age       = p_max_age,
    surv_movers   = p_surv_movers,
    condition     = p_condition,
    occupancy     = p_occupancy,
    raster        = p_raster,
    spawn_traits  = p_spawn_traits,
    spawn_timing  = p_spawn_timing
  ))
}


plot_population <- function(ibm_long, habitat_df, exp_start_date = NULL) {

  # ── Shared palette, labels, and derived objects ──────────────────────────────
  strategy_pal <- c(
    resident_warm = "tomato",
    resident_cold = "steelblue",
    optimal_mover = "darkgreen"
  )
  strategy_labs <- c(
    resident_warm = "Resident (warm)",
    resident_cold = "Resident (cold)",
    optimal_mover = "Optimal mover"
  )
  age_pal <- c(
    "Age-0" = "#d4e6f1", "Age-1" = "#5dade2", "Age-2" = "#1a5276",
    "Age-3" = "#a9dfbf", "Age-4" = "#1e8449", "Age-5" = "#784212",
    "Age-6" = "#f0b27a", "Age-7" = "#e67e22", "Age-8+" = "#922b21"
  )

  ibm_long  <- ibm_long |> mutate(age_class = floor(age))
  vline_yr  <- if (!is.null(exp_start_date))
    geom_vline(xintercept = year(exp_start_date),  linetype = "dashed") else NULL
  vline_dt  <- if (!is.null(exp_start_date))
    geom_vline(xintercept = exp_start_date,         linetype = "dashed") else NULL

  summer_shading <- data.frame(
    xmin = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-06-01")),
    xmax = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-08-31")),
    ymin = -Inf, ymax = Inf
  )

  # October 1 annual census
  census <- ibm_long |>
    filter(month(date) == 10, day(date) == 1, survived == 1) |>
    mutate(year = year(date))

  census_summary <- census |>
    group_by(year, strategy) |>
    summarise(
      n_total   = n(),
      biomass_g = sum(weight),
      mean_wt   = mean(weight),
      .groups   = "drop"
    ) |>
    mutate(strategy_lbl = recode(strategy, !!!strategy_labs))

  census_age <- census |>
    mutate(
      age_class_f  = factor(
        case_when(age_class >= 8L ~ "Age-8+",
                  TRUE            ~ paste0("Age-", age_class)),
        levels = c(paste0("Age-", 0:7), "Age-8+")
      ),
      strategy_lbl = recode(strategy, !!!strategy_labs)
    ) |>
    group_by(year, strategy, strategy_lbl, age_class_f) |>
    summarise(n = n(), .groups = "drop")

  lambda_df <- census_summary |>
    arrange(strategy, year) |>
    group_by(strategy, strategy_lbl) |>
    mutate(lambda = n_total / lag(n_total)) |>
    filter(!is.na(lambda))

  cohort_census <- census |>
    group_by(strategy, cohort, age_class) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(strategy, cohort, age_class)

  apparent_survival <- cohort_census |>
    group_by(strategy, cohort) |>
    mutate(surv_rate = lead(n) / n) |>
    filter(!is.na(surv_rate), age_class <= 8) |>
    mutate(strategy_lbl = recode(strategy, !!!strategy_labs))

  cohort_max_age <- ibm_long |>
    group_by(strategy, cohort, pid) |>
    summarise(max_age = max(age), .groups = "drop") |>
    group_by(strategy, cohort) |>
    mutate(cohort_n0 = n()) |>
    ungroup()

  lx_df <- cohort_max_age |>
    crossing(age_class = 0:12) |>
    group_by(strategy, cohort, age_class) |>
    summarise(lx = sum(max_age > age_class) / first(cohort_n0), .groups = "drop") |>
    filter(lx > 0) |>
    mutate(strategy_lbl = recode(strategy, !!!strategy_labs))

  ibm_summary_pop <- ibm_long |>
    filter(survived == 1) |>
    group_by(strategy, dayofsim, date) |>
    summarise(
      n_alive      = n(),
      mean_weight  = mean(weight, na.rm = TRUE),
      prop_warm    = mean(patch == "warm", na.rm = TRUE),
      .groups      = "drop"
    ) |>
    mutate(strategy_lbl = recode(strategy, !!!strategy_labs))


  # ── 1. DAILY TIME SERIES ─────────────────────────────────────────────────────
  p_daily_abund <- ibm_summary_pop |>
    ggplot(aes(x = date, y = n_alive, color = strategy_lbl)) +
    geom_rect(data = summer_shading,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
    vline_dt +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = NULL, y = "Number alive", title = "Daily abundance by strategy") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "top")

  p_daily_biomass <- ibm_summary_pop |>
    mutate(biomass_kg = (n_alive * mean_weight) / 1000) |>
    ggplot(aes(x = date, y = biomass_kg, color = strategy_lbl)) +
    geom_rect(data = summer_shading,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
    vline_dt +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Date", y = "Biomass (kg)", title = "Daily biomass by strategy") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")

  p_daily_census <- p_daily_abund / p_daily_biomass

  # ── 2. ANNUAL CENSUS ─────────────────────────────────────────────────────────
  p_ann_abund <- ggplot(
      census_summary |> filter(strategy != "resident_warm"),
      aes(x = year, y = n_total, color = strategy_lbl, group = strategy_lbl)
    ) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    vline_yr +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    labs(x = NULL, y = "Count", title = "Annual abundance (October 1 census)") +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "top")

  p_ann_bio <- ggplot(
      census_summary |> filter(strategy != "resident_warm"),
      aes(x = year, y = biomass_g / 1000, color = strategy_lbl, group = strategy_lbl)
    ) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    vline_yr +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    labs(x = "Year", y = "Biomass (kg)", title = "Annual biomass (October 1 census)") +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

  p_ann_census <- p_ann_abund / p_ann_bio

  p_ann_age_struct <- census_age |>
    filter(strategy != "resident_warm") |>
    ggplot(aes(x = year, y = n, fill = age_class_f)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = age_pal, name = "Age class") +
    vline_yr +
    facet_wrap(~strategy_lbl, ncol = 1, scales = "free_y") +
    labs(x = "Year", y = "Count",
         title = "Annual age-structured abundance (October 1)") +
    theme_bw() +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))


  # ── 3. WEIGHT AT AGE ─────────────────────────────────────────────────────────
  waa <- ibm_long |>
    filter(day(date) == 1) |>
    group_by(strategy, cohort, age_class) |>
    summarise(
      mean_wt = mean(weight),
      sd_wt   = sd(weight, na.rm = TRUE),
      age_mid = mean(age),
      .groups = "drop"
    ) |>
    mutate(strategy_lbl = recode(strategy, !!!strategy_labs),
           cohort_f     = factor(cohort))

  p_waa <- ggplot(waa |> filter(strategy != "resident_warm"),
                  aes(x = age_mid, y = mean_wt,
                      group = cohort_f, color = cohort_f)) +
    geom_ribbon(aes(ymin = mean_wt - sd_wt, ymax = mean_wt + sd_wt,
                    fill = cohort_f),
                alpha = 0.08, color = NA) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    scale_color_viridis_d(option = "turbo", name = "Cohort") +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    facet_wrap(~strategy_lbl) +
    labs(x = "Age (years)", y = "Mean weight (g)",
         title = "Weight-at-age by cohort",
         subtitle = "Ribbons = ±1 SD; each line = one cohort") +
    theme_bw() + theme(panel.grid = element_blank())

  waa_annual <- census |>
    filter(strategy != "resident_warm") |>
    mutate(strategy_lbl = recode(strategy, !!!strategy_labs)) |>
    group_by(year, strategy, strategy_lbl, age_class) |>
    summarise(
      mean_wt = mean(weight), sd_wt = sd(weight), n = n(), .groups = "drop"
    ) |>
    filter(n >= 3, age_class <= 7) |>
    left_join(census_summary |> select(year, strategy, n_total),
              by = c("year", "strategy"))

  p_waa_annual <- waa_annual |>
    mutate(age_class_f = paste0("Age-", age_class)) |>
    ggplot(aes(x = year, y = mean_wt,
               color = strategy_lbl, group = strategy_lbl)) +
    geom_ribbon(aes(ymin = mean_wt - sd_wt, ymax = mean_wt + sd_wt,
                    fill = strategy_lbl),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(aes(size = n), alpha = 0.8) +
    vline_yr +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    scale_fill_manual(values  = setNames(strategy_pal, strategy_labs), guide = "none") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    facet_wrap(~age_class_f, scales = "free_y", ncol = 4) +
    labs(x = "Census year", y = "Mean weight (g)",
         title = "Mean weight at age over time (October 1 census)",
         subtitle = "Ribbon = ±1 SD  |  point size = sample size") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"), legend.position = "top")

  p_density_growth <- waa_annual |>
    filter(strategy == "optimal_mover", age_class <= 4) |>
    filter(year > min(year) + 9) |>
    mutate(age_class_f = paste0("Age-", age_class)) |>
    ggplot(aes(x = n_total, y = mean_wt)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linewidth = 0.7, alpha = 0.2) +
    geom_point(aes(color = year), size = 3, alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", name = "Year") +
    facet_wrap(~age_class_f, scales = "free", ncol = 2) +
    labs(x = "Total population size (October 1)",
         y = "Mean weight (g)",
         title = "Density–growth relationship by age class (optimal mover)",
         subtitle = "Point colour = year  |  line = OLS fit") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))


  # ── 4. ANNUAL SURVIVAL ───────────────────────────────────────────────────────
  # Each cohort contributes exactly one survival rate per age class per year
  # (year = cohort + age_class). Facet by age class, colored by strategy.
  surv_time <- apparent_survival |>
    filter(strategy != "resident_warm", age_class <= 6) |>
    mutate(
      year         = cohort + age_class,
      age_class_f  = factor(paste0("Age-", age_class),
                            levels = paste0("Age-", 0:6))
    )

  p_ann_surv <- ggplot(surv_time,
                        aes(x = year, y = surv_rate,
                            color = strategy_lbl, group = strategy_lbl)) +
    vline_yr +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.9) +
    facet_wrap(~age_class_f, ncol = 4, scales = "free_y") +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Census year", y = "Apparent annual survival",
         title = "Apparent survival by age class over time (October 1 census)",
         subtitle = "Survival from age X to X+1") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"), legend.position = "top")


  # ── 6. LAMBDA ────────────────────────────────────────────────────────────────
  p_lambda <- ggplot(
      lambda_df |> filter(strategy != "resident_warm"),
      aes(x = year, y = lambda, color = strategy_lbl, group = strategy_lbl)
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.9) +
    vline_yr +
    scale_color_manual(values = setNames(strategy_pal, strategy_labs), name = NULL) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
    labs(x = "Year", y = "\u03bb (log scale)",
         title = "Annual finite rate of increase (October 1 census)",
         subtitle = "Dashed line = \u03bb = 1 (stable population)") +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "top")


  # ── 7. COHORT SURVIVORSHIP ───────────────────────────────────────────────────
  p_lx <- lx_df |>
    filter(strategy != "resident_warm") |>
    mutate(cohort_f = factor(cohort)) |>
    ggplot(aes(x = age_class, y = lx, group = cohort_f, color = cohort_f)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    scale_color_viridis_d(option = "turbo", name = "Cohort") +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_x_continuous(breaks = 0:12) +
    facet_wrap(~strategy_lbl) +
    labs(x = "Age class", y = "Fraction of cohort surviving",
         title = "Cohort survivorship curves (lx)",
         subtitle = "Each line = one birth cohort") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))


  # ── Print ────────────────────────────────────────────────────────────────────
  print(p_daily_census)
  print(p_ann_census)
  print(p_ann_age_struct)
  print(p_waa)
  print(p_waa_annual)
  print(p_density_growth)
  print(p_ann_surv)
  print(p_lambda)
  print(p_lx)

  invisible(list(
    daily_census   = p_daily_census,
    ann_census     = p_ann_census,
    ann_age_struct = p_ann_age_struct,
    waa            = p_waa,
    waa_annual     = p_waa_annual,
    density_growth = p_density_growth,
    ann_surv       = p_ann_surv,
    lambda         = p_lambda,
    lx             = p_lx
  ))
}


plot_patchvalue <- function(ibm_long, habitat_df, exp_start_date = NULL) {

  # ── Shared setup ─────────────────────────────────────────────────────────────
  patch_pal <- c(warm = "tomato", cold = "steelblue")

  vline_yr <- if (!is.null(exp_start_date))
    geom_vline(xintercept = year(exp_start_date), linetype = "dashed", color = "grey50")
  else NULL

  hab_value <- ibm_long |>
    filter(!is.na(ggd)) |>
    mutate(cal_year = year(date), daily_growth_g = weight * ggd)

  annual_area <- habitat_df |>
    group_by(cal_year = year(date)) |>
    summarise(A_warm = mean(A_warm), A_cold = mean(A_cold), .groups = "drop")

  # Helper: add season label + season_year (meteorological winter convention)
  add_season <- function(df) {
    df |>
      mutate(
        season = factor(
          case_when(
            month(date) %in% 3:5  ~ "Spring (Mar–May)",
            month(date) %in% 6:8  ~ "Summer (Jun–Aug)",
            month(date) %in% 9:11 ~ "Autumn (Sep–Nov)",
            TRUE                  ~ "Winter (Dec–Feb)"
          ),
          levels = c("Spring (Mar–May)", "Summer (Jun–Aug)",
                     "Autumn (Sep–Nov)", "Winter (Dec–Feb)")
        ),
        season_year = if_else(month(date) %in% 1:2, cal_year - 1L, cal_year)
      )
  }

  seasonal_area <- habitat_df |>
    mutate(cal_year = year(date)) |>
    add_season() |>
    group_by(season_year, season) |>
    summarise(A_warm = mean(A_warm), A_cold = mean(A_cold), .groups = "drop")

  # Shared equivalence-curve annotation layers (above/below zero shading)
  equiv_bands <- list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin =  0, ymax =  Inf,
             fill = "steelblue", alpha = 0.05),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = "tomato",    alpha = 0.05),
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
  )

  # ── Single-patch guard ───────────────────────────────────────────────────────
  # When warm area is 0 throughout (e.g. null_cold scenario), there is no second
  # patch to compare against. Equivalence curves collapse to a single degenerate
  # point and per-area warm biomass produces NaN (0/0). Return early rather than
  # produce misleading output.
  if (all(habitat_df$A_warm == 0)) {
    message(
      "plot_patchvalue(): warm habitat area is 0 for all time steps. ",
      "Habitat-value and equivalence plots require two patches — skipping."
    )
    return(invisible(NULL))
  }

  om <- hab_value |> filter(strategy == "optimal_mover")


  # ── 1. HABITAT SELECTION ─────────────────────────────────────────────────────
  # Selection ratio = prop time in warm / prop area that is warm.
  # Values > 1 indicate warm preference; < 1 indicate cold preference; 1 = no preference.
  # Computed at the individual level then averaged, so the SD reflects among-individual
  # variation in selection, not variation in area.
  sel_ann <- om |>
    group_by(pid, cal_year) |>
    summarise(prop_warm = mean(patch == "warm"), .groups = "drop") |>
    left_join(annual_area |> mutate(prop_warm_area = A_warm / (A_warm + A_cold)),
              by = "cal_year") |>
    mutate(sel_ratio = prop_warm / prop_warm_area) |>
    group_by(cal_year) |>
    summarise(mean_sel = mean(sel_ratio), sd_sel = sd(sel_ratio), .groups = "drop")

  p_sel_ann <- ggplot(sel_ann, aes(x = cal_year, y = mean_sel)) +
    vline_yr +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = mean_sel - sd_sel, ymax = mean_sel + sd_sel),
                alpha = 0.2, fill = "tomato") +
    geom_line(linewidth = 0.8, color = "tomato") +
    geom_point(size = 2, color = "tomato") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
    labs(x = "Year", y = "Warm habitat selection ratio (use / availability)",
         title   = "Annual warm habitat selection (optimal movers)",
         subtitle = "Dashed line = proportional use (no preference)",
         caption = "Band = ±1 SD") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  sel_seas <- om |>
    add_season() |>
    group_by(pid, season_year, season) |>
    summarise(prop_warm = mean(patch == "warm"), .groups = "drop") |>
    left_join(seasonal_area |> mutate(prop_warm_area = A_warm / (A_warm + A_cold)),
              by = c("season_year", "season")) |>
    mutate(sel_ratio = prop_warm / prop_warm_area) |>
    group_by(season_year, season) |>
    summarise(mean_sel = mean(sel_ratio), sd_sel = sd(sel_ratio), .groups = "drop")

  p_sel_seas <- ggplot(sel_seas, aes(x = season_year, y = mean_sel)) +
    vline_yr +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = mean_sel - sd_sel, ymax = mean_sel + sd_sel),
                alpha = 0.2, fill = "tomato") +
    geom_line(linewidth = 0.8, color = "tomato") +
    geom_point(size = 2, color = "tomato") +
    facet_wrap(~season, ncol = 2) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Year", y = "Warm habitat selection ratio (use / availability)",
         title   = "Seasonal warm habitat selection (optimal movers)",
         subtitle = "Dashed line = proportional use (no preference)",
         caption = "Band = ±1 SD") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))


  # ── 2. TOTAL BIOMASS ─────────────────────────────────────────────────────────
  biomass_ann <- om |>
    group_by(patch, cal_year) |>
    summarise(total_kg = sum(daily_growth_g, na.rm = TRUE) / 1000, .groups = "drop")

  p_biomass_ann <- ggplot(biomass_ann, aes(x = cal_year, y = total_kg, color = patch)) +
    vline_yr +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
    labs(x = "Year", y = "Total biomass accrued (kg)",
         title = "Annual biomass accrued by habitat patch (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  biomass_seas <- om |>
    add_season() |>
    group_by(patch, season_year, season) |>
    summarise(total_kg = sum(daily_growth_g, na.rm = TRUE) / 1000, .groups = "drop")

  p_biomass_seas <- ggplot(biomass_seas,
                            aes(x = season_year, y = total_kg, color = patch)) +
    vline_yr +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~season, ncol = 2, scales = "free_y") +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Year", y = "Total biomass accrued (kg)",
         title = "Seasonal biomass accrued by habitat patch (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))

  biomass_equiv <- biomass_ann |>
    pivot_wider(names_from = patch, values_from = total_kg, values_fill = 0) |>
    mutate(biomass_diff = cold - warm) |>
    left_join(annual_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = "cal_year")

  p_biomass_equiv <- ggplot(biomass_equiv, aes(x = prop_cold, y = biomass_diff)) +
    equiv_bands +
    geom_path(aes(color = cal_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = cal_year), size = 2.5) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Biomass accrued (cold) \u2212 biomass accrued (warm)  (kg)",
         title = "Total biomass equivalence curve (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  biomass_equiv_seas <- biomass_seas |>
    pivot_wider(names_from = patch, values_from = total_kg, values_fill = 0) |>
    mutate(biomass_diff = cold - warm) |>
    left_join(seasonal_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = c("season_year", "season"))

  p_biomass_equiv_seas <- ggplot(biomass_equiv_seas,
                                  aes(x = prop_cold, y = biomass_diff)) +
    equiv_bands +
    geom_path(aes(color = season_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = season_year), size = 2) +
    facet_wrap(~season, ncol = 2) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Biomass accrued (cold) \u2212 biomass accrued (warm)  (kg)",
         title = "Total biomass equivalence curve by season (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))


  # ── 3. PER-AREA BIOMASS ──────────────────────────────────────────────────────
  per_area_ann <- biomass_ann |>
    pivot_wider(names_from = patch, values_from = total_kg) |>
    left_join(annual_area, by = "cal_year") |>
    mutate(warm = warm / A_warm, cold = cold / A_cold) |>
    pivot_longer(c(warm, cold), names_to = "patch", values_to = "kg_per_area")

  p_per_area_ann <- ggplot(per_area_ann,
                            aes(x = cal_year, y = kg_per_area, color = patch)) +
    vline_yr +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
    labs(x = "Year", y = "Biomass accrued per unit area (kg)",
         title = "Annual per-area biomass: warm vs. cold habitat (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  per_area_seas <- biomass_seas |>
    pivot_wider(names_from = patch, values_from = total_kg, values_fill = 0) |>
    left_join(seasonal_area, by = c("season_year", "season")) |>
    mutate(warm = warm / A_warm, cold = cold / A_cold) |>
    pivot_longer(c(warm, cold), names_to = "patch", values_to = "kg_per_area")

  p_per_area_seas <- ggplot(per_area_seas,
                             aes(x = season_year, y = kg_per_area, color = patch)) +
    vline_yr +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~season, ncol = 2, scales = "free_y") +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Year", y = "Biomass accrued per unit area (kg)",
         title = "Seasonal per-area biomass: warm vs. cold habitat (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))

  per_area_equiv <- biomass_ann |>
    pivot_wider(names_from = patch, values_from = total_kg) |>
    left_join(annual_area, by = "cal_year") |>
    mutate(diff_per_area = cold / A_cold - warm / A_warm,
           prop_cold     = A_cold / (A_cold + A_warm))

  p_per_area_equiv <- ggplot(per_area_equiv, aes(x = prop_cold, y = diff_per_area)) +
    equiv_bands +
    geom_path(aes(color = cal_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = cal_year), size = 2.5) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Biomass/area (cold) \u2212 biomass/area (warm)  (kg)",
         title = "Per-area biomass equivalence curve (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  per_area_equiv_seas <- biomass_seas |>
    pivot_wider(names_from = patch, values_from = total_kg, values_fill = 0) |>
    left_join(seasonal_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = c("season_year", "season")) |>
    mutate(diff_per_area = cold / A_cold - warm / A_warm)

  p_per_area_equiv_seas <- ggplot(per_area_equiv_seas,
                                   aes(x = prop_cold, y = diff_per_area)) +
    equiv_bands +
    geom_path(aes(color = season_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = season_year), size = 2) +
    facet_wrap(~season, ncol = 2) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Biomass/area (cold) \u2212 biomass/area (warm)  (kg)",
         title = "Per-area biomass equivalence curve by season (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))


  # ── 4. GROWTH RATE ───────────────────────────────────────────────────────────
  ggd_ann <- om |>
    group_by(patch, cal_year) |>
    summarise(mean_ggd = mean(ggd, na.rm = TRUE),
              sd_ggd   = sd(ggd,   na.rm = TRUE), .groups = "drop")

  p_ggd_ann <- ggplot(ggd_ann, aes(x = cal_year, y = mean_ggd,
                                    color = patch, fill = patch)) +
    vline_yr +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_ribbon(aes(ymin = mean_ggd - sd_ggd, ymax = mean_ggd + sd_ggd),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_fill_manual(values  = patch_pal, name = "Patch") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
    labs(x = "Year", y = "Mean specific growth rate (g/g/d)",
         title   = "Annual realized growth rate by habitat patch (optimal movers)",
         caption = "Band = ±1 SD") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  ggd_seas <- om |>
    add_season() |>
    group_by(patch, season_year, season) |>
    summarise(mean_ggd = mean(ggd, na.rm = TRUE),
              sd_ggd   = sd(ggd,   na.rm = TRUE), .groups = "drop")

  p_ggd_seas <- ggplot(ggd_seas,
                        aes(x = season_year, y = mean_ggd,
                            color = patch, fill = patch)) +
    vline_yr +
    geom_ribbon(aes(ymin = mean_ggd - sd_ggd, ymax = mean_ggd + sd_ggd),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~season, ncol = 2, scales = "free_y") +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_fill_manual(values  = patch_pal, name = "Patch") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Year", y = "Mean specific growth rate (g/g/d)",
         title   = "Seasonal realized growth rate by habitat patch (optimal movers)",
         caption = "Band = ±1 SD") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))

  ggd_equiv <- ggd_ann |>
    select(-sd_ggd) |>
    pivot_wider(names_from = patch, values_from = mean_ggd) |>
    mutate(ggd_diff = cold - warm) |>
    left_join(annual_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = "cal_year")

  p_ggd_equiv <- ggplot(ggd_equiv, aes(x = prop_cold, y = ggd_diff)) +
    equiv_bands +
    geom_path(aes(color = cal_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = cal_year), size = 2.5) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Mean ggd (cold) \u2212 mean ggd (warm)  (g/g/d)",
         title = "Growth equivalence curve (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  ggd_equiv_seas <- ggd_seas |>
    select(-sd_ggd) |>
    pivot_wider(names_from = patch, values_from = mean_ggd) |>
    mutate(ggd_diff = cold - warm) |>
    left_join(seasonal_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = c("season_year", "season"))

  p_ggd_equiv_seas <- ggplot(ggd_equiv_seas, aes(x = prop_cold, y = ggd_diff)) +
    equiv_bands +
    geom_path(aes(color = season_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = season_year), size = 2) +
    facet_wrap(~season, ncol = 2) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Mean ggd (cold) \u2212 mean ggd (warm)  (g/g/d)",
         title = "Growth equivalence curve by season (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))


  # ── 5. SURVIVAL ──────────────────────────────────────────────────────────────
  surv_ann <- om |>
    group_by(patch, cal_year) |>
    summarise(mean_psurv = mean(p_survive, na.rm = TRUE),
              sd_psurv   = sd(p_survive,   na.rm = TRUE), .groups = "drop")

  p_surv_ann <- ggplot(surv_ann, aes(x = cal_year, y = mean_psurv,
                                      color = patch, fill = patch)) +
    vline_yr +
    geom_ribbon(aes(ymin = mean_psurv - sd_psurv, ymax = mean_psurv + sd_psurv),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_fill_manual(values  = patch_pal, name = "Patch") +
    scale_y_continuous(labels = scales::label_percent(accuracy = 0.01)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
    labs(x = "Year", y = "Mean daily survival probability",
         title   = "Annual mean daily survival by habitat patch (optimal movers)",
         caption = "Band = ±1 SD") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  surv_seas <- om |>
    add_season() |>
    group_by(patch, season_year, season) |>
    summarise(mean_psurv = mean(p_survive, na.rm = TRUE),
              sd_psurv   = sd(p_survive,   na.rm = TRUE), .groups = "drop")

  p_surv_seas <- ggplot(surv_seas,
                         aes(x = season_year, y = mean_psurv,
                             color = patch, fill = patch)) +
    vline_yr +
    geom_ribbon(aes(ymin = mean_psurv - sd_psurv, ymax = mean_psurv + sd_psurv),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
    facet_wrap(~season, ncol = 2, scales = "free_y") +
    scale_color_manual(values = patch_pal, name = "Patch") +
    scale_fill_manual(values  = patch_pal, name = "Patch") +
    scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Year", y = "Mean daily survival probability",
         title   = "Seasonal mean daily survival by habitat patch (optimal movers)",
         caption = "Band = ±1 SD") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))

  surv_equiv <- surv_ann |>
    select(-sd_psurv) |>
    pivot_wider(names_from = patch, values_from = mean_psurv) |>
    mutate(psurv_diff = cold - warm) |>
    left_join(annual_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = "cal_year")

  p_surv_equiv <- ggplot(surv_equiv, aes(x = prop_cold, y = psurv_diff)) +
    equiv_bands +
    geom_path(aes(color = cal_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = cal_year), size = 2.5) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Mean p_survive (cold) \u2212 mean p_survive (warm)",
         title = "Survival equivalence curve (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  surv_equiv_seas <- surv_seas |>
    select(-sd_psurv) |>
    pivot_wider(names_from = patch, values_from = mean_psurv) |>
    mutate(psurv_diff = cold - warm) |>
    left_join(seasonal_area |> mutate(prop_cold = A_cold / (A_cold + A_warm)),
              by = c("season_year", "season"))

  p_surv_equiv_seas <- ggplot(surv_equiv_seas, aes(x = prop_cold, y = psurv_diff)) +
    equiv_bands +
    geom_path(aes(color = season_year), linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = season_year), size = 2) +
    facet_wrap(~season, ncol = 2) +
    scale_color_viridis_c(name = "Year", option = "plasma") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    labs(x = "Cold habitat as proportion of total",
         y = "Mean p_survive (cold) \u2212 mean p_survive (warm)",
         title = "Survival equivalence curve by season (optimal movers)") +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       strip.text = element_text(face = "bold"))


  # ── Combine annual and seasonal panels side by side ──────────────────────────
  p_sel          <- p_sel_ann | p_sel_seas
  p_biomass      <- (p_biomass_ann      | p_biomass_seas)      + plot_layout(guides = "collect")
  p_biomass_equiv  <- (p_biomass_equiv  | p_biomass_equiv_seas)  + plot_layout(guides = "collect")
  p_per_area     <- (p_per_area_ann     | p_per_area_seas)     + plot_layout(guides = "collect")
  p_per_area_equiv <- (p_per_area_equiv | p_per_area_equiv_seas) + plot_layout(guides = "collect")
  p_ggd          <- (p_ggd_ann          | p_ggd_seas)          + plot_layout(guides = "collect")
  p_ggd_equiv      <- (p_ggd_equiv      | p_ggd_equiv_seas)      + plot_layout(guides = "collect")
  p_surv         <- (p_surv_ann         | p_surv_seas)         + plot_layout(guides = "collect")
  p_surv_equiv     <- (p_surv_equiv     | p_surv_equiv_seas)     + plot_layout(guides = "collect")

  # ── Print ────────────────────────────────────────────────────────────────────
  print(p_sel)
  print(p_biomass);      print(p_biomass_equiv)
  print(p_per_area);     print(p_per_area_equiv)
  print(p_ggd);          print(p_ggd_equiv)
  print(p_surv);         print(p_surv_equiv)

  invisible(list(
    sel            = p_sel,
    biomass        = p_biomass,
    biomass_equiv  = p_biomass_equiv,
    per_area       = p_per_area,
    per_area_equiv = p_per_area_equiv,
    ggd            = p_ggd,
    ggd_equiv      = p_ggd_equiv,
    surv           = p_surv,
    surv_equiv     = p_surv_equiv
  ))
}


