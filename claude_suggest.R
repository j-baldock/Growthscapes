# =============================================================================
# Two-Patch Individual-Based Model: Migratory Strategy & Population Dynamics
# =============================================================================
# Framework:
#   Patch 1 (warm/productive): high temperature, high food availability
#   Patch 2 (cold/unproductive): low temperature, low food availability
#
# Three migratory strategies are compared:
#   "resident_warm"  – remains in warm patch year-round
#   "resident_cold"  – remains in cold patch year-round
#   "migrant"        – moves to warm patch in summer, cold patch in winter
#
# Growth follows a simplified bioenergetics model:
#   delta_W = consumption(T, F, W) - respiration(T, W)
# Fecundity scales with body mass. Survival is size- and condition-dependent.
# =============================================================================

library(tidyverse)
library(ggplot2)

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Global parameters
# -----------------------------------------------------------------------------

params <- list(
  
  # Simulation
  n_individuals  = 200,        # individuals per strategy
  n_timesteps    = 52,        # weekly time steps (~10 years)
  dt             = 1,          # time step (weeks)
  
  # Patch environment --------------------------------------------------------
  # Warm patch: seasonal temperature oscillation around a warm mean
  T_warm_mean    = 18,         # deg C, mean annual temperature
  T_warm_amp     = 6,          # seasonal amplitude
  F_warm_mean    = 0.8,        # food availability [0,1], mean
  F_warm_sd      = 0.05,       # stochastic noise on food
  
  # Cold patch
  T_cold_mean    = 8,
  T_cold_amp     = 4,
  F_cold_mean    = 0.3,
  F_cold_sd      = 0.05,
  
  # Bioenergetics ------------------------------------------------------------
  # Consumption: C = Cmax(T) * F * W^alpha
  Cmax_ref       = 0.06,       # max mass-specific consumption at T_opt (per week)
  T_opt          = 16,         # optimum temperature for consumption (deg C)
  T_sigma        = 6,          # thermal breadth (deg C)
  alpha          = 0.80,       # allometric exponent for consumption
  # Respiration: R = r0 * exp(r1 * T) * W^beta
  r0             = 0.012,      # baseline respiration coefficient
  r1             = 0.07,       # temperature sensitivity of respiration
  beta           = 0.75,       # allometric exponent for respiration

  # Movement cost ------------------------------------------------------------
  # Deducted from expected net growth in the alternative patch before
  # the patch comparison is made. The fish only moves if:
  #   net_growth(alt_patch) - move_cost*W  >  net_growth(current_patch)
  move_cost      = 0.02,       # proportion of body mass lost per move event
  
  # Life history -------------------------------------------------------------
  W_init_mean    = 10,         # initial body mass (g), mean
  W_init_sd      = 2,          # variation in initial mass
  W_maturity     = 80,         # mass threshold for maturation (g)
  fecundity_a    = 0.8,        # fecundity = fecundity_a * W^fecundity_b
  fecundity_b    = 1.2,
  # Survival: logistic function of condition factor (K = W / W_expected(age))
  surv_base      = 0.995,      # baseline weekly survival
  surv_size_coef = 0.002,      # additional survival benefit per gram above threshold
  W_surv_thresh  = 30,         # mass threshold below which size penalty applies
  # Density dependence: applied to recruitment (Beverton-Holt)
  BH_alpha       = 0.8,        # max per-capita recruitment rate
  BH_beta        = 0.002,      # density-dependence scaling
  
  # Migration phenology (week of year) ---------------------------------------
  # Migrants move to warm patch at week 14 (spring) and back at week 40 (autumn)
  migrate_to_warm  = 14,
  migrate_to_cold  = 40
)

# -----------------------------------------------------------------------------
# 2. Environment functions
# -----------------------------------------------------------------------------

# Patch temperature at time t (week)
patch_temperature <- function(t, patch, p) {
  week_of_year <- (t - 1) %% 52 + 1
  phase <- 2 * pi * (week_of_year - 26) / 52  # peak in summer (~week 26)
  if (patch == "warm") {
    p$T_warm_mean + p$T_warm_amp * sin(phase)
  } else {
    p$T_cold_mean + p$T_cold_amp * sin(phase)
  }
}

# plot warm/cold patch temperatures through time
plot(patch_temperature(c(1:52), "warm", params) ~ c(1:52), type = "l", col = "red", xlab = "week", ylab = "temperature", ylim = c(0,24))
lines(patch_temperature(c(1:52), "cold", params) ~ c(1:52), type = "l", col = "blue")


# Patch food availability at time t (stochastic)
patch_food <- function(patch, p) {
  if (patch == "warm") {
    max(0, min(1, rnorm(1, p$F_warm_mean, p$F_warm_sd)))
  } else {
    max(0, min(1, rnorm(1, p$F_cold_mean, p$F_cold_sd)))
  }
}

# plot warm/cold patch food availability through time
plot(replicate(n = 52, patch_food("warm", params)) ~ c(1:52), type = "l", col = "red", xlab = "week", ylab = "temperature", ylim = c(0,1))
lines(replicate(n = 52, patch_food("cold", params)) ~ c(1:52), type = "l", col = "blue")


# -----------------------------------------------------------------------------
# 3. Bioenergetics functions
# -----------------------------------------------------------------------------

# Gaussian thermal performance curve for consumption
thermal_scalar <- function(T, T_opt, T_sigma) {
  exp(-0.5 * ((T - T_opt) / T_sigma)^2)
}
plot(thermal_scalar(seq(from = 0, to = 30, by = 1), params$T_opt, params$T_sigma) ~ seq(from = 0, to = 30, by = 1), type = "l", xlab = "temperature", ylab = "thermal performance (consumption)")


# Net weekly growth in grams
net_growth <- function(W, T, F_avail, p) {
  Cmax  <- p$Cmax_ref * thermal_scalar(T, p$T_opt, p$T_sigma)
  C     <- Cmax * F_avail * (W ^ p$alpha)        # gross consumption (g/week)
  R     <- p$r0 * exp(p$r1 * T) * (W ^ p$beta)   # respiration (g/week)
  C - R                                            # net somatic growth
}



# -----------------------------------------------------------------------------
# 4. Energetics-based patch choice
# -----------------------------------------------------------------------------
# For "optimal_mover" fish: evaluate expected net growth in both patches,
# penalise the *alternative* patch by the movement cost, then choose the
# patch with higher post-cost expected growth.
#
# Arguments:
#   current_patch  – the patch the fish currently occupies ("warm"/"cold")
#   W              – current body mass of the individual (g)
#   T_warm, T_cold – this time-step temperatures for each patch
#   F_warm, F_cold – this time-step food levels for each patch
#   p              – parameter list
#
# Returns:
#   list(patch = chosen patch, moved = logical)

choose_patch <- function(current_patch, W, T_warm, T_cold, F_warm, F_cold, p) {
  
  g_warm <- net_growth(W, T_warm, F_warm, p)
  g_cold <- net_growth(W, T_cold, F_cold, p)
  
  # Apply movement cost to whichever patch is *not* the current one
  if (current_patch == "warm") {
    g_cold_net <- g_cold - p$move_cost * W   # cost of leaving warm
    g_warm_net <- g_warm                      # staying: no cost
    if (g_cold_net > g_warm_net) {
      return(list(patch = "cold", moved = TRUE))
    } else {
      return(list(patch = "warm", moved = FALSE))
    }
  } else {
    g_warm_net <- g_warm - p$move_cost * W   # cost of leaving cold
    g_cold_net <- g_cold
    if (g_warm_net > g_cold_net) {
      return(list(patch = "warm", moved = TRUE))
    } else {
      return(list(patch = "cold", moved = FALSE))
    }
  }
}



# -----------------------------------------------------------------------------
# 5. IBM: simulate one cohort under a given strategy
# -----------------------------------------------------------------------------
# The loop is now individual-level for "optimal_mover" fish because each
# fish may independently be in a different patch at any time step.
# Residents are vectorised as before (all fish share the same patch).

simulate_strategy <- function(strategy, p) {
  
  n  <- p$n_individuals
  dt <- p$dt
  
  # Initialize individuals
  fish <- data.frame(
    id      = 1:n,
    W       = pmax(1, rnorm(n, p$W_init_mean, p$W_init_sd)),
    age     = 0,
    alive   = TRUE,
    mature  = FALSE,
    patch   = ifelse(strategy == "resident_cold", "cold", "warm")
  )
  
  results <- vector("list", p$n_timesteps)
  
  for (t in seq_len(p$n_timesteps)) {
    
    alive_idx    <- which(fish$alive)
    if (length(alive_idx) == 0) break
    week_of_year <- (t - 1) %% 52 + 1
    
    # --- Sample patch environments for this time step -----------------------
    T_warm <- patch_temperature(t, "warm", p)
    T_cold <- patch_temperature(t, "cold", p)
    F_warm <- patch_food("warm", p)
    F_cold <- patch_food("cold", p)
    
    # --- Patch assignment & growth update -----------------------------------
    
    if (strategy == "resident_warm") {
      # All fish stay in warm patch; vectorised update
      fish$patch[alive_idx] <- "warm"
      dW <- net_growth(fish$W[alive_idx], T_warm, F_warm, p)
      fish$W[alive_idx] <- pmax(0.1, fish$W[alive_idx] + dW * dt)
      n_moved   <- 0
      prop_warm <- 1
      
    } else if (strategy == "resident_cold") {
      fish$patch[alive_idx] <- "cold"
      dW <- net_growth(fish$W[alive_idx], T_cold, F_cold, p)
      fish$W[alive_idx] <- pmax(0.1, fish$W[alive_idx] + dW * dt)
      n_moved   <- 0
      prop_warm <- 0
      
    } else {
      # optimal_mover: each fish independently chooses its patch --------------
      # Pre-compute expected growth in each patch for every fish mass
      # (vectorised over W; movement cost applied inside choose_patch loop)
      moved <- logical(length(alive_idx))
      
      for (i in seq_along(alive_idx)) {
        idx    <- alive_idx[i]
        choice <- choose_patch(
          current_patch = fish$patch[idx],
          W       = fish$W[idx],
          T_warm  = T_warm, T_cold = T_cold,
          F_warm  = F_warm, F_cold = F_cold,
          p       = p
        )
        moved[i]         <- choice$moved
        fish$patch[idx]  <- choice$patch
      }
      
      # Apply growth conditional on chosen patch
      in_warm <- alive_idx[fish$patch[alive_idx] == "warm"]
      in_cold <- alive_idx[fish$patch[alive_idx] == "cold"]
      
      if (length(in_warm) > 0) {
        dW_warm <- net_growth(fish$W[in_warm], T_warm, F_warm, p)
        # Subtract movement cost for fish that just moved into warm patch
        moved_to_warm <- alive_idx[moved & fish$patch[alive_idx] == "warm"]
        dW_warm[in_warm %in% moved_to_warm] <-
          dW_warm[in_warm %in% moved_to_warm] - p$move_cost * fish$W[moved_to_warm]
        fish$W[in_warm] <- pmax(0.1, fish$W[in_warm] + dW_warm * dt)
      }
      if (length(in_cold) > 0) {
        dW_cold <- net_growth(fish$W[in_cold], T_cold, F_cold, p)
        moved_to_cold <- alive_idx[moved & fish$patch[alive_idx] == "cold"]
        dW_cold[in_cold %in% moved_to_cold] <-
          dW_cold[in_cold %in% moved_to_cold] - p$move_cost * fish$W[moved_to_cold]
        fish$W[in_cold] <- pmax(0.1, fish$W[in_cold] + dW_cold * dt)
      }
      
      n_moved   <- sum(moved)
      prop_warm <- mean(fish$patch[alive_idx] == "warm")
    }
    
    # --- Maturation ---------------------------------------------------------
    fish$mature[alive_idx] <- fish$W[alive_idx] >= p$W_maturity
    
    # --- Survival -----------------------------------------------------------
    surv_prob <- p$surv_base +
      p$surv_size_coef * pmax(0, fish$W[alive_idx] - p$W_surv_thresh)
    surv_prob <- pmin(surv_prob, 0.9999)
    survived  <- runif(length(alive_idx)) < surv_prob
    fish$alive[alive_idx[!survived]] <- FALSE
    
    # --- Age ----------------------------------------------------------------
    fish$age[alive_idx] <- fish$age[alive_idx] + dt
    
    # --- Reproduction (week 26 each year) -----------------------------------
    total_recruits <- 0
    if (week_of_year == 26) {
      mature_alive <- which(fish$alive & fish$mature)
      if (length(mature_alive) > 0) {
        eggs_per_fish  <- p$fecundity_a * (fish$W[mature_alive] ^ p$fecundity_b)
        total_eggs     <- sum(eggs_per_fish)
        total_recruits <- round(
          (p$BH_alpha * total_eggs) / (1 + p$BH_beta * total_eggs)
        )
      }
    }
    
    # --- Summary record -----------------------------------------------------
    alive_fish   <- fish[fish$alive, ]
    n_alive      <- nrow(alive_fish)
    results[[t]] <- data.frame(
      t              = t,
      week_of_year   = week_of_year,
      strategy       = strategy,
      n_alive        = n_alive,
      mean_W         = if (n_alive > 0) mean(alive_fish$W) else NA,
      sd_W           = if (n_alive > 0) sd(alive_fish$W)   else NA,
      mean_age       = if (n_alive > 0) mean(alive_fish$age) else NA,
      prop_mature    = if (n_alive > 0) mean(alive_fish$mature) else NA,
      total_recruits = total_recruits,
      T_warm         = T_warm,
      T_cold         = T_cold,
      F_warm         = F_warm,
      F_cold         = F_cold,
      prop_in_warm   = if (strategy == "optimal_mover") prop_warm else
        as.numeric(strategy == "resident_warm"),
      n_moved        = n_moved
    )
  }
  
  bind_rows(results)
}


# -----------------------------------------------------------------------------
# 6. Run simulation across all three strategies
# -----------------------------------------------------------------------------

cat("Running simulations...\n")
strategies  <- c("optimal_mover")
sim_results <- map_dfr(strategies, ~ {
  cat("  Strategy:", .x, "\n")
  simulate_strategy(.x, params)
})

# Rolling mean for smoother body-mass plot
sim_results <- sim_results %>%
  group_by(strategy) %>%
  arrange(t) %>%
  mutate(mean_W_smooth = zoo::rollmean(mean_W, k = 4, fill = NA,
                                       align = "right")) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 7. Derived population metrics
# -----------------------------------------------------------------------------

annual_census <- sim_results %>%
  filter(week_of_year == 52) %>%
  mutate(year = ceiling(t / 52))

r_estimates <- sim_results %>%
  filter(n_alive > 0) %>%
  group_by(strategy) %>%
  summarise(
    r_estimate         = coef(lm(log(n_alive) ~ t))[2],
    mean_final_W       = mean(tail(mean_W[!is.na(mean_W)], 52)),
    mean_prop_mature   = mean(prop_mature, na.rm = TRUE),
    total_recruitment  = sum(total_recruits, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Population-level summary ---\n")
print(r_estimates)

# -----------------------------------------------------------------------------
# 8. Visualization
# -----------------------------------------------------------------------------

theme_fish <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey50")
  )

strategy_colors <- c(
  "resident_warm" = "#D85A30",
  "resident_cold" = "#185FA5",
  "optimal_mover" = "#1D9E75"
)

strategy_labels <- c(
  "resident_warm" = "Resident (warm)",
  "resident_cold" = "Resident (cold)",
  "optimal_mover" = "Optimal mover"
)

# --- Plot 1: Mean body mass over time ---------------------------------------
p1 <- ggplot(sim_results,
             aes(x = t / 52, y = mean_W_smooth, color = strategy)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_vline(xintercept = 1:10, linetype = "dashed",
             color = "grey80", linewidth = 0.3) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  labs(title    = "Mean body mass over time",
       subtitle = "Rolling 4-week average",
       x = "Year", y = "Mean body mass (g)", color = NULL) +
  theme_fish

# --- Plot 2: Population abundance -------------------------------------------
p2 <- ggplot(sim_results, aes(x = t / 52, y = n_alive, color = strategy)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  labs(title    = "Population abundance over time",
       subtitle = "Individuals alive per strategy",
       x = "Year", y = "N individuals", color = NULL) +
  theme_fish

# --- Plot 3: Annual census body mass ----------------------------------------
p3 <- ggplot(annual_census,
             aes(x = year, y = mean_W, color = strategy)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  labs(title    = "Annual mean body mass at census",
       subtitle = "Recorded at week 52 of each year",
       x = "Year", y = "Mean body mass (g)", color = NULL) +
  theme_fish

# --- Plot 4: Emergent proportion of optimal movers in warm patch ------------
# This replaces the fixed-schedule plot: shows when the warm patch is
# energetically preferred across the year, averaged across years.

mover_patch_use <- sim_results %>%
  filter(strategy == "optimal_mover") %>%
  group_by(week_of_year) %>%
  summarise(
    mean_prop_warm = mean(prop_in_warm, na.rm = TRUE),
    mean_T_warm    = mean(T_warm),
    mean_T_cold    = mean(T_cold),
    .groups = "drop"
  )

p4 <- ggplot(mover_patch_use, aes(x = week_of_year)) +
  geom_ribbon(aes(ymin = mean_T_cold, ymax = mean_T_warm),
              fill = "grey92") +
  geom_line(aes(y = mean_T_warm, color = "Warm patch temp"), linewidth = 0.9) +
  geom_line(aes(y = mean_T_cold, color = "Cold patch temp"), linewidth = 0.9) +
  geom_line(aes(y = mean_prop_warm * 30, color = "Prop. in warm patch"),
            linewidth = 1, linetype = "solid") +
  scale_y_continuous(
    name     = "Temperature (°C)",
    sec.axis = sec_axis(~ . / 30, name = "Proportion in warm patch")
  ) +
  scale_color_manual(values = c(
    "Warm patch temp"      = "#D85A30",
    "Cold patch temp"      = "#185FA5",
    "Prop. in warm patch"  = "#1D9E75"
  )) +
  labs(title    = "Emergent patch use by optimal movers",
       subtitle = "Proportion in warm patch tracks temperature seasonality",
       x = "Week of year", color = NULL) +
  theme_fish

# --- Combine and save -------------------------------------------------------
combined_plot <- (p4 | p1) / (p3 | p2) +
  plot_annotation(
    title   = "Two-patch IBM: Energetics-based migration",
    caption = paste0(
      "Warm patch: T_mean=", params$T_warm_mean, "°C, food=", params$F_warm_mean,
      " | Cold patch: T_mean=", params$T_cold_mean, "°C, food=", params$F_cold_mean,
      " | move_cost=", params$move_cost, " | n=", params$n_individuals, " per strategy"
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("/mnt/user-data/outputs/two_patch_IBM_energetic_migration.png",
       combined_plot, width = 14, height = 9, dpi = 150)
cat("\nPlot saved.\n")

# -----------------------------------------------------------------------------
# 9. Strategy comparison table
# -----------------------------------------------------------------------------

comparison_table <- r_estimates %>%
  mutate(
    strategy          = strategy_labels[strategy],
    r_estimate        = round(r_estimate, 5),
    mean_final_W      = round(mean_final_W, 1),
    mean_prop_mature  = round(mean_prop_mature, 3),
    total_recruitment = round(total_recruitment)
  ) %>%
  rename(
    Strategy             = strategy,
    `r (intrinsic rate)` = r_estimate,
    `Mean body mass (g)` = mean_final_W,
    `Prop. mature`       = mean_prop_mature,
    `Total recruits`     = total_recruitment
  )

cat("\n--- Strategy comparison table ---\n")
print(comparison_table)
