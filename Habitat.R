#' ---
#' title: Simulate Habitat
#' ---
#' 


library(tidyverse)
library(ggpubr)
library(egg)


nyears_burnin <- 20
nyears_experiment <- 20

mindate <- as.Date("2001-05-01")
maxdate <- mindate + years(nyears_burnin) + years(nyears_experiment)

dates <- seq(mindate, maxdate, by="day")
n <- length(dates)
doy <- as.numeric(format(dates, "%j")) # Day of year


# Set seed for reproducibility
set.seed(123)

# 2. Simulate Stream Temperature (Sine wave + Noise)
# Parameters: Midpoint, Amplitude, Offset
base_temp <- 11        # Average annual temperature
amplitude <- 14         # Seasonal fluctuation amplitude
# Peak temperature occurs late in the year (July/Aug)
stream_temp <- base_temp + amplitude * sin(2 * pi * (doy - 100) / 365)

# 3. Add random noise (daily weather variation)
noise <- rnorm(n, mean = 0, sd = 0) # no variation for simplicity
simulated_temp <- stream_temp + noise

# 4. Force values near 0 to small postive, max = 25
simulated_temp[simulated_temp < 0.5] <- 0.5
simulated_temp[simulated_temp > 25] <- 25

# 4. Simulate cold regime using a multiplier
simulated_temp_cold <- simulated_temp*0.6

# 4. Create Data Frame
habitat_df <- tibble(date = dates, 
                     doy = doy, 
                     dayofsim = 1:length(dates),
                     temp_warm = round(simulated_temp, digits = 1), 
                     temp_cold = round(simulated_temp_cold, digits = 1))

# 5. Plot the simulation
summer_shading <- data.frame(
  xmin = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-06-01")),
  xmax = as.Date(paste0(unique(format(habitat_df$date, "%Y")), "-08-31")),
  ymin = -Inf, ymax = Inf
)

habitat_df %>% ggplot() +
  geom_rect(data = summer_shading,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(aes(x = date, y = temp_warm), color = "red") +
  geom_line(aes(x = date, y = temp_cold), color = "blue") +
  theme_bw() + 
  xlab("Date") + 
  ylab("Simulated daily stream temperature (deg C)") + 
  labs(title = "Daily stream temperature with Jun-Aug shaded")


# function to retrieve water temperature in patch = patch and at doy/time = t
get_patchtemp <- function(t, patch) {
  if (patch == "warm") {
    habitat_df$temp_warm[habitat_df$dayofsim == t]
  } else {
    habitat_df$temp_cold[habitat_df$dayofsim == t]
  }
}


# --- Maximum proportion of maximum consumption ---
pcmax_warm_val <- 0.5
pcmax_cold_val <- 0.5

# --- Patch area (arbitrary units; keep at 1.0 to preserve default behaviour) ---
A_warm_val <- 1.0
A_cold_val <- 1.0

# --- Half-saturation density for hyperbolic ration reduction ---
# Competitor density (fish per unit area) at which ration is halved
K_warm_val <- 500
K_cold_val <- 500

# --- Maximum daily survival probability (size-based survival ceiling) ---
S_max_warm_val <- 0.9994
S_max_cold_val <- 0.9994

# Append as time-invariant columns (replace with time-varying vectors to simulate change)
habitat_df <- habitat_df %>%
  mutate(
    pcmax_warm = pcmax_warm_val,
    pcmax_cold = pcmax_cold_val,
    A_warm     = A_warm_val,
    A_cold     = A_cold_val,
    K_warm     = K_warm_val,
    K_cold     = K_cold_val,
    S_max_warm = S_max_warm_val,
    S_max_cold = S_max_cold_val
  )


# First row of the experiment period
exp_start_date <- mindate + lubridate::years(nyears_burnin)
d_exp_start    <- which(habitat_df$date == exp_start_date)

# --- Target values for each parameter (NA = no change) ---
pcmax_warm_target  <- NA
pcmax_cold_target  <- NA
A_warm_target      <- 2
A_cold_target      <- 0
K_warm_target      <- NA
K_cold_target      <- NA
S_max_warm_target  <- NA
S_max_cold_target  <- NA

# --- Transition: 0 = step change; > 0 = linear ramp over this many years ---
ramp_years <- nyears_experiment

# Linearly ramps col from its end-of-burn-in value to target_val over ramp_yrs years,
# then holds target_val for the remainder of the experiment period.
apply_ramp <- function(col, d_start, target_val, ramp_yrs) {
  baseline  <- col[d_start - 1]
  n_ramp    <- max(round(ramp_yrs * 365.25), 1)  # min 1 to avoid /0
  exp_idx   <- d_start:length(col)
  ramp_frac <- pmin(seq_along(exp_idx) / n_ramp, 1)
  col[exp_idx] <- baseline + (target_val - baseline) * ramp_frac
  col
}

params_to_change <- list(
  pcmax_warm = pcmax_warm_target,
  pcmax_cold = pcmax_cold_target,
  A_warm     = A_warm_target,
  A_cold     = A_cold_target,
  K_warm     = K_warm_target,
  K_cold     = K_cold_target,
  S_max_warm = S_max_warm_target,
  S_max_cold = S_max_cold_target
)

for (param in names(params_to_change)) {
  target <- params_to_change[[param]]
  if (!is.na(target)) {
    habitat_df[[param]] <- apply_ramp(habitat_df[[param]], d_exp_start, target, ramp_years)
  }
}


p1 <- habitat_df %>% ggplot() +
  geom_vline(xintercept = exp_start_date, linetype = "dashed", color = "grey40") +
  geom_line(aes(x = date, y = pcmax_warm), color = "red") +
  geom_line(aes(x = date, y = pcmax_cold), color = "blue") +
  theme_bw() +
  xlab("Date") +
  ylab("Maximum P_Cmax")

p2 <- habitat_df %>% ggplot() +
  geom_vline(xintercept = exp_start_date, linetype = "dashed", color = "grey40") +
  geom_line(aes(x = date, y = A_warm), color = "red") +
  geom_line(aes(x = date, y = A_cold), color = "blue") +
  theme_bw() +
  xlab("Date") +
  ylab("Patch area (unitless)")

p3 <- habitat_df %>% ggplot() +
  geom_vline(xintercept = exp_start_date, linetype = "dashed", color = "grey40") +
  geom_line(aes(x = date, y = K_warm), color = "red") +
  geom_line(aes(x = date, y = K_cold), color = "blue") +
  theme_bw() +
  xlab("Date") +
  ylab("Strength of density-dependence\n(half-saturation density)")

p4 <- habitat_df %>% ggplot() +
  geom_vline(xintercept = exp_start_date, linetype = "dashed", color = "grey40") +
  geom_line(aes(x = date, y = S_max_warm), color = "red") +
  geom_line(aes(x = date, y = S_max_cold), color = "blue") +
  theme_bw() +
  xlab("Date") +
  ylab("Maximum probability of survival")

egg::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)


str(habitat_df)
head(habitat_df)


#quarto::qmd_to_r_script("Habitat.qmd")

