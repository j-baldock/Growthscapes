#' ---
#' title: Simulate Habitat
#' ---
#' 


library(tidyverse)
library(ggpubr)
library(egg)


# Set seed for reproducibility
set.seed(123)

# 1. Create a date range (1 year)
dates <- seq(as.Date("2001-05-01"), as.Date("2020-12-31"), by="day")
n <- length(dates)
doy <- as.numeric(format(dates, "%j")) # Day of year

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


# define habitat specific ration size
ration_cold <- 0.10
ration_warm <- 0.10

# define random variation in ration size
ration_sd <- 0.00 # no variation for simplicity

# function to generate patch food availability at time t (stochastic)
patch_food <- function(patch) {
  if (patch == "warm") {
    max(0, min(1, rnorm(1, ration_warm, ration_sd)))
  } else {
    max(0, min(1, rnorm(1, ration_cold, ration_sd)))
  }
}

# define food in habitat tibble
habitat_df <- habitat_df %>% 
  mutate(ration_warm = replicate(n = dim(.)[1], patch_food("warm")),
         ration_cold = replicate(n = dim(.)[1], patch_food("cold")))

# plot ration size
habitat_df %>% ggplot() + 
  geom_line(aes(x = date, y = ration_warm), color = "red") +
  geom_line(aes(x = date, y = ration_cold), color = "blue") +
  theme_bw() + 
  xlab("Date") + 
  ylab("Simulated daily ration size")


# function to retrieve ration in patch = patch and at doy/time = t
get_patchration <- function(t, patch) {
  if (patch == "warm") {
    habitat_df$ration_warm[habitat_df$dayofsim == t]
  } else {
    habitat_df$ration_cold[habitat_df$dayofsim == t]
  }
}


# define habitat specific ration size
pcmax_cold <- 0.5
pcmax_warm <- 0.5

# define random variation in ration size
pcmax_sd <- 0.00 # no variation for simplicity

# function to generate patch food availability at time t (stochastic)
patch_pcmax <- function(patch) {
  if (patch == "warm") {
    max(0, min(1, rnorm(1, pcmax_warm, ration_sd)))
  } else {
    max(0, min(1, rnorm(1, pcmax_cold, ration_sd)))
  }
}

# define food in habitat tibble
habitat_df <- habitat_df %>% 
  mutate(pcmax_warm = replicate(n = dim(.)[1], patch_pcmax("warm")),
         pcmax_cold = replicate(n = dim(.)[1], patch_pcmax("cold")))

# plot ration size
habitat_df %>% ggplot() + 
  geom_line(aes(x = date, y = pcmax_warm), color = "red") +
  geom_line(aes(x = date, y = pcmax_cold), color = "blue") +
  theme_bw() + 
  xlab("Date") + 
  ylab("Simulated daily P_Cmax")


# function to retrieve ration in patch = patch and at doy/time = t
get_patchpcmax <- function(t, patch) {
  if (patch == "warm") {
    habitat_df$pcmax_warm[habitat_df$dayofsim == t]
  } else {
    habitat_df$pcmax_cold[habitat_df$dayofsim == t]
  }
}


str(habitat_df)
head(habitat_df)


#quarto::qmd_to_r_script("Habitat.qmd")

