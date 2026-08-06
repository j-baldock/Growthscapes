#' ---
#' title: IBM Simulation Loop
#' ---
#' 


library(tidyverse)
library(lubridate)


run_simulation <- function(habitat_df, params = list(), wt_growth) {

  # ── Parameters ──────────────────────────────────────────────────────────────
  # Configuration options:
  move_stochastic   <- params$move_stochastic   %||% "prob"           # individual movement stochasticity: options = "none", "prob", "indiv"
  food_densdepen    <- params$food_densdepen     %||% "hyperbolic"     # density dependence: options = "fullerton", "hyperbolic"
  sense_environment <- params$sense_environment %||% "density_current" # how omniscient are fish to temperature and density-dependent growth across the habitat patches?
                                                                        #   - "abiotic_only" = Omniscient to abiotic (temp-based) GP only
                                                                        #   - "density_current" = Sense density dependent GP in current habitat, and compare to abiotic GP in the alternative habitat (perhaps the most biologically realistic)
                                                                        #   - "density_all" = Omniscient to density dependent GP in all habitats

  # misc. parameters
  n_fish_per_strat  <- params$n_fish_per_strat  %||% 50        # fish per strategy
  start_wt          <- params$start_wt          %||% 0.5       # initial weight (g)
  wt_sd             <- params$wt_sd             %||% 0.05      # SD of starting weight
  MaxDensity4Growth <- params$MaxDensity4Growth %||% 50        # growth is not depressed further above this density (currently this is arbitrary)
  movecost_c        <- params$movecost_c        %||% 0.1       # movement cost allometric function: `c` = scaling constant
  movecost_b        <- params$movecost_b        %||% 0.3       # movement cost allometric function `b` = decay exponent (higher = steeper drop-off)
  sigma_bold        <- params$sigma_bold        %||% 0.002     # SD of boldness distribution, for dispersal propensity (g/g/d units)
  tau               <- params$tau               %||% 0.001     # sensitivity to growth difference (for probabilistic/softmax variation in movement)
  s_min             <- params$s_min             %||% 0.96      # minimum size-based daily survival probability
  s_w0              <- params$s_w0              %||% 7         # inflection point of size-based survival logistic curve
  s_k               <- params$s_k               %||% 1         # steepness of size-based survival logistic curve
  T1_mort           <- params$T1_mort           %||% 30        # temperature (°C) at which daily survival from thermal stress = 0.1
  T9_mort           <- params$T9_mort           %||% 25.8      # temperature (°C) at which daily survival from thermal stress = 0.9
  K9_starv          <- params$K9_starv          %||% 0.55      # relative condition (W_current/W_peak) at which starvation survival = 0.9
  K1_starv          <- params$K1_starv          %||% 0.45      # relative condition at which starvation survival = 0.1
  dominance_beta    <- params$dominance_beta    %||% 1         # size-based dominance exponent for effective density (0 = pure scramble, 1 = linear dominance, Inf = pure contest)
  egg_wt            <- params$egg_wt            %||% 0.07      # weight of a single egg (g); energetic cost per offspring
  repro_cost        <- params$repro_cost        %||% 0.2       # energetic cost of reproduction in percentage of body mass
  sigma_fecund      <- params$sigma_fecund      %||% 0.085     # random variation in weight-fecundity relationship
  egg_surv                   <- params$egg_surv                   %||% 0.1    # egg-to-fry (hatching) survival probability ("pass-through" filter for fecundity)
  age_structured_competition <- params$age_structured_competition %||% TRUE   # if TRUE, fry (age-0) compete only with fry; age-1+ compete only with age-1+
  crit_pcmax_lo              <- params$crit_pcmax_lo              %||% 0.30   # pcmax_adjusted_dd at which ~1% critical-period survival is achieved
  crit_pcmax_hi              <- params$crit_pcmax_hi              %||% 0.40   # pcmax_adjusted_dd at which ~99% critical-period survival is achieved (negligible cost)
  crit_period_days           <- params$crit_period_days           %||% 60L    # length of the critical period (days post-hatch)
  sim_seed                   <- params$sim_seed                   %||% 7843

  # Competition grouping vector: used in group_by(across(all_of(eff_grp))) in Steps 2 and 3.
  # When age_structured_competition = TRUE, eff_density is computed within patch × age class.
  # When FALSE, all fish in the same patch compete (original behaviour).
  eff_grp <- if (age_structured_competition) c("patch", "age_class") else "patch"

  # A_warm, A_cold, K_warm, K_cold, S_max_warm, S_max_cold are read from habitat_df each day (see Habitat.qmd)

  # Fullerton has movement cost coded slightly differently, where movement cost is a function of
  # the instantaneous growth rate, distance moved, and a constant...rather than simply a fixed cost
  # as is the case here. But because we only have two patches that are spatially implicit, we can
  # ignore this more complex parameterization...at least for now.

  # ── Local helpers ────────────────────────────────────────────────────────────
  get_patchtemp <- function(t, patch) {
    if (patch == "warm") habitat_df$temp_warm[habitat_df$dayofsim == t]
    else                 habitat_df$temp_cold[habitat_df$dayofsim == t]
  }

  # Lookup sequences — must match the axes of the wt_growth array
  wt_seq_ibm <- seq(0.1,   25,  0.1)    # water temperature (250 values)
  ra_seq_ibm <- seq(0.001, 0.4, 0.001)  # ration (400 values)
  min_wt     <- 0.25                     # absolute minimum fish weight (g)

  # ── Initialise fish ─────────────────────────────────────────────────────────
  set.seed(sim_seed)

  # optimal movers are randomly placed into either warm or cold patches
  fish_pop <- tibble(
    strategy = "optimal_mover",
    patch    = sample(c("warm", "cold"), size = n_fish_per_strat, replace = TRUE),
    weight   = pmax(min_wt, rnorm(n_fish_per_strat, start_wt, wt_sd))
  ) |>
    mutate(
      pid               = row_number(),
      cmax_allometric   = NA_real_,
      pcmax_baseline    = NA_real_,
      pcmax_adjusted    = NA_real_,
      pcmax_adjusted_dd = NA_real_,
      func_temp         = NA_real_,
      ration            = NA_real_,
      # Individual movement threshold — drawn once, fixed for life
      # Positive = reluctant to move (needs a clear advantage)
      # Negative = bold/prone to move (moves even when slightly disadvantaged)
      move_threshold    = rnorm(n_fish_per_strat, mean = 0, sd = sigma_bold),
      peak_weight       = weight,      # track historical max weight for condition factor
      age_days          = 0L,          # age in days, incremented each iteration independent of d
      spawned_this_year = FALSE,       # has this fish already spawned in the current calendar year?
      parent_pid        = NA_integer_, # pid of parent fish (NA for founder cohort)
      cohort            = year(habitat_df$date[1]),  # birth year (cohort)
      prob_surv         = 1,
      survive           = 1
    ) |>
    group_by(patch) |>
    mutate(density = n() / if_else(first(patch) == "warm",
                                   habitat_df$A_warm[1], habitat_df$A_cold[1])) |>
    ungroup()

  # ── Pre-allocate storage ─────────────────────────────────────────────────────
  n_days <- nrow(habitat_df)
  n_fish <- nrow(fish_pop)

  # fish_registry: grows throughout the simulation as offspring are born.
  # Used in place of fish_pop_init in summaries once reproduction is active.
  # One row per fish (founders + all offspring), ordered by pid.
  fish_registry <- fish_pop |>
    select(pid, strategy, parent_pid, cohort) |>
    mutate(birth_dayofsim = 0L)

  # Pre-allocate per-day record list; one data.frame per loop iteration.
  # Avoids large sparse matrices — only alive (fish × day) combos are stored.
  ibm_records <- vector("list", n_days)
  switches    <- integer(n_fish)  # count patch switches per fish (indexed by pid)
  next_pid    <- n_fish + 1L      # next available pid (grows as offspring are born)

  # Spawn log: one row per spawning event, recording parent-offspring relationships
  spawn_log <- tibble(
    parent_pid          = integer(),  # pid of the spawning parent
    dayofsim            = integer(),  # simulation day of spawning
    n_offspring         = integer(),  # number of offspring produced
    weight              = numeric(),  # pre-spawn body weight (g)
    condition           = numeric(),  # pre-spawn relative condition (weight / peak_weight)
    offspring_pid_start = integer(),  # first pid assigned to offspring
    offspring_pid_end   = integer()   # last pid assigned to offspring
  )

  # ── Simulation loop ──────────────────────────────────────────────────────────
  set.seed(sim_seed)

  for (d in seq_len(n_days)) {
    t            <- habitat_df$dayofsim[d]
    current_year <- year(habitat_df$date[d])     # calendar year of this simulation day
    fish_pop$age_days <- fish_pop$age_days + 1L  # increment age for all living fish

    # 1. READ TIME-VARYING HABITAT PARAMETERS ─────────────────────────────────
    A_warm     <- habitat_df$A_warm[d]
    A_cold     <- habitat_df$A_cold[d]
    K_warm     <- habitat_df$K_warm[d]
    K_cold     <- habitat_df$K_cold[d]
    S_max_warm <- habitat_df$S_max_warm[d]
    S_max_cold <- habitat_df$S_max_cold[d]
    pcmax_warm <- habitat_df$pcmax_warm[d]
    pcmax_cold <- habitat_df$pcmax_cold[d]

    # Temperature indices for this day — shared by all fish
    T_warm      <- get_patchtemp(t, "warm")             # get warm patch temp on day t
    T_cold      <- get_patchtemp(t, "cold")             # get cold patch temp on day t
    wt_idx_warm <- which.min(abs(wt_seq_ibm - T_warm))  # water temp index for warm patch on day t
    wt_idx_cold <- which.min(abs(wt_seq_ibm - T_cold))  # water temp index for cold patch on day t

    # If a patch has collapsed to zero area, force all fish out of it before any
    # density calculations (A = 0 would cause 0/0 = NaN in the density step).
    if (A_cold == 0 && any(fish_pop$patch == "cold")) {
      fish_pop$patch[fish_pop$patch == "cold"] <- "warm"
    }
    if (A_warm == 0 && any(fish_pop$patch == "warm")) {
      fish_pop$patch[fish_pop$patch == "warm"] <- "cold"
    }

    # If both patches have zero area, there is no valid habitat: exit immediately.
    # (When A_cold = A_warm = 0, the guard above cycles fish cold→warm→cold, leaving
    # them in cold with A_cold = 0. Any 1-fish age-class group then computes
    # (n-1)/A = 0/0 = NaN, which propagates through eff_density → ration → error.)
    if (A_cold == 0 && A_warm == 0) break

    # Calculate habitat, temperature, and body size dependent rations following convo. with Jonny
    fish_pop <- fish_pop |>
      mutate(
        func_temp        = ifelse(patch == "warm", fncTempDepend(T_warm), fncTempDepend(T_cold)),
        pcmax_baseline   = ifelse(patch == "warm", pcmax_warm, pcmax_cold),
        pcmax_adjusted   = ifelse(func_temp < pcmax_baseline, func_temp, pcmax_baseline),
        cmax_allometric  = fncAllomCmax(weight)
      )

    # 2. CHOOSE PATCH: OPTIMAL MOVERS ─────────────────────────────────────────
    # Mirrors fncGrowthPossible() logic, vectorized across all movers at once.
    # Patch choice senses baseline temperature and ration, but not density.
    mover_rows <- which(fish_pop$strategy == "optimal_mover")  # get row indices for optimal movers

    if (length(mover_rows) > 0) {
      # Mass indices for each mover (clamped to array bounds)
      ma_idx_vec    <- pmax(1L, pmin(4500L, round(fish_pop$weight[mover_rows])))  # get row index for current mass
      ## 2.A.
      in_warm       <- fish_pop$patch[mover_rows] == "warm"  # is the fish in warm patch? T/F
      move_cost_vec <- fncMoveCost_allometric(fish_pop$weight[mover_rows],
                                              c = movecost_c, b = movecost_b)  # calculate size-dependent cost of movement

      # Pre-compute each mover's hypothetical ration in BOTH patches for patch choice sensing.
      # This initial step is done without accounting for the effect of density on p_cmax/ration.
      # pcmax scalars are day-level constants; cmax_allometric is already in fish_pop from step 1.
      pcmax_warm_adj <- min(fncTempDepend(T_warm), pcmax_warm)
      pcmax_cold_adj <- min(fncTempDepend(T_cold), pcmax_cold)
      ra_warm_vec    <- fish_pop$cmax_allometric[mover_rows] * pcmax_warm_adj
      ra_cold_vec    <- fish_pop$cmax_allometric[mover_rows] * pcmax_cold_adj
      ra_idx_warm_vec <- pmax(1L, pmin(400L, map_int(ra_warm_vec, ~which.min(abs(ra_seq_ibm - .x)))))
      ra_idx_cold_vec <- pmax(1L, pmin(400L, map_int(ra_cold_vec, ~which.min(abs(ra_seq_ibm - .x)))))

      if (sense_environment == "abiotic_only") {
        ### 2.A.1. Individuals sense only temperature-based growth potential across all habitats
        g_warm_vec <- wt_growth[cbind(wt_idx_warm, ra_idx_warm_vec, ma_idx_vec)]
        g_cold_vec <- wt_growth[cbind(wt_idx_cold, ra_idx_cold_vec, ma_idx_vec)]
        g_stay     <- ifelse(in_warm, g_warm_vec, g_cold_vec)  # growth if staying in current patch
        g_move_net <- ifelse(in_warm, g_cold_vec - move_cost_vec, g_warm_vec - move_cost_vec)  # growth if moving to alternate patch

      } else if (sense_environment == "density_current") {
        ### 2.A.2. Individuals sense density-dependent GP in current habitat, abiotic GP in alternative habitat
        # Each mover uses its own dominance-weighted effective competitor density (fish per unit
        # area, self excluded) in its current patch.
        # Effective competitor density: grouped by patch × age_class when
        # age_structured_competition = TRUE, by patch only when FALSE.
        eff_density_all      <- fish_pop |>
          mutate(age_class = if_else(cohort == current_year, "age0", "age1plus")) |>
          group_by(across(all_of(eff_grp))) |>
          mutate(
            A_patch     = if_else(first(patch) == "warm", A_warm, A_cold),
            eff_density = (fncEffDensity(weight, beta = dominance_beta) - 1) / A_patch
          ) |>
          ungroup() |>
          pull(eff_density)
        eff_dens_current_vec <- eff_density_all[mover_rows]

        # DD-adjusted ration uses each mover's per-area effective competitor density
        ra_warm_dd_vec     <- fish_pop$cmax_allometric[mover_rows] * pcmax_warm_adj *
                              (K_warm / (K_warm + eff_dens_current_vec))
        ra_cold_dd_vec     <- fish_pop$cmax_allometric[mover_rows] * pcmax_cold_adj *
                              (K_cold / (K_cold + eff_dens_current_vec))
        ra_idx_warm_dd_vec <- pmax(1L, pmin(400L, map_int(ra_warm_dd_vec, ~which.min(abs(ra_seq_ibm - .x)))))
        ra_idx_cold_dd_vec <- pmax(1L, pmin(400L, map_int(ra_cold_dd_vec, ~which.min(abs(ra_seq_ibm - .x)))))

        # Warm fish sense dd-adjusted warm / abiotic cold; cold fish sense abiotic warm / dd-adjusted cold
        g_warm_vec <- wt_growth[cbind(wt_idx_warm,
                                      ifelse(in_warm, ra_idx_warm_dd_vec, ra_idx_warm_vec),
                                      ma_idx_vec)]
        g_cold_vec <- wt_growth[cbind(wt_idx_cold,
                                      ifelse(in_warm, ra_idx_cold_vec, ra_idx_cold_dd_vec),
                                      ma_idx_vec)]
        g_stay     <- ifelse(in_warm, g_warm_vec, g_cold_vec)
        g_move_net <- ifelse(in_warm, g_cold_vec - move_cost_vec, g_warm_vec - move_cost_vec)

      } else if (sense_environment == "density_all") {
        ### 2.A.3. Individuals sense density-dependent GP across all habitats
        # Each mover evaluates both patches using its own dominance-weighted effective density:
        #   Current patch: effective density among current co-occupants (fish is already there).
        #   Alternative patch: counterfactual effective density — how the fish would rank if it
        #     joined the alternative patch's current occupants as a new entrant.
        mover_weights_vec <- fish_pop$weight[mover_rows]
        age_class_all     <- if_else(fish_pop$cohort == current_year, "age0", "age1plus")
        mover_age_class   <- age_class_all[mover_rows]

        # Current-patch effective competitor density (per unit area, self excluded).
        # Respects age_structured_competition toggle via eff_grp.
        eff_density_all      <- fish_pop |>
          mutate(age_class = age_class_all) |>
          group_by(across(all_of(eff_grp))) |>
          mutate(
            A_patch     = if_else(first(patch) == "warm", A_warm, A_cold),
            eff_density = (fncEffDensity(weight, beta = dominance_beta) - 1) / A_patch
          ) |>
          ungroup() |>
          pull(eff_density)
        eff_dens_current_vec <- eff_density_all[mover_rows]

        # Helper: vectorized outer-product effective density (per unit area)
        eff_join <- function(focal_wts, comp_wts, area) {
          if (length(focal_wts) == 0) return(numeric(0))
          if (length(comp_wts)  == 0) return(rep(0, length(focal_wts)))
          rowSums(outer(focal_wts, comp_wts,
                        FUN = function(mi, mj) ifelse(mj >= mi, 1, (mj/mi)^dominance_beta))) / area
        }

        # Counterfactual effective competitor density in the alternative patch if the mover joins.
        if (age_structured_competition) {
          # Competitors restricted to the same age class as the focal mover.
          warm_wts_age0  <- fish_pop$weight[fish_pop$patch == "warm" & age_class_all == "age0"]
          warm_wts_age1p <- fish_pop$weight[fish_pop$patch == "warm" & age_class_all == "age1plus"]
          cold_wts_age0  <- fish_pop$weight[fish_pop$patch == "cold" & age_class_all == "age0"]
          cold_wts_age1p <- fish_pop$weight[fish_pop$patch == "cold" & age_class_all == "age1plus"]

          idx_age0  <- which(mover_age_class == "age0")
          idx_age1p <- which(mover_age_class == "age1plus")

          eff_dens_if_join_warm <- numeric(length(mover_rows))
          eff_dens_if_join_cold <- numeric(length(mover_rows))

          if (length(idx_age0)  > 0) {
            eff_dens_if_join_warm[idx_age0]  <- eff_join(mover_weights_vec[idx_age0],  warm_wts_age0,  A_warm)
            eff_dens_if_join_cold[idx_age0]  <- eff_join(mover_weights_vec[idx_age0],  cold_wts_age0,  A_cold)
          }
          if (length(idx_age1p) > 0) {
            eff_dens_if_join_warm[idx_age1p] <- eff_join(mover_weights_vec[idx_age1p], warm_wts_age1p, A_warm)
            eff_dens_if_join_cold[idx_age1p] <- eff_join(mover_weights_vec[idx_age1p], cold_wts_age1p, A_cold)
          }
        } else {
          # All fish in the alternative patch are competitors (original behaviour).
          warm_weights          <- fish_pop$weight[fish_pop$patch == "warm"]
          cold_weights          <- fish_pop$weight[fish_pop$patch == "cold"]
          eff_dens_if_join_warm <- eff_join(mover_weights_vec, warm_weights, A_warm)
          eff_dens_if_join_cold <- eff_join(mover_weights_vec, cold_weights, A_cold)
        }

        # Assign per-area effective competitor density per patch per mover:
        #   current patch → eff_dens_current_vec (fish is already there, self excluded)
        #   alternative   → eff_dens_if_join_*   (existing residents only, self excluded)
        eff_dens_warm_vec <- ifelse(in_warm, eff_dens_current_vec, eff_dens_if_join_warm)
        eff_dens_cold_vec <- ifelse(in_warm, eff_dens_if_join_cold, eff_dens_current_vec)

        ra_warm_dd_vec     <- fish_pop$cmax_allometric[mover_rows] * pcmax_warm_adj *
                              (K_warm / (K_warm + eff_dens_warm_vec))
        ra_cold_dd_vec     <- fish_pop$cmax_allometric[mover_rows] * pcmax_cold_adj *
                              (K_cold / (K_cold + eff_dens_cold_vec))
        ra_idx_warm_dd_vec <- pmax(1L, pmin(400L, map_int(ra_warm_dd_vec, ~which.min(abs(ra_seq_ibm - .x)))))
        ra_idx_cold_dd_vec <- pmax(1L, pmin(400L, map_int(ra_cold_dd_vec, ~which.min(abs(ra_seq_ibm - .x)))))

        # Both patches sensed with dominance-adjusted individual rations
        g_warm_vec <- wt_growth[cbind(wt_idx_warm, ra_idx_warm_dd_vec, ma_idx_vec)]
        g_cold_vec <- wt_growth[cbind(wt_idx_cold, ra_idx_cold_dd_vec, ma_idx_vec)]
        g_stay     <- ifelse(in_warm, g_warm_vec, g_cold_vec)
        g_move_net <- ifelse(in_warm, g_cold_vec - move_cost_vec, g_warm_vec - move_cost_vec)
      }

      prev_patch <- fish_pop$patch[mover_rows]

      ## 2.B. Impose stochasticity in movement
      if (move_stochastic == "none") {
        ### 2.B.1. Deterministic: All movers do the same thing
        # Update patch based on whether moving is beneficial.
        # 1. If fish is in warm patch and growth potential of moving exceeds growth potential of staying, move to cold.
        # 2. If fish is in cold patch and growth potential of moving exceeds growth potential of staying, move to warm.
        fish_pop$patch[mover_rows] <- case_when(
          in_warm  & g_move_net > g_stay ~ "cold",
          !in_warm & g_move_net > g_stay ~ "warm",
          TRUE                           ~ fish_pop$patch[mover_rows]
        )
      } else if (move_stochastic == "prob") {
        ### 2.B.2. Probabilistic (softmax) rule: Instead of a hard threshold, fish switch with a
        ### probability that increases with the growth advantage. One parameter: τ (sensitivity;
        ### small = nearly deterministic, large = nearly random).
        # Use cost-adjusted net growth rates: staying is free, moving incurs move_cost_vec.
        # g_stay     = growth in current patch (no cost)
        # g_move_net = growth in alternate patch - move_cost_vec
        # Re-map to warm/cold perspective for fncMoveSoftmax:
        gwarm_eff  <- ifelse(in_warm, g_stay, g_move_net)  # effective growth if occupying warm patch
        gcold_eff  <- ifelse(in_warm, g_move_net, g_stay)  # effective growth if occupying cold patch
        p_warm     <- fncMoveSoftmax(gwarm = gwarm_eff, gcold = gcold_eff, tau = tau)  # calculate probability of selecting warm habitat
        choose_warm <- runif(length(mover_rows)) < p_warm  # choose warm if p_warm exceeds a random number between 0 and 1
        fish_pop$patch[mover_rows] <- if_else(choose_warm, "warm", "cold")  # update patch
      } else if (move_stochastic == "indiv") {
        ### 2.B.3. Individual-level movement threshold: each fish has a unique threshold for
        ### movement, representing individual variation in boldness/dispersal propensity
        fish_pop$patch[mover_rows] <- case_when(
          in_warm  & (g_move_net - g_stay) > fish_pop$move_threshold[mover_rows] ~ "cold",
          !in_warm & (g_move_net - g_stay) > fish_pop$move_threshold[mover_rows] ~ "warm",
          TRUE                                                                    ~ prev_patch
        )
      }

      # Tally switches (indexed by pid so counts remain correct as fish die)
      switched <- fish_pop$patch[mover_rows] != prev_patch
      switches[fish_pop$pid[mover_rows]] <- switches[fish_pop$pid[mover_rows]] +
                                            as.integer(switched)
    } # end patch choice

    # Recompute func_temp and pcmax_adjusted based on post-choice patch.
    # Required because some fish may have switched patches during step 2;
    # step 1 values reflect the pre-choice patch and would be stale for switchers.
    fish_pop <- fish_pop |>
      mutate(
        func_temp      = ifelse(patch == "warm", fncTempDepend(T_warm), fncTempDepend(T_cold)),
        pcmax_baseline = ifelse(patch == "warm", pcmax_warm, pcmax_cold),
        pcmax_adjusted = ifelse(func_temp < pcmax_baseline, func_temp, pcmax_baseline)
      )

    # 3. UPDATE HABITAT QUALITY (DENSITY-DEPENDENT PCMAX / RATION) ─────────────
    # Re-apply zero-area guard after patch choice: movers may have re-selected a
    # collapsed patch (A = 0) during step 2. Force them back to the surviving patch
    # before density calculations to prevent division-by-zero / NaN rations.
    if (A_cold == 0 && any(fish_pop$patch == "cold")) {
      fish_pop$patch[fish_pop$patch == "cold"] <- "warm"
    }
    if (A_warm == 0 && any(fish_pop$patch == "warm")) {
      fish_pop$patch[fish_pop$patch == "warm"] <- "cold"
    }

    # Update density (fish per unit area) and dominance-weighted effective competitor density
    # after fish have moved. density counts all fish in the patch regardless of age class.
    # eff_density is computed within the grouping defined by eff_grp:
    #   age_structured_competition = TRUE  → group by patch × age_class (fry vs. age-1+)
    #   age_structured_competition = FALSE → group by patch only (all fish compete together)
    # This reflects ontogenetic habitat segregation in salmonids when enabled. Self is
    # excluded via the (-1) correction from fncEffDensity.
    fish_pop <- fish_pop |>
      mutate(age_class = if_else(cohort == current_year, "age0", "age1plus")) |>
      group_by(patch) |>
      mutate(
        A_patch = if_else(first(patch) == "warm", A_warm, A_cold),
        density = n() / A_patch
      ) |>
      group_by(across(all_of(eff_grp))) |>
      mutate(
        eff_density = (fncEffDensity(weight, beta = dominance_beta) - 1) / A_patch
      ) |>
      ungroup() |>
      select(-A_patch, -age_class)

    if (food_densdepen == "fullerton") {
      ### 3a. Fullerton approach (uses raw patch density; dominance not yet integrated here)
      fdens <- fish_pop$density
      fdens[fdens > MaxDensity4Growth] <- MaxDensity4Growth
      density_effect <- fncRescale((1 - c(fdens, 0.01, MaxDensity4Growth)), c(0.5, 1))
      density_effect <- density_effect[-c(length(density_effect), length(density_effect) - 1)]
      fish_pop <- fish_pop |> mutate(ration = ration * density_effect)
    } else if (food_densdepen == "hyperbolic") {
      ### 3b. Hyperbolic reduction — each fish's P_Cmax is penalised by its own effective density
      # (eff_density, in fish per unit area, self excluded). Large dominant fish face a lower
      # eff_density → less suppression; small subordinate fish face a higher eff_density → more.
      # No -1 correction needed: self is already excluded from eff_density.
      fish_pop <- fish_pop |>
        mutate(
          k                 = if_else(patch == "warm", K_warm, K_cold),
          pcmax_adjusted_dd = pcmax_adjusted * (k / (k + eff_density)),
          ration            = cmax_allometric * pcmax_adjusted_dd
        ) |>
        select(-k)
    }

    # 4. GROW FISH ─────────────────────────────────────────────────────────────
    # Growth lookup for all fish via fncGrowthFish. Growth in g/g/d.
    growth_df <- fncGrowthFish(NA, as.data.frame(fish_pop), T_warm, T_cold)

    # Apply movement cost: deduct from growth rate for fish that actually switched patches.
    # move_cost_vec is in g/g/d and is indexed over mover_rows; switched is the logical
    # subset of those rows that changed patch this step.
    if (length(mover_rows) > 0 && any(switched)) {
      switched_rows <- mover_rows[switched]
      growth_df$growth[switched_rows] <- growth_df$growth[switched_rows] - move_cost_vec[switched]
    }

    # Update weight: W_{t+1} = W_t + (growth_rate * W_t)
    fish_pop$weight      <- fish_pop$weight + (growth_df$growth * fish_pop$weight)
    # Update peak weight: ratchet up only, never down
    fish_pop$peak_weight <- pmax(fish_pop$peak_weight, fish_pop$weight)

    # 5. SPAWNING AND REPRODUCTION ─────────────────────────────────────────────
    # Reset spawned_this_year flag at the calendar year boundary
    if (d > 1) {
      prev_year <- year(habitat_df$date[d - 1])
      if (current_year != prev_year) fish_pop$spawned_this_year <- FALSE
    }

    # Day-of-year and relative condition needed for spawning probability
    doy_d         <- yday(habitat_df$date[d])
    condition_spw <- fish_pop$weight / fish_pop$peak_weight

    # Combined daily spawning probability: size × condition × date
    p_spawn <- fncMaturitySize(fish_pop$weight) *
               fncMaturityCondition(condition_spw) *
               fncMaturityDate(doy_d)

    # Fish that have already spawned this year are ineligible
    p_spawn[fish_pop$spawned_this_year] <- 0

    # Bernoulli trial: which fish spawn today?
    spawns <- as.logical(rbinom(nrow(fish_pop), size = 1, prob = p_spawn))

    # For spawning fish: calculate fecundity, impose energetic cost, mark as spawned
    if (any(spawns)) {
      spawner_idx <- which(spawns)

      # Fecundity: number of eggs as a function of weight, with log-scale noise
      n_eggs      <- round(fncFecundBromage(fish_pop$weight[spawner_idx],
                                            sigma = sigma_fecund, survival = egg_surv))

      # Capture pre-cost weight for logging (spawn decision was made at this weight)
      weight_at_spawn  <- fish_pop$weight[spawner_idx]

      # Energetic cost of reproduction: deduct from parent weight
      # Weight floor at start_wt — spawning cannot reduce a fish below hatch weight
      reproduction_cost <- weight_at_spawn * repro_cost
      fish_pop$weight[spawner_idx]          <- pmax(start_wt, weight_at_spawn - reproduction_cost)

      # Flag these fish as having spawned; they are ineligible to spawn again this year
      fish_pop$spawned_this_year[spawner_idx] <- TRUE

      # Create offspring rows and add to population --------------------------------
      new_fish_list <- vector("list", length(spawner_idx))

      for (i in seq_along(spawner_idx)) {
        si      <- spawner_idx[i]
        parent  <- fish_pop[si, ]
        parent$weight <- weight_at_spawn[i]  # restore pre-cost weight for logging
        n_off   <- n_eggs[i]
        new_pids <- seq(next_pid, next_pid + n_off - 1L)

        # Initial offspring weights drawn from same distribution as founders
        off_wt <- pmax(min_wt, rnorm(n_off, start_wt, wt_sd))

        # Offspring inherit strategy and patch from parent; movement threshold drawn fresh
        off_thresh <- if (parent$strategy == "optimal_mover") {
          rnorm(n_off, mean = 0, sd = sigma_bold)
        } else {
          rep(Inf, n_off)
        }

        new_fish_list[[i]] <- tibble(
          strategy          = parent$strategy,
          patch             = parent$patch,
          weight            = off_wt,
          pid               = new_pids,
          cmax_allometric   = NA_real_,
          pcmax_baseline    = NA_real_,
          pcmax_adjusted    = NA_real_,
          pcmax_adjusted_dd = NA_real_,
          func_temp         = NA_real_,
          ration            = NA_real_,
          move_threshold    = off_thresh,
          peak_weight       = off_wt,
          age_days          = 0L,
          spawned_this_year = FALSE,
          parent_pid        = parent$pid,
          cohort            = current_year,
          prob_surv         = 1,
          survive           = 1
        )

        # Log the spawning event (weight and condition are pre-spawn values)
        spawn_log <- bind_rows(spawn_log, tibble(
          parent_pid          = parent$pid,
          dayofsim            = d,
          n_offspring         = n_off,
          weight              = parent$weight,
          condition           = parent$weight / parent$peak_weight,
          offspring_pid_start = next_pid,
          offspring_pid_end   = next_pid + n_off - 1L
        ))

        next_pid <- next_pid + n_off
      }

      # Bind offspring and add them to the population BEFORE survival.
      # Size-based and density-driven starvation mortality will thin the cohort naturally.
      new_fish <- bind_rows(new_fish_list)
      n_new    <- nrow(new_fish)
      switches <- c(switches, integer(n_new))  # extend switch counter for new fish

      # Extend growth_df with stub rows for offspring so survival indexing stays aligned.
      # WT.actual is set to their patch temperature; growth is NA (they didn't feed today).
      off_wt_actual <- ifelse(new_fish$patch == "warm", T_warm, T_cold)
      growth_df <- bind_rows(growth_df,
                             data.frame(WT.actual = off_wt_actual, growth = NA_real_))

      # Register and add offspring to live population
      fish_registry <- bind_rows(
        fish_registry,
        new_fish |>
          select(pid, strategy, parent_pid, cohort) |>
          mutate(birth_dayofsim = d)
      )
      fish_pop <- bind_rows(fish_pop, new_fish)
    }

    # 6. SURVIVAL ─────────────────────────────────────────────────────────────
    s_max_vec <- ifelse(fish_pop$patch == "warm", S_max_warm, S_max_cold)

    # Size-dependent survival: logistic sigmoid in weight, rising from minprob floor to habitat-specific maxprob
    p_sg      <- fncSurviveSize(fish_pop$weight, minprob = s_min, maxprob = s_max_vec, w0 = s_w0, k = s_k)[[1]]

    # Temperature-dependent survival: multiplied with size/growth probability so both act simultaneously
    p_temp    <- fncSurviveTemp(growth_df$WT.actual, T1 = T1_mort, T9 = T9_mort)

    # Condition-based (starvation) survival: relative condition = current / peak weight
    condition <- fish_pop$weight / fish_pop$peak_weight
    p_starv   <- fncSurviveStarve(condition, K9 = K9_starv, K1 = K1_starv)

    # Age-based survival (senescence): survival probability declines past 5 years old,
    # fish older than 8 have ~0 chance of survival.
    # p_age <- fncSurviveAge(fish_pop$age_days / 365)

    # Critical-period consumption-based survival (Elliott 1989) — disabled.
    # Replaced by the sigmoid shape of fncSurviveSize, which sustains strong size-dependent
    # mortality below w0 (~7g) and achieves the same growth-mediated density dependence
    # without requiring a separate consumption threshold module.
    # p_crit <- fncSurviveConsumption(
    #   pcmax_dd         = fish_pop$pcmax_adjusted_dd,
    #   age_days         = fish_pop$age_days,
    #   crit_pcmax_lo    = crit_pcmax_lo,
    #   crit_pcmax_hi    = crit_pcmax_hi,
    #   crit_period_days = crit_period_days
    # )

    # Calculate combined daily survival rate
    prb.srv   <- pmin(p_sg * p_temp * p_starv, 1)

    # Minimum weight floor: fish that drop below hatch weight die immediately (backstop)
    prb.srv[fish_pop$weight < start_wt] <- 0

    survivors <- rbinom(nrow(growth_df), size = 1, prob = prb.srv)
    growth_df <- growth_df |>
      mutate(prob_surv = prb.srv, survive = survivors)

    # 7. STORE RESULTS AND REMOVE NON-SURVIVORS ────────────────────────────────
    ibm_records[[d]] <- data.frame(
      pid       = fish_pop$pid,
      dayofsim  = d,
      strategy  = fish_pop$strategy,
      weight    = fish_pop$weight,
      patch     = fish_pop$patch,
      temp      = growth_df$WT.actual,   # includes offspring stub (their patch temp)
      ggd       = growth_df$growth,      # NA for offspring on birth day
      survived  = growth_df$survive,
      condition = condition,
      age       = fish_pop$age_days / 365,
      p_survive = prb.srv                # worth keeping — used in later analyses
    )

    # Remove non-survivors — only living fish carry forward to next iteration
    fish_pop <- fish_pop[growth_df$survive == 1, ]

    # exit loop if all fish have died
    if (nrow(fish_pop) == 0) break
  }

  # ── Collate ibm_long ─────────────────────────────────────────────────────────
  ibm_long <- bind_rows(ibm_records) |>
    left_join(fish_registry |> select(pid, parent_pid, cohort, birth_dayofsim), by = "pid") |>
    left_join(habitat_df |> select(dayofsim, date), by = "dayofsim")

  # ── Standard daily summary by strategy ───────────────────────────────────────
  ibm_summary <- ibm_long |>
    group_by(strategy, dayofsim, date) |>
    summarise(
      mean_weight    = mean(weight,    na.rm = TRUE),
      sd_weight      = sd(weight,      na.rm = TRUE),
      mean_temp      = mean(temp,      na.rm = TRUE),
      sd_temp        = sd(temp,        na.rm = TRUE),
      mean_ggd       = mean(ggd,       na.rm = TRUE),
      sd_ggd         = sd(ggd,         na.rm = TRUE),
      prop_warm      = mean(patch == "warm", na.rm = TRUE),
      n_alive        = sum(survived,   na.rm = TRUE),
      mean_condition = mean(condition, na.rm = TRUE),
      sd_condition   = sd(condition,   na.rm = TRUE),
      .groups        = "drop"
    )

  list(
    ibm_records   = ibm_records,
    ibm_long      = ibm_long,
    ibm_summary   = ibm_summary,
    fish_registry = fish_registry,
    spawn_log     = spawn_log,
    switches      = switches,
    params        = params,
    habitat_df    = habitat_df
  )
}


