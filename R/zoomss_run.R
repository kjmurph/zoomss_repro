#' Run ZooMSS Model Forward in Time
#'
#' @title Execute the main ZooMSS simulation loop with dynamic environmental forcing
#' @description Runs the ZooMSS model forward in time, updating environmental conditions
#'   and population dynamics at each time step using the McKendrick-von Foerster framework.
#' @details This is the core simulation engine of ZooMSS that:
#'
#'   **Environmental Dynamics:**
#'   - Updates phytoplankton abundance spectrum based on chlorophyll time series
#'   - Applies temperature effects on zooplankton and fish metabolism
#'   - Recalculates feeding kernels with current environmental conditions
#'
#'   **Population Dynamics:**
#'   - Solves McKendrick-von Foerster equation for size-structured growth
#'   - Updates feeding interactions between all size classes and groups
#'   - Calculates mortality from predation, senescence, and fishing
#'   - Handles recruitment and boundary conditions for each functional group
#'
#'   **Time Integration:**
#'   - Processes model through all time steps with adaptive environmental forcing
#'   - Saves output at specified intervals for memory efficiency
#'   - Maintains mass balance and numerical stability throughout simulation
#'
#'   Unlike static models, this version dynamically updates phytoplankton spectra
#'   and temperature effects at each time step based on provided environmental data.
#'
#' @param model Model object created by zoomss_setup containing:
#'   - param: Complete parameter list with environmental time series
#'   - Feeding kernels and biological rate parameters
#'   - Initial conditions and model structure
#'
#' @return List containing complete model output:
#'   \itemize{
#'     \item param: Model parameters used in simulation
#'     \item N: Abundance time series (time x groups x size classes)
#'     \item gg: Growth rate time series
#'     \item diet: Diet composition time series
#'     \item Z: Mortality rate time series
#'     \item time: Time values corresponding to saved results (accounting for isave)
#'     \item w: Size class weights (g)
#'     \item Additional time series data and model results
#'   }
#'
#' @examples
#' \dontrun{
#' # Set up model parameters and structure
#' params <- zoomss_params(Groups, input_params)
#' model <- zoomss_setup(params)
#'
#' # Run the simulation
#' results <- zoomss_run(model)
#'
#' # Access final abundances
#' final_abundances <- results$N[dim(results$N)[1],,]
#' }
#'
#' @noRd
#'
zoomss_run <- function(model){

  # Pull out some useful parameters - just a shortcut
  param <- model$param
  dt <- model$param$dt[1]  # Use scalar value instead of vector
  dx <- model$param$dx     # This should already be scalar
  ngrps <- model$param$ngrps
  ngrid <- model$param$ngrid
  w <- model$param$w
  W0 <- param$Groups$W0
  Wmax <- param$Groups$Wmax
  fish_grps <- model$param$fish_grps
  zoo_grps <- model$param$zoo_grps
  num_fish <- param$num_fish

  # Extract effort-driven fishing parameters
  effort_fishing <- param$effort_fishing
  if (effort_fishing) {
    effort <- param$effort          # ntime x num_fish matrix
    q_fish <- param$q               # catchability per fish group (length num_fish)
    selectivity <- model$selectivity # ngrps x ngrid matrix
    fish_mort_dynamic <- matrix(0, nrow = ngrps, ncol = ngrid)  # working matrix
  }

  # Extract dual-pathway parameters
  assim_eff_gge <- model$assim_eff_gge  # GGE*Carbon matrix for zoo pathway (ngrps x ngrid)
  assim_by_prey <- model$assim_by_prey  # (1 - def) by pred-prey pair, fish only
  K_growth <- model$K_growth            # Growth fraction of assimilated (fish only)
  R_frac <- model$R_frac                # Reproduction fraction (0 for zoo, >0 for fish)
  mat_ogive <- model$mat_ogive          # Maturity ogive (ngrps x ngrid)
  repro_on <- model$repro_on            # Reproduction enabled flag
  repro_eff <- model$repro_eff          # Reproductive efficiency

  dynam_dietkernel <- model$dynam_dietkernel
  dynam_growthkernel <- model$dynam_growthkernel
  dynam_mortkernel <- model$dynam_mortkernel
  dynam_diffkernel <- model$dynam_diffkernel

  # Get pre-calculated phytoplankton and temperature time series from param
  phyto_int <- param$phyto_int
  phyto_slope <- param$phyto_slope
  temp_eff_zoo <- param$temp_eff_zoo
  temp_eff_fish <- param$temp_eff_fish

  curr_min_size <- vector()
  curr_max_size <- vector()
  for (i in 1:ngrps){
    curr_min_size[i] <- which(round(log10(w), digits = 2) == W0[i])
    curr_max_size[i] <- which(round(log10(w), digits = 2) == Wmax[i])
  }

  # Store min/max size indices for reproduction
  min_size_idx <- param$min_size_idx
  max_size_idx <- param$max_size_idx

  idx_iter <- 2:ngrid
  idx <- 2:(ngrid-1)
  itimemax  <- param$itimemax  # Number of time points to simulate

  if(length(zoo_grps) > 1){ # If there's only one zoo group, then you do not need w0idx
    w0idx <- which(W0 > min(W0) & is.na(param$Groups$Prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    props_z <- param$Groups$Prop[w0idx] # Zooplankton proportions

    for(i in seq_along(w0idx)){
      w0mins[i] <- which(round(log10(w), digits = 2) == W0[w0idx[i]])
    }
  }

  # Matrices for MvF and MvF-D numeric solution
  A_iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  C_iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  S_iter <- matrix(0, nrow = ngrps, ncol = ngrid)

  A <- matrix(0, nrow = ngrps, ncol = ngrid)
  B <- matrix(0, nrow = ngrps, ncol = ngrid)
  C <- matrix(0, nrow = ngrps, ncol = ngrid)
  S <- matrix(0, nrow = ngrps, ncol = ngrid)

  # Temporary matrices for energy budget calculations
  repro_rate <- matrix(0, nrow = ngrps, ncol = ngrid)  # Reproductive rate per group/size

  # Temporary Matrices that get updated each time step
  N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid)

  pb <- progress::progress_bar$new(
    format = "ZooMSS Time Loop [:bar] :percent eta: :eta elapsed: :elapsed",
    total = itimemax,
    width = 60,
    show_after = 0
  )

  # BIG TIME LOOP
  for (itime in 1:itimemax){

    pb$tick() # Update progress bar

    # DYNAMIC ENVIRONMENTAL FORCING - Update for current time step
    current_phyto_int <- phyto_int[itime]
    current_phyto_slope <- phyto_slope[itime]

    # Update phytoplankton abundance spectrum with current environmental parameters
    model$nPP <- 10^(current_phyto_int)*(param$w_phyto^(current_phyto_slope))

    # Update temperature effects matrix for all groups based on current SST
    if(param$num_zoo > 0) {
      for(i in 1:param$num_zoo) {
        zoo_grp_idx <- zoo_grps[i]
        model$temp_eff[zoo_grp_idx, ] <- temp_eff_zoo[itime, i]
      }
    }

    if(num_fish > 0) {
      for(i in 1:num_fish) {
        fish_grp_idx <- fish_grps[i]
        model$temp_eff[fish_grp_idx, ] <- temp_eff_fish[itime, i]
      }
    }

    # ==========================================================================
    # PHYTOPLANKTON FEEDING
    # ==========================================================================
    # Phyto kernels are pre-multiplied in setup with per-group assimilation:
    #   Zoo groups: GrossGEscale * cc_phyto (final GGE, no further partitioning)
    #   Fish groups: (1 - def_phyto) (pure assimilation, K_growth applied below)
    current_ingested_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2))
    current_diff_phyto <- model$temp_eff^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2))

    # Update senescence mortality with current temperature effects
    model$M_sb <- model$M_sb_base * model$temp_eff

    # ==========================================================================
    # GROWTH FROM DYNAMIC SPECTRUM — DUAL PATHWAY
    # ==========================================================================

    # --- ZOOPLANKTON PATHWAY: Original GGE approach (vectorised) ---
    # growth_multiplier[w] = sum_j N[j,w] * GrossGEscale[j] * Carbon[j]
    # This is the prey-weighted growth efficiency summed over all prey groups.
    growth_multiplier_gge <- colSums(N * assim_eff_gge)  # 1 x ngrid

    # Apply temperature effects to growth kernel
    temp_growth_kernel_3d <- sweep(dynam_growthkernel, c(1,2), model$temp_eff, '*')

    # Vectorised kernel multiplication (computes growth for ALL predators using GGE)
    temp_gk_2d <- temp_growth_kernel_3d
    dim(temp_gk_2d) <- c(ngrps * ngrid, ngrid)
    cs <- .colSums(growth_multiplier_gge * t(temp_gk_2d), m = ngrid, n = ngrps * ngrid)
    dim(cs) <- c(ngrps, ngrid)

    gg_dynam <- cs  # Growth from dynamic spectrum (GGE-based for all rows)

    # --- FISH PATHWAY: Explicit energy budget (overwrite fish rows) ---
    # For fish predators, replace with prey-specific defecation calculation.
    # assim_by_prey[pred, prey] = (1 - def[pred, prey]) — pure assimilation.
    if (num_fish > 0) {
      for (f in 1:num_fish) {
        fg <- fish_grps[f]

        # Compute prey-weighted assimilation for this fish predator
        fish_growth_mult <- numeric(ngrid)
        for (prey in 1:ngrps) {
          fish_growth_mult <- fish_growth_mult + N[prey, ] * assim_by_prey[fg, prey]
        }

        # Matrix-vector product: kernel[pred_size, prey_size] %*% growth_mult[prey_size]
        gg_dynam[fg, ] <- as.vector(temp_growth_kernel_3d[fg, , ] %*% fish_growth_mult)
      }
    }

    # ==========================================================================
    # TOTAL GROWTH — COMBINE PATHWAYS
    # ==========================================================================
    # For zooplankton: phyto (GGE-based) + dynamic (GGE-based) = FINAL growth rate
    # For fish: phyto ((1-def)-based) + dynamic ((1-def)-based) = total ASSIMILATED energy
    gg <- current_ingested_phyto + gg_dynam

    # ==========================================================================
    # FISH ENERGY PARTITIONING
    # ==========================================================================
    # For fish only: partition assimilated energy into K_growth and R_frac.
    # Zooplankton gg values are ALREADY the final growth rate — no partitioning.

    gg_total <- gg  # Store total assimilated (needed for fish reproduction)

    # Reproduction rate (fish only)
    repro_rate <- matrix(0, nrow = ngrps, ncol = ngrid)

    if (num_fish > 0) {
      for (f in 1:num_fish) {
        fg <- fish_grps[f]

        if (repro_on[fg] == 1 && R_frac[fg] > 0) {
          # Reproductive rate: R_frac portion of assimilated energy, scaled by maturity
          repro_rate[fg, ] <- R_frac[fg] * gg_total[fg, ] * mat_ogive[fg, ]
        }

        # Maturity-dependent growth for fish:
        # Below Wmat (mat_ogive ~ 0): R_frac redirected to somatic growth
        # Above Wmat (mat_ogive ~ 1): R_frac goes to reproduction
        # gg_fish = K_growth * gg_total + R_frac * gg_total * (1 - mat_ogive)
        gg[fg, ] <- K_growth[fg] * gg_total[fg, ] +
          R_frac[fg] * gg_total[fg, ] * (1 - mat_ogive[fg, ])
      }
    }
    # After this block:
    #   gg[zoo_grps, ] = GGE-based growth rate (unchanged, final)
    #   gg[fish_grps, ] = energy-budget-partitioned growth rate

    # ==========================================================================
    # MORTALITY
    # ==========================================================================
    predation_multiplier <- N * model$temp_eff

    sw2 <- sweep(dynam_mortkernel, c(2,3), predation_multiplier, '*')
    ap2 <- aperm(sw2, c(2,3,1))
    M2 <- .colSums(colSums(ap2), ngrid, ngrid)

    # Fishing mortality: effort-driven (dynamic) or static (fallback)
    if (effort_fishing) {
      fish_mort_dynamic[] <- 0
      for (f in 1:num_fish) {
        fg <- fish_grps[f]
        # F(w,t) = Effort(t) * q * Selectivity(w)
        fish_mort_dynamic[fg, ] <- effort[itime, f] * q_fish[f] * selectivity[fg, ]
      }
      current_fish_mort <- fish_mort_dynamic
    } else {
      current_fish_mort <- model$fish_mort
    }

    Z <- sweep(model$M_sb + current_fish_mort, 2, M2, '+')
    rm(sw2, ap2)

    # ==========================================================================
    # DIFFUSION — DUAL PATHWAY
    # ==========================================================================

    # --- ZOOPLANKTON PATHWAY: Original GGE approach (vectorised) ---
    diffusion_multiplier_gge <- colSums(N * (assim_eff_gge^2))  # 1 x ngrid

    temp_diff_kernel_3d <- sweep(dynam_diffkernel, c(1,2), model$temp_eff^2, '*')

    temp_dk_2d <- temp_diff_kernel_3d
    dim(temp_dk_2d) <- c(ngrps * ngrid, ngrid)
    cs_diff <- .colSums(diffusion_multiplier_gge * t(temp_dk_2d), m = ngrid, n = ngrps * ngrid)
    dim(cs_diff) <- c(ngrps, ngrid)

    diff_dynam <- cs_diff  # Diffusion from dynamic spectrum (GGE-based for all rows)

    # --- FISH PATHWAY: Overwrite fish rows ---
    if (num_fish > 0) {
      for (f in 1:num_fish) {
        fg <- fish_grps[f]

        fish_diff_mult <- numeric(ngrid)
        for (prey in 1:ngrps) {
          fish_diff_mult <- fish_diff_mult + N[prey, ] * (assim_by_prey[fg, prey]^2)
        }

        diff_dynam[fg, ] <- as.vector(temp_diff_kernel_3d[fg, , ] %*% fish_diff_mult)
      }
    }

    diff <- current_diff_phyto + diff_dynam

    # ==========================================================================
    # FISH DIFFUSION PARTITIONING
    # ==========================================================================
    # For zooplankton: diffusion is already based on GGE^2 — no further scaling.
    # For fish: diffusion must scale with effective growth fraction squared.
    #   eff_K = K_growth + R_frac * (1 - mat_ogive)
    #   (immature fish redirect R_frac to growth, so eff_K > K_growth)

    if (num_fish > 0) {
      for (f in 1:num_fish) {
        fg <- fish_grps[f]
        eff_K <- K_growth[fg] + R_frac[fg] * (1 - mat_ogive[fg, ])
        diff[fg, ] <- diff[fg, ] * (eff_K^2)
      }
    }

    # ==========================================================================
    # McKendrick-von Foerster NUMERICAL SOLUTION
    # ==========================================================================
    A_iter[,idx_iter] <- dt/dx * gg[,idx_iter-1]
    C_iter[,idx_iter] <- 1 + dt * Z[,idx_iter] + dt/dx * gg[,idx_iter]
    S_iter[,idx_iter] <- N[,idx_iter]
    N_iter <- N
    N_iter[1,1] <- N[1,1]

    A[,idx] <- dt/dx * (gg[,idx-1] + diff[,idx-1] * (log(10)/2+1/(2*dx)))
    B[,idx] <- diff[,idx+1] * dt/(2*dx^2)
    C[,idx] <- 1 + dt * Z[,idx] + dt/dx*(gg[,idx] + diff[,idx] * (log(10)/2+1/dx))
    S[,idx] <- N[,idx]

    # Solve MvF equation
    N <- zoomss_mvf(ngrps, curr_min_size, curr_max_size,
                    A_iter, C_iter, N_iter, S_iter,
                    A, B, C, N, S)

    # ==========================================================================
    # BOUNDARY CONDITIONS
    # ==========================================================================
    # Zooplankton: maintain current proportional closure
    if(length(zoo_grps) > 1){
      for(i in seq_along(w0idx)){
        w_min_curr <- w0mins[i]
        exclude_mins <- w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] <- props_z[i] * sum(N[-exclude_mins, w_min_curr])
      }
    }

    # ==========================================================================
    # FISH RECRUITMENT FROM REPRODUCTION
    # ==========================================================================
    # For fish with repro_on = 1: calculate recruitment from reproductive investment
    # For fish with repro_on = 0: use original boundary condition

    fish_mins <- unlist(lapply(W0[fish_grps],
                               function(x){which(round(log10(w), digits = 2) == x)}))

    for (f in 1:num_fish) {
      fg <- fish_grps[f]

      if (repro_on[fg] == 1) {
        # Calculate total reproductive output (biomass-weighted reproductive rate)
        # R_total = sum over mature sizes: repro_rate * N * w * dx
        min_idx <- min_size_idx[fg]
        max_idx <- max_size_idx[fg]

        R_total <- sum(repro_rate[fg, min_idx:max_idx] *
                         N[fg, min_idx:max_idx] *
                         w[min_idx:max_idx]) * dx

        # Recruitment flux (following DBPM equation)
        # recruits = (R_total * repro_eff) / (w_min * dx)
        recruitment_flux <- (R_total * repro_eff[fg]) / (w[min_idx] * dx)

        # Semi-implicit boundary condition for recruitment
        # Balances recruitment input against mortality (Z) and growth-out (gg/dx)
        # at the minimum size class, consistent with the interior MvF scheme.
        # Without the denominator, N[min_idx] accumulates without bound because
        # the MvF solver (which starts at min_idx+1) never applies mortality or
        # advection at the boundary.
        N[fg, min_idx] <- (N[fg, min_idx] + recruitment_flux * dt) / (1 + dt * Z[fg, min_idx] + (dt / dx) * gg[fg, min_idx])

      } else {
        # Original boundary condition for fish without reproduction
        if(length(fish_grps) > 1 && length(zoo_grps) > 1){
          N[fg, fish_mins[f]] <- (1/num_fish) * sum(N[-fish_grps, fish_mins[f]])
        } else {
          N[fg, fish_mins[f]] <- (1/num_fish) * sum(N[-fish_grps, fish_mins[f]])
        }
      }
    }

    # ==========================================================================
    # SAVE OUTPUT
    # ==========================================================================
    save_this_step <- (itime %% param$isave) == 0

    if(save_this_step){
      isav <- itime/param$isave

      ## Phytoplankton diet
      current_pico_diet <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$param$w_phyto) < -11.5), "*"), dims = 2))
      current_nano_diet <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$param$w_phyto) >= -11.5 & log10(model$param$w_phyto) < -8.5), "*"), dims = 2))
      current_micro_diet <- model$temp_eff*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$param$w_phyto) >= -8.5), "*"), dims = 2))

      pico_phyto_diet <- rowSums(current_pico_diet*N)
      nano_phyto_diet <- rowSums(current_nano_diet*N)
      micro_phyto_diet <- rowSums(current_micro_diet*N)

      ## Functional group diet
      N_array_temp <- aperm(replicate(ngrid, N), c(3,1,2))
      N_array <- aperm(replicate(ngrps, N_array_temp), c(4,1,2,3))
      temp_dynam_dietkernel <- sweep(dynam_dietkernel, c(1,2), model$temp_eff, "*")
      dynam_diet <- rowSums(aperm(rowSums(sweep(temp_dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)

      # Save standard outputs
      model$time[isav] <- param$time[itime]
      model$N[isav,,] <- N
      model$Z[isav,,] <- Z
      model$gg[isav,,] <- gg
      model$gg_total[isav,,] <- gg_total  # Total assimilated energy before partitioning
      model$diet[isav,,1:3] <- cbind(pico_phyto_diet, nano_phyto_diet, micro_phyto_diet)
      model$diet[isav,,c(4:(dim(param$Groups)[1]+3))] <- dynam_diet

      # Save reproduction outputs
      model$repro_rate[isav,,] <- repro_rate

      # Calculate and save fish-specific reproduction metrics
      for (f in 1:num_fish) {
        fg <- fish_grps[f]
        min_idx <- min_size_idx[fg]
        max_idx <- max_size_idx[fg]

        # SSB: Spawning Stock Biomass (biomass of mature individuals)
        model$SSB[isav, f] <- sum(mat_ogive[fg, min_idx:max_idx] *
                                    N[fg, min_idx:max_idx] *
                                    w[min_idx:max_idx]) * dx

        # Total reproductive output
        R_total <- sum(repro_rate[fg, min_idx:max_idx] *
                         N[fg, min_idx:max_idx] *
                         w[min_idx:max_idx]) * dx
        model$total_repro_output[isav, f] <- R_total

        # Recruitment
        model$recruitment[isav, f] <- (R_total * repro_eff[fg]) / (w[min_idx] * dx)
      }

      # Save effort-driven fishing outputs
      if (effort_fishing) {
        model$Fmort_ts[isav,,] <- current_fish_mort
        for (f in 1:num_fish) {
          fg <- fish_grps[f]
          min_idx <- min_size_idx[fg]
          max_idx <- max_size_idx[fg]
          # Catch = integral of F(w) * N(w) * w dw over fished size range
          model$catch[isav, f] <- sum(current_fish_mort[fg, min_idx:max_idx] *
                                        N[fg, min_idx:max_idx] *
                                        w[min_idx:max_idx]) * dx
        }
      }
    }
  } # End of time loop

  return(model)
}
