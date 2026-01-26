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

  # Extract energy budget parameters
  assim_by_prey <- model$assim_by_prey  # Pre-calculated: (1 - def) * K_growth by predator-prey pair
  K_growth <- model$K_growth            # Growth fraction of assimilated
  R_frac <- model$R_frac                # Reproduction fraction of assimilated
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
    # PHYTOPLANKTON FEEDING (pre-multiplied with assim_phyto in setup)
    # ==========================================================================
    # These kernels already incorporate (1 - def_phyto) * K_growth
    current_ingested_phyto <- model$temp_eff*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2))
    current_diff_phyto <- model$temp_eff^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2))

    # Update senescence mortality with current temperature effects
    model$M_sb <- model$M_sb_base * model$temp_eff

    # ==========================================================================
    # DYNAMIC SPECTRUM FEEDING WITH EXPLICIT ENERGY BUDGET
    # ==========================================================================
    # Calculate growth from dynamic spectrum with prey-specific defecation
    # assim_by_prey[pred, prey] = (1 - def[pred,prey]) * K_growth[pred]

    # Growth multiplier now uses prey-specific assimilation efficiency
    # For each predator, sum over all prey: N[prey] * assim_by_prey[pred, prey]
    # This creates a matrix: rows = predators, cols = prey size classes
    # Then sum over prey groups

    growth_multiplier_by_prey <- matrix(0, nrow = ngrps, ncol = ngrid)
    for (pred in 1:ngrps) {
      # For this predator, calculate growth contribution from each prey group
      for (prey in 1:ngrps) {
        growth_multiplier_by_prey[pred, ] <- growth_multiplier_by_prey[pred, ] +
          N[prey, ] * assim_by_prey[pred, prey]
      }
    }

    # Apply temperature effects and sum using growth kernels
    temp_growth_kernel <- sweep(dynam_growthkernel, c(1,2), model$temp_eff, '*')

    # Calculate growth from dynamic spectrum for each predator
    gg_dynam <- matrix(0, nrow = ngrps, ncol = ngrid)
    for (pred in 1:ngrps) {
      for (prey in 1:ngrps) {
        # Growth from this prey group
        kernel_contribution <- temp_growth_kernel[pred, , ] %*%
          (N[prey, ] * assim_by_prey[pred, prey] * dx)
        gg_dynam[pred, ] <- gg_dynam[pred, ] + kernel_contribution
      }
    }

    # Total growth = phytoplankton + dynamic spectrum
    gg <- current_ingested_phyto + gg_dynam

    # ==========================================================================
    # REPRODUCTION CALCULATION (Fish only)
    # ==========================================================================
    # For fish with repro_on = 1:
    # - Immature individuals (mat_ogive ~ 0): all energy to growth (gg includes R_frac)
    # - Mature individuals (mat_ogive ~ 1): R_frac goes to reproduction, not growth

    # Calculate reproductive investment rate (only for fish with repro_on = 1)
    repro_rate <- matrix(0, nrow = ngrps, ncol = ngrid)

    for (f in 1:num_fish) {
      fg <- fish_grps[f]
      if (repro_on[fg] == 1) {
        # Calculate the additional growth rate that would go to immature individuals
        # This is the R_frac portion of assimilated energy
        # R_frac / K_growth gives the ratio of reproduction to growth allocation

        # For mature fish, R_frac portion goes to reproduction instead of growth
        # repro_rate = (R_frac / K_growth) * gg * mat_ogive
        # This represents the reproductive investment per unit biomass

        if (K_growth[fg] > 0) {
          repro_ratio <- R_frac[fg] / K_growth[fg]
          repro_rate[fg, ] <- repro_ratio * gg[fg, ] * mat_ogive[fg, ]

          # Adjust growth: mature individuals lose R_frac portion
          # gg_adjusted = gg * (1 - mat_ogive * R_frac/(K_growth + R_frac))
          # Simplifies to: gg * (K_growth + R_frac*(1-mat_ogive)) / (K_growth + R_frac)
          total_prod <- K_growth[fg] + R_frac[fg]
          growth_frac <- (K_growth[fg] + R_frac[fg] * (1 - mat_ogive[fg, ])) / total_prod
          gg[fg, ] <- gg[fg, ] * growth_frac
        }
      }
    }

    # ==========================================================================
    # MORTALITY
    # ==========================================================================
    predation_multiplier <- N * model$temp_eff

    sw2 <- sweep(dynam_mortkernel, c(2,3), predation_multiplier, '*')
    ap2 <- aperm(sw2, c(2,3,1))
    M2 <- .colSums(colSums(ap2), ngrid, ngrid)
    Z <- sweep(model$M_sb + model$fish_mort, 2, M2, '+')
    rm(sw2, ap2)

    # ==========================================================================
    # DIFFUSION
    # ==========================================================================
    # Diffusion multiplier with prey-specific assimilation (squared)
    diffusion_multiplier_by_prey <- matrix(0, nrow = ngrps, ncol = ngrid)
    for (pred in 1:ngrps) {
      for (prey in 1:ngrps) {
        diffusion_multiplier_by_prey[pred, ] <- diffusion_multiplier_by_prey[pred, ] +
          N[prey, ] * (assim_by_prey[pred, prey]^2)
      }
    }

    temp_diff_kernel <- sweep(dynam_diffkernel, c(1,2), model$temp_eff^2, '*')

    diff_dynam <- matrix(0, nrow = ngrps, ncol = ngrid)
    for (pred in 1:ngrps) {
      for (prey in 1:ngrps) {
        kernel_contribution <- temp_diff_kernel[pred, , ] %*%
          (N[prey, ] * (assim_by_prey[pred, prey]^2) * dx)
        diff_dynam[pred, ] <- diff_dynam[pred, ] + kernel_contribution
      }
    }

    diff <- current_diff_phyto + diff_dynam

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

        # Add recruitment to smallest size class
        N[fg, min_idx] <- N[fg, min_idx] + recruitment_flux * dt

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
    }
  } # End of time loop

  return(model)
}
