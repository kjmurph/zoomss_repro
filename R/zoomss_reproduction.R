# zoomss_reproduction.R
# Fish reproduction functions for ZooMSS model
# Maintains implicit metabolism approach - energy already net of metabolic costs

#' Calculate maturity ogive for fish groups
#' 
#' @param w Weight vector (g)
#' @param Wmat Log10 weight at 50% maturity
#' @param steepness Steepness of maturity curve (default 1)
#' @return Vector of maturity proportions (0-1)
#' @export
calculate_maturity_ogive <- function(w, Wmat, steepness = 1) {
  # Logistic maturity ogive
  # 50% mature at Wmat
  ogive <- 1 / (1 + exp(-steepness * (log10(w) - Wmat)))
  return(ogive)
}

#' Calculate size-dependent reproductive investment
#' 
#' @param w Weight vector (g)
#' @param ReproInvest Base reproductive investment at maturity
#' @param ReproExp Scaling exponent (typically negative)
#' @param Wmat Log10 weight at maturity
#' @return Vector of reproductive investment fractions
#' @export
calculate_reproductive_investment <- function(w, ReproInvest, ReproExp, Wmat) {
  # Handle missing values
  if(is.na(ReproInvest) || is.na(ReproExp)) {
    return(rep(0, length(w)))
  }
  
  # Only invest in reproduction above maturity size
  w_mat <- 10^Wmat
  
  # Size-dependent reproductive investment
  # Following Gunderson 1997, investment decreases with size
  repro_invest <- rep(0, length(w))
  mature_idx <- which(w >= w_mat)
  
  if(length(mature_idx) > 0) {
    # Scale investment relative to maturity size
    repro_invest[mature_idx] <- ReproInvest * (w[mature_idx]/w_mat)^ReproExp
    
    # Cap at maximum reasonable investment (e.g., 30% of net production)
    repro_invest[repro_invest > 0.3] <- 0.3
  }
  
  return(repro_invest)
}

#' Calculate spawning stock biomass (SSB) from reproductive output
#' 
#' @param N Abundance matrix (groups x size classes)
#' @param repro_energy Reproductive energy matrix (groups x size classes)
#' @param w Weight vector (g)
#' @param fish_grps Indices of fish groups
#' @return Vector of SSB for each fish group
#' @export
calculate_spawning_stock_biomass <- function(N, repro_energy, w, fish_grps) {
  SSB <- rep(0, length(fish_grps))
  
  for(i in seq_along(fish_grps)) {
    grp <- fish_grps[i]
    # SSB is total reproductive output weighted by abundance and body mass
    # Note: repro_energy already accounts for energy allocation
    SSB[i] <- sum(repro_energy[grp,] * N[grp,] * w, na.rm = TRUE)
  }
  
  return(SSB)
}

#' Calculate fish recruitment using Beverton-Holt function
#' 
#' @param SSB Current spawning stock biomass
#' @param params List with recruitment parameters (R0, h, SSB0)
#' @param env_factor Environmental scaling factor (default 1)
#' @return Recruitment abundance
#' @export
calculate_beverton_holt_recruitment <- function(SSB, params, env_factor = 1) {
  # Beverton-Holt stock-recruitment relationship
  # R0: virgin recruitment
  # h: steepness (proportion of R0 when SSB = 0.2*SSB0)
  # SSB0: virgin spawning stock biomass
  
  R0 <- params$R0
  h <- params$h
  SSB0 <- params$SSB0
  
  # Prevent division by zero
  if(SSB0 == 0) SSB0 <- 1e-10
  
  # Calculate Beverton-Holt parameters
  alpha <- 4 * h * R0 / (SSB0 * (1 - h))
  beta <- (5 * h - 1) / (SSB0 * (1 - h))
  
  # Calculate recruitment
  recruitment <- alpha * SSB / (1 + beta * SSB)
  
  # Apply environmental factor
  recruitment <- recruitment * env_factor
  
  # Ensure non-negative
  recruitment[recruitment < 0] <- 0
  
  return(recruitment)
}

#' Calculate fish recruitment using Ricker function
#' 
#' @param SSB Current spawning stock biomass
#' @param params List with recruitment parameters (alpha, beta)
#' @param env_factor Environmental scaling factor (default 1)
#' @return Recruitment abundance
#' @export
calculate_ricker_recruitment <- function(SSB, params, env_factor = 1) {
  # Ricker stock-recruitment relationship
  # Allows for overcompensation at high SSB
  
  alpha <- params$alpha
  beta <- params$beta
  
  # Calculate recruitment
  recruitment <- alpha * SSB * exp(-beta * SSB)
  
  # Apply environmental factor
  recruitment <- recruitment * env_factor
  
  # Ensure non-negative
  recruitment[recruitment < 0] <- 0
  
  return(recruitment)
}

#' Update fish boundary conditions with dynamic recruitment
#' 
#' @param N Current abundance matrix
#' @param SSB Vector of spawning stock biomass for fish groups
#' @param Groups Groups dataframe with parameters
#' @param fish_grps Indices of fish groups
#' @param w Weight vector
#' @param recruit_params List of recruitment parameters for each fish group
#' @param method Recruitment function ("beverton_holt" or "ricker")
#' @param env_factor Optional environmental scaling factor
#' @return Updated abundance matrix with new recruitment
#' @export
update_fish_recruitment <- function(N, SSB, Groups, fish_grps, w, 
                                   recruit_params, method = "beverton_holt",
                                   env_factor = 1) {
  
  for(i in seq_along(fish_grps)) {
    grp <- fish_grps[i]
    
    # Find index of egg size (W0)
    egg_idx <- which.min(abs(log10(w) - Groups$W0[grp]))
    
    # Calculate recruitment based on method
    if(method == "beverton_holt") {
      recruitment <- calculate_beverton_holt_recruitment(SSB[i], recruit_params[[i]], env_factor)
    } else if(method == "ricker") {
      recruitment <- calculate_ricker_recruitment(SSB[i], recruit_params[[i]], env_factor)
    } else {
      stop("Unknown recruitment method: ", method)
    }
    
    # Ensure recruitment is not infinite or NA
    if(!is.finite(recruitment)) recruitment <- 0
    
    # Update abundance at smallest size class
    N[grp, egg_idx] <- recruitment
  }
  
  return(N)
}

#' Initialize recruitment parameters for fish groups
#' 
#' @param Groups Groups dataframe
#' @param fish_grps Indices of fish groups
#' @param N Initial abundance matrix
#' @param w Weight vector
#' @param h Steepness parameter (default 0.75)
#' @return List of recruitment parameters for each fish group
#' @export
initialize_recruitment_params <- function(Groups, fish_grps, N, w, h = 0.75) {
  recruit_params <- list()
  
  for(i in seq_along(fish_grps)) {
    grp <- fish_grps[i]
    
    # Get initial recruitment (current abundance at W0)
    egg_idx <- which.min(abs(log10(w) - Groups$W0[grp]))
    R0 <- N[grp, egg_idx]
    
    # Ensure R0 is positive
    if(R0 <= 0) R0 <- 1
    
    # Estimate virgin SSB (assuming current state is ~60% of virgin)
    # This is a simplification - in reality would need spin-up
    mature_idx <- which(log10(w) >= Groups$Wmat[grp])
    
    if(length(mature_idx) > 0) {
      # Assume 10% of mature biomass goes to reproduction initially
      current_SSB <- sum(N[grp, mature_idx] * w[mature_idx]) * 0.1
      SSB0 <- current_SSB / 0.6  # Assume current is 60% of virgin
    } else {
      # Fallback if no mature individuals
      SSB0 <- R0 * 100  # Rough approximation
    }
    
    # Ensure SSB0 is positive
    if(SSB0 <= 0) SSB0 <- 1
    
    recruit_params[[i]] <- list(
      R0 = R0,
      h = h,  # Steepness (typically 0.6-0.9 for fish)
      SSB0 = SSB0,
      # For Ricker alternative
      alpha = R0 / SSB0 * exp(1),  # Maximum recruitment at SSB0/e
      beta = 1 / SSB0
    )
  }
  
  return(recruit_params)
}

#' Diagnostic function to check reproduction is working
#' 
#' @param model Model output object
#' @param fish_grps Fish group indices  
#' @return List of diagnostic metrics
#' @export
check_reproduction_diagnostics <- function(model, fish_grps) {
  
  diagnostics <- list()
  
  # Check if reproduction arrays exist and are properly structured
  if(!is.null(model$repro_energy)) {
    # Check if this is a 3D array (time x groups x size) or 2D array (groups x size)
    if(length(dim(model$repro_energy)) == 3) {
      # 3D array: time x groups x size
      final_time <- dim(model$repro_energy)[1]  # Last time step
      diagnostics$repro_energy_total <- rowSums(model$repro_energy[final_time, fish_grps, , drop = FALSE])
      diagnostics$mean_repro_investment <- mean(model$repro_energy[final_time, fish_grps, ] / 
                                                (model$gg[final_time, fish_grps, ] + model$repro_energy[final_time, fish_grps, ] + 1e-10))
    } else {
      # 2D array: groups x size  
      diagnostics$repro_energy_total <- rowSums(model$repro_energy[fish_grps, , drop = FALSE])
      diagnostics$mean_repro_investment <- mean(model$repro_energy[fish_grps, ] / 
                                                (model$gg[fish_grps, ] + model$repro_energy[fish_grps, ] + 1e-10))
    }
  }
  
  # Check SSB time series
  if(!is.null(model$SSB)) {
    if(length(dim(model$SSB)) == 2) {
      # SSB is fish_groups x time
      diagnostics$SSB_mean <- rowMeans(model$SSB, na.rm = TRUE)
      diagnostics$SSB_cv <- apply(model$SSB, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
    } else {
      # SSB is a vector
      diagnostics$SSB_mean <- mean(model$SSB, na.rm = TRUE)
      diagnostics$SSB_cv <- sd(model$SSB, na.rm = TRUE)/mean(model$SSB, na.rm = TRUE)
    }
  }
  
  # Check recruitment time series
  if(!is.null(model$N)) {
    W0_indices <- model$param$Groups$W0[fish_grps]
    
    # N is time x groups x size
    n_time <- dim(model$N)[1]
    recruitment <- matrix(NA, length(fish_grps), n_time)
    
    for(i in seq_along(fish_grps)) {
      grp <- fish_grps[i]
      egg_idx <- which.min(abs(log10(model$param$w) - W0_indices[i]))
      recruitment[i, ] <- model$N[, grp, egg_idx]
    }
    
    diagnostics$recruitment_mean <- rowMeans(recruitment, na.rm = TRUE)
    diagnostics$recruitment_cv <- apply(recruitment, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
  }
  
  return(diagnostics)
}

