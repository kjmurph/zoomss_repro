test_that("backward compatibility: static fishing still works", {
  skip_if_not_installed("zoomss")
  
  # Create simple environmental data
  n_years <- 5
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  
  # Get default groups
  Groups <- getGroups(source = "default")
  
  # Set some static fishing mortality
  Groups$Fmort[1] <- 0.2
  Groups$Fmort_W0[1] <- Groups$W0[1]
  Groups$Fmort_Wmax[1] <- Groups$Wmax[1]
  
  # Should run without errors (no fishing_params = static fishing)
  expect_no_error({
    results <- zoomss_model(input_params, Groups, isave = 10)
  })
  
  # Check that model ran and produced results
  expect_true(!is.null(results))
  expect_true("abundance" %in% names(results))
})


test_that("dynamic fishing params validation works correctly", {
  # Create simple environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Invalid: missing elements
  fishing_bad1 <- list(effort = rep(100, n_steps))
  expect_error(
    zoomss_model(input_params, Groups, fishing_params = fishing_bad1),
    "must contain 'effort', 'catchability', and 'selectivity'"
  )
  
  # Invalid: effort wrong length
  fishing_bad2 <- list(
    effort = rep(100, 10),  # Wrong length
    catchability = rep(0.01, ngrps),
    selectivity = vector("list", ngrps)
  )
  expect_error(
    zoomss_model(input_params, Groups, fishing_params = fishing_bad2),
    "does not match input_params rows"
  )
  
  # Invalid: catchability wrong length
  fishing_bad3 <- list(
    effort = rep(100, n_steps),
    catchability = rep(0.01, 2),  # Wrong length
    selectivity = vector("list", ngrps)
  )
  expect_error(
    zoomss_model(input_params, Groups, fishing_params = fishing_bad3),
    "does not match number of groups"
  )
})


test_that("dynamic fishing runs successfully with knife-edge selectivity", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 3
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters with knife-edge selectivity
  effort <- seq(0, 200, length.out = n_steps)  # Increasing effort
  catchability <- rep(0.01, ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = (Groups$W0[g] + Groups$Wmax[g]) / 2),  # Midpoint
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Should run without errors
  expect_no_error({
    results <- zoomss_model(input_params, Groups, isave = 6, fishing_params = fishing_params)
  })
  
  # Check that results exist
  expect_true(!is.null(results))
  expect_true("abundance" %in% names(results))
})


test_that("dynamic fishing runs successfully with logistic selectivity", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters with logistic selectivity
  effort <- rep(150, n_steps)
  catchability <- rep(0.015, ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    w_mid <- (Groups$W0[g] + Groups$Wmax[g]) / 2
    selectivity[[g]] <- list(
      type = "logistic",
      params = list(
        L50 = w_mid,
        L95 = w_mid + 0.5
      ),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Should run without errors
  expect_no_error({
    results <- zoomss_model(input_params, Groups, isave = 4, fishing_params = fishing_params)
  })
  
  expect_true(!is.null(results))
})


test_that("dynamic fishing effort time series affects results", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Setup common selectivity
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = Groups$W0[g]),
      log_scale = TRUE
    )
  }
  
  catchability <- rep(0.01, ngrps)
  
  # Run 1: Zero effort (no fishing)
  fishing_zero <- list(
    effort = rep(0, n_steps),
    catchability = catchability,
    selectivity = selectivity
  )
  
  results_zero <- zoomss_model(input_params, Groups, isave = 6, fishing_params = fishing_zero)
  
  # Run 2: High effort
  fishing_high <- list(
    effort = rep(500, n_steps),
    catchability = catchability,
    selectivity = selectivity
  )
  
  results_high <- zoomss_model(input_params, Groups, isave = 6, fishing_params = fishing_high)
  
  # Final abundance should be lower with high fishing effort
  # (assuming fishing reduces abundance)
  final_biomass_zero <- sum(results_zero$biomass[nrow(results_zero$biomass), , ])
  final_biomass_high <- sum(results_high$biomass[nrow(results_high$biomass), , ])
  
  # With fishing, biomass should generally be lower (though ecosystem complexity may vary)
  # At minimum, they should be different
  expect_false(isTRUE(all.equal(final_biomass_zero, final_biomass_high)))
})


test_that("fishing parameters are stored correctly in model output", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 1
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters
  effort_values <- seq(100, 200, length.out = n_steps)
  catchability <- seq(0.01, 0.02, length.out = ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = Groups$W0[g] + 1),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort_values,
    catchability = catchability,
    selectivity = selectivity
  )
  
  results <- zoomss_model(input_params, Groups, isave = 2, fishing_params = fishing_params)
  
  # Check that fishing parameters are in model output
  expect_true(results$param$dynamic_fishing)
  expect_equal(results$param$effort_ts, effort_values)
  expect_equal(results$param$catchability, catchability)
  expect_equal(length(results$param$selectivity_params), ngrps)
})


test_that("zero fishing effort gives same result as static zero fishing", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  
  # Get groups with zero static fishing
  Groups <- getGroups(source = "default")
  Groups$Fmort <- 0  # Set static fishing to zero
  
  # Run 1: Static zero fishing (no fishing_params)
  results_static <- zoomss_model(input_params, Groups, isave = 6)
  
  # Run 2: Dynamic zero fishing
  ngrps <- nrow(Groups)
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = Groups$W0[g]),
      log_scale = TRUE
    )
  }
  
  fishing_dynamic <- list(
    effort = rep(0, n_steps),
    catchability = rep(0.01, ngrps),  # Doesn't matter, effort is zero
    selectivity = selectivity
  )
  
  results_dynamic <- zoomss_model(input_params, Groups, isave = 6, fishing_params = fishing_dynamic)
  
  # Final abundances should be very similar (allowing for numerical precision)
  expect_equal(
    results_static$abundance[nrow(results_static$abundance), , ],
    results_dynamic$abundance[nrow(results_dynamic$abundance), , ],
    tolerance = 1e-6
  )
})


test_that("selectivity matrix is calculated correctly for each group", {
  # Create parameters
  n_years <- 1
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters with different selectivity for each group
  effort <- rep(100, n_steps)
  catchability <- rep(0.01, ngrps)
  
  w_min_values <- seq(-10, 0, length.out = ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = w_min_values[g]),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run model
  results <- zoomss_model(input_params, Groups, isave = 2, fishing_params = fishing_params)
  
  # Results should contain parameter information
  expect_true(results$param$dynamic_fishing)
  expect_equal(length(results$param$selectivity_params), ngrps)
  
  # Each group should have different selectivity threshold
  for (g in 1:ngrps) {
    expect_equal(results$param$selectivity_params[[g]]$params$w_min, w_min_values[g])
  }
})


test_that("fishing mortality output has correct dimensions and values", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters
  effort <- seq(0, 200, length.out = n_steps)
  catchability <- rep(0.01, ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = Groups$W0[g]),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run model
  results <- zoomss_model(input_params, Groups, isave = 4, fishing_params = fishing_params)
  
  # Check that fishing_mortality exists in output
  expect_true("fishing_mortality" %in% names(results))
  
  # Check dimensions: (nsave x ngrps x ngrid)
  nsave <- floor(n_steps / 4)
  ngrid <- length(results$param$w)
  
  expect_equal(dim(results$fishing_mortality), c(nsave, ngrps, ngrid))
  
  # Check that F values are non-negative
  expect_true(all(results$fishing_mortality >= 0, na.rm = TRUE))
  
  # Check that F increases over time (effort increases)
  # Compare first and last timestep
  F_first <- results$fishing_mortality[1, , ]
  F_last <- results$fishing_mortality[nsave, , ]
  
  # Last timestep should have higher F than first (effort increases)
  expect_true(mean(F_last, na.rm = TRUE) > mean(F_first, na.rm = TRUE))
})


test_that("catch output has correct dimensions and values (total only)", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters
  effort <- rep(150, n_steps)
  catchability <- rep(0.015, ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "logistic",
      params = list(
        L50 = (Groups$W0[g] + Groups$Wmax[g]) / 2,
        L95 = (Groups$W0[g] + Groups$Wmax[g]) / 2 + 0.5
      ),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run model with save_catch_by_size = FALSE (default)
  results <- zoomss_model(input_params, Groups, isave = 3, 
                          fishing_params = fishing_params, 
                          save_catch_by_size = FALSE)
  
  # Check that catch exists in output
  expect_true("catch" %in% names(results))
  
  # Check that catch_by_size does NOT exist
  expect_false("catch_by_size" %in% names(results))
  
  # Check dimensions: (nsave x ngrps)
  nsave <- floor(n_steps / 3)
  expect_equal(dim(results$catch), c(nsave, ngrps))
  
  # Check that catch values are non-negative
  expect_true(all(results$catch >= 0, na.rm = TRUE))
  
  # Check that catch is positive (non-zero fishing)
  expect_true(any(results$catch > 0, na.rm = TRUE))
})


test_that("catch_by_size output has correct dimensions and values", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing parameters with knife-edge selectivity
  effort <- rep(150, n_steps)
  catchability <- rep(0.015, ngrps)
  
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = (Groups$W0[g] + Groups$Wmax[g]) / 2),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run model with save_catch_by_size = TRUE
  results <- zoomss_model(input_params, Groups, isave = 3, 
                          fishing_params = fishing_params, 
                          save_catch_by_size = TRUE)
  
  # Check that both catch outputs exist
  expect_true("catch" %in% names(results))
  expect_true("catch_by_size" %in% names(results))
  
  # Check dimensions
  nsave <- floor(n_steps / 3)
  ngrid <- length(results$param$w)
  
  expect_equal(dim(results$catch), c(nsave, ngrps))
  expect_equal(dim(results$catch_by_size), c(nsave, ngrps, ngrid))
  
  # Check that catch_by_size values are non-negative
  expect_true(all(results$catch_by_size >= 0, na.rm = TRUE))
  
  # Check that catch_by_size is positive (non-zero fishing)
  expect_true(any(results$catch_by_size > 0, na.rm = TRUE))
  
  # Verify that total catch equals sum of catch_by_size
  for (t in 1:nsave) {
    for (g in 1:ngrps) {
      total_from_size <- sum(results$catch_by_size[t, g, ], na.rm = TRUE)
      total_direct <- results$catch[t, g]
      expect_equal(total_from_size, total_direct, tolerance = 1e-10)
    }
  }
})


test_that("catch_by_size respects selectivity patterns", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 1
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Create fishing with strong knife-edge selectivity
  effort <- rep(200, n_steps)
  catchability <- rep(0.02, ngrps)
  
  # Set selectivity threshold at midpoint of size range
  selectivity <- vector("list", ngrps)
  w_min_threshold <- (Groups$W0[1] + Groups$Wmax[1]) / 2
  
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = w_min_threshold),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run model with size-resolved catch
  results <- zoomss_model(input_params, Groups, isave = 2, 
                          fishing_params = fishing_params, 
                          save_catch_by_size = TRUE)
  
  # Check that catch is zero for small sizes, positive for large sizes
  w_log10 <- results$param$w_log10
  small_sizes <- w_log10 < w_min_threshold
  large_sizes <- w_log10 >= w_min_threshold
  
  # For at least one group and timestep, verify selectivity pattern
  final_timestep <- nrow(results$catch_by_size)
  group_1_catch <- results$catch_by_size[final_timestep, 1, ]
  
  # Catch should be zero (or very small) for sizes below threshold
  if (any(small_sizes)) {
    expect_true(all(group_1_catch[small_sizes] < 1e-6 | is.na(group_1_catch[small_sizes])))
  }
  
  # Catch should be positive for at least some sizes above threshold
  # (May not be positive for all due to low abundance at large sizes)
  if (any(large_sizes)) {
    expect_true(any(group_1_catch[large_sizes] > 0, na.rm = TRUE))
  }
})


test_that("zero fishing produces zero catch", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 1
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  Groups <- getGroups(source = "default")
  ngrps <- nrow(Groups)
  
  # Zero effort
  selectivity <- vector("list", ngrps)
  for (g in 1:ngrps) {
    selectivity[[g]] <- list(
      type = "knife_edge",
      params = list(w_min = Groups$W0[g]),
      log_scale = TRUE
    )
  }
  
  fishing_params <- list(
    effort = rep(0, n_steps),
    catchability = rep(0.01, ngrps),
    selectivity = selectivity
  )
  
  # Run model
  results <- zoomss_model(input_params, Groups, isave = 2, fishing_params = fishing_params)
  
  # All F should be zero
  expect_equal(results$fishing_mortality, array(0, dim = dim(results$fishing_mortality)))
  
  # All catch should be zero
  expect_equal(results$catch, matrix(0, nrow = nrow(results$catch), ncol = ncol(results$catch)))
})


test_that("static fishing still produces outputs", {
  skip_if_not_installed("zoomss")
  
  # Create environmental data
  n_years <- 2
  n_steps <- n_years * 12
  time <- seq(0, n_years, length.out = n_steps)
  sst <- rep(15, n_steps)
  chl <- rep(0.5, n_steps)
  
  input_params <- data.frame(time = time, sst = sst, chl = chl)
  
  # Use static fishing
  Groups <- getGroups(source = "default")
  Groups$Fmort[1] <- 0.3  # Add some fishing mortality
  Groups$Fmort_W0[1] <- Groups$W0[1]
  Groups$Fmort_Wmax[1] <- Groups$Wmax[1]
  
  # Run model without fishing_params (static fishing)
  results <- zoomss_model(input_params, Groups, isave = 4)
  
  # Should still have F and catch outputs
  expect_true("fishing_mortality" %in% names(results))
  expect_true("catch" %in% names(results))
  
  # F should be constant over time for static fishing
  nsave <- nrow(results$fishing_mortality)
  F_first <- results$fishing_mortality[1, 1, ]
  F_last <- results$fishing_mortality[nsave, 1, ]
  expect_equal(F_first, F_last)
  
  # Catch should be positive for group with fishing
  expect_true(any(results$catch[, 1] > 0, na.rm = TRUE))
})
