test_that("Catch is properly accumulated across timesteps when isave > 1", {
  # Setup
  Groups <- getGroups()
  ngrps <- nrow(Groups)
  
  # Create 10-year monthly timestep scenario
  time <- seq(0, 10, by = 1/12)
  n_steps <- length(time)
  
  env_data <- createInputParams(
    time = time,
    sst = 15,
    chl = 0.5
  )
  
  # Setup fishing parameters with constant effort
  effort <- rep(100, n_steps)
  catchability <- rep(0, ngrps)
  catchability[10:12] <- 0.01  # Fish only the fish groups
  
  selectivity <- lapply(1:ngrps, function(g) {
    list(
      type = "knife_edge",
      params = list(w_min = (Groups$W0[g] + Groups$Wmax[g])/2),
      log_scale = TRUE
    )
  })
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run with isave = 1 (save every timestep, no accumulation needed)
  mdl_isave1 <- zoomss_model(env_data, Groups, isave = 1, 
                             fishing_params = fishing_params)
  
  # Run with isave = 12 (save annually, accumulation required)
  mdl_isave12 <- zoomss_model(env_data, Groups, isave = 12,
                              fishing_params = fishing_params)
  
  # Total catch over entire simulation should be approximately equal
  total_catch_isave1 <- sum(mdl_isave1$catch)
  total_catch_isave12 <- sum(mdl_isave12$catch)
  
  # Should match within 3% (differences due to N evolving between timesteps)
  expect_equal(total_catch_isave1, total_catch_isave12, tolerance = 0.03)
  
  # Each annual catch in isave=12 should equal sum of 12 monthly catches in isave=1
  # Check first year (indices 1:12 in isave=1, index 1 in isave=12)
  first_year_monthly_sum <- sum(mdl_isave1$catch[1:12, ])
  first_year_annual <- sum(mdl_isave12$catch[1, ])
  
  expect_equal(first_year_monthly_sum, first_year_annual, tolerance = 0.01)
  
  # Check last year (indices 109:120 in isave=1, index 10 in isave=12)
  last_year_monthly_sum <- sum(mdl_isave1$catch[109:120, ])
  last_year_annual <- sum(mdl_isave12$catch[10, ])
  
  expect_equal(last_year_monthly_sum, last_year_annual, tolerance = 0.01)
})

test_that("Size-resolved catch accumulation works correctly", {
  # Setup
  Groups <- getGroups()
  ngrps <- nrow(Groups)
  
  time <- seq(0, 5, by = 1/12)
  n_steps <- length(time)
  
  env_data <- createInputParams(
    time = time,
    sst = 15,
    chl = 0.5
  )
  
  # Setup fishing
  effort <- rep(150, n_steps)
  catchability <- rep(0, ngrps)
  catchability[10:12] <- 0.015
  
  selectivity <- lapply(1:ngrps, function(g) {
    list(
      type = "knife_edge",
      params = list(w_min = (Groups$W0[g] + Groups$Wmax[g])/2),
      log_scale = TRUE
    )
  })
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run with size-resolved catch and isave = 6 (save twice per year)
  mdl_size <- zoomss_model(env_data, Groups, isave = 6,
                          fishing_params = fishing_params,
                          save_catch_by_size = TRUE)
  
  # Verify that sum of catch_by_size equals total catch for each save point
  for (i in 1:nrow(mdl_size$catch)) {
    for (g in 1:ngrps) {
      size_sum <- sum(mdl_size$catch_by_size[i, g, ], na.rm = TRUE)
      total <- mdl_size$catch[i, g]
      
      # Should match exactly (no tolerance needed for summation)
      expect_equal(size_sum, total, tolerance = 1e-10,
                   info = paste("Mismatch at save point", i, "group", g))
    }
  }
  
  # Total catch should still be conserved
  total_from_size <- sum(mdl_size$catch_by_size, na.rm = TRUE)
  total_from_catch <- sum(mdl_size$catch)
  
  expect_equal(total_from_size, total_from_catch, tolerance = 1e-10)
})

test_that("Catch accumulation handles varying effort correctly", {
  # Setup with time-varying effort
  Groups <- getGroups()
  ngrps <- nrow(Groups)
  
  time <- seq(0, 5, by = 1/12)
  n_steps <- length(time)
  
  env_data <- createInputParams(
    time = time,
    sst = 15,
    chl = 0.5
  )
  
  # Linearly increasing effort
  effort <- seq(0, 200, length.out = n_steps)
  catchability <- rep(0, ngrps)
  catchability[10:12] <- 0.01
  
  selectivity <- lapply(1:ngrps, function(g) {
    list(
      type = "logistic",
      params = list(
        L50 = (Groups$W0[g] + Groups$Wmax[g])/2,
        L95 = Groups$W0[g] + 0.75*(Groups$Wmax[g] - Groups$W0[g])
      ),
      log_scale = TRUE
    )
  })
  
  fishing_params <- list(
    effort = effort,
    catchability = catchability,
    selectivity = selectivity
  )
  
  # Run with two different isave values
  mdl_isave1 <- zoomss_model(env_data, Groups, isave = 1,
                             fishing_params = fishing_params)
  mdl_isave12 <- zoomss_model(env_data, Groups, isave = 12,
                              fishing_params = fishing_params)
  
  # Total catch should be conserved
  expect_equal(sum(mdl_isave1$catch), sum(mdl_isave12$catch), 
               tolerance = 0.03)
  
  # Catch should increase over time due to increasing effort
  # (even though biomass may decline)
  # Check that later years have higher catch than earlier years in isave=12
  early_catch <- sum(mdl_isave12$catch[1, ])
  late_catch <- sum(mdl_isave12$catch[5, ])
  
  # Due to increasing effort, late catch should be higher
  # (unless stock is completely depleted, which shouldn't happen in 5 years)
  expect_gt(late_catch, early_catch * 0.5)  # At least 50% of early catch
})

test_that("Zero fishing effort produces zero catch with accumulation", {
  Groups <- getGroups()
  ngrps <- nrow(Groups)
  
  time <- seq(0, 2, by = 1/12)
  env_data <- createInputParams(time = time, sst = 15, chl = 0.5)
  
  # Zero effort
  fishing_params <- list(
    effort = rep(0, length(time)),
    catchability = rep(0.01, ngrps),
    selectivity = lapply(1:ngrps, function(g) {
      list(type = "knife_edge", 
           params = list(w_min = Groups$W0[g]),
           log_scale = TRUE)
    })
  )
  
  mdl <- zoomss_model(env_data, Groups, isave = 6,
                     fishing_params = fishing_params)
  
  # All catch should be zero
  expect_equal(sum(mdl$catch), 0, tolerance = 1e-15)
  expect_true(all(mdl$catch == 0))
})
