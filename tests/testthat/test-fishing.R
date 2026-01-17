test_that("calc_fishing_mortality validates inputs correctly", {
  # Valid inputs
  ngrps <- 3
  ngrid <- 50
  q <- c(0.01, 0.02, 0.015)
  selectivity <- matrix(runif(ngrps * ngrid, 0, 1), nrow = ngrps, ncol = ngrid)
  
  expect_no_error(calc_fishing_mortality(100, q, selectivity))
  expect_no_error(calc_fishing_mortality(0, q, selectivity))  # Zero effort valid
  
  # Invalid effort
  expect_error(calc_fishing_mortality(c(100, 200), q, selectivity))  # Vector
  expect_error(calc_fishing_mortality(-10, q, selectivity))  # Negative
  expect_error(calc_fishing_mortality(NULL, q, selectivity))
  
  # Invalid q
  expect_error(calc_fishing_mortality(100, numeric(0), selectivity))  # Empty
  expect_error(calc_fishing_mortality(100, c(-0.01, 0.02), selectivity[1:2, ]))  # Negative
  expect_error(calc_fishing_mortality(100, "not_numeric", selectivity))
  
  # Invalid selectivity
  expect_error(calc_fishing_mortality(100, q, as.vector(selectivity)))  # Not matrix
  expect_error(calc_fishing_mortality(100, q, selectivity[1:2, ]))  # Wrong dimensions
  
  # Selectivity out of range (warning, not error)
  sel_bad <- selectivity
  sel_bad[1, 1] <- 1.5
  expect_warning(calc_fishing_mortality(100, q, sel_bad), "outside")
})


test_that("calc_fishing_mortality calculates F correctly", {
  ngrps <- 3
  ngrid <- 50
  
  # Simple case: uniform selectivity
  q <- c(0.01, 0.02, 0.03)
  selectivity <- matrix(1, nrow = ngrps, ncol = ngrid)  # All sizes vulnerable
  effort <- 100
  
  F_matrix <- calc_fishing_mortality(effort, q, selectivity)
  
  # Check dimensions
  expect_equal(dim(F_matrix), c(ngrps, ngrid))
  
  # Check values: F should equal effort * q for each group
  expect_equal(F_matrix[1, 1], 100 * 0.01)
  expect_equal(F_matrix[2, 1], 100 * 0.02)
  expect_equal(F_matrix[3, 1], 100 * 0.03)
  
  # All size classes should have same F for uniform selectivity
  expect_true(all(F_matrix[1, ] == F_matrix[1, 1]))
  expect_true(all(F_matrix[2, ] == F_matrix[2, 1]))
})


test_that("calc_fishing_mortality handles zero effort", {
  ngrps <- 2
  ngrid <- 20
  q <- c(0.01, 0.02)
  selectivity <- matrix(runif(ngrps * ngrid, 0, 1), nrow = ngrps, ncol = ngrid)
  
  F_matrix <- calc_fishing_mortality(0, q, selectivity)
  
  # All F should be zero when effort is zero
  expect_equal(F_matrix, matrix(0, nrow = ngrps, ncol = ngrid))
})


test_that("calc_fishing_mortality respects selectivity patterns", {
  ngrps <- 2
  ngrid <- 50
  w_log10 <- seq(-2, 3, length.out = ngrid)
  
  q <- c(0.01, 0.01)  # Same catchability
  effort <- 100
  
  # Different selectivity patterns
  selectivity <- matrix(0, nrow = ngrps, ncol = ngrid)
  selectivity[1, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = -2))  # All sizes vulnerable
  selectivity[2, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 2))   # Only large sizes
  
  F_matrix <- calc_fishing_mortality(effort, q, selectivity)
  
  # Group 1 should have F > 0 everywhere (w_min = -2, at lower bound)
  expect_true(all(F_matrix[1, ] > 0))
  
  # Group 2 should have F = 0 for small sizes, F > 0 for large sizes
  small_sizes <- w_log10 < 2
  large_sizes <- w_log10 >= 2
  
  expect_equal(F_matrix[2, small_sizes], rep(0, sum(small_sizes)))
  expect_true(all(F_matrix[2, large_sizes] > 0))
  
  # For large sizes, both groups should have same F (same q and effort)
  expect_equal(F_matrix[1, large_sizes], F_matrix[2, large_sizes])
})


test_that("calc_catch validates inputs correctly", {
  ngrps <- 3
  ngrid <- 50
  F_w <- matrix(runif(ngrps * ngrid, 0, 0.5), nrow = ngrps, ncol = ngrid)
  N_w <- matrix(runif(ngrps * ngrid, 0, 1e6), nrow = ngrps, ncol = ngrid)
  w <- 10^seq(-2, 3, length.out = ngrid)
  
  # Valid inputs
  expect_no_error(calc_catch(F_w, N_w, w))
  expect_no_error(calc_catch(F_w, N_w, w, by_size = TRUE))
  
  # Invalid F_w
  expect_error(calc_catch(as.vector(F_w), N_w, w))  # Not matrix
  expect_error(calc_catch(F_w[1:2, ], N_w, w))  # Wrong dimensions
  
  # Invalid N_w
  expect_error(calc_catch(F_w, as.vector(N_w), w))  # Not matrix
  expect_error(calc_catch(F_w, N_w[, 1:40], w))  # Wrong dimensions
  
  # Invalid w
  expect_error(calc_catch(F_w, N_w, w[1:40]))  # Wrong length
  expect_error(calc_catch(F_w, N_w, c(-1, w[-1])))  # Negative mass
  expect_error(calc_catch(F_w, N_w, "not_numeric"))
  
  # Invalid dw and dt
  expect_error(calc_catch(F_w, N_w, w, dw = -0.1))
  expect_error(calc_catch(F_w, N_w, w, dt = 0))
  expect_error(calc_catch(F_w, N_w, w, dw = c(0.1, 0.2)))
  
  # Negative F or N (warning, not error)
  F_bad <- F_w
  F_bad[1, 1] <- -0.1
  expect_warning(calc_catch(F_bad, N_w, w), "Negative fishing mortality")
  
  N_bad <- N_w
  N_bad[1, 1] <- -100
  expect_warning(calc_catch(F_w, N_bad, w), "Negative abundance")
})


test_that("calc_catch calculates catch correctly", {
  ngrps <- 2
  ngrid <- 3
  
  # Simple case with known values
  F_w <- matrix(c(0.1, 0.2,
                  0.1, 0.2,
                  0.1, 0.2), nrow = ngrps, ncol = ngrid, byrow = FALSE)
  
  N_w <- matrix(c(1000, 2000,
                  500, 1000,
                  100, 200), nrow = ngrps, ncol = ngrid, byrow = FALSE)
  
  w <- c(1, 10, 100)  # Simple masses
  dw <- 0.1
  dt <- 1.0
  
  # Calculate expected catch for first group, first size
  # C = F * N * w * log(10) * dw * dt
  expected_g1_w1 <- 0.1 * 1000 * 1 * log(10) * 0.1 * 1.0
  
  catch_by_size <- calc_catch(F_w, N_w, w, dw, dt, by_size = TRUE)
  expect_equal(catch_by_size[1, 1], expected_g1_w1, tolerance = 1e-10)
  
  # Calculate total catch
  catch_total <- calc_catch(F_w, N_w, w, dw, dt, by_size = FALSE)
  
  # Total should equal sum across sizes
  expect_equal(catch_total[1], sum(catch_by_size[1, ]))
  expect_equal(catch_total[2], sum(catch_by_size[2, ]))
  expect_equal(length(catch_total), ngrps)
})


test_that("calc_catch handles edge cases", {
  ngrps <- 2
  ngrid <- 10
  F_w <- matrix(0.1, nrow = ngrps, ncol = ngrid)
  N_w <- matrix(1000, nrow = ngrps, ncol = ngrid)
  w <- 10^seq(-2, 1, length.out = ngrid)
  
  # Zero fishing mortality
  F_zero <- matrix(0, nrow = ngrps, ncol = ngrid)
  catch_zero <- calc_catch(F_zero, N_w, w)
  expect_equal(catch_zero, c(0, 0))
  
  # Zero abundance
  N_zero <- matrix(0, nrow = ngrps, ncol = ngrid)
  catch_extinct <- calc_catch(F_w, N_zero, w)
  expect_equal(catch_extinct, c(0, 0))
  
  # Both zero
  catch_none <- calc_catch(F_zero, N_zero, w)
  expect_equal(catch_none, c(0, 0))
})


test_that("calc_catch scales correctly with time step", {
  ngrps <- 2
  ngrid <- 10
  F_w <- matrix(0.1, nrow = ngrps, ncol = ngrid)
  N_w <- matrix(1000, nrow = ngrps, ncol = ngrid)
  w <- 10^seq(-2, 1, length.out = ngrid)
  dw <- 0.1
  
  # Annual catch
  catch_annual <- calc_catch(F_w, N_w, w, dw, dt = 1.0)
  
  # Monthly catch (dt = 1/12)
  catch_monthly <- calc_catch(F_w, N_w, w, dw, dt = 1/12)
  
  # Monthly should be 1/12 of annual
  expect_equal(catch_monthly, catch_annual / 12, tolerance = 1e-10)
})


test_that("calc_catch gives realistic size-based patterns", {
  ngrps <- 1
  ngrid <- 50
  w <- 10^seq(-2, 3, length.out = ngrid)
  w_log10 <- log10(w)
  
  # Power-law size spectrum (common in ecology)
  N_w <- matrix(1e6 * w^(-2), nrow = ngrps, ncol = ngrid, byrow = TRUE)
  
  # Knife-edge selectivity
  selectivity <- matrix(calc_selectivity(w_log10, "knife_edge", list(w_min = 1)),
                        nrow = ngrps, ncol = ngrid, byrow = TRUE)
  
  # Calculate F and catch
  q <- 0.01
  effort <- 100
  F_w <- calc_fishing_mortality(effort, q, selectivity)
  
  catch_by_size <- calc_catch(F_w, N_w, w, dw = 0.1, by_size = TRUE)
  
  # Catch should be zero for small sizes (below selectivity threshold)
  small_sizes <- w_log10 < 1
  expect_equal(catch_by_size[1, small_sizes], rep(0, sum(small_sizes)))
  
  # Catch should be > 0 for vulnerable sizes
  vulnerable_sizes <- w_log10 >= 1
  expect_true(all(catch_by_size[1, vulnerable_sizes] > 0))
  
  # Total catch should be positive
  catch_total <- calc_catch(F_w, N_w, w, dw = 0.1)
  expect_true(catch_total[1] > 0)
})


test_that("fishing functions work together correctly", {
  # Integration test: selectivity -> F -> catch
  ngrps <- 3
  ngrid <- 50
  w <- 10^seq(-2, 3, length.out = ngrid)
  w_log10 <- log10(w)
  
  # Create selectivity for each group (different thresholds)
  selectivity <- matrix(0, nrow = ngrps, ncol = ngrid)
  selectivity[1, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 0.5))
  selectivity[2, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 1.0))
  selectivity[3, ] <- calc_selectivity(w_log10, "knife_edge", list(w_min = 1.5))
  
  # Create abundance (power-law)
  N_w <- matrix(0, nrow = ngrps, ncol = ngrid)
  for (g in 1:ngrps) {
    N_w[g, ] <- 1e6 * w^(-2)
  }
  
  # Calculate fishing mortality
  q <- c(0.01, 0.015, 0.02)
  effort <- 100
  F_w <- calc_fishing_mortality(effort, q, selectivity)
  
  # Calculate catch
  catch_total <- calc_catch(F_w, N_w, w, dw = 0.1)
  
  # Basic sanity checks
  expect_true(all(catch_total > 0))  # All groups should have some catch
  expect_equal(length(catch_total), ngrps)
  
  # Group with higher q should have higher catch (all else equal)
  # (This may not always hold due to selectivity differences, but worth checking)
  
  # Check that zero effort gives zero catch
  F_zero <- calc_fishing_mortality(0, q, selectivity)
  catch_zero <- calc_catch(F_zero, N_w, w, dw = 0.1)
  expect_equal(catch_zero, rep(0, ngrps))
})


test_that("calc_catch handles NA values appropriately", {
  ngrps <- 2
  ngrid <- 10
  F_w <- matrix(0.1, nrow = ngrps, ncol = ngrid)
  N_w <- matrix(1000, nrow = ngrps, ncol = ngrid)
  w <- 10^seq(-2, 1, length.out = ngrid)
  
  # Introduce NA
  F_w[1, 5] <- NA
  N_w[2, 3] <- NA
  
  # Should still calculate (rowSums with na.rm = TRUE)
  catch_total <- calc_catch(F_w, N_w, w)
  expect_type(catch_total, "double")
  expect_length(catch_total, ngrps)
  
  # NA contributions should be ignored (treated as zero)
  expect_true(!is.na(catch_total[1]))
  expect_true(!is.na(catch_total[2]))
})
