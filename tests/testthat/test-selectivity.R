test_that("calc_selectivity validates inputs correctly", {
  w <- seq(-2, 3, by = 0.1)
  
  # Valid inputs should not throw errors
  expect_no_error(calc_selectivity(w, "knife_edge", list(w_min = 1)))
  expect_no_error(calc_selectivity(w, "logistic", list(L50 = 1.5, L95 = 2.0)))
  
  # Invalid w
  expect_error(calc_selectivity(NULL, "knife_edge", list(w_min = 1)))
  expect_error(calc_selectivity(character(0), "knife_edge", list(w_min = 1)))
  expect_error(calc_selectivity(numeric(0), "knife_edge", list(w_min = 1)))
  
  # Invalid type
  expect_error(calc_selectivity(w, c("knife_edge", "logistic"), list(w_min = 1)))
  expect_error(calc_selectivity(w, 123, list(w_min = 1)))
  expect_error(calc_selectivity(w, "unknown_type", list(w_min = 1)))
  
  # Invalid params
  expect_error(calc_selectivity(w, "knife_edge", "not_a_list"))
  expect_error(calc_selectivity(w, "knife_edge", list()))  # Missing w_min
  expect_error(calc_selectivity(w, "logistic", list(L50 = 1.5)))  # Missing L95
  expect_error(calc_selectivity(w, "logistic", list(L50 = 2.0, L95 = 1.5)))  # L95 <= L50
})


test_that("knife-edge selectivity works correctly", {
  w <- seq(-2, 3, by = 0.1)
  
  # Test with w_min = 1.0 (log10 scale, so 10g in linear scale)
  s <- calc_selectivity(w, type = "knife_edge", params = list(w_min = 1.0))
  
  # Check output properties
  expect_type(s, "double")
  expect_length(s, length(w))
  expect_true(all(s %in% c(0, 1)))
  
  # Check threshold behavior
  expect_equal(s[w < 1.0], rep(0, sum(w < 1.0)))
  expect_equal(s[w >= 1.0], rep(1, sum(w >= 1.0)))
  
  # Test edge cases
  expect_equal(calc_selectivity(1.0, "knife_edge", list(w_min = 1.0)), 1)
  expect_equal(calc_selectivity(0.9, "knife_edge", list(w_min = 1.0)), 0)
  
  # Test with different thresholds
  s_low <- calc_selectivity(w, "knife_edge", list(w_min = -1.0))
  expect_true(sum(s_low == 1) > sum(s == 1))
  
  s_high <- calc_selectivity(w, "knife_edge", list(w_min = 2.5))
  expect_true(sum(s_high == 1) < sum(s == 1))
})


test_that("logistic selectivity works correctly", {
  w <- seq(-2, 3, by = 0.01)
  
  # Test with L50 = 1.5, L95 = 2.0
  s <- calc_selectivity(w, type = "logistic", params = list(L50 = 1.5, L95 = 2.0))
  
  # Check output properties
  expect_type(s, "double")
  expect_length(s, length(w))
  expect_true(all(s >= 0 & s <= 1))
  
  # Check that curve passes through specified points (within tolerance)
  s_at_L50 <- calc_selectivity(1.5, "logistic", list(L50 = 1.5, L95 = 2.0))
  expect_equal(s_at_L50, 0.5, tolerance = 1e-10)
  
  s_at_L95 <- calc_selectivity(2.0, "logistic", list(L50 = 1.5, L95 = 2.0))
  expect_equal(s_at_L95, 0.95, tolerance = 1e-10)
  
  # Check monotonicity (selectivity increases with size)
  expect_true(all(diff(s) >= 0))
  
  # Check asymptotic behavior
  s_low <- calc_selectivity(-10, "logistic", list(L50 = 1.5, L95 = 2.0))
  expect_true(s_low < 0.01)
  
  s_high <- calc_selectivity(10, "logistic", list(L50 = 1.5, L95 = 2.0))
  expect_true(s_high > 0.99)
})


test_that("custom selectivity works correctly", {
  w <- seq(-2, 3, by = 0.1)
  
  # Test with simple custom function (constant selectivity)
  constant_sel <- function(w) rep(0.5, length(w))
  s <- calc_selectivity(w, type = "custom", params = list(fun = constant_sel))
  expect_equal(s, rep(0.5, length(w)))
  
  # Test with dome-shaped selectivity
  dome_sel <- function(w) exp(-((w - 1.5)^2) / 0.5)
  s_dome <- calc_selectivity(w, type = "custom", params = list(fun = dome_sel))
  
  expect_type(s_dome, "double")
  expect_length(s_dome, length(w))
  expect_true(all(s_dome >= 0 & s_dome <= 1))
  
  # Peak should be near w = 1.5
  peak_idx <- which.max(s_dome)
  expect_equal(w[peak_idx], 1.5, tolerance = 0.1)
  
  # Test error handling for invalid custom function
  expect_error(calc_selectivity(w, "custom", list(fun = "not_a_function")))
  expect_error(calc_selectivity(w, "custom", list()))  # Missing fun
})


test_that("selectivity handles edge cases", {
  # Single value
  expect_equal(calc_selectivity(1.5, "knife_edge", list(w_min = 1.0)), 1)
  expect_equal(calc_selectivity(0.5, "knife_edge", list(w_min = 1.0)), 0)
  
  # Vector with single value
  expect_equal(calc_selectivity(c(1.5), "knife_edge", list(w_min = 1.0)), 1)
  
  # Negative sizes (valid in log10 scale)
  w_neg <- seq(-5, 0, by = 0.5)
  s_neg <- calc_selectivity(w_neg, "knife_edge", list(w_min = -2.0))
  expect_equal(sum(s_neg), sum(w_neg >= -2.0))
})


test_that("selectivity with log_scale parameter works", {
  # The log_scale parameter is documented for future use but currently
  # doesn't change behavior (both linear and log10 scales use same comparison)
  
  w_log <- seq(-2, 3, by = 0.1)  # log10 scale
  w_lin <- 10^w_log              # linear scale
  
  # Both should work (though with different w_min values)
  expect_no_error(calc_selectivity(w_log, "knife_edge", list(w_min = 1.0), log_scale = TRUE))
  expect_no_error(calc_selectivity(w_lin, "knife_edge", list(w_min = 10), log_scale = FALSE))
  
  # Results should be identical if w_min is appropriately scaled
  s_log <- calc_selectivity(w_log, "knife_edge", list(w_min = 1.0), log_scale = TRUE)
  s_lin <- calc_selectivity(w_lin, "knife_edge", list(w_min = 10), log_scale = FALSE)
  expect_equal(s_log, s_lin)
})


test_that("selectivity handles boundary cases correctly", {
  w <- c(-2, -1, 0, 1, 2, 3)
  
  # Knife-edge at exact boundary
  s <- calc_selectivity(w, "knife_edge", list(w_min = 1))
  expect_equal(s, c(0, 0, 0, 1, 1, 1))
  
  # Logistic curve properties
  s_logistic <- calc_selectivity(w, "logistic", list(L50 = 0.5, L95 = 1.5))
  
  # Should be monotonically increasing
  expect_true(all(diff(s_logistic) > 0))
  
  # Values should be strictly between 0 and 1 (except at extremes)
  expect_true(all(s_logistic[2:5] > 0 & s_logistic[2:5] < 1))
})


test_that("selectivity internal functions work correctly", {
  # Test internal knife_edge function directly
  w <- seq(-2, 3, by = 0.1)
  s_knife <- selectivity_knife_edge(w, w_min = 1.0)
  
  expect_type(s_knife, "double")
  expect_true(all(s_knife %in% c(0, 1)))
  expect_equal(sum(s_knife), sum(w >= 1.0))
  
  # Test internal logistic function directly
  s_logistic <- selectivity_logistic(w, L50 = 1.5, L95 = 2.0)
  
  expect_type(s_logistic, "double")
  expect_true(all(s_logistic >= 0 & s_logistic <= 1))
  expect_equal(selectivity_logistic(1.5, 1.5, 2.0), 0.5, tolerance = 1e-10)
  expect_equal(selectivity_logistic(2.0, 1.5, 2.0), 0.95, tolerance = 1e-10)
})


test_that("selectivity warns about out-of-range values from custom functions", {
  w <- seq(-2, 3, by = 0.1)
  
  # Custom function returning values outside [0, 1]
  bad_custom <- function(w) w  # Returns values outside [0, 1]
  
  expect_warning(
    calc_selectivity(w, "custom", list(fun = bad_custom)),
    "Selectivity values outside"
  )
  
  # Check that values are clamped
  s <- suppressWarnings(calc_selectivity(w, "custom", list(fun = bad_custom)))
  expect_true(all(s >= 0 & s <= 1))
})
