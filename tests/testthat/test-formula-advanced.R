test_that("checkGoldenRules validates formulas correctly", {
    # Valid natural molecules
    g_glc <- checkGoldenRules("C6H12O6")
    expect_true(g_glc$valid)
    expect_true(g_glc$senior_rules)
    expect_true(g_glc$hc_pass)
    expect_equal(g_glc$hc_ratio, 2)
    expect_equal(g_glc$dbe, 1)

    g_caf <- checkGoldenRules("C8H10N4O2")
    expect_true(g_caf$valid)

    # Chemically impossible formula with extreme H/C and negative DBE
    g_bad <- checkGoldenRules("C3H32N12O3")
    expect_false(g_bad$valid)
    expect_false(g_bad$hc_pass)
    expect_false(g_bad$dbe_pass)

    # Formula violating Senior rules
    g_rad <- checkGoldenRules("CH3", z = 0) # radical (odd valence sum = 7)
    expect_false(g_rad$senior_rules)
    expect_false(g_rad$valid)
})

test_that("scoreIsotopes calculates cosine and log-likelihood", {
    obs_m <- c(180.0634, 181.0668, 182.0712)
    obs_i <- c(100.0, 6.7, 1.2)
    theo_m <- c(180.06339, 181.06674, 182.0718)
    theo_i <- c(1.0, 0.066, 0.012)

    sc <- scoreIsotopes(obs_m, obs_i, theo_m, theo_i, tolerance = 0.005)
    expect_gt(sc$cosine, 0.99)
    expect_gt(sc$weighted_cosine, 0.99)
    expect_equal(sc$matched_peaks, 3)
    expect_true(is.numeric(sc$log_likelihood))
})

test_that("getformula supports vectorization and golden rules", {
    mzs <- c(180.0634, 194.0804)
    res <- getformula(mzs, golden_rules = TRUE, nthreads = 1)
    expect_length(res, 2)
    expect_true("C6H12O6" %in% res[[1]])
    expect_true("C8H10N4O2" %in% res[[2]])
})
