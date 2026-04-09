context("getmass")

test_that("getmass returns correct mass for simple formula", {
    mass <- getmass("H2O")
    expect_equal(mass, 18.010565, tolerance = 1e-4)
})

test_that("getmass returns correct mass for single element", {
    mass <- getmass("C")
    expect_equal(mass, 12.0, tolerance = 1e-4)
})

test_that("getmass returns correct mass for reaction (subtraction)", {
    # C2H4 - H2 = ethylene minus hydrogen
    mass <- getmass("C2H4-H2")
    # C2H4 ~ 28.0313, H2 ~ 2.01565, difference ~ 26.01565
    expect_equal(mass, 26.01565, tolerance = 1e-3)
})

test_that("getmass returns numeric value", {
    mass <- getmass("CH4")
    expect_true(is.numeric(mass))
    expect_equal(length(mass), 1)
})
