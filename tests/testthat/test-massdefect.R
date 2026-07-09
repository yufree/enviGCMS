context("Mass defect helpers")

test_that("getrmd returns correct relative mass defect", {
    # getrmd(mz) = round((round(mz) - mz) / mz * 10^6)
    # For mz = 28.0313: round(28) = 28, (28 - 28.0313)/28.0313*1e6 = -1117
    rmd <- getrmd(28.0313)
    expected <- round((round(28.0313) - 28.0313) / 28.0313 * 10^6)
    expect_equal(rmd, expected)
})

test_that("getrmd works with vectors", {
    rmd <- getrmd(c(100.05, 200.10))
    expect_equal(length(rmd), 2)
    expect_true(is.numeric(rmd))
})

test_that("getmdr returns correct raw mass defect", {
    # getmdr(mz) = round((round(mz) - mz) * 10^3)
    md <- getmdr(200.3)
    expected <- round((round(200.3) - 200.3) * 10^3)
    expect_equal(md, expected)
})

test_that("getmdr works with vectors", {
    md <- getmdr(c(100.05, 200.10))
    expect_equal(length(md), 2)
    expect_true(is.numeric(md))
})

test_that("getmdh returns data.frame with correct columns for single cus", {
    mz <- c(28.0313, 44.0262)
    result <- getmdh(mz, cus = 'CH2', method = 'round')
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)
    expect_true("mz" %in% colnames(result))
    expect_true("MD1" %in% colnames(result))
})

test_that("getmdh returns data.frame with MD1 and MD2 for two cus", {
    mz <- c(28.0313)
    result <- getmdh(mz, cus = 'CH2,H2', method = 'round')
    expect_true(is.data.frame(result))
    expect_true(all(c("mz", "MD1", "MD2") %in% colnames(result)))
})

test_that("getmdh works with floor method", {
    mz <- c(28.0313)
    result <- getmdh(mz, cus = 'CH2', method = 'floor')
    expect_true(is.data.frame(result))
    expect_true("MD1" %in% colnames(result))
})

test_that("getmdh works with ceiling method", {
    mz <- c(28.0313)
    result <- getmdh(mz, cus = 'CH2', method = 'ceiling')
    expect_true(is.data.frame(result))
    expect_true("MD1" %in% colnames(result))
})

test_that("getmdh returns MD1/MD2/MD3 for three cus across all methods", {
    mz <- c(28.0313)
    for (m in c('round', 'floor', 'ceiling')) {
        result <- getmdh(mz, cus = 'CH2,H2,O', method = m)
        expect_true(all(c("mz", "MD1", "MD2", "MD3") %in% colnames(result)),
                    info = paste("method =", m))
    }
})

test_that("getmdh handles more than three cus with floor method", {
    # regression: floor branch previously assigned MD1_3 instead of MD3
    mz <- c(28.0313)
    expect_silent(res <- suppressMessages(
        getmdh(mz, cus = 'CH2,H2,O,NH', method = 'floor')))
    expect_true(all(c("mz", "MD1", "MD2", "MD3") %in% colnames(res)))
})

test_that("getmassdefect returns data.frame with correct dimensions", {
    mass <- c(100.1022, 245.2122, 267.3144)
    sf <- 0.9988
    mf <- getmassdefect(mass, sf)
    expect_equal(nrow(mf), 3)
    expect_equal(ncol(mf), 3)
    expect_true(all(c("mass", "sm", "sd") %in% colnames(mf)))
})
