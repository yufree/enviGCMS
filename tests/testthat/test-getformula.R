context("getformula")

test_that("getformula finds H2O for water mass", {
    # Use wider window since decomposeMass needs sufficient tolerance
    result <- getformula(18.010565, window = 0.01,
                         elements = list(H = c(0, 10), O = c(0, 5),
                                         C = c(0, 5), N = c(0, 5)))
    expect_true(is.list(result))
    expect_true(length(result) >= 1)
    expect_true(length(result[[1]]) > 0)
    # H2O should be among the valid formulae
    expect_true("H2O" %in% result[[1]])
})

test_that("getformula returns list with one entry per mass", {
    result <- getformula(c(18.010565, 28.031300), window = 0.01)
    expect_true(is.list(result))
    expect_equal(length(result), 2)
})

test_that("getformula returns character vector entries", {
    result <- getformula(18.010565, window = 0.01,
                         elements = list(H = c(0, 10), O = c(0, 5),
                                         C = c(0, 5), N = c(0, 5)))
    expect_true(is.character(result[[1]]))
})
