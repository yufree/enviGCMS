context("getalign")

test_that("getalign finds overlapping m/z values without rt", {
    mz1 <- c(100.0001, 200.0002, 300.0003)
    mz2 <- c(100.0002, 400.0004)
    result <- getalign(mz1, mz2)
    expect_true(is.data.frame(result))
    # mz1[1] and mz2[1] should match within 10 ppm
    expect_true(nrow(result) >= 1)
    expect_true("mz1" %in% colnames(result))
    expect_true("mz2" %in% colnames(result))
})

test_that("getalign finds overlapping with rt filtering", {
    mz1 <- c(100.0001, 200.0002)
    mz2 <- c(100.0002, 200.0003)
    rt1 <- c(100, 200)
    rt2 <- c(105, 300)
    # deltart = 10, so first pair matches (rt diff = 5), second does not (rt diff = 100)
    result <- getalign(mz1, mz2, rt1, rt2, deltart = 10)
    expect_true(is.data.frame(result))
    expect_true(nrow(result) >= 1)
    expect_true("rt1" %in% colnames(result))
})

test_that("getalign returns message when no overlap", {
    mz1 <- c(100.0)
    mz2 <- c(200.0)
    expect_message(getalign(mz1, mz2), "No result")
})
