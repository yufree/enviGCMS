context("Mass Defect")

test_that("Length of Mass defect dataframe", {
    mass <- c(100.1022, 245.2122, 267.3144, 400.1222, 707.2294)
    sf <- 0.9988
    mf <- getmassdefect(mass, sf)
    expect_equal(nrow(mf), 5)
    
})
