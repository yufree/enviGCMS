context("Peak integration")

# Synthetic gaussian chromatogram with a single peak near RT 8.65 and a flat
# baseline, used to exercise the RT-window subsetting.
make_chrom <- function() {
    rt <- seq(8.0, 9.5, by = 0.01)
    int <- 1000 * exp(-((rt - 8.65)^2) / (2 * 0.05^2)) + 50
    data.frame(rt = rt, intensity = int)
}

test_that("integration returns a positive finite area for an ascending RT range", {
    # regression: subset used `> rt[2] & < rt[1]` (impossible), so the peak
    # window was always empty and the function errored / returned NA.
    d <- make_chrom()
    area <- integration(d, rt = c(8.3, 9), brt = c(8.0, 8.2), smoothit = FALSE)
    expect_true(is.finite(area))
    expect_gt(area, 0)
})

test_that("getintegration returns a positive area on the same data", {
    d <- make_chrom()
    res <- getintegration(d, rt = c(8.3, 9))
    expect_true(is.finite(res$area))
    expect_gt(res$area, 0)
})
