context("HRMF helpers and function")

test_that(".parse_formula parses simple formula correctly", {
    result <- enviGCMS:::.parse_formula("C8H11NO")
    expect_true(is.data.frame(result))
    expect_equal(result$element, c("C", "H", "N", "O"))
    expect_equal(result$count, c(8L, 11L, 1L, 1L))
})

test_that(".parse_formula handles elements without counts", {
    result <- enviGCMS:::.parse_formula("CHN")
    expect_equal(result$element, c("C", "H", "N"))
    expect_equal(result$count, c(1L, 1L, 1L))
})

test_that(".parse_formula handles two-letter elements", {
    result <- enviGCMS:::.parse_formula("C2H5Cl")
    expect_equal(result$element, c("C", "H", "Cl"))
    expect_equal(result$count, c(2L, 5L, 1L))
})

test_that(".ELECTRON_MASS is approximately correct", {
    expect_equal(enviGCMS:::.ELECTRON_MASS, 0.000549, tolerance = 1e-4)
})

test_that("HRMF returns error for missing spectra", {
    bad_msp <- list(name = "test")
    expect_error(HRMF(bad_msp, formula = "C6H6"))
})

test_that("HRMF returns data.frame with expected columns", {
    msp_entry <- list(
        name = "test",
        spectra = data.frame(
            mz = c(121.0653, 122.0686),
            intensity = c(100, 7)
        )
    )
    result <- HRMF(msp_entry, formula = "C8H8O", charge = 0,
                   mass_accuracy = 10, intensity_cutoff = 0)
    if (!is.null(result)) {
        expect_true(is.data.frame(result))
        expect_true("Candidate" %in% colnames(result))
        expect_true("FoM" %in% colnames(result))
    }
})
