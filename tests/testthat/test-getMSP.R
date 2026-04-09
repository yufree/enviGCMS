context("getMSP")

test_that("getMSP parses a simple MSP file correctly", {
    # Create a temporary MSP file
    msp_content <- c(
        "BEGIN IONS",
        "Name: TestCompound",
        "Formula: C6H12O6",
        "Num Peaks: 3",
        "60 500",
        "73 999",
        "147 300",
        "END IONS"
    )
    tmp <- tempfile(fileext = ".msp")
    writeLines(msp_content, tmp)
    on.exit(unlink(tmp))

    result <- getMSP(tmp)
    expect_true(is.list(result))
    expect_equal(length(result), 1)

    comp <- result[[1]]
    expect_equal(comp$name, "TestCompound")
    expect_equal(comp$formula, "C6H12O6")
})

test_that("getMSP extracts spectra with correct m/z values", {
    msp_content <- c(
        "BEGIN IONS",
        "Name: TestCompound",
        "Num Peaks: 3",
        "60 500",
        "73 999",
        "147 300",
        "END IONS"
    )
    tmp <- tempfile(fileext = ".msp")
    writeLines(msp_content, tmp)
    on.exit(unlink(tmp))

    result <- getMSP(tmp)
    comp <- result[[1]]

    expect_true(!is.null(comp$spectra))
    expect_true(is.data.frame(comp$spectra))
    expect_true(all(c("mz", "intensity") %in% colnames(comp$spectra)))
    expect_equal(comp$spectra$mz, c(60, 73, 147))
    # Intensities are normalized to max = 100
    expect_equal(max(comp$spectra$intensity), 100)
})

test_that("getMSP handles multiple compounds", {
    msp_content <- c(
        "BEGIN IONS",
        "Name: Compound1",
        "Num Peaks: 2",
        "100 500",
        "200 1000",
        "END IONS",
        "BEGIN IONS",
        "Name: Compound2",
        "Num Peaks: 2",
        "150 800",
        "250 400",
        "END IONS"
    )
    tmp <- tempfile(fileext = ".msp")
    writeLines(msp_content, tmp)
    on.exit(unlink(tmp))

    result <- getMSP(tmp)
    expect_equal(length(result), 2)
    expect_equal(result[[1]]$name, "Compound1")
    expect_equal(result[[2]]$name, "Compound2")
})
