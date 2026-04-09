context("writeMSP and round-trip I/O")

test_that("writeMSP creates an MSP file", {
    tmp_dir <- tempdir()
    tmp_name <- file.path(tmp_dir, "testwrite")

    intensity <- c(10000, 20000, 10000, 30000, 5000)
    mz <- c(101, 143, 189, 221, 234)
    test_list <- list(list(
        name = "TestCompound",
        formula = "C10H12",
        spectra = cbind.data.frame(mz = mz, intensity = intensity)
    ))

    writeMSP(test_list, name = tmp_name)
    msp_file <- paste0(tmp_name, ".msp")
    on.exit(unlink(msp_file))

    expect_true(file.exists(msp_file))
})

test_that("writeMSP round-trip preserves m/z values", {
    tmp_dir <- tempdir()
    tmp_name <- file.path(tmp_dir, "testroundtrip")

    intensity <- c(10000, 20000, 30000)
    mz <- c(101, 143, 189)
    test_list <- list(list(
        name = "RoundTrip",
        formula = "C8H10",
        spectra = cbind.data.frame(mz = mz, intensity = intensity)
    ))

    writeMSP(test_list, name = tmp_name)
    msp_file <- paste0(tmp_name, ".msp")
    on.exit(unlink(msp_file))

    # Read back
    result <- getMSP(msp_file)
    expect_true(is.list(result))
    expect_equal(length(result), 1)
    expect_equal(result[[1]]$name, "RoundTrip")
    expect_equal(result[[1]]$spectra$mz, mz)
})

test_that("writeMSP round-trip preserves formula", {
    tmp_dir <- tempdir()
    tmp_name <- file.path(tmp_dir, "testformula")

    intensity <- c(5000, 10000)
    mz <- c(100, 200)
    test_list <- list(list(
        name = "FormulaTest",
        formula = "C6H12O6",
        spectra = cbind.data.frame(mz = mz, intensity = intensity)
    ))

    writeMSP(test_list, name = tmp_name)
    msp_file <- paste0(tmp_name, ".msp")
    on.exit(unlink(msp_file))

    result <- getMSP(msp_file)
    expect_equal(result[[1]]$formula, "C6H12O6")
})
