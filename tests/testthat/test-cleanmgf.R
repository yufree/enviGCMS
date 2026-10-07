context("cleanMGF")

mgf_path <- local({
    tmp <- tempfile(fileext = ".mgf")
    writeLines(c(
        "BEGIN IONS",
        "TITLE=caffeine_test",
        "PEPMASS=195.087652",
        "CHARGE=1+",
        "195.087652 100.0",
        "140.081838 60.0",
        "124.086924 35.0",
        "97.063449 20.0",
        "301.234567 10.0",
        "END IONS",
        "",
        "BEGIN IONS",
        "TITLE=no_pepmass",
        "100.0 50.0",
        "200.0 25.0",
        "END IONS",
        "",
        "BEGIN IONS",
        "TITLE=impossible_precursor",
        "PEPMASS=0.5",
        "CHARGE=1+",
        "150.0 100.0",
        "250.0 40.0",
        "END IONS"
    ), tmp)
    tmp
})

test_that(".parse_mgf parses blocks, headers and peaks", {
    sp <- enviGCMS:::.parse_mgf(mgf_path)
    expect_length(sp, 3)
    expect_true(any(grepl("^PEPMASS=", sp[[1]]$headers)))
    expect_equal(sp[[1]]$peaks$mz,
                 c(195.087652, 140.081838, 124.086924, 97.063449, 301.234567))
    expect_equal(sp[[2]]$peaks$intensity, c(50, 25))
})

test_that(".mgf_charge parses common formats", {
    expect_equal(enviGCMS:::.mgf_charge("1+"), 1L)
    expect_equal(enviGCMS:::.mgf_charge("2+"), 2L)
    expect_equal(enviGCMS:::.mgf_charge("+2"), 2L)
    expect_equal(enviGCMS:::.mgf_charge("2-"), -2L)
    expect_equal(enviGCMS:::.mgf_charge(NA_character_), 1L)
})

test_that(".explain_peaks flags only sub-formula peaks", {
    mz <- c(195.087652, 140.081838, 97.063449, 301.234567)
    mask <- enviGCMS:::.explain_peaks(mz, "C8H10N4O2", 1.007276, 1, 5)
    expect_equal(mask, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("cleanMGF keeps explainable peaks and writes output", {
    out <- tempfile(fileext = ".mgf")
    res <- cleanMGF(mgf_path, out_file = out, adduct = "[M+H]+",
                    min_peaks = 2)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 3)

    # Caffeine spectrum: 4 real sub-formula peaks kept, noise removed
    expect_equal(res$n_kept[1], 4)
    expect_equal(res$n_before[1], 5)
    expect_equal(res$status[1], "cleaned")
    expect_false(is.na(res$formula[1]))

    # Spectrum without PEPMASS passes through unchanged
    expect_equal(res$status[2], "kept (no pepmass)")
    expect_equal(res$n_kept[2], 2)

    # Spectrum with impossible precursor mass has no candidates
    expect_equal(res$status[3], "kept (no formula candidates)")

    # Output file contains FORMULA header and dropped noise peak is gone
    txt <- readLines(out)
    expect_true(any(grepl("^FORMULA=", txt)))
    expect_false(any(grepl("301\\.234567", txt)))
    expect_true(any(grepl("^140\\.081838", txt)))
    # PEPMASS-less spectrum kept in output
    expect_true(any(grepl("100\\.0 50\\.0", txt)))
    # Impossible-precursor spectrum kept in output
    expect_true(any(grepl("^150\\.0 100\\.0", txt)))
})

test_that("cleanMGF drops spectra below min_peaks", {
    out <- tempfile(fileext = ".mgf")
    tmp <- tempfile(fileext = ".mgf")
    writeLines(c(
        "BEGIN IONS",
        "TITLE=only_noise",
        "PEPMASS=195.087652",
        "CHARGE=1+",
        "301.234567 100.0",
        "402.345678 50.0",
        "END IONS"
    ), tmp)
    res <- cleanMGF(tmp, out_file = out, adduct = "[M+H]+", min_peaks = 2)
    expect_equal(res$status[1], "dropped")
    txt <- readLines(out)
    expect_false(any(grepl("only_noise", txt)))
})

test_that("cleanMGF rejects multiple adducts", {
    expect_error(cleanMGF(mgf_path, adduct = c("[M+H]+", "[M+Na]+")),
                 "single adduct")
})
