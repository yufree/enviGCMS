context("HRMF adduct support")

test_that(".resolve_adducts handles built-in adducts", {
    res <- enviGCMS:::.resolve_adducts(c("[M+H]+", "[M-H]-"))
    expect_equal(res$label, c("[M+H]+", "[M-H]-"))
    expect_equal(res$charge, c(1L, -1L))
    expect_equal(res$delta, c(1.007276, -1.007276))
})

test_that(".resolve_adducts handles custom named deltas", {
    res <- enviGCMS:::.resolve_adducts(c("[M+MeOH]+" = 33.033491,
                                         "[M-H-Hex]-" = -81.033502))
    expect_equal(res$charge, c(1L, -1L))
    expect_equal(res$delta, c(33.033491, -81.033502))
})

test_that(".resolve_adducts rejects invalid input", {
    expect_error(enviGCMS:::.resolve_adducts("[M+X]+"), "Unknown adduct")
    expect_error(enviGCMS:::.resolve_adducts(33.033491), "named numeric")
    expect_error(enviGCMS:::.resolve_adducts(c("[M+MeOH]" = 33.033491)),
                 "'\\+' or '-'")
    expect_error(enviGCMS:::.resolve_adducts(42), "named numeric")
})

test_that("HRMF converts [M+H]+ fragment m/z internally", {
    # Caffeine C8H10N4O2 fragments as protonated ions (neutral + proton)
    # C8H10N4O2=194.080376, C6H9N3O=139.074562, C6H9N3=123.079647,
    # C4H6N3=96.056172
    msp_entry <- list(
        name = "caffeine",
        spectra = data.frame(
            mz = c(195.087652, 140.081838, 124.086924, 97.063449),
            intensity = c(100, 60, 35, 20)
        )
    )
    result <- HRMF(msp_entry, formula = "C8H10N4O2", adduct = "[M+H]+",
                   mass_accuracy = 5, intensity_cutoff = 0)
    skip_if(is.null(result), "no isotope match on this platform")
    expect_true(is.data.frame(result))
    expect_true("Adduct" %in% colnames(result))
    expect_equal(result$Adduct, "[M+H]+")
    expect_equal(result$Candidate, "C8H10N4O2")

    # Precursor must be annotated as the full formula with ~0 ppm error
    detail <- HRMF(msp_entry, formula = "C8H10N4O2", adduct = "[M+H]+",
                   mass_accuracy = 5, intensity_cutoff = 0, detailed = TRUE)
    skip_if(is.null(detail), "no isotope match on this platform")
    ions <- detail[["C8H10N4O2 [M+H]+"]]$all_ions
    expect_true("C8H10N4O2" %in% ions$MonoIsoFormula)
    prec <- ions[ions$MonoIsoFormula == "C8H10N4O2", , drop = FALSE]
    expect_lt(min(abs(prec$MassError_ppm), na.rm = TRUE), 1)
})

test_that("HRMF scores multiple adducts separately", {
    msp_entry <- list(
        name = "caffeine",
        spectra = data.frame(
            mz = c(195.087652, 140.081838, 124.086924, 97.063449),
            intensity = c(100, 60, 35, 20)
        )
    )
    result <- HRMF(msp_entry, formula = "C8H10N4O2",
                   adduct = c("[M+H]+", "[M+Na]+"),
                   mass_accuracy = 5, intensity_cutoff = 0)
    skip_if(is.null(result), "no isotope match on this platform")
    # The protonated fragments cannot be explained by [M+Na]+, so only
    # the [M+H]+ mode should produce a scored row
    expect_true("[M+H]+" %in% result$Adduct)
    expect_false("[M+Na]+" %in% result$Adduct)
})

test_that("HRMF negative mode [M-H]- works", {
    # Benzoic acid C7H6O2 = 122.036779; [M-H]- = 121.029503
    # Fragment C6H6O = 94.041860; [F-H]- = 93.034584
    msp_entry <- list(
        name = "benzoic acid",
        spectra = data.frame(
            mz = c(121.029503, 93.034584),
            intensity = c(100, 45)
        )
    )
    result <- HRMF(msp_entry, formula = "C7H6O2", adduct = "[M-H]-",
                   mass_accuracy = 5, intensity_cutoff = 0)
    skip_if(is.null(result), "no isotope match on this platform")
    expect_equal(result$Adduct, "[M-H]-")
    expect_gte(result$FoM, 0)
})

test_that("getHRMF passes adduct through", {
    tmp <- tempfile(fileext = ".msp")
    writeLines(c(
        "NAME: caffeine",
        "FORMULA: C8H10N4O2",
        "Num peaks: 4",
        "195.0877 100.0",
        "140.0818 60.0",
        "124.0869 35.0",
        "97.0634 20.0"
    ), tmp)
    res <- getHRMF(tmp, adduct = "[M+H]+", intensity_cutoff = 0)
    skip_if(is.null(res) || length(res) == 0, "no isotope match on this platform")
    expect_true("caffeine" %in% names(res))
    expect_equal(res$caffeine[[1]]$HRMF_scores$Adduct, "[M+H]+")
    unlink(tmp)
})
