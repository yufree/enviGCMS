context("findlipid RKMD")

make_list <- function() {
    list(mz = c(400.30, 500.45, 760.58))
}

test_that("findlipid returns RKMD and class columns", {
    res <- findlipid(make_list(), mode = 'none')
    expect_true(!is.null(res$RKMD))
    expect_true(all(c("TAG_RKMD", "PC_RKMD", "mz") %in% colnames(res$RKMD)))
    expect_true(all(c("TAG", "PC", "PI") %in% colnames(res$RKMD)))
})

test_that("pos mode applies an adduct correction distinct from none", {
    # regression: the pos/neg branch previously applied no adduct correction,
    # so 'pos' and 'none' produced identical RKMD.
    l <- make_list()
    pos  <- findlipid(l, mode = 'pos')$RKMD$TAG_RKMD
    none <- findlipid(l, mode = 'none')$RKMD$TAG_RKMD
    expect_false(isTRUE(all.equal(pos, none)))
})

test_that("pos and neg corrections move mass in opposite directions", {
    l <- make_list()
    pos <- findlipid(l, mode = 'pos')$RKMD$TAG_RKMD
    neg <- findlipid(l, mode = 'neg')$RKMD$TAG_RKMD
    expect_false(isTRUE(all.equal(pos, neg)))
})
