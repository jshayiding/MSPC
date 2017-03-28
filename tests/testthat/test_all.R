#' unit test for MSPC Package
#'

library(MSPC)
library(rtracklayer)
library(GenomicRanges)

## prepare test data set
beds <- getPeakFile()[1:3]
grs <- readPeakFiles(beds, pvalueBase = 1L)
total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04)
options(scipen = 0)

#' check file extension
#'
context("Testthat context ...")
test_that("getPeakFile", {
    files <- getPeakFile()
    expect_equal(length(files),11)
    lapply(files, function(x) {
        expect_match(grep("\\.bed$", x, value = TRUE), ".bed")
    })
})

test_that("check MSPC pipeline", {
    expect_true(class(grs[[1]])=="GRanges")
    expect_equal(length(mcols(grs[[1]])), 3)
    len <- lapply(total.ERs, length)
    expect_equivalent(len, list(283, 459, 297))
    hit <- peakOverlapping(peakset = total.ERs, FUN = which.max)
    expect_true(class(hit[[1]])=="CompressedIntegerList")
    keepList <- filterByOverlapHit(.ovHit = hit, peakset = total.ERs,
                                   replicate.type = "Biological", isSuffOverlap = TRUE)
    expect_true(class(keepList[[1]])=="CompressedIntegerList")
    initDiscERs <- filterByOverlapHit(.ovHit = hit, peakset = total.ERs,
                                      replicate.type = "Biological", isSuffOverlap = FALSE)
    expect_true(class(initDiscERs[[1]])=="GRanges")
    confirmedERs <- filterBycombStringency(total.ERs, .hitList = keepList,
                                           cmbstrgThreshold = 1.0E-08, isFisherPass = TRUE)
    len_confER <- lapply(confirmedERs, length)
    expect_equivalent(len_confER, list(241, 465,324))
    # ERs that failed from minimum overlapping peak requirement
    initDiscERs <- filterByOverlapHit(.ovHit = hit, peakset = total.ERs,
                                      replicate.type = "Biological", isSuffOverlap = FALSE)
    ## ERs that failed from fisher method
    fisherDiscERs <- filterBycombStringency(total.ERs, .hitList = keepList,
                                            cmbstrgThreshold = 1.0E-08, isFisherPass = FALSE)
    discardedERs <- mapply(c, initDiscERs, fisherDiscERs)
    len_disc <- lapply(discardedERs, length)
    expect_equivalent(len_disc, list(42, 195,107))
})
