#' unit test for MSPC Package
#'

library(MSPC)
library(rtracklayer)
library(GenomicRanges)

## prepare test data set
files <- getPeakFile()[1:3]
grs <- readPeakFiles(files, pvalueBase=1L)
options(scipen = 0)

## Exclude background noise
total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
                         overwrite = TRUE)

## explore set of confirmed, discarde peaks
confirmedERs <- runMSPC(peakset = total.ERs, whichType = "max",
                        cmbStrgThreshold = 1.0E-08, isConfirmed = TRUE)
discardedERs <- runMSPC(peakset = total.ERs, whichType = "max",
                        cmbStrgThreshold = 1.0E-08, isConfirmed = FALSE)

#' check file extension
#'
context("Testthat context ...")
test_that("getPeakFile", {
  files <- getPeakFile()
  expect_equal(length(files),11)
  lapply(files, function(x) {
    expect_match(grep("\\.bed$", x, value = TRUE), ".bed")
  })
  expect_true(class(grs)=="GRangesList")
  len <- lapply(total.ERs, length)
  expect_equivalent(len, list(283, 459,297))
  len_confER <- lapply(confirmedERs, length)
  expect_equivalent(len_confER, list(243, 267,193))
  len_disc <- lapply(discardedERs, length)
  expect_equivalent(len_disc, list(40, 192,104))
})

