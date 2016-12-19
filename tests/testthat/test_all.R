#' unit test for MSPC Package
#'

library(MSPC)

#' check file extension
#'
context("check context ...")
test_that("getPeakFile", {
  files <- getPeakFile()
  expect_equal(length(files),8)
  lapply(files, function(x) {
    expect_match(grep("\\.bed$", x, value = TRUE), ".bed")
  })
})



