## Bioconductor Package for Multiple Sample Peak Calling

library(MSPC)

## check all ERs are stored in GRanges object

context("check Chip-seq replicates imported and all ERs are stored as GRanges")
test_that("readPeakFile", {
  files <- getPeakFile()
  rdBED <- readPeakFiles(files)
  lapply(rdBED, function(x) {
    expect_true(is(x, "GRanges"))
  })
})




