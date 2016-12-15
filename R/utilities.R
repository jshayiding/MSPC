##' load input peak bed files
##'
##' Input peak bed files can be detected by system.file(),
##' and let package example data available for use.
##'
##' @return list of Chip-seq replicates' name
##' @export
##' @importFrom IRanges gsub
##' @author Julaiti Shayiding

getPeakFile <- function() {
  dir <- system.file("extdata", package="MSPC")
  files <- list.files(dir)
  ER <- gsub(pattern='wgEncode\\d+_(\\w+_\\w+)_.*', replacement='\\1',files)
  patt <- sub("_Chip.+", "", ER)
  res <- paste(dir, files, sep="/")
  res <- as.list(res)
  names(res) <- ER
  return(res)
}


