#' load input peak bed files
#'
#' Input peak bed files can be detected by system.file(),
#' and let package example data available for use.
#'
#' @return list of Chip-seq replicates' name
#' @export
#' @importFrom IRanges gsub
#' @author Jurat Shahidin
#'
#' @examples
#' ## get bed file
#' files <- getPeakFile()

getPeakFile <- function() {
    dir <- system.file("extdata", package="MSPC")
    files <- list.files(dir, recursive = TRUE, full.names = TRUE)
    # ChipPeak <- gsub(pattern='wgEncode\\d+_(\\w+_\\w+)_.*',
    #                  replacement='\\1',files)
    # ChipPeak <- sub("_Chip.+", "", ChipPeak)
    # res <- paste(dir, files, sep="/")
    # res <- as.list(res)
    # names(res) <- ChipPeak
    setNames(files, basename(files))
}
