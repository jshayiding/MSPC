#' Import Chip-seq replicates and all ERs are stored in GRanges objects.
#'
#' readPeakFiles can read peak file in standard BED format using \link[rtracklayer]{import.bed}
#' and stored in GRanges object, where several peak files can be read simultaneously using lapply
#'
#' Passed to \link{denoise_ERs}
#'
#' @param peakFolder set of ChipSeq peak replicate in BED format file
#' @param pvalueBase user has options to choose the p-value format
#' @param verbose logical that control whether the output is printed or not
#' @return GRanges object
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom methods as
#' @importFrom stats setNames
#' @author  Julaiti Shayiding

readPeakFiles <- function(peakFolder, pvalueBase = 10L, verbose=FALSE) {
  # input param checking
  if (missing(peakFolder)) {
    stop("Missing required argument peakFolder, please choose input Chip-seq replicates in BED file!")
  }
  if (verbose)
    cat(">> reading all peakfiles from given folder...\t\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  stopifnot(length(peakFolder)>0)
  stopifnot(is.numeric(pvalueBase))
  f.read <- lapply(peakFolder, function(ele_) {
    .gr <- as(import.bed(ele_), "GRanges")
    if(is.null(.gr$p.value)) {
      .gr <- pvalueConversion(.gr)
    }
  })
  return(f.read)
}
