##' Import Chip-seq replicates and all ERs are stored in GRanges objects.
##'
##' readPeakFiles can read peak file in standard BED format and stored in GRanges object,
##' where several peak files (A.K.A, Chip-seq replicates) can be read simultaneously using lapply
##'
##' @title readPeakFile
##'
##' @description
##' read Chip-seq replicated and all enriched regiions are stored in GRanges object.
##'
##' @param peakFolder set of ChipSeq peak replicate in BED format
##' @param pvalueBase user has options to choose the p-value format(- log(p-value), -10 log(p-value), -100 log(p-value)). By default, use -10 log(p-value) ;
##' @param verbose logical that control whether the output is printed or not
##' @return GRanges object
##' @export
##' @importFrom rtracklayer import.bed
##' @importFrom stats setNames
##' @importFrom tools file_path_sans_ext
##' @author  Julaiti Shayiding

readPeakFiles <- function(peakFolder, pvalueBase = 1L, verbose=FALSE) {
  # input param checking
  if (missing(peakFolder)) {
    stop("Missing required argument peakFolder, please choose input Chip-seq replicates in BED file!")
  }
  if (verbose)
    cat(">> reading all peakfiles from given folder...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  stopifnot(length(peakFolder)>0)
  stopifnot(is.numeric(pvalueBase))
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(peakFolder, function(ele_) {
      .gr <- as(import.bed(ele_), "GRanges")
      if(is.null(.gr$p.value)) {
        .gr <- pvalueConversion(.gr, 1L)
      }
    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}
