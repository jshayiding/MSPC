##' pre-processing set of Chip-seq replicate & exclude all noise peak from samples
##'
##' we set up threshold for signal' significant value of each enriched region,
##' where extremely weakly enriched regions won't be processed, so we did purification on original input dataset;
##' Extremenly weakly enriched regions (A.K.A, background noise signals) are exported as standard BED file for the sake of giving clear biological evidence.
##' The purpose of exposting noise peak as BED file, we'll evaluate each Chip-seq replicates that bearing different output set with clear biological evidence.
##'
##' @title denoise_ERs
##' @param peakGRs list of Chip-seq replicate imported and all enriched regions stored in GRanges objects
##' @param tau.w threshold value for weakly enriched peak, all enrichred regions' pvalue above this thrshold, are considered background signal (A.K.A, noise peak)
##' @param .fileName user has option to name background signal by their own preference.
##' @param outDir user control where the exported BED files goes
##' @param verbose control whether the output is printed or not
##' @return GRangesList
##' @export
##' @importFrom rtracklayer export.bed
##' @importFrom stats setNames
##' @author Julaiti Shayiding

denoise_ERs <- function(peakGRs, tau.w= 1.0E-04, .fileName="", outDir=getwd(), verbose=FALSE, ...) {
  # check input param
  if (missing(peakGRs)) {
    stop("Missing required argument peakGRs, please choose imported Chip-seq replicates!")
  }
  stopifnot(inherits(peakGRs[[1L]], "GRanges"))
  stopifnot(length(peakGRs)>0)
  stopifnot(is.numeric(tau.w))
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  if(!dir.exists(outDir)) {
    dir.create(file.path(outDir))
    setwd(file.path(outDir))
  }
  res <- lapply(seq_along(peakGRs), function(x) {
    .gr <- peakGRs[[x]]
    .grNM <- names(peakGRs)[x]
    .drop <- .gr[.gr$p.value > tau.w]
    export.bed(.drop, sprintf("%s/%s.%s.bed", outDir, .fileName, .grNM))
    .keep <- .gr[.gr$p.value <= tau.w]
    return(.keep)
  })
  rslt <- setNames(res, names(peakGRs))
  return(rslt)
}

##' @example
##' total.ERs <- .denoise.ERs(myData, tau.w = 1.0E-04, .fileName = "noisePeak", outDir = "test/")
