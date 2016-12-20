#' Exclude all background signal (a.k.a, noise) from input Chip-seq samples
#'
#' Give the output of \link{readPeakFiles}, this function
#' exclude all bacckground noise from input Chip-seq replicates.
#' using permissive threshold tau.w for signal' significant
#' value of each enriched region. Extremely weakly enriched regions
#' won't be processed and excluded from input Chip-seq replicates.
#' Background signal (a.k.a, noise) also exported as standard BED file
#' by using \link[rtracklayer]{export.bed} for the sake of
#' evaluate each Chip-seq replicate that bearing
#' different output set with clear biological evidence.
#'
#' Passed to \link{peakOverlapping}
#'
#' @param peakGRs list of Chip-seq replicate imported and
#' all enriched regions stored in \code{GRanges} objects
#' @param tau.w threshold value for weakly enriched peak,
#' all enrichred regions' pvalue above this thrshold,
#' are considered background signal (a.k.a, noise)
#' @param fileName user has option to name background signal
#' by their own preference.
#' @param outDir user control where the exported BED files goes
#' @param verbose control whether the output is printed or not
#'
#' @return list of enriched regions stored in \link[GenomicRanges]{GRanges}
#' @export
#' @importFrom rtracklayer export.bed
#' @importFrom stats setNames
#' @author Julaiti Shayiding
#'
#' @examples
#' require(rtracklayer)
#' require(GenomicRanges)
#'
#' ## example peaks in GRanges object
#' bar=GRanges(
#'     seqnames=Rle("chr1", 3),
#'     ranges=IRanges(c(12,21,37), c(14,29,45)), strand=Rle(c("*"),3),
#'     rangeName=c("a1", "a2", "a3"), score=c(22, 6,13)
#' )
#'
#' cat=GRanges(
#'     seqnames=Rle("chr1", 6),
#'     ranges=IRanges(c(5,12,16,21,37,78), c(9,14,19,29,45,84)),
#'     strand=Rle(c("*"),6), rangeName=c("b1", "b2","b3", "b4", "b6", "b7"),
#'     score=c(12, 5, 11, 8, 4, 3)
#' )
#'
#' ## Add p.value as metadata column
#' grs <- GRangesList("bar"=bar, "cat"=cat)
#' grs <- lapply(grs, pvalueConversion)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                         fileName = "noise", outDir = getwd())
#' ## Explore all stringent and weak enriched regions
#' total.ERs

denoise_ERs <- function(peakGRs=NULL,tau.w = 1.0E-04,fileName ="",
                        outDir = getwd(),verbose = FALSE) {
    # check input param
    if (missing(peakGRs)) {
        stop("Missing required argument peakGRs,
             please choose imported Chip-seq replicates!")
    }
    stopifnot(inherits(peakGRs[[1L]], "GRanges"))
    stopifnot(length(peakGRs)>0)
    stopifnot(is.numeric(tau.w))
    if (verbose) {
        cat(">> filter out all background noise peaks whose pvalue
            above threshold \t\t", format(Sys.time(),
                                          "%Y-%m-%d %X"), "\n")
    }
    if(!dir.exists(outDir)) {
        dir.create(file.path(outDir))
        setwd(file.path(outDir))
    }
    res <- lapply(seq_along(peakGRs), function(x) {
        gr <- peakGRs[[x]]
        grNM <- names(peakGRs)[x]
        drop <- gr[gr$p.value > tau.w]
        export.bed(drop, sprintf("%s/%s.%s.bed", outDir, grNM, fileName))
        keep <- gr[gr$p.value <= tau.w]
        return(keep)
    })
    rslt <- setNames(res, names(peakGRs))
    return(rslt)
}
