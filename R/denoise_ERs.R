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
#'
#' @param tau.w threshold value for weakly enriched peak,
#' all enrichred regions' pvalue above this thrshold,
#' are considered background signal (a.k.a, noise)
#'
#' @param fileName user has option to name background signal
#' by their own preference. by default name it as "noise"
#'
#' @param outDir only noise can be exported as standard BED file;
#' user has an option to control where the exported files goes
#'
#' @param verbose control whether the output is printed or not
#'
#' @return list of enriched regions stored in \link[GenomicRanges]{GRanges}
#' @export
#' @importFrom rtracklayer export.bed
#' @importFrom stats setNames
#' @importFrom methods hasArg
#' @author Julaiti Shayiding
#'
#' @examples
#' require(rtracklayer)
#' require(GenomicRanges)
#'
#' ## read peak file as GRanges object
#' files <- getPeakFile()[7:8]
#' grs <- readPeakFiles(files, pvalueBase=1L)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                         fileName = "noise", outDir = getwd())
#' ## Explore all stringent and weak enriched regions
#' total.ERs

denoise_ERs <- function(peakGRs=NULL,tau.w = 1.0E-04,fileName ="noise",
                        outDir = getwd(),verbose = FALSE) {
    # check input param
    if (!hasArg(peakGRs)) {
        stop("required argument peakGRs is missing,
             please choose imported Chip-seq replicates!")
    }
    if(!hasArg(tau.w)) {
        stop("permissive threshold value must be specified")
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
