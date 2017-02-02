#' Exclude all background signal (a.k.a, noise) from input Chip-seq samples
#'
#' Given the output of \link{readPeakFiles}, this function
#' exclude all bacckground noise from input Chip-seq replicates.
#' using permissive threshold tau.w for signal' significant
#' value of each enriched region. Extremely weakly enriched regions
#' won't be processed and excluded from all sample.
#' Background signal (a.k.a, noise) also exported as standard BED file
#' by using \link[rtracklayer]{export.bed} for the sake of
#' evaluate each Chip-seq replicate that bearing
#' different output set with clear biological evidence.
#'
#' Passed to \link{runMSPC}
#'
#' @param peakGRs list of Chip-seq replicate imported and
#' all enriched regions stored in \code{GRanges} objects
#'
#' @param tau.w threshold value for weakly enriched peak,
#' all enrichred regions' pvalue above this thrshold,
#' are considered background signal (a.k.a, noise)
#'
#' @param noiLab user has option to name background signal
#' by their own preference. by default name it as "noise"
#'
#' @param outDir only noise can be exported as standard BED file;
#' user has an option to control where the exported files goes
#'
#' @param overwrite logical whether noise files are existed or not
#'
#' @param verbose control whether the output is printed or not
#'
#' @return list of enriched regions stored in \link[GenomicRanges]{GRanges}
#' @export
#' @importFrom rtracklayer export.bed
#' @importFrom methods hasArg
#' @author Jurat Shahidin
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
#'                         noiLab = "noise", outDir = getwd(), overwrite=FALSE)
#' ## Explore all stringent and weak enriched regions
#' total.ERs


denoise_ERs <- function(peakGRs,
                        tau.w = 1.0E-04,
                        noiLab ="noise",
                        overwrite=FALSE,
                        outDir = tempdir()) {
    # check input param
    if (class(peakGRs) != "GRangesList") {
        stop("input must be a GRangesList Object")
    }
    stopifnot(length(peakGRs)>0)
    stopifnot(is.numeric(tau.w))
    if(!dir.exists(outDir)) {
        dir.create(file.path(outDir))
        setwd(file.path(outDir))
    }
    gr <- unlist(peakGRs, use.names = FALSE)
    idx <- factor(rep(names(peakGRs), lengths(peakGRs)),
                  levels = names(peakGRs))
    Drop <- gr$p.value > tau.w
    noise <- split(gr[Drop], idx[Drop])
    if (!overwrite)
        stop("noise peak files are already exist")
    for(i in seq_along(noise)) {
        filename <- paste(names(noise)[i], ".bed", sep = "")
        export.bed(noise[[i]], sprintf("%s/%s.%s", outDir, noiLab,filename))
    }
    res <- split(gr[!Drop], idx[!Drop])
    return(res)
}
