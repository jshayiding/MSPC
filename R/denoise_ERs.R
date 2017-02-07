#' Data cleaning (exclude noise from input sample)
#'
#' Given the output of \link{readPeakFiles}, We set up
#' permissive p-value threshold \code{tau.w} for weakly
#' enriched regions associated with score and p-value.
#' We are interested in moderately enriched regions
#' (or weak peak), but extremely weakly enriched regions
#' must be removed from all candidate sample because it
#' won’t give any biological insight. Note that p-value
#' threshold can be tuned by users, so using different
#' threshold value will result in different filtering results.
#' Background signal (A.K.A, noise) also exported as standard
#' BED file by using \link[rtracklayer]{export.bed} for the sake of
#' evaluate each ChIP-Seq replicate that bearing different
#' output set with clear biological evidence.
#'
#' Passed to \link{runMSPC}
#'
#' @param peakGRs Chip-seq replicate imported and
#' all enriched regions stored in \code{GRanges} objects
#'
#' @param tau.w permissive p-value threshold for weakly
#' enriched peak, all enrichred regions' pvalue above
#' this thrshold, are considered background noise
#'
#' @param nmtab Character string that assigned as name to the noise
#'
#' @param dest.dir dirctory that exported files can be placed.
#' only noise can be exported as standard BED file;
#'
#' @param overwrite default FALSE indicating whether
#' existing files with identical name should be over-written.
#'
#' @return Peaks without background noise is
#' stoted in \link[GenomicRanges]{GRangesList}
#'
#' @export
#' @importFrom rtracklayer export.bed
#' @importFrom methods hasArg
#' @importFrom GenomicRanges GRangesList
#' @author Jurat Shahidin
#'
#' @examples
#' require(rtracklayer)
#' require(GenomicRanges)
#'
#' ## read peak file as GRanges object
#' files <- getPeakFile()[1:3]
#' grs <- readPeakFiles(files, pvalueBase=1L)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                         nmtab = "noise", dest.dir = getwd(), overwrite = TRUE)
#' ## Explore all stringent and weak enriched regions
#' total.ERs

denoise_ERs <- function(peakGRs,
                        tau.w = 1.0E-04,
                        nmtab ="noise",
                        dest.dir = tempdir(),
                        overwrite=FALSE) {
    # check input param
    if (class(peakGRs) != "GRangesList") {
        stop("input must be a GRangesList Object")
    }
    stopifnot(length(peakGRs)>0)
    stopifnot(is.numeric(tau.w))
    if (!dir.exists(dest.dir)) {
        if (!file.exists(dest.dir))
            dir.create(dest.dir)
        else
            stop("'dest.dir' exists but is not a directory")
    }
    gr <- unlist(peakGRs, use.names = FALSE)
    idx <- factor(rep(names(peakGRs), lengths(peakGRs)),
                  levels = names(peakGRs))
    Drop <- gr$p.value > tau.w
    noise <- split(gr[Drop], idx[Drop])
    if (overwrite) {
        for(i in seq_along(noise)) {
            filename <- paste(names(noise)[i], ".bed", sep = "")
            export.bed(noise[[i]], sprintf("%s/%s.%s", dest.dir, nmtab, filename))
        }
    } else {
        stop(paste("prevent overwriting files; please delete existing file:", dest.dir))
    }
    res <- GRangesList(split(gr[!Drop], idx[!Drop]))
    return(res)
}
