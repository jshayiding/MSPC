#' Asess peak overlapping across multiple Chip-seq replicate
#'
#' We assess the presense of overlapping enriched regions
#' across multiple Chip-seq replicates in parallel.
#' Each enriched regions in current replicates
#' (A.K.A, chosen replicates) is evaluated with
#' the support of rest of Chip-seq replicates for
#' overlapping enriched regions by element-wise.
#' Finding set of overlapping enriched regions
#' across multiple replicates can give rise to ambiguities.
#' Therefore, global approach can depends on
#' order of Chip-seq replicates when proceed finding overlap.
#' Due to processing each genomic region (all peak intervals
#' are stored in GRanges objects) by element-wise
#' to find overlap is quite inefficient, \code{peakOverlapping}
#' function efficiently vectorize overlapping peaks
#' and retrieve overlapped regions from multiple Chip-seq replicates
#' as list-like vector. Using List to get correct gemotry of
#' overlap-hit index. Using \link[GenomicRanges]{findOverlaps} efficiently
#' identify overlapping genomic regions, where overlap position hit index
#' gives us correct geomety of peak overlapping.
#' It may happen that an enriched region(A.K.A, peak) from current
#' Chip-seq replicate overlap multiple ERs from other supported replicates.
#' peakOverlapping retrieve only one overlapping peak from each supported
#' Chip-seq replicate when multiple overlapping regions were detected.
#'
#' Passed to \link{filterByOverlapHit}
#'
#' @param peakset the output of \link{denoise_ERs}.
#' set of Chip-seq replicate imported and all peaks are stored in
#' GRanges object, where all background signal were excluded.
#'
#' @param FUN user has options to keep most stringent(with loest p-value)
#' or least stringent(highest p-value) peakif multiple overlapping ERs
#' were detected. By default, keep most stringent overlaped ER.
#'
#' @return overlap-hit list in \link[IRanges]{IntegerList} object
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom XVector extractList
#' @importFrom rtracklayer score
#' @importFrom IRanges which.max
#' @importFrom IRanges which.min
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges as.matrix
#' @importFrom methods as
#' @importFrom methods hasArg
#' @author Julaiti Shayiding
#'
#' @examples
#'
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(S4Vectors)
#' require(XVector)
#'
#' ## read peak file as GRanges object
#' files <- getPeakFile()[7:8]
#' grs <- readPeakFiles(files, pvalueBase=1L)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                         fileName = "noise", outDir = getwd())
#' ## find peak overlapping
#' hit <- peakOverlapping(total.ERs, FUN=which.max)
#' ## explore overlap-hit index table
#' print(hit)

#' @note Thanks to Martin Morgan's revision on this implementation

peakOverlapping <- function(peakset, whichType=c("max","min")) {
    # sanity input param checking
    if (class(peakset) != "GRangesList") {
        stop("input must be a GRangesList Object")
    }
    whichType = match.arg(whichType)
    stopifnot(length(peakset)>0)
    # flatten out input peak list
    #allPeaks <- unlist(GRangesList(peakset))
    allPeaks <- unlist(peakset, use.names = TRUE)
    message("find peak overlapping across multiple sample simultanously")
    hitIdx <- findOverlaps(allPeaks)
    new_mcols <- DataFrame(
        query = queryHits(hitIdx),
        subject = subjectHits(hitIdx),
        score = allPeaks$score[subjectHits(hitIdx)],
        p.value = allPeaks$p.value[subjectHits(hitIdx)]
    )
    new_mcols <- new_mcols[order(new_mcols$score),]
    len_vect <- rep(seq_along(peakset), lengths(peakset))
    bindHit <- cbind(new_mcols$query, len_vect[new_mcols$subject])
    if(whichType == "max")
        maxHit <- !duplicated(bindHit, fromLast=TRUE)
    maxHit <- !duplicated(bindHit, fromLast=FALSE)
    res <- new_mcols[maxHit,]
    return(res)
}
