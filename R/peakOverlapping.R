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
#' @author Julaiti Shayiding
#'
#' @examples
#'
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(S4Vectors)
#' require(XVector)
#'
#' ## example peak interval in GRanges objects
#' bar=GRanges(
#'     seqnames=Rle("chr1", 3),ranges=IRanges(c(12,21,37), c(14,29,45)),
#'     strand=Rle(c("*"),3), rangeName=c("a1", "a2", "a3"), score=c(22, 6,13)
#' )
#'
#' cat=GRanges(
#'     seqnames=Rle("chr1", 6),ranges=IRanges(c(5,12,16,21,37,78),
#'     c(9,14,19,29,45,84)),
#'     strand=Rle(c("*"),6), rangeName=c("b1", "b2","b3", "b4", "b6", "b7"),
#'     score=c(12, 5, 11, 8, 4, 3)
#' )
#'
#' foo=GRanges(
#'     seqnames=Rle("chr1", 7),ranges=IRanges(c(2,8,18,35, 42,59,81),
#'     c(6,13,27,40,46,63,114)), strand=Rle(c("*"),7),
#'     rangeName=c("c1", "c2", "c3", "c4","c5","c8","c11"),
#'     score=c(2.1, 3, 5.1, 3.5, 7, 12, 10)
#' )
#' ## create GRangesList
#' grs <- GRangesList("bar"=bar, "cat"=cat, "foo"=foo)
#'
#' ## add p.value as metadata column
#' grs <- lapply(grs, pvalueConversion)
#'
#' ## Exclude all background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                          fileName = "noise", outDir = getwd())
#' ## find peak overlapping
#' hit <- peakOverlapping(total.ERs, FUN=which.max)
#'

peakOverlapping <- function(peakset, FUN=which.max) {
    # input param checking
    if (missing(peakset)) {
        stop("Missing required argument peakset,
             please choose the set of pre-processed peaks!")
    }
    if (missing(FUN)) {
        stop("Missing required argument FUN, please specify the
             peak type to be selected from multiple overlapping ERs!")
    }
    stopifnot(inherits(peakset[[1]], "GRanges"))
    res <- list()
    for(i in seq_along(peakset)) {
        que <- peakset[[i]]
        queHit <- as(findOverlaps(que), "List")
        supHit <- lapply(peakset[- i], function(ele_) {
            ans <- as(findOverlaps(que, ele_), "List")
            out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
            out.idx0 <- out.idx0[!is.na(out.idx0),]
            ans <- ans[out.idx0]
        })
        res[[i]] = DataFrame(c(list(que=queHit), sup=supHit))
        names(res[[i]]) = c(names(peakset[i]),names(peakset[- i]))
    }
    rslt <- lapply(res, function(x) as.matrix(x[names(res[[1]])]))
    rslt <- DataFrame(rbind(rslt[[1]],
                            unique(do.call("rbind", rslt[2: length(rslt)]))))
    rslt <- lapply(rslt, function(x) as(x, "CompressedIntegerList"))
    return(rslt)
}
