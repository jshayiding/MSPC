#' Merge discarded ERs from different level permissive threshold
#'
#' @param discER_A discarded ERs that failed from minimum
#' overlapping peak requirement
#'
#' @param discER_B discarded ERs that failed from combined
#' stringency test
#'
#' @return list of discarded ERs that failed from minimum
#' overlapping peak requirement and Fisher's combine test,
#' all discarded ERs are stored in \link[GenomicRanges]{GRanges} object
#'
#' @export
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#'
#' ## example peak
#'
#' initDiscERs <- GRangesList(
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(7,19,31), c(13,28,43)),
#'         strand = Rle(c("*"),3), rangeName=c("b3","b6","b7"),
#'         score=c(14,9,17)),
#'     bar = GRanges(
#'         seqnames=Rle("chr1",2),ranges=IRanges(c(1,16), c(4,23)),
#'         strand = Rle(c("*"),2), rangeName=c("a1","a3"),
#'         score=c(22,13))
#' )
#'
#' fisherDiscERs <- GRangesList(
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 2),ranges=IRanges(c(15,47), c(18,55)),
#'         strand = Rle(c("*"),2), rangeName=c("b4","b9"), score=c(1,6)),
#'     bar = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(6,25,40), c(12,33,49)),
#'         strand = Rle(c("*"),3), rangeName=c("a2","a5","a8"), score=c(3,2,4))
#' )
#'
#' ## merge
#' DiscardedERs <- mergeDiscERs(initDiscERs,fisherDiscERs)

mergeDiscERs <- function(discER_A, discER_B) {
    # sanity check
    if(missing(discER_A) | missing(discER_B)) {
        stop("please select discarded ERs that required to be merged")
    }
    res <- suppressWarnings(
        mapply(c, discER_A, discER_B)
    )
    return(res)
}
