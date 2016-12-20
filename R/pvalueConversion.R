#' Data conversion for significance of peak signal
#'
#' Some data sources provides Chip-seq enriched regions (A.K.A, peaks)
#' without p-value under specific conditions.
#' \code{pvalueConversion} is the utility function to do data conversion of
#' peak' score as p-value. In standard BED file,
#' significant value of peak signal defined as Score column,
#' Thus we need to convert it
#' as p.value (- log(p.value), -10 log(p.value), -100 log(p.value))
#'
#' Passed to \link{readPeakFiles}
#'
#' @param x GRanges objects All Enriched regions (a.k.a, peaks)
#' are stored in GRanges object.
#'
#' @param pvalueBase User has option to select
#' p-value format (- log(p.value), -10 log(p.value), -100 log(p.value))
#'
#' @return GRanges
#' @export
#' @importFrom rtracklayer score
#  @importFrom rtracklayer mcols
#  @importFrom IRanges colnames
#' @author Julaiti Shayiding
#'
#' @examples
#' library(GenomicRanges)
#'
#' ## example GRanges object
#' grs <- GRangesList(
#'     bar=GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(1,6,16), c(4,12,23)),
#'         strand=Rle(c("*"),3), rangeName=c("a1", "a2", "a3"),
#'         score=c(22,6,13)),
#'     cat=GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(7,19,31), c(13,28,43)),
#'         strand=Rle(c("*"),3), rangeName=c("b3","b6","b7"),
#'         score=c(14,9,17))
#' )
#'
#' ## pvalue conversion
#' grs <- lapply(grs, pvalueConversion)
#'

pvalueConversion <- function(x, pvalueBase = 10L) {
    # input param checking
    stopifnot(class(x) == "GRanges")
    stopifnot(is.numeric(pvalueBase))
    # explore score of all features
    if(is.null(x$pvalue)){
        x$p.value <- 10^(score(x)/(- pvalueBase))
        #colnames(mcols(x))[3] <- "p.value"
    } else {
        x
    }
    res <- x
    return(res)
}
