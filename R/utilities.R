#' load input peak bed files
#'
#' Input peak bed files can be detected by system.file(),
#' and let package example data available for use.
#'
#' @return list of Chip-seq replicates' name
#' @export
#' @importFrom IRanges gsub
#' @author Julaiti Shayiding
#'
#' @examples
#' ## get bed file
#' files <- getPeakFile()

getPeakFile <- function() {
    dir <- system.file("extdata", package="MSPC")
    files <- list.files(dir)
    ChipPeak <- gsub(pattern='wgEncode\\d+_(\\w+_\\w+)_.*',
                     replacement='\\1',files)
    ChipPeak <- sub("_Chip.+", "", ChipPeak)
    res <- paste(dir, files, sep="/")
    res <- as.list(res)
    names(res) <- ChipPeak
    return(res)
}

#' helper function to get combine p-value
#'
#' @param p pvalue vector
#' @return combined pvalue
#' @export
#' @importFrom stats pchisq
#' @importFrom stats na.omit
#' @author Julaiti
#'
#' @examples
#' ## pvalue vector
#' df <- data.frame(v1=c(1e-02,1e-11, 1e-09),
#'                    v2=c(1e-18,1e-07, 1e-13))
#' ## combined pvalue
#' cmbP <- apply(df,1,fisherCmbp)
#'

fisherCmbp <- function (p) {
    # sanity check
    if (any(is.na(p))) {
        warning("some p-values are missing")
        p <- na.omit(p)
    }
    rng.pv <- (p > 0) & (p <= 1)
    lg.pv <- log(p[rng.pv])
    chisq <- (-2) * sum(lg.pv)
    df <- 2 * length(lg.pv)
    res <- pchisq(chisq, df, lower.tail = FALSE)
    return(res)
}


#' retrieve pvalue of ERs that comply minimum overlapping peak requirement
#'
#' @param hit overlap hit
#' @param gr GRanges object
#' @return pvalue
#' @export
#' @importFrom XVector extractList
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(XVector)
#'
#' ## example peak interval
#' grs <- GRangesList(
#'     bar = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(6,25,40), c(12,33,49)),
#'         strand = Rle(c("*"),3), rangeName=c("a2","a5","a8"),
#'         score=c(3,12,4),p.value=c(1e-03,1e-12,1e-04)),
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(15,19,47), c(18,28,55)),
#'         strand = Rle(c("*"),3), rangeName=c("b4","b6","b9"),
#'         score=c(11,3,6),p.value=c(1e-11,1e-03,1e-06))
#' )
#'
#' ## find overlap
#' hit <- peakOverlapping(grs, which.max)
#'
#' ## explore pvalue
#' pv <- Map(get.pvalue, hit, grs)


get.pvalue <- function(hit, gr) {
    # input param checking
    stopifnot(class(gr)=="GRanges")
    stopifnot(inherits(hit, "CompressedIntegerList"))
    res <- extractList(gr$p.value, hit)
    return(res)
}
