#' Retrieve Fisher's global score through Fisher Method
#'
#' We assess the presence of overlapping enrichred regions across
#' multiple Chip-seq replicates. Therefore, the significance of
#' overlapping enriched regions are rigorously combined
#' with Fisher's method to obtain global Fisher score.
#' Using \link[metap]{sumlog} to get global Fisher'score.
#'
#' Passed to \link{filterBycombStringency}
#'
#' @param hitList output of \link{filterByOverlapHit} take \code{isSuffOverlap}
#' as \code{TRUE}. Only Enriched regions that comply with
#' minimum overlapping peak requirement can be further evaluated. We kept ERs
#' that fulfill minimum overlapping requirement in \link[IRanges]{IntegerList}.
#'
#' @param peakset output of \link{denoise_ERs}. Set of Chip-seq replicate
#' imported and all peaks are stored in GRanges object, where all
#' background noise won't be further processed.
#'
#' @return numeric vector, Global Fisher's score are in vector
#' @export
#' @importFrom metap sumlog
#' @importFrom XVector extractList
#' @importFrom rtracklayer as.data.frame
#' @importFrom IRanges as.matrix
#'
#' @author Julaiti Shayiding
#'
#' @examples
#' library(GenomicRanges)
#' library(rtracklayer)
#' library(XVector)
#' library(metap)
#'
#' ## Example Peaks in GRanges
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
#' ##Add p.value as metadata column
#' grs <- GRangesList("bar"=bar, "cat"=cat)
#' grs <- lapply(grs, pvalueConversion)
#'
#' ##Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                          fileName = "noise", outDir = getwd())
#'
#' ##peak overlapping
#' hit <- peakOverlapping(total.ERs, FUN=which.max)
#'
#' ##Retrieve pvalue of ERs that comply minimum overlapping peak requirement
#' keepList <- filterByOverlapHit(hit, peakset = total.ERs,
#'                                replicate.type = "Biological",
#'                                isSuffOverlap=TRUE)
#' \dontrun{
#' ## Global Fisher' score
#' comb.p <- Fisher_stats(hitList = keepList , peakset = total.ERs)
#' comb.p
#' }


Fisher_stats <- function(hitList, peakset) {
    # input param checking
    if (missing(peakset)) {
        stop("Missing required argument peakset, please
             choose the set of pre-processed peaks!")
    }
    if (missing(hitList)) {
        stop("Missing required argument hitList, please choose overlap
             hit list that comply minimum overlapping peak requirement!")
    }
    message("retrieve pvalue of peaks")
    .get.pvalue <- function(hit, obj) {
        # input param checking
        res <- extractList(obj$p.value, hit)
        return(res)
    }
    pval_List <- mapply(.get.pvalue, hitList, peakset)
    .helper.PV <- function(p.list) {
        res <- sapply(p.list, function(x) {
            out <- ifelse(length(x)>0,
                          x, 0.000000e+00)
        })
    }
    pval.TB <- as.data.frame(mapply(.helper.PV, pval_List))
    # pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
    # Npresent <- rowSums( !is.na(pval.TB) )
    # comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
    comb.pval <- suppressWarnings(
        .globSC <- apply(pval.TB[, 1:length(pval.TB)], 1, function(ro) {
            res <- sumlog(ro)$p
        })
    )
    comb.pval <- as.matrix(comb.pval)
    return(comb.pval)
}
