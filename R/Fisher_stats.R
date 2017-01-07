#' Retrieve Fisher's global score through Fisher Method
#'
#' We assess the presence of overlapping enrichred regions across
#' multiple Chip-seq replicates. Therefore, the significance of
#' overlapping enriched regions are rigorously combined
#' with Fisher's method to obtain global Fisher score.
#' Using \link{fisherCmbp} to get global Fisher'score.
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
#' @importFrom rtracklayer as.data.frame
#' @importFrom IRanges as.matrix
#'
#' @author Julaiti Shayiding
#'
#' @examples
#' library(GenomicRanges)
#' library(rtracklayer)
#'
#' # read peak file as GRanges object
#' files <- getPeakFile()[7:8]
#' grs <- readPeakFiles(files, pvalueBase=1L)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                         fileName = "noise", outDir = getwd())
#' ## find peak overlapping
#' hit <- peakOverlapping(total.ERs, FUN=which.max)
#'
#' ## check whether each ERs comply with minimum overlapping peak requirement
#' keepList <- filterByOverlapHit(hit, peakset = total.ERs,
#'                                replicate.type = "Biological",
#'                                isSuffOverlap=TRUE)
#'
#' ## Global Fisher' score
#' comb.p <- Fisher_stats(hitList = keepList , peakset = total.ERs)
#'

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
    pval_List <- mapply(get.pvalue, hitList, peakset)
    .helper.PV <- function(p.list) {
        res <- sapply(p.list, function(x) {
            out <- ifelse(length(x)>0,
                          x, 0.000000e+00)
        })
    }
    pval.TB <- as.data.frame(mapply(.helper.PV, pval_List))
    comb.pval <- suppressWarnings(
        .globSC <- apply(pval.TB[, 1:length(pval.TB)],1, function(ro) {
            fisherCmbp(ro)
        })
    )
    comb.pval <- as.matrix(comb.pval)
    return(comb.pval)
}
