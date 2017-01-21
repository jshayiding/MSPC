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
#' @param ovHit output of \link{filterByOverlapHit} take \code{isSuffOverlap}
#' as \code{TRUE}. Only Enriched regions that comply with
#' minimum overlapping peak requirement can be further evaluated.
#' We used \link[S4Vectors]{DataFrame} to hold all metadata that we need to
#' evaluate for successive downstream analysis.
#'
#' @return DataFrame with combined pvalue
#' @export
#' @importFrom methods hasArg
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
#' comb.p <- Fisher_stats(ovHit = keepList)
#'

Fisher_stats <- function(ovHit) {
    # sanity input param checking
    if(!c("subject", "query", "p.value") %in% colnames(ovHit)) {
        stop("required columns are missing")
    }
    DF <- DataFrame(split(ovHit, ovHit$query))
    pvList <- as.matrix(DF$p.value)
    pvList[is.na(pvList)] <- 0
    message("getting combined pvalue")
    cmbp <- apply(pvList,1,fisherCmbp)
    DF$comb.pv <- cmbp
    return(DF)
}
