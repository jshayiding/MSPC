#' Combined stringency test
#'
#' Main idea behind our method is that repeated evidence across
#' Chip-seq replicates can compensate for a lower significance in
#' a single sample, which is implemented though the Fisher method.
#' The significace of ovelapping enriched regions is rigorously combined
#' with the Fisher's method to obtain global Fisher score.
#' However, given the output of \link{filterByOverlapHit} take
#' \code{isSuffOverlap} as \code{TRUE}, we need to check whether
#' each ERs fulfill requirement of combined stringency test.
#' Global Fisher' socre of enriched regions below combined stringency
#' threshold, are rescued, we call this as confirmed peak;
#' Meanwhile Peaks are discarded, due to failing for Fisher's combined test.
#'
#' Passed to \link{FDR_stats}
#' ERs that failed from Fisher method also passed to \link{mergeDiscERs}
#'
#' @param ERs output of \link{denoise_ERs},
#' set of enriched regions stored in GRanges,
#' where background signal won't be further evaluated.
#'
#' @param .hitList output of \link{filterByOverlapHit},
#' take \code{isSuffOverlap} as \code{TRUE}, only all enriched regions
#' that comply minimum overlapping peak requirements are further evaluated.
#'
#' @param cmbstrgThreshold combined stringency threshold
#' against all enriched regions p-value, so we could identify
#' whether ERs fulfill stringency test,
#' and result can be set of confirmed ERs, and discarded ERs respectively.
#'
#' @param isFisherPass logical vector that check whether
#' all enriched regins are passed in Fisher' combined method.
#'
#' \code{TRUE}: return all enriched regions that fulfill Fisher's combined test,
#' classified as \code{confirmedERs} in \link[GenomicRanges]{GRangesList}
#'
#' \code{FALSE} : return all ERs that failing for combined stringency test,
#' classified as \code{discardedERs} in \link[GenomicRanges]{GRangesList}
#'
#' @return list of ERs in \link[GenomicRanges]{GRanges} object.
#' \code{isFisherPass} is \code{TRUE},
#' output is list of confirmed ERs thourgh Fisher's combined test.
#'
#' \code{isFisherPass} is \code{FALSE},
#' output is list of discarded ERs due to failing for combined stringency test.
#' ERs set that failed from Fisher method must be merged with ERs that
#' failed from minimum overlapping peak requirement, to do so,
#' using \link{mergeDiscERs}
#'
#' @export
#' @importFrom XVector extractList
#' @importFrom rtracklayer as.data.frame
#' @importFrom stats setNames
#' @importFrom BiocGenerics Map
#' @importFrom data.table setDT
#' @importFrom data.table .N
#' @importFrom data.table .I
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(XVector)
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
#' ## check whether ERs fulfill combined stringency test
#' Confirmed.ERs <- filterBycombStringency(total.ERs, keepList,
#' cmbstrgThreshold=1.0E-08, isFisherPass = TRUE)
#' fisherDiscERs <- filterBycombStringency(total.ERs, keepList,
#' cmbstrgThreshold=1.0E-08, isFisherPass = FALSE)
#'

filterBycombStringency <- function(ERs,.hitList,
                                   cmbstrgThreshold=1.0E-08,
                                   isFisherPass=TRUE){
    # input param checking
    if (missing(ERs)) {
        stop("Missing required argument ERs,
             please choose the set of pre-processed peaks!")
    }
    if (missing(.hitList)) {
        stop("please choose the set of overlap hit index
             that comply minimum overlapping peak requirement!")
    }
    stopifnot(is.numeric(cmbstrgThreshold))
    comb.p <- Fisher_stats(.hitList, ERs)
    if(isFisherPass==TRUE) {
        .filtHelper <- function(ele_) {
            keepMe <- sapply(comb.p, function(x) x<=cmbstrgThreshold)
            res <- ele_[keepMe]
        }
    } else if(isFisherPass==FALSE) {
        .filtHelper <- function(ele_) {
            drop_ <- sapply(comb.p, function(x) x > cmbstrgThreshold)
            res <- ele_[drop_]
        }
    }
    .hitIdx <- lapply(.hitList, .filtHelper)
    .expandAsGR <- Map(unlist,
                       mapply(extractList, ERs, .hitIdx))
    DF <- lapply(.expandAsGR, as.data.frame)
    # conditional duplicate remobal for each sample
    # res <- Map(function(x,y)
    #     setDT(x)[x[,
    #                .I[
    #                    if(.N > y) seq_len(pmax(y-1, 1))
    #                    else seq_len(.N)],
    #                .(name, score, p.value)]$V1],DF, 1:length(DF))
    res <- Map(function(x,y)
        setDT(x)[x[, .I[(1:.N) <= y] ,
                   .(name, score, p.value)]$V1], DF, 1:length(DF))
    asGR <- lapply(res, function(x) as(x, "GRanges"))
    rslt <- setNames(asGR, names(ERs))
    return (rslt)
}
