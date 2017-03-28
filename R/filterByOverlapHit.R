#' Check minimum overlapping peak requirement over overlap-hit list
#'
#' Given output of \link{peakOverlapping}, we need to check whether
#' each enriched regions comply minimum overlapping peak
#' requiremnet or not. \link{peakOverlapping} function give us
#' correct geometry of overlap peaks as list-like object.
#' Getting cardinality of overlapping ERs in parallel.
#' \code{min.c} is minimum overlapping peak requirement,
#' that determined by type of Chip-seq replicates, and number of
#' samples are chosen. Based on \code{min.c} parameter, we could
#' identify whether all enriched regions comply minimum overlapping
#' peak requirement. Where only ERs that comply the minimum
#' requirement, we could proceed our analysis For Fisher' combined test.
#' For the sake of efficiency, all ERs fullfill the requirement
#' kept in \link[IRanges]{IntegerList} objects.
#'
#' Passed to \link{filterBycombStringency}
#'
#' ERs set that failed from minimum overlapping peak requirement also
#' passed to \link{mergeDiscERs}
#'
#' @param .ovHit output of \link{peakOverlapping},
#' where all overlap-hit pair is in \link[IRanges]{IntegerList}
#'
#' @param peakset output of \link{denoise_ERs}. Set of Chip-seq
#' replicate imported and all peaks are stored in GRanges object,
#' where all background noise were excluded from
#' \link{denoise_ERs} function
#'
#' @param replicate.type A charcter vector used to select type of
#' Chip-seq replicate ( Biological / Technical replicate)
#'
#' @param isSuffOverlap logical vector that check whether ERs
#' comply minimum overlapping peak requirement.
#'
#' \code{TRUE}: return overlap hit list where all enriched regions
#' comply minimum overlapping peak requirement.
#'
#' \code{FALSE}: return all discarded enriched regions as GRange object
#'  due to failing to comply minimum overlapping peak requirement.
#'
#' @return
#' \code{isSuffOverlap} is \code{True}:
#' return list of ERs in \link[IRanges]{IntegerList} that comply with
#'  minimum overlapping peak requirement
#'
#' \code{isSuffOverlap} is \code{False} :
#' return list of \link[GenomicRanges]{GRanges} objects that
#' all discarded peaks due to failing for minimum requirement.
#'
#' @export
#' @importFrom S4Vectors Reduce
#' @importFrom S4Vectors lengths
#' @importFrom XVector extractList
#' @importFrom BiocGenerics Map
#' @importFrom methods hasArg
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(XVector)
#'
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
#' initDisc.ERs <- filterByOverlapHit(hit,
#'                                    peakset = total.ERs,
#'                                    replicate.type = "Biological",
#'                                    isSuffOverlap=FALSE)
#' # explore ERs set
#' lapply(keepList, head)
#' lapply(initDisc.ERs, head)

filterByOverlapHit <- function(.ovHit, peakset,
                               replicate.type=c("Biological", "Technical"),
                               isSuffOverlap= TRUE) {
    # check input param
    stopifnot(length(peakset)>0)
    stopifnot(inherits(peakset[[1]], "GRanges"))
    replicate.type = match.arg(replicate.type)
    if (!hasArg(peakset)) {
        stop("Missing required argument peakset, please
             choose the set of ERs without noise!")
    }
    if (!hasArg(.ovHit)) {
        stop("please choose the set of list of overlap hit index for ERs
             that comply minimum overlapping peak requirement!")
    }
    if(!hasArg(replicate.type)) {
        stop("replicate type must be specified")
    }
    min.c <- ifelse(replicate.type=="Biological",
                    length(peakset)-1,
                    length(peakset))
    cnt.ovHit <- as.matrix(Reduce('+', lapply(.ovHit, lengths)))
    if(isSuffOverlap==TRUE) {
        keepHit <- lapply(.ovHit, function(ele_) {
            keepMe <- sapply(cnt.ovHit, function(x) x >= min.c)
            res <- ele_[keepMe]
        })
    } else if(isSuffOverlap==FALSE){
        dropHit <- lapply(.ovHit, function(ele_) {
            droped <- sapply(cnt.ovHit, function(x) x < min.c)
            res <- ele_[droped]
        })
        rslt <- Map(unlist,
                    mapply(extractList, peakset, dropHit))
    }
}
