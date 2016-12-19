#' Check minimum overlapping peak requirement over overlap-hit list
#'
#' Given output of \link{peakOverlapping}, we need to check whether
#' each enriched regions comply minimum overlapping peak requiremnet or not.
#' \link{peakOverlapping} function give us correct geometry of overlap peaks as
#' list-like object. Getting cardinality of overlapping ERs in parallel.
#' \code{min.c} is minimum overlapping peak requirement, that determined by
#' type of Chip-seq replicates, and number of samples are chosen.
#' Based on \code{min.c} parameter, we could identify whether all enriched regions
#' comply minimum overlapping peak requirement. Where only ERs that
#' comply the minimum requirement, we could proceed our analysis For
#' Fisher' combined test. For the sake of efficiency, all ERs
#' fullfill the requirement kept in \link[IRanges]{IntegerList} objects.
#'
#' Passed to \link{filterBycombStringency}
#'
#' @param .ovHit output of \link{peakOverlapping},
#' where all overlap-hit pair is in \link[IRanges]{IntegerList}
#'
#' @param peakset output of \link{denoise_ERs}. Set of Chip-seq replicate imported
#' and all peaks are stored in GRanges object,
#' where all background noise were excluded from \link{denoise_ERs} function
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
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(XVector)
#'
#' ## example peaks in GRanges object
#' bar=GRanges(
#'   seqnames=Rle("chr1", 3),ranges=IRanges(c(12,21,37), c(14,29,45)),
#'   strand=Rle(c("*"),3), rangeName=c("a1", "a2", "a3"), score=c(22, 6,13))
#'
#' cat=GRanges(
#'   seqnames=Rle("chr1", 6),ranges=IRanges(c(5,12,16,21,37,78),
#'   c(9,14,19,29,45,84)),
#'   strand=Rle(c("*"),6), rangeName=c("b1", "b2","b3", "b4", "b6", "b7"),
#'   score=c(12, 5, 11, 8, 4, 3))
#'
#' ## Add p.value as metadata column
#' grs <- GRangesList("bar"=bar, "cat"=cat)
#' grs <- lapply(grs, pvalueConversion)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                          fileName = "noise", outDir = getwd())
#' ## find overlapping
#' hit <- peakOverlapping(total.ERs, FUN=which.max)
#'
#' ## check whether each ERs comply with minimum overlapping peak requirement
#' keepList <- filterByOverlapHit(hit, peakset = total.ERs,
#'                                replicate.type = "Biological",
#'                                isSuffOverlap=TRUE)
#' initDisc.ERs <- filterByOverlapHit(hit, peakset = total.ERs,
#'                                    replicate.type = "Biological",
#'                                    isSuffOverlap=FALSE)
#'

filterByOverlapHit <- function(.ovHit, peakset, replicate.type=c("Biological", "Technical"), isSuffOverlap= TRUE) {
  # check input param
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  if (missing(peakset)) {
    stop("Missing required argument peakset, please choose the set of pre-processed peaks!")
  }
  if (missing(.ovHit)) {
    stop("Missing required argument .ovHit, please choose overlapHit index from previous workflow!")
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
