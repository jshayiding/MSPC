#' Combined stringency test
#'
#' Main idea behind our method is that repeated evidence across Chip-seq replicates can compensate for a
#' lower significance in a single sample, which is implemented though the Fisher method. The significace of
#' ovelapping enriched regions is rigorously combined with the Fisher's method to obtain global Fisher score (A.K.A, combined p-value)
#' However, we need to check whether each ERs fulfill requirement of combined stringency test.
#' Global Fisher' socre of enriched regions below combined stringency threshold, are rescued, we call this as confirmed peak.
#' Peaks are discarded due to failing for Fisher's combined test.
#'
#' Passed to \link{FDR_stats}
#'
#' @param ERs output of \link{denoise_ERs}, set of enriched regions stored in GRanges, where background signal won't be further evaluated.
#' @param .hitList output of \link{filterByOverlapHit}, take \code{isSuffOverlap} as \code{TRUE},
#' only all enriched regions that comply minimum overlapping peak requirements are further evaluated.
#' @param cmbstrgThreshold combined stringency threshold against all enriched regions p-value,
#' we could identify whether ERs fulfill stringency test, and result can be set of confirmed ERs, and discarded ERs respectively.
#' @param isFisherPass logical vector that check whether all enriched regins are passed in Fisher' combined method.
#' \code{TRUE}: return all enriched regions that fulfill Fisher's combined test, classified as \code{confirmedERs} in \link[GenomicRanges]{GRangesList}
#' \code{FALSE} : return all ERs that failing for combined stringency test, classified as \code{discardedERs} in \link[GenomicRanges]{GRangesList}
#'
#' @export
#' @importFrom XVector extractList
#' @importFrom stats setNames
#' @importFrom BiocGenerics Map
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#'
#' ## Prepare example peak Interval in GRanges objects
#' bar = GRanges(seqnames=Rle("chr1", 3), ranges=IRanges(c(12,21,37), c(14,29,45)),
#'               strand = Rle(c("*"),3), rangeName=c("a1", "a2", "a3"), score=c(22, 6,13))
#' cat = GRanges(seqnames=Rle("chr1", 6), ranges=IRanges(c(5,12,16,21,37,78), c(9,14,19,29,45,84)),
#'               strand = Rle(c("*"),6),
#'               rangeName=c("b1", "b2","b3","b4", "b6", "b7"), score=c(12, 5, 11, 8, 4, 3))
#' ## Add p.value as metadata column
#' grs <- GRangesList("bar"=bar, "cat"=cat)
#' grs <- lapply(grs, pvalueConversion)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                          fileName = "noise", outDir = getwd())
#'
#' ## peak overlapping
#' hit <- peakOverlapping(total.ERs, FUN=which.max)
#'
#' ## Retrieve pvalue of ERs that comply minimum overlapping peak requirement
#' keepList <- filterByOverlapHit(hit, peakset = total.ERs,
#'                                replicate.type = "Biological",
#'                                isSuffOverlap=TRUE)
#' \dontrun{
#' ## Get global Fisher's score
#' comb.p <- Fisher_stats(hitList = keepList, peakset = total.ERs)
#'
#' ## check whether ERs fulfill combined stringency test
#' Confirmed.ERs <- filterBycombStringency(total.ERs, keepList,
#' cmbstrgThreshold=1.0E-08, comb.p, isFisherPass = TRUE)
#' Discarded.ERs <- filterBycombStringency(total.ERs, keepList,
#' cmbstrgThreshold=1.0E-08, comb.p, isFisherPass = FALSE)
#' }
#'

filterBycombStringency <- function(ERs, .hitList, cmbstrgThreshold=1.0E-08 ,isFisherPass=TRUE) {
  # input param checking
  if (missing(ERs)) {
    stop("Missing required argument ERs, please choose the set of pre-processed peaks!")
  }
  if (missing(.hitList)) {
    stop("please choose the set of overlap hit index that comply minimum overlapping peak requirement!")
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
  .expandAsGR[[1L]] <- unique(.expandAsGR[[1L]])
  .expandAsGR <- setNames(.expandAsGR, names(ERs))
  return (.expandAsGR)
}
