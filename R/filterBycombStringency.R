##' Filter ERs by combined stringency test & create output set
##'
##' Main idea behind our method is that repeated evidence across Chip-seq replicates can compensate for a
##' lower significance in a single sample, which is implemented though the Fisher method. The significace of
##' ovelapping enriched regions is rigorously combined with the Fisher's method to obtain global Fisher score (A.K.A, combined p-value)
##' to evaluate combined stringency of all enriched regions and classified each peak as confirmed or discarded.
##'
##' @title filterBycombStringency
##'
##' @description perform combined stringency test
##'
##' @param ERs set of pre-processed enriched regions stored in GRanges, where background signal (A.K.A, noise peak) were pre-processed and won't involve further workflow
##' @param .hitLiist list of overlap hit that set of enriched regions comply minimum overlapping peak requirement.
##' @param cmbstrgThreshold combined stringency threshold against all enriched regions p-value, and each peaks are classified as confirmed or discarded peak accordingly.
##' @param isFisherPass logical vector that check whether all enriched regins are passed in Fisher' combined method. TRUE : return GRangesList of all enriched regions are classified as confirmed. FALSE: return GRangesList of all discarded peaks. both is required
##' @export
##' @importFrom XVector extractList
##' @importFrom stats setNames
##' @importFrom BiocGenerics Map
##' @author Julaiti Shayiding

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

#' @example
# .Confirmed.ERs <- filterBycombStringency(total.ERs, keepList, cmbstrgThreshold=1.0E-08, comb.p, isFisherPass = TRUE)
# .Discarded.ERs <- filterBycombStringency(total.ERs, keepList, cmbstrgThreshold=1.0E-08, comb.p, isFisherPass = FALSE)
