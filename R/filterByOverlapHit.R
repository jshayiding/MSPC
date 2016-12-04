##' Filter ERs by overlap-hit index list
##'
##' Once we obtained list of overlap-hit index where all possible overlap pair is included, we need to
##' further investigate hit index with our defined method. In previous workflow, `peakOverlapping` function give us correct geometry of overlap peaks as list-like object by parallel.
##' we need to do vector sum for getting overall overlaped peak numbers and proceed our evaluation. Where all enriched regions comply minimum overlap peak requirement,
##' are returned as overlap hit list and used for next workflow, while all enriched regions fail to comply the minimum overlap peak requirement, are returned as GRanges objects.
##' Due to further need of evaluating all enriched regions in different level in order to give clear Biological evidence, we decide to expand those overlap hit index as GRanges objects.
##'
##' @title filterByOverlapHit
##' @param .ovHit list of overlap hit index. Once we obtained overlap hit list by calling `peakOverlapping` function in previous workflow, resuled hit-list ready to be used.
##' @param peakset set of Chip-seq replicate imported and all peaks are stored in GRanges object.
##' @param replicate.type A charcter vector used to select type of Chip-seq replicate ( Biological / Technical replicate)
##' @param isSuffOverlap logical vector that check whether comply minimum overlapping peak requirement. TRUE: return overlap hit list where all enriched regions comply minimum overlapping peak requirement. FALSE : return all discarded enriched regions as GRange object due to failing to comply minimum overlapping peak requirement. both is required
##' @return list-like object
##' @export
##' @importFrom S4Vectors Reduce
##' @importFrom S4Vectors lengths
##' @importFrom XVector extractList
##' @importFrom BiocGenerics Map
##' @author Julaiti Shayiding

filterByOverlapHit <- function(.ovHit, peakset, replicate.type=c("Biological", "Technical"),
                               isSuffOverlap= c(TRUE, FALSE), verbose=FALSE, ...) {
  # check input param
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(missing(isSuffOverlap))
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
  cnt.ovHit <- as.matrix(Reduce('+', lapply(hit.List, lengths)))
  if(isSuffOverlap==TRUE) {
    keepHit <- lapply(hit.List, function(ele_) {
      keepMe <- sapply(cnt.ovHit, function(x) x >= min.c)
      res <- ele_[keepMe]
    })
  } else if(isSuffOverlap==FALSE){
    dropHit <- lapply(hit.List, function(ele_) {
      droped <- sapply(cnt.ovHit, function(x) x < min.c)
      res <- ele_[droped]
    })
    rslt <- Map(unlist,
                mapply(extractList, peakset, dropHit))
  }
}

##' @example
## keepList <- .filterByOverlapHit(Hit, peakset = total.ERs, replicate.type = "Biological", isSuffOverlap=TRUE)
## initDisc.ERs <- filterByOverlapHit(Hit, peakset = total.ERs, replicate.type = "Biological", isSuffOverlap=FALSE)
