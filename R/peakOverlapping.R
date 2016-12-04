##' Asess peak overlapping across multiple Chip-seq replicate
##'
##' we assess the presense of overlapping enriched regions across multiple Chip-seq replicates.
##' Each enriched regions in current replicates (A.K.A, chosen replicates) is evaluated
##' with the support of rest of Chip-seq replicates in input dataset for overlapping enriched regions by element-wise.
##' Due to processing each genomic region (all peak intervals are stored in GRanges objects) by element-wise to find overlap is quite inefficient,
##' peakOverlapping function efficiently vectorize overlapping peaks and retrieving overlapped regions from multiple Chip-seq replicates as list-like vector
##' where overlap position hit index gives us correct geomety of peak overlapping. It may happen that an enriched region(A.K.A, peak or ER) from current Chip-seq replicate overlap multiple ERs from other supported replicates.
##' peakOverlapping function retrieve only one overlapping peak from each supported Chip-seq replicate when multiple overlapping regions were detected.
##'
##' @title peakOverlapping
##'
##' @description
##' asess peak overlapping across multiple sample and capture correct geometry of overlaping ERs
##'
##' @param peakset set of Chip-seq replicate imported and all peaks are stored in GRanges object, where all background signal were pre-processed and won't involve in further downstream analysis.
##' @param FUN user has options to keep most stringent(with loest p-value) or least stringent(highest p-value) peak if multiple overlapping were detected. By default, keep most stringent overlaped peak.
##' @return list-like object
##' @export
##' @importFrom GenomicRanges findOverlaps
##' @importFrom XVector extractList
##' @importfrom rtracklayer score
##' @importFrom IRanges which.max
##' @importFrom IRanges which.min
##' @importFrom S4Vectors DataFrame
##' @importFrom IRanges as.matrix
##' @author Julaiti Shayiding

peakOverlapping <- function(peakset, FUN=which.max, ...) {
  # input param checking
  if (missing(peakset)) {
    stop("Missing required argument peakset, please choose the set of pre-processed peaks!")
  }
  if (missing(FUN)) {
    stop("Missing required argument FUN, please specify the peak type to be selected from multiple overlapping ERs!")
  }
  stopifnot(inherits(peakset[[1]], "GRanges"))
  res <- list()
  for(i in seq_along(peakset)) {
    que <- peakset[[i]]
    queHit <- as(findOverlaps(que), "List")
    supHit <- lapply(peakset[- i], function(ele_) {
      ans <- as(findOverlaps(que, ele_), "List")
      out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
      out.idx0 <- out.idx0[!is.na(out.idx0),]
      ans <- ans[out.idx0]
      ans
    })
    res[[i]] = DataFrame(c(list(que=queHit), sup=supHit))
    names(res[[i]]) = c(names(peakset[i]),names(peakset[- i]))
  }
  rslt <- lapply(res, function(x) as.matrix(x[names(res[[1]])]))
  rslt <- DataFrame(rbind(rslt[[1]],
                          unique(do.call("rbind", rslt[2: length(rslt)]))))
  rslt <- lapply(rslt, function(x) as(x, "CompressedIntegerList"))
  return(rslt)
}
