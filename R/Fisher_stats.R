##' Combine stringecy test to rigorously combined overlapping peaks
##'
##' we assess the presence of overlapping enrichred regions across multiple Chip-seq replicates.
##' Therefore, the significance of overlapping regions is rigorously combined with Fisher's method to obtain global Fisher score.
##' However, in next workflow, we are using Fisher combined p-value against combined stringency threshold
##' to evaluate combined stringency to all enriched regions in second level classification for all enriched regions.
##' We set up combined stringency threshold for all enriched regions and getting confirmed / discarded peak set accordingly.
##'
##' @title Fisher_stats
##'
##' @description
##' retrieve pvalue of ERs that comply minimum overlapping peak requirement and perfrom Fisher method
##'
##' @param .hitList overlap hit index for all enrichred regions that comply minimum overlapping peak requirement.
##' @param peakset set of Chip-seq replicate imported and all peaks are stored in GRanges object, where all background noise pre-processed and won't involve in further downstream analysis.
##' @return numeric vector
##' @export
##' @importFrom metap sumlog
##' @importFrom XVector extractList
##' @importFrom rtracklayer as.data.frame
##' @importFrom IRanges as.matrix
##'
##' @author Julaiti Shayiding

Fisher_stats <- function(.hitList, peakset) {
  # input param checking
  if (missing(peakset)) {
    stop("Missing required argument peakset, please choose the set of pre-processed peaks!")
  }
  if (missing(.hitList)) {
    stop("Missing required argument .hitList, please choose overlap hit list that comply minimum overlapping peak requirement!")
  }
  stopifnot(inherits(peakset[[1]], "GRanges"))
  message("retrieve pvalue of peaks")
  pval_List <- mapply(.get.pvalue, hitTB, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0.000000e+00)
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

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}
