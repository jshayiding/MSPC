##' Intermediate set purification and to create output set through BH correction test
##'
##' We distinguish Technical and Biological Chip-seq replicates. Technical replicates aim at controlling
##' the variability of the experimental procedure used to obtain the data should yield exactly the same result in
##' the absence of experimental noise. If user choose replicate type as `Technical`,
##' then enriched regions were passed in one test (A.K.A, comply minimum overlapping peak requirement) but failed the another (A.K.A, combined stringency test with Fisher's method),
##' these peaks won't be appeared in final output set. We call this procedure as intermediate set purification. I believe using comparing list-like data structure is intuitive and fast,
##' using `anti_join` method from dplyr packages made this job easily. Then using the Benjamini-Hochberg multiple testing correction with user-specified false discovery rate to create output set accordingly.
##' confirmed peak in output set are exported standard BED file.
##'
##' @title FDR_stats
##'
##' @description proceed set purification for set of ERs and create output set
##'
##' @param peakList_A set of all confirmed enriched regions. `filterByFisherMethod` return set of confirmed peaks in combined stringency method
##' @param peakList_B set of all discarded enriched regions. `filterByFisherMethod` return set of discarded peaks in combined stringency method.
##' @param pAdjustMethod adjusted pvalue for multiple comparison
##' @param .fdr parameter for false discovery rate. User can tune false discovery rate.#' @param replicate.type A charcter vector used to select type of Chip-seq replicate ( Biological / Technical replicate)
##' @param outDir user can control where the exported BED file goes
##'
##' @return BED file
##' @export
##' @importFrom stats p.adjust
##' @importFrom dplyr anti_join
##' @importFrom rtracklayer export.bed
##' @importFrom rtracklayer as.data.frame
##' @author Julaiti Shayiding

FDR_stats <- function(peakList_A, peakList_B, pAdjustMethod="BH", .fdr = 0.05
                      , replicate.type=c("Biological", "Technical"), outDir=getwd(), ...) {
  # check input param
  if (missing(peakList_A)) {
    stop("Missing required argument peakList_A, please choose the list of all confirmed enriched regions in previous workflow!")
  }
  if (missing(peakList_B)) {
    stop("Missing required argument peakList_B, please choose the list of all discarded enriched regions in previous workflow!")
  }
  pAdjustMethod = match.arg(pAdjustMethod)
  replicate.type = match.arg(replicate.type)
  stopifnot(is.numeric(.fdr))
  message("set purification on set of confirmed, discarded peaks")
  # peakList_A, peakList_B must be casted to data.frame


  peakList_A <- lapply(peakList_A, as.data.frame)
  peakList_B <- lapply(peakList_B, as.data.frame)

  if (!dir.exists(outDir)) {
    dir.create(file.path(outDir))
  }
  if(replicate.type=="Biological") {
    .setPurf <- peakList_A
  } else {
    .setPurf <- Map(anti_join, peakList_A, peakList_B)
  }
  .setPurf <- lapply(.setPurf, function(x) as(x, "GRanges"))
  res <- lapply(.setPurf, function(ele_) {
    if(is.null(ele_$p.value)) {
      stop("p.value is required")
    } else {
      p <- ele_$p.value
      ele_$p.adj <- p.adjust(p, method = "BH")
      .filt <- split(ele_,
                     ifelse(ele_$p.adj <= .fdr,
                            "BH_Pass", "BH_Failed"))
    }
  })
  res <- lapply(res, unique)
  rslt <- lapply(names(res), function(ele_) {
    mapply(export.bed,
           res[[ele_]],
           paste0(outDir, ele_, ".", names(res[[ele_]]), ".bed"))
  })
  return(rslt)
}

##' @example
## fdr.rslt <- .FDR.stats(.Confirmed.ERs, .Discarded.ERs, pAdjustMethod = "BH",
##                       .fdr = 0.05, replicate.type = "Bio", outDir = "test/")
