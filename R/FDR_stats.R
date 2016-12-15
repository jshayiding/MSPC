#' Multiple Testing Correction
#'
#' We distinguish Technical and Biological Chip-seq replicates.
#' Technical replicates aim at controlling the variability of
#' the experimental procedure used to obtain the data should
#' yield exactly the same result in the absence of experimental noise.
#' If user choose replicate type as `Technical`, then enriched regions
#' were passed in one test (A.K.A, comply minimum overlapping peak requirement)
#' but failed in the another (A.K.A, combined stringency test with Fisher's method),
#' these peaks won't be appeared in final output.
#' Given the output of \link{filterBycombStringency}, we decided to use
#' intermediate set purification to further check any ERs both exist in
#' confirmed, discarded ERs set. We kept all confirmed ERs, discarded ERs in the list,
#' so using `anti_join` method from dplyr packages make set purification easy.
#' However, we need to correct the p.value of ERs using the Benjamini-Hochberg
#' multiple testing correction with user-specified false discovery rate,
#' using \link[stats]{p.adjust} to do BH correction method,  which yields
#' multiple-testing confirmed or discarded ERs.
#'
#' @param peakList_A output of \link{filterBycombStringency}, is set of all confirmed ERs in \link[GenomicRanges]{GRanges} objects.
#' @param peakList_B output of \link{filterBycombStringency}, is set of all discarded ERs in \link[GenomicRanges]{GRanges} objects.
#' @param fdr parameter for false discovery rate. User can tune false discovery rate.
#' @param pAdjustMethod pvalue adjustment method
#' @param replicate.type A charcter vector used to select type of Chip-seq replicate ( Biological / Technical replicate).
#'
#' @return output of multiple testing correction can be exported as BED file by using \link[rtracklayer]{export.bed} method.
#' @export
#' @importFrom stats p.adjust
#' @importFrom dplyr anti_join
#' @importFrom utils write.csv
#' @importFrom rtracklayer as.data.frame
#' @author Julaiti Shayiding

FDR_stats <- function(peakList_A, peakList_B, pAdjustMethod="BH", fdr = 0.05
                      , replicate.type=c("Biological", "Technical")) {
  # check input param
  if (missing(peakList_A)) {
    stop("Missing required argument peakList_A, please choose the list of all confirmed enriched regions in previous workflow!")
  }
  if (missing(peakList_B)) {
    stop("Missing required argument peakList_B, please choose the list of all discarded enriched regions in previous workflow!")
  }
  pAdjustMethod = match.arg(pAdjustMethod)
  replicate.type = match.arg(replicate.type)
  stopifnot(is.numeric(fdr))
  message("set purification on set of confirmed, discarded peaks")
  ## peakList_A, peakList_B must be casted to data.frame
  peakList_A <- lapply(peakList_A, as.data.frame)
  peakList_B <- lapply(peakList_B, as.data.frame)

  # if (!dir.exists(outDir)) {
  #   dir.create(file.path(outDir))
  # }
  if(replicate.type=="Biological") {
    setPurf <- peakList_A
  } else {
    setPurf <- Map(anti_join, peakList_A, peakList_B)
  }
  setPurf <- lapply(setPurf, function(x) as(x, "GRanges"))
  byFDR <- lapply(setPurf, function(ele_) {
    if(is.null(ele_$p.value)) {
      stop("p.value is required")
    } else {
      p <- ele_$p.value
      ele_$p.adj <- p.adjust(p, method = "BH")
      ele_ <- split(ele_,
                     ifelse(ele_$p.adj <= fdr,
                            "BH_Pass", "BH_Failed"))
      ele_
    }
  })
  byFDR <- lapply(byFDR, unique)
  rslt <- lapply(names(byFDR), function(elm) {
    mapply(write.csv,
           byFDR[[elm]],
           paste0(elm, ".", names(byFDR[[elm]]), ".csv"))
  })
  return(rslt)
}
