#' Multiple Testing Correction
#'
#' We distinguish Technical and Biological Chip-seq replicates.
#' Technical replicates aim at controlling the variability of
#' the experimental procedure used to obtain the data should
#' yield exactly the same result in the absence of experimental noise.
#' If user choose replicate type as `Technical`, then enriched regions
#' were passed in one test
#' (A.K.A, comply minimum overlapping peak requirement) but failed
#' in the another (A.K.A, combined stringency test with Fisher's method),
#' these peaks won't be appeared in final output.
#' Given the output of \link{filterBycombStringency},
#' we decided to use intermediate set purification
#' to further check any ERs both exist in confirmed, discarded ERs set.
#' We kept all confirmed ERs, discarded ERs in the list,
#' so using `anti_join` method from dplyr packages
#' make set purification easy. However, we need to correct
#' the p.value of ERs using the Benjamini-Hochberg
#' multiple testing correction with user-specified false discovery rate,
#' using \link[stats]{p.adjust} to do BH correction method,  which yields
#' multiple-testing confirmed or discarded ERs.
#'
#' @param peakList_A output of \link{filterBycombStringency},
#' is set of all confirmed ERs
#' in \link[GenomicRanges]{GRanges} objects.
#'
#' @param peakList_B output of \link{filterBycombStringency},
#' is set of all discarded ERs
#' in \link[GenomicRanges]{GRanges} objects.
#'
#' @param fdr parameter for false discovery rate;
#'
#' @param pAdjustMethod pvalue adjustment method
#'
#' @param replicate.type A charcter vector used to select
#' type of Chip-seq replicate ( Biological / Technical replicate).
#'
#' @return output of multiple testing correction can be exported as
#' BED file by using \link[rtracklayer]{export.bed} method.
#'
#' @export
#' @importFrom stats p.adjust
#' @importFrom dplyr anti_join
#' @importFrom rtracklayer export.bed
#' @importFrom rtracklayer as.data.frame
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(XVector)
#'
#' ## prepare list of confirmedERs, discardedERs
#' confirmedERs <- GRangesList(
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(7,19,31), c(13,28,43)),
#'         strand = Rle(c("*"),3), rangeName=c("b3","b6","b7"),
#'         score=c(14,9,17),p.value=c(1e-14,1e-09,1e-17)),
#'     bar = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(1,6,16), c(4,12,23)),
#'         strand = Rle(c("*"),3), rangeName=c("a1", "a2", "a3"),
#'         score=c(22,6,13),p.value=c(1e-22,1e-06,1e-13))
#' )
#'
#' discardedERs <- GRangesList(
#'     bar = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(6,25,40), c(12,33,49)),
#'         strand = Rle(c("*"),3), rangeName=c("a2","a5","a8"),
#'         score=c(3,12,4),p.value=c(1e-03,1e-12,1e-04)),
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(15,19,47), c(18,28,55)),
#'         strand = Rle(c("*"),3), rangeName=c("b4","b6","b9"),
#'         score=c(11,3,6),p.value=c(1e-11,1e-03,1e-06))
#' )
#'
#' ## call FDR_stats to create output set
#' BH_output <- FDR_stats(
#'     peakList_A = confirmedERs,peakList_B = discardedERs,
#'     pAdjustMethod = "BH", fdr = 0.05,
#'     replicate.type = "Technical")
#'
#' ## Explore output of multiple testting correction
#' BH_output


FDR_stats <- function(peakList_A, peakList_B, pAdjustMethod="BH",
                      fdr = 0.05,
                      replicate.type=c("Biological", "Technical")) {
    # check input param
    if (missing(peakList_A)) {
        stop("Missing required argument peakList_A, please choose the
             list of all confirmed enriched regions in previous workflow!")
    }
    if (missing(peakList_B)) {
        stop("Missing required argument peakList_B, please choose the
             list of all discarded enriched regions in previous workflow!")
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
        mapply(export.bed,
               byFDR[[elm]],
               paste0(elm, ".", names(byFDR[[elm]]), ".bed"))
        })
    return(rslt)
}
