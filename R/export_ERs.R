#' Identify and Export stringent/ weak ERs
#'
#' Chip-seq detects genome-wide DNA protein interation,
#' returing enriched regions which associated with significance score.
#' Using permissive threshold tau.s for signal' significant
#' value of stringent enriched region, we could identify
#' set of stringent, weakly enriched regions by the output of
#' \link{filterBycombStringency}. All ERs in different group can be
#' exported as BED file by using \link[rtracklayer]{export.bed}.
#'
#' @param peakList_A output of \link{filterBycombStringency},
#' is set of ERs that fulfill combined stringency test,
#' and rescued by Fisher's method, also known as confirmed ERs.
#'
#' @param peakList_B output of \link{filterBycombStringency},
#' is set of discarded ERs that failing for combined stringency test
#'
#' @param tau.s permissive threshold value for stringent ERs,
#' all ERs' pvalue below this value, are considered as stringent ERs.
#'
#' @return stringent/weak ERs can be exported as BED file by using
#' \link[rtracklayer]{export.bed}
#'
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom methods as
#' @importFrom rtracklayer export.bed
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(tidyr)
#' require(tidyverse)
#'
#' ## list of confirmedERs, discardedERs
#' confirmedERs <- GRangesList(
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(7,19,31), c(13,28,43)),
#'         strand = Rle(c("*"),3), rangeName=c("b3","b6","b7"),
#'         score=c(14,9,17)),
#'     bar = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(1,6,16), c(4,12,23)),
#'         strand = Rle(c("*"),3), rangeName=c("a1", "a2", "a3"),
#'         score=c(22,6,13))
#' )
#'
#' discardedERs <- GRangesList(
#'     bar = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(6,25,40), c(12,33,49)),
#'         strand = Rle(c("*"),3), rangeName=c("a2","a5","a8"), score=c(3,2,4)),
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(15,19,47), c(18,28,55)),
#'         strand = Rle(c("*"),3), rangeName=c("b4","b6","b9"), score=c(1,3,6))
#' )
#'
#' ## add pvalue
#' confirmedERs <- lapply(confirmedERs, pvalueConversion)
#' discardedERs <- lapply(discardedERs, pvalueConversion)
#'
#' ## cast GRangesList to data.frame list
#' confirmedDF <- lapply(confirmedERs, as.data.frame)
#' discardedDF <- lapply(discardedERs, as.data.frame)
#'
#' ## Identify and Export Stringent/Weak ERs
#'
#' output <- export_ERs(
#'     peakList_A=confirmedDF,peakList_B=discardedDF, tau.s=1.0E-08)

#' ## Explore ERs in different group
#' print(output)
#'

export_ERs <- function(peakList_A, peakList_B, tau.s=1.0E-08) {
    # input param checking
    if(missing(peakList_A) | missing(peakList_B)) {
        stop("required arguments is missing,
             please choose set of confirmed, discarded ERs")
    }
    peakList_A <- lapply(peakList_A, as.data.frame)
    peakList_B <- lapply(peakList_B, as.data.frame)
    stopifnot(is.numeric(tau.s))
    allERs <-
        bind_rows(c(confirmed = peakList_A,
                    discarded = peakList_B), .id = "id") %>%
        separate(id, c("isConfirmed", "Replicate")) %>%
        mutate(peakStringency = ifelse(p.value <= tau.s,
                                       "Stringent", "Weak")) %>%
        arrange(Replicate, isConfirmed, desc(peakStringency))
    res <- allERs %>% split(list(.$Replicate,
                                 .$peakStringency, .$isConfirmed))
    DF <- lapply(res, function(x) {
        x[ , -which(names(x) %in% c("isConfirmed",
                                    "peakStringency", "Replicate"))]
    })
    asGRs <- lapply(DF, function(x) as(x, "GRanges"))
    rslt <- mapply(
        export.bed, asGRs, paste0(names(asGRs), ".bed"))
    return(rslt)
}
