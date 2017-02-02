# MSPC - an R/Bioconductor package for Multiple Sample Peak Calling
#
#' Assess overlapping enriched regions across multiple sample
#'
#'
#' @param peakset the output of \link{denoise_ERs}.
#' set of Chip-seq replicate imported and all peaks are stored in
#' GRanges object, where all background signal were excluded.
#'
#' @param whichType user has options to keep most stringent(with loest p-value)
#' or least stringent(highest p-value) peak if multiple overlapping peaks
#' were detected. By default, keep most stringent one.
#'
#' @param replicate.type A charcter vector used to select type of input
#' Chip-seq replicate ( Biological / Technical replicate).
#'
#' @param cmbStrgThreshold combined stringency threshold
#' against all enriched regions p-value, so we could identify
#' whether ERs fulfill stringency test,
#' and result can be set of confirmed/discarded peaks respectively.
#'
#' @param isConfirmed logical vector that check whether ERs
#' comply combined stringency test.
#'
#' @return
#' \code{isConfirmed} is \code{True}:
#' return set of enriched regions passed from combined stringency test
#'
#' \code{isConfirmed} is \code{False} :
#' return set of enriched regions that failed from Fisher method
#'
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom IRanges splitAsList
#' @importFrom utils relist
#' @importFrom stats pchisq
#' @author Jurat Shahidin
#'
#' @note Special thanks to Martin Morgan's contribution on this revision

runMSPC <- function(peakset,
                    whichType=c("max","min"),
                    replicate.type=c("Biological", "Technical"),
                    cmbStrgThreshold=1.0E-08,
                    isConfirmed=TRUE) {
    # sanity input param checking
    if (class(peakset) != "GRangesList") {
        stop("input must be a GRangesList Object")
    }
    whichType = match.arg(whichType)
    replicate.type = match.arg(replicate.type)
    stopifnot(length(peakset)>0)
    stopifnot(is.numeric(cmbStrgThreshold))
    min.c <- ifelse(replicate.type=="Biological",
                    length(peakset)-1,
                    length(peakset))
    # flatten out input peak list
    allPeaks <- unlist(peakset, use.names = FALSE)
    hitIdx <- findOverlaps(allPeaks)
    new_mcols <- DataFrame(
        query = factor(queryHits(hitIdx), levels=seq_along(allPeaks)),
        subject = factor(subjectHits(hitIdx), levels=seq_along(allPeaks)),
        score = allPeaks$score[subjectHits(hitIdx)],
        p.value = allPeaks$p.value[subjectHits(hitIdx)]
    )
    new_mcols <- new_mcols[order(new_mcols$score),]
    len_vect <- rep(seq_along(peakset), lengths(peakset))
    bindHit <- cbind(new_mcols$query, len_vect[new_mcols$subject])
    if(whichType == "max") {
        message("use most stringent peaks")
        multiOv <- !duplicated(bindHit, fromLast=TRUE)
    }
    multiOv <- !duplicated(bindHit, fromLast=FALSE)
    res <- new_mcols[multiOv,]
    # check cardinality for minimum overlapping peak requirement
    K <- tabulate(res$query) >= min.c
    keepIdx <- K[res$query]
    res <- res[keepIdx,]
    pvlist <- splitAsList(res$p.value, res$query)
    allPeaks$comb.pv <- pchisq(-2 * sum(log(pvlist)), df=lengths(pvlist), lower.tail=FALSE)
    Keep <- allPeaks$comb.pv <= cmbStrgThreshold
    if(!isConfirmed)
        Keep <- !(Keep)
    gr <- relist(allPeaks, peakset)[relist(Keep, peakset)]
    # Make sure `gr` returned as `data.frame`
    # Find out better way to produce data.frame list, instead of using `lapply`
    res <- lapply(gr, as.data.frame)
    return(res)
}
