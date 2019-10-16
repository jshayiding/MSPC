#' Import Chip-seq replicates and all ERs are stored in GRanges objects.
#'
#' This function can read peak file in standard BED file
#' using \link[rtracklayer]{import.bed} or data.frame like object
#' using \link[GenomicRanges]{makeGRangesFromDataFrame}
#' and stored in \link[GenomicRanges]{GRanges} object,
#' where several peak files can be read simultaneously
#' using lapply. Note that some data sources provide ChIP-Seq
#' enriched regions without p-value under specific
#' conditions. In our workflow, we are going to evaluate each
#' peak interval with p-value threshold, so represent peak's score
#' as p-value is needed.
#'
#' Passed to \link{denoise_ERs}
#'
#' @param peakFolder ChipSeq peak replicate either in BED format or data.frame
#' @param verbose logical that control whether the output is printed or not
#' @param pvalueBase User has option to select
#' p-value format (- log(p.value), -10 log(p.value), -100 log(p.value))
#'
#' @return GRanges object
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom rtracklayer score
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges GRangesList
#' @importFrom methods as
#' @importFrom methods hasArg
#' @importFrom stats setNames
#' @importFrom tools file_ext
#' @author  Jurat Shahidin
#'
#' @examples
#' require(rtracklayer)
#' require(GenomicRanges)
#'
#' ## get bed file
#' files <- getPeakFile()[1:3]
#' readPeakFiles(files, pvalueBase=1L)

readPeakFiles <- function(peakFolder,
                          pvalueBase=1L,
                          verbose=FALSE) {
    # input param checking
    if (missing(peakFolder)) {
        stop("Missing required argument peakFolder, please
             choose input Chip-seq replicates in BED file!")
    }
    if(missing(pvalueBase)) {
        stop("please specify pvalue convention")
    }
    stopifnot(length(peakFolder)>0)
    stopifnot(is.numeric(pvalueBase))
    if (verbose)
        cat(">> reading all peakfiles from given folder...\t\t\t",
            format(Sys.time(), "%Y-%M-%D %X"), "\n")
    f.read <- lapply(peakFolder, function(ele_) {
        if(tools::file_ext(ele_)=="data.frame")
            ele_ <- makeGRangesFromDataFrame(ele_)
        if(tools::file_ext(ele_)=="bed")
            ele_ <- as(import.bed(ele_), "GRanges")
        if(is.null(ele_$p.value)){
            ele_$p.value <- 10^(as.numeric(score(ele_)/(- pvalueBase)))
            format(ele_$p.value, scientific = TRUE)
        }
        ele_
    })
    return(GRangesList(f.read))
}
