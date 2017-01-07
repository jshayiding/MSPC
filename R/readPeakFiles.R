#' Import Chip-seq replicates and all ERs are stored in GRanges objects.
#'
#' readPeakFiles can read peak file in standard BED format
#' using \link[rtracklayer]{import.bed} and stored in
#' \link[GenomicRanges]{GRanges} object, where several peak files
#' can be read simultaneously using lapply. Note that some data
#' sources provides Chip-seq enriched regions (A.K.A, peaks)
#' without p-value under specific conditions.
#'
#' Passed to \link{denoise_ERs}
#'
#' @param peakFolder set of ChipSeq peak replicate in BED format file
#' @param verbose logical that control whether the output is printed or not
#' @param pvalueBase User has option to select
#' p-value format (- log(p.value), -10 log(p.value), -100 log(p.value))
#' @return GRanges object
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom rtracklayer score
#' @importFrom methods as
#' @importFrom stats setNames
#' @author  Julaiti Shayiding
#'
#' @examples
#' require(rtracklayer)
#' require(GenomicRanges)
#' ## get bed file
#' files <- getPeakFile()[1:3]
#' rd <- readPeakFiles(files, pvalueBase=1L)

readPeakFiles <- function(peakFolder, pvalueBase=1L, verbose=FALSE) {
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
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    f.read <- lapply(peakFolder, function(ele_) {
        .gr <- as(import.bed(ele_), "GRanges")
        if(is.null(.gr$p.value)) {
            .gr$p.value <- 10^(score(.gr)/(- pvalueBase))
        } else {
            .gr
        }
        return(.gr)
    })
}
