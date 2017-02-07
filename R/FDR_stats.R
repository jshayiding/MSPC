# MSPC - an R/Bioconductor Package for Multiple Sample Peak Calling
#'
#' Multiple testing correction
#'
#' Given the output of \link{runMSPC}, we'll get set peaks which comply with
#' combined stringency test by Fisher method, and classified as confirmed peaks.
#' However, we need to further evaluate set of confirmed peaks
#' with multiple testing corrections procedure, and produce final output set.
#'
#' @param peakList set of confirmed peaks through combined stringency test.
#' @param pAdjustMethod pvalue adjustment method
#' @param fdr parameter for false discovery rate
#' @param asPlot logical whether produce graphical plot or not
#' @return set of enriched regions in BED format file
#' @export
# @importFrom rtracklayer export.bed
#' @importFrom utils write.csv
#' @importFrom stats p.adjust
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate
#' @importFrom dplyr mutate
#' @importFrom dplyr count
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @author Jurat Shahidin
#'
#' @examples
#' # set up
#' library(GenomicRanges)
#' library(rtracklayer)
#'
#' # load peak files
#' files <- getPeakFile()[1:3]
#' grs <- readPeakFiles(files, pvalueBase=1L)
#'
#' ## Exclude background noise
#' total.ERs <- denoise_ERs(peakGRs = grs, tau.w = 1.0E-04,
#'                         overwrite = TRUE)
#'
#' ## explore set of confirmed, discarde peaks
#' confirmedERs <- runMSPC(peakset = total.ERs, whichType = "max",
#'                         cmbStrgThreshold = 1.0E-08, isConfirmed = TRUE)
#'
#' # multiple testing correction
#' FDR_stats(peakList = confirmedERs,
#'           pAdjustMethod = "BH",
#'           fdr = 0.05, asPlot =FALSE )
#'

FDR_stats <- function(peakList,
                      pAdjustMethod="BH",
                      fdr = 0.05,
                      asPlot=TRUE) {
    # sanity check for input param
    if(!hasArg(peakList)) {
        stop("required arguments is missing,
             please choose set of all confirmed ERs")
    }
    pAdjustMethod = match.arg(pAdjustMethod)
    stopifnot(is.numeric(fdr))
    peakList <- lapply(peakList, data.frame)
    res <- bind_rows(peakList, .id = "id") %>%
        separate(id, "sample") %>%
        mutate(adjust.pvalue = p.adjust(p.value, method = pAdjustMethod),
               Output = ifelse(adjust.pvalue <= fdr, "Output_mtc", "Output_mtd")) %>%
        split(list(.$sample, .$Output))
    if(asPlot) {
        plotDat <- bind_rows(res, .id = "Catg") %>% count(sample, Catg)
        ggplot(plotDat, aes(x = sample, y = n, fill = Output)) +
            geom_col(aes(fill = Catg), position = "dodge")
    } else {
        DF <- res %>% lapply(select, -c(sample, Output))
        ## Give option to export them as `csv` or `bed`
        #res <- lapply(DF, function(x) as(x, "GRanges"))
        rslt <- mapply(write.csv, DF, paste0(names(DF), ".csv"))
        return(rslt)
    }
}
