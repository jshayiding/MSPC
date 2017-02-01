# MSPC - an R/Bioconductor Package for Multiple Sample Peak Calling
#'
#'
#' @param peakList set of confirmed peaks through combined stringency test.
#' @param pAdjustMethod pvalue adjustment method
#' @param fdr parameter for false discovery rate
#' @param asPlot logical whether produce graphical plot or not
#' @return
#' @export
#'
#' @importFrom rtracklayer export.bed
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
#'
#' @importFrom
#' @author Jurat Shahidin

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
        res <- lapply(DF, function(x) as(x, "GRanges"))
        rslt <- mapply(export.bed, res, paste0(names(res), ".bed"))
        return(rslt)
    }
}
