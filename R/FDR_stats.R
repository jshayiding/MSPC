# MSPC - an R/Bioconductor Package for Multiple Sample Peak Calling
#'
#'
#' @param peakList
#' @param pAdjustMethod
#' @param fdr
#' @param asPlot
#' @return
#' @export
#' @importFrom
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
