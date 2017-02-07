#' Identify and Export stringent/ weak ERs
#'
#' ChIP-Seq detects genome-wide DNA protein interaction,
#' returning enriched regions which associated with significance score.
#' Using permissive threshold tau.s for signal significant
#' value of stringent enriched region, we could identify
#' set of stringent, weakly enriched regions by the output of
#' \link{runMSPC}. All ERs in different group can be
#' exported as either BED format file by using \link[rtracklayer]{export.bed},
#' or csv.
#'
#' @param peakList_A output of \link{runMSPC}, is set of ERs
#' that fulfill combined stringency test, rescued by Fisher's method,
#' also known as confirmed ERs.
#'
#' @param peakList_B output of \link{runMSPC}, is set of discarded ERs
#' that failing from combined stringency test and minimum overlapping peak requirement
#'
#' @param tau.s permissive threshold value for stringent ERs,
#' all ERs' pvalue below this value, are considered as stringent ERs.
#'
#' @param exportFormat user has an option to export the result either BED format file or csv.
#'
#' @return stringent/weak peak set
#'
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom purrr map
#' @importFrom purrr walk
#' @importFrom purrr walk2
#' @importFrom methods as
#' @importFrom rtracklayer export.bed
#' @importFrom utils write.csv
#' @importFrom methods hasArg
#'
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
#' discardedERs <- runMSPC(peakset = total.ERs, whichType = "max",
#'                         cmbStrgThreshold = 1.0E-08, isConfirmed = FALSE)
#'
#' ## explore the peaks in different output set for file bar
#' export_ERs(peakList_A = confirmedERs,
#'            peakList_B = discardedERs,
#'            tau.s = 1.0E-08, exportFormat = "csv")
#'

export_ERs <- function(peakList_A,
                       peakList_B,
                       tau.s=1.0E-08,
                       exportFormat=c("bed","csv")) {
    # sanity check for input param
    if(!hasArg(peakList_A)) {
        stop("required arguments is missing,
             please choose set of all confirmed peaks")
    }
    if(!hasArg(peakList_B)) {
        stop("required arguments is missing,
             please choose set of all discarded peaks")
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
        arrange(Replicate, isConfirmed, desc(peakStringency)) %>%
        split(list(.$Replicate, .$peakStringency, .$isConfirmed)) %>%
        map(~.[setdiff(names(.), c("isConfirmed",
                                   "peakStringency","Replicate"))])
    if(exportFormat=="csv") {
        # put exported csv files into individual folder
        # purrr::walk2 works for data.frame like object
        dir_names <- gsub("\\..+","",names(allERs)) %>%
            unique() %>% walk(dir.create)
        purrr::walk2(allERs,names(allERs), function(x,y)
            write.csv(x,paste0(gsub("\\..+","",y),"/",y,".csv")))
    } else{
        asGRs <- lapply(allERs, function(x) as(x, "GRanges"))
        mapply(export.bed, asGRs, paste0(names(asGRs), ".bed"))
    }
}
