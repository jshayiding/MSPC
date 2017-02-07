#' Graphical view of different ERs set for each Chip-seq replicates.
#'
#' This function is served as graphical version of \link{export_ERs}.
#' To help user gaining deeper insight and biological evaluation of
#' analysis result, using \link[ggplot2]{ggplot} to generate
#' stack bar plot for each Chip-seq replicates can be done.
#'
#' @param peakList_A output of \link{runMSPC},
#' is set of all confirmed ERs in \link[GenomicRanges]{GRanges} objects.
#'
#' @param peakList_B output of \link{runMSPC},
#' is set of all discarded ERs in \link[GenomicRanges]{GRanges} objects.
#'
#' @param tau.s permissive threshold for stringent enriched regions,
#' all enriched regions below this threshold, are considered stringent ERs
#'
#' @return using \link[ggplot2]{ggplot} to generate stack bar plot for file bar
#' @export
#' @importFrom methods hasArg
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
#' @importFrom dplyr setdiff
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom dplyr tally
#' @importFrom dplyr ungroup
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 label_wrap_gen
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 position_stack
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
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
#' # Visualize the output set for file bar
#' getPlot(peakList_A = confirmedERs,
#'         peakList_B = discardedERs, tau.s = 1.0E-08)

getPlot <- function(peakList_A,
                    peakList_B,
                    tau.s=1.0E-08) {
    # sanity check for input param
    if(!hasArg(peakList_A)) {
        stop("required arguments is missing,
             please choose set of all confirmed ERs")
    }
    if(!hasArg(peakList_B)) {
        stop("required arguments is missing,
             please choose set of all discarded ERs")
    }
    peakList_A <- lapply(peakList_A, data.frame)
    peakList_B <- lapply(peakList_B, data.frame)
    stopifnot(is.numeric(tau.s))
    plotDat <- bind_rows(c(confirmed = peakList_A,
                          discarded = peakList_B), .id = "id") %>%
        separate(id, c("isConfirmed", "Sample")) %>%
        mutate(peakStringency = ifelse(p.value <= tau.s,
                                       "Stringent", "Weak")) %>%
        arrange(Sample, isConfirmed, desc(peakStringency))%>%
        split(list(.$Sample, .$peakStringency, .$isConfirmed))
    res <- plotDat %>% bind_rows %>%
        group_by(peakStringency, isConfirmed, Sample) %>%
        tally %>% ungroup %>%
        setNames(c("Replicate", "Output", "Sample", "n")) %>%
        {
            bind_rows(., setNames(., c("Output", "Replicate", "Sample", "n")))
        } %>%
        ggplot(aes(x=Replicate, y=n, fill=Output)) + geom_col() +
        facet_wrap(~paste0(substr(Sample, 1, 20), "\n",
                           substr(Sample, 21, nchar(Sample))),
                   labeller = label_wrap_gen(25),scales = "free_x")+
        geom_text(aes(label=n), position=position_stack(vjust = 0.85))+
        labs(x = "Replicate", y = "count") +
        theme(strip.text = element_text(size = 18),
              axis.text.x = element_text(size=16, angle=70, vjust=0.6),
              axis.title.x = element_text(size=18,angle=0,hjust=.5,vjust=0),
              axis.title.y = element_text(size=18,angle=90,hjust=.5,vjust=.5))
    return(res)
}

