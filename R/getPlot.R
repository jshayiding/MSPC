#' Graphical view of different ERs set for each Chip-seq replicates.
#'
#' This function is served as graphical version of \link{export_ERs}.
#' To help user gaining deeper insight and biological evaluation of
#' analysis result, using \link[ggplot2]{ggplot} to generate
#' stack bar plot for each Chip-seq replicates can be done.
#'
#' @param peakList_A output of \link{filterBycombStringency},
#' is set of all confirmed ERs in \link[GenomicRanges]{GRanges} objects.
#'
#' @param peakList_B output of \link{mergeDiscERs},
#' is set of all discarded ERs in \link[GenomicRanges]{GRanges} objects.
#'
#' @param tau.s permissive threshold for stringent enriched regions,
#' all enriched regions below this threshold, are considered stringent ERs
#'
#' @return using \link[ggplot2]{ggplot} to generate stack bar plot for file bar
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate_
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr tally
#' @importFrom dplyr ungroup
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 position_stack
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 label_wrap_gen
#' @importFrom ggplot2 element_text
#' @author Julaiti Shayiding
#'
#' @examples
#' require(GenomicRanges)
#' require(rtracklayer)
#' require(XVector)
#' require(ggplot2)
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
#'         score=c(13,2,4),p.value=c(1e-13,1e-02,1e-04)),
#'     cat = GRanges(
#'         seqnames=Rle("chr1", 3),ranges=IRanges(c(15,19,47), c(18,28,55)),
#'         strand = Rle(c("*"),3), rangeName=c("b4","b6","b9"),
#'         score=c(11,3,6),p.value=c(1e-11,1e-03,1e-06))
#' )
#'
#' ## visualize output
#' output <- getPlot(
#'   peakList_A=confirmedERs, peakList_B=discardedERs, tau.s=1.0E-08)
#'


getPlot <- function(peakList_A, peakList_B, tau.s=1.0E-08) {
    # input param checking
    if (missing(peakList_A)) {
        stop("Missing required argument peakList_A, please choose
             the list of all confirmed enriched regions in previous workflow!")
    }
    if (missing(peakList_B)) {
        stop("Missing required argument peakList_B, please choose
             the list of all discarded enriched regions in previous workflow!")
    }
    stopifnot(is.numeric(tau.s))
    # cast peakList_A, peakList_B in data.frame
    peakList_A <- lapply(peakList_A, as.data.frame)
    peakList_B <- lapply(peakList_B, as.data.frame)
    names(peakList_A) <- paste("Confirmed", names(peakList_A), sep = ".")
    names(peakList_B) <- paste("Discarded", names(peakList_B), sep = ".")
    combDF <- do.call(rbind, c(peakList_A, peakList_B))
    combDF %<>% rownames_to_column(var = "cn")
    combDF %<>% separate_("cn",
                          c("original_list", "letters", "seq"), sep = "\\.")
    combDF %<>% mutate(peakStringency = ifelse(p.value <= tau.s ,
                                               "Stringent", "Weak"))
    res <- combDF %>% split(list(.$letters, .$peakStringency, .$original_list))
    res %>%
        bind_rows %>%
        group_by(peakStringency, original_list, letters) %>%
        tally %>%
        ungroup %>%
        setNames(c("Replicate", "output", "letters", "n")) %>%
        {
            bind_rows(., setNames(., c("output", "Replicate", "letters", "n")))
        } %>%
        ggplot(aes(x=Replicate, y=n, fill=output)) + geom_col() +
        facet_wrap(~paste0(substr(letters, 1, 20), "\n",
                           substr(letters, 21, nchar(letters))),
                   labeller = label_wrap_gen(15),scales = "free_x")+
        geom_text(aes(label=n), position=position_stack(vjust = 0.85))+
        theme(strip.text = element_text(size = 10),
              axis.text.x = element_text(angle=70, vjust=0.6))
}
