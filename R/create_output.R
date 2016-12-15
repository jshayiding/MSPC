#' Graphical view of different output for each Chip-seq replicates.
#'
#' Chip-seq detects genome-wide DNA protein interation,
#' returing enriched regions which associated with significance score.
#' Using permissive threshold tau.s for signal' significant
#' value of stringent enriched region, we could identify
#' set of stringent, weakly enriched regions by manipulating the output of
#' \link{filterBycombStringency} function.  To help user gaining
#' deeper insight and biological evaluation of analysis result, using \link[ggplot2]{ggplot}
#' to generate stack bar plot for each Chip-seq replicates can be done.
#'
#' @param peakList_A output of \link{filterBycombStringency},
#' is set of all confirmed ERs in \link[GenomicRanges]{GRanges} objects.
#'
#' @param peakList_B output of \link{filterBycombStringency},
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
#' @importFrom tidyr separate
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
#'
#' @author Julaiti Shayiding

create_output <- function(peakList_A, peakList_B, tau.s=1.0E-08) {
  # input param checking
  if (missing(peakList_A)) {
    stop("Missing required argument peakList_A, please choose the list of all confirmed enriched regions in previous workflow!")
  }
  if (missing(peakList_B)) {
    stop("Missing required argument peakList_B, please choose the list of all discarded enriched regions in previous workflow!")
  }
  stopifnot(is.numeric(tau.s))
  # if (!dir.exists(output_path)) {
  #   dir.create(file.path(output_path))
  # }
  names(peakList_A) <- paste("Confirmed", names(peakList_A), sep = ".")
  names(peakList_B) <- paste("Discarded", names(peakList_B), sep = ".")
  combDF <- do.call(rbind, c(peakList_A, peakList_B))
  combDF %<>% rownames_to_column(var = "cn")
  combDF %<>% separate(cn, c("original_list", "letters", "seq"), sep = "\\.")
  combDF %<>% mutate(peakStringency = ifelse(p.value <= tau.s , "Stringent", "Weak"))
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
    facet_wrap(~letters)+ geom_text(aes(label=n), position=position_stack(vjust = 0.85))
}
