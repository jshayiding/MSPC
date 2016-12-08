##' Filter ERs by permissive stringency threshold
##'
##' Chip-seq detects genome-wide DNA protein interation, returing enriched regions which associated with significance score.
##' In our methodological framework, we design permissive threshold for all enriched regions' score (A.K.A, peak signal's significance),
##' to classify each peak as stringent or weakly enriched by its p-value.Once getting set of confirmed, discarded enriched regions for all Chip-seq replicates,
##' it is intuitive to generate stringent, weak peak set accordingly. To help user gaining deeper insight and biological evaluation of analysis result,
##' graphical visualization of corresponding peak set also created by handy.
##'
##' @title create_output
##' @description
##' create stringnet/weak ERs by permissive stringent threshold
##'
##' @param peaklist_A set of confirmed enriched regions for all Chip-seq replicates
##' @param peaklist_B set of discarded enriched regions for all Chip-seq replicates
##' @param tau.s permissive threshold for stringent enriched regions, all enriched peaks below this threshold, are considered stringent peaks
##' @param outDir user can control where exported BED file goes
##' @return BED file
##' @export
##' @importFrom rtracklayer export.bed
##' @importFrom magrittr %>%
##' @importFrom magrittr %<>%
##' @importFrom tibble rownames_to_column
##' @importFrom tidyr separate
##' @importFrom dplyr mutate
##' @importFrom dplyr bind_rows
##' @importFrom dplyr group_by
##' @importFrom dplyr tally
##' @importFrom dplyr ungroup
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_col
##' @importFrom ggplot2 facet_wrap
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 position_stack
##' @importFrom ggplot2 geom_text
##'
##' @author Julaiti Shayiding
## create_output(peakList_A = confirmedERs, peakList_B = discardedERs, tau.s = 1.0E-08, output_path ="test/")


create_output <- function(peakList_A, peakList_B , tau.s=1.0E-08, output_path=getwd()) {
  # input param checking
  if (missing(peakList_A)) {
    stop("Missing required argument peakList_A, please choose the list of all confirmed enriched regions in previous workflow!")
  }
  if (missing(peakList_B)) {
    stop("Missing required argument peakList_B, please choose the list of all discarded enriched regions in previous workflow!")
  }
  stopifnot(is.numeric(tau.s))
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path))
  }
  peakList_A <- lapply(peakList_A, as.data.frame)
  peakList_B <- lapply(peakList_B, as.data.frame)

  names(peakList_A) <- paste("confirmed", names(peakList_A), sep = ".")
  names(peakList_B) <- paste("discarded", names(peakList_B), sep = ".")
  combDF <- do.call(rbind, c(peakList_A, peakList_B))
  combDF %<>% rownames_to_column(var = "cn")
  combDF %<>% separate(cn, c("original_list", "Replicate", "seq"), sep = "\\.")
  combDF %<>% mutate(peakStringency = ifelse(p.value <= tau.s , "Stringent", "Weak"))
  res <- combDF %>% split(list(.$Replicate, .$peakStringency, .$original_list))
  res %>%
    bind_rows %>%
    group_by(peakStringency, original_list, Replicate) %>%
    tally %>%
    ungroup %>%
    setNames(c("var", "val", "Replicate", "n")) %>%
    {
      bind_rows(., setNames(., c("val", "var", "Replicate", "n")))
    } %>%
    ggplot(aes(x=var, y=n, fill=val)) + geom_col() +
    facet_wrap(~Replicate)+ geom_text(aes(label=n), position=position_stack(vjust = 0.5))
  out_names <- paste0(output_path, names(combDF), ".csv")
  return(mapply(write.csv, combDF, out_names))
}

##' @example
## create_output(peakList_A = confirmedERs, peakList_B = discardedERs, tau.s = 1.0E-08, output_path ="test/")
