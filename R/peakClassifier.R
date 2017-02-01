# MSPC Package - an R/Bioconductor Package for Multiple Sample Peak Calling
#'
#'
#' @param peakList_A
#' @param peakList_B
#' @param tau.s
#' @param graphical.output
#' @return
#' @export
#' @importFrom methods hasArg
#' @importFrom rtracklayer export.bed
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
#' @importFrom ggplot2 position_stack
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggsave
#'
#' @author Jurat Shahidin

peakClassifier <- function(peakList_A,
                           peakList_B,
                           tau.s=1.0E-08,
                           graphical.output=TRUE,
                           exportType=c("bed", "csv")) {
    # sanity check for input param
    if(!hasArg(peakList_A)) {
        stop("required arguments is missing,
             please choose set of all confirmed ERs")
    }
    if(!hasArg(peakList_B)) {
        stop("required arguments is missing,
             please choose set of all discarded ERs")
    }
    stopifnot(is.numeric(tau.s))
    allERs <- bind_rows(c(confirmed = peakList_A, discarded = peakList_B), .id = "id") %>%
        separate(id, c("isConfirmed", "Sample")) %>%
        mutate(peakStringency = ifelse(p.value <= tau.s,
                                       "Stringent", "Weak")) %>%
        arrange(Sample, isConfirmed, desc(peakStringency))%>%
        split(list(.$Sample, .$peakStringency, .$isConfirmed))
    if(graphical.output) {
        allERs %>%
            bind_rows %>%
            group_by(peakStringency, isConfirmed, Sample) %>%
            tally %>%
            ungroup %>%
            setNames(c("Replicate", "output", "Sample", "n")) %>%
            {
                bind_rows(., setNames(., c("output", "Replicate", "Sample", "n")))
            } %>%
            ggplot(aes(x=Replicate, y=n, fill=output)) + geom_col() +
            facet_wrap(~paste0(substr(Sample, 1, 20), "\n",
                               substr(Sample, 21, nchar(Sample))),
                       labeller = label_wrap_gen(25),scales = "free_x")+
            geom_text(aes(label=n), position=position_stack(vjust = 0.85))+
            theme(strip.text = element_text(size = 15),
                  axis.text.x = element_text(angle=70, vjust=0.6))
        ggsave("graphical_output.png", width = 14, height = 10)
    } else {
        DF <- allERs %>% map(~.[setdiff(names(.), c("isConfirmed","peakStringency","Sample"))])
        if(exportType=="bed"){
            asGRs <- lapply(DF, function(x) as(x, "GRanges"))
            mapply(export.bed, asGRs, paste0(names(asGRs), ".bed"))
        } else {
            dir_names <- gsub("\\..+","",names(DF)) %>%
                unique() %>% walk(dir.create)
            purrr::walk2(DF,names(DF), function(x,y)
                write.csv(x,paste0(gsub("\\..+","",y),"/",y,".csv")))
        }
    }
}

