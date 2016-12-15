#' Data conversion for significance of peak signal
#'
#' Some data sources provides Chip-seq enriched regions (A.K.A, peaks) without p-value under specific conditions.
#' pvalueConversion is the utility function to do data conversion of peak' score as p-value
#' In standard BED file, significant value of peak signal defined as Score column,
#' Thus we need to convert it as p.value (- log(p.value), -10 log(p.value), -100 log(p.value))
#'
#' Passed to \link{readPeakFiles}
#'
#' @title pvalueCoversion
#' @param x GRanges objects All Enriched regions (a.k.a, peaks) are stored in GRanges object.
#' @param pvalueBase User has option to set up p-value format (- log(p.value), -10 log(p.value), -100 log(p.value))
#' @return GRanges
#' @export
#' @importFrom rtracklayer score
#' @importFrom rtracklayer mcols
#' @importFrom IRanges colnames
#' @author Julaiti Shayiding

pvalueConversion <- function(x, pvalueBase = 10L) {
  # input param checking
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$pvalue)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  res <- x
  return(res)
}
