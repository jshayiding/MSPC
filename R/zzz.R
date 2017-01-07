#' R CMD check NOTE avoidance
#'
#' @importFrom utils globalVariables

.onLoad <- function(libname = find.package("MSPC"), pkgname = "MSPC"){
    # R CMD check Note avoidance
    if(getRversion() >= "2.15.1")
        utils::globalVariables(
            c(".","p.value","peakStringency","original_list",
              "id","Replicate", "isConfirmed","n","output","name")
        )
    invisible()
}
