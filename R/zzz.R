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

.onAttach = function(libname, pkgname) {
    msg = "Welcome to 'MSPC'. Here is the initialization of MSPC package"
    msg = strwrap(msg, exdent=4, indent=4)
    packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}
