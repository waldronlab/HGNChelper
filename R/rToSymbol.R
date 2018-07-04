#' Title function to reverse the conversion made by symbolToR
#' @description This function reverses the actions of the symbolToR function.
#' @param x the character vector returned by the symbolToR function.
#' @return a character vector of HGNC gene symbols, which are not in general valid R names.
#' @seealso \code{\link{symbolToR}}
#' @export
#'
rToSymbol <- function(x){
    output <- gsub("hyphen", "-", x, fixed=TRUE)
    output <- gsub("ampersand", "@", output, fixed=TRUE)
    output <- gsub("forwardslash", "/", output, fixed=TRUE)
    output <- sub("symbol.", "", output, fixed=TRUE)
    return(output)
}
