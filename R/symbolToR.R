#' Title function to *reversibly* convert HGNC gene symbols to valid R names.
#' @description 
#' This function reversibly converts HGNC gene symbols to valid R
#' names by prepending "symbol.", and making the following
#' substitutions: "-" to "hyphen", "@" to "ampersand", and "/" to
#' "forwardslash".
#' @param x vector of HGNC symbols
#' @return
#' a vector of valid R names, of the same length as x, which can be
#' converted to the same HGNC symbols using the rToSymbol function.
#' @seealso \code{\link{rToSymbol}}
#' @export
#'
symbolToR <- function(x){
    output <- gsub("-", "hyphen", x, fixed=TRUE)
    output <- gsub("@", "ampersand", output, fixed=TRUE)
    output <- gsub("/", "forwardslash", output, fixed=TRUE)
    output <- paste("symbol.", output, sep="")
    return(output)
}
