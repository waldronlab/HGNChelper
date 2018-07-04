#' Title function to convert Affymetrix probeset identifiers to valid R names
#' @description 
#' This function simply prepends "affy." to the probeset IDs to create
#' valid R names.  Reverse operation is done by the \link{rToAffy}
#' function.
#' @param x 
#' vector of Affymetrix probeset identifiers, or any identifier which may with a digit.
#' @return
#' a character vector that is simply x with "affy." prepended to each value.
#' @export
affyToR <- function(x){
    output <- paste("affy.", x, sep="")
    return(output)
}
