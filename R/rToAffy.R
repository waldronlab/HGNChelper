#' Title function to convert the output of affyToR back to the original Affymetrix probeset identifiers.
#' @description This function simply strips the "affy." added by the \link{affyToR} function.
#' @param x the character vector returned by the affyToR function.
#' @return a character vector of Affymetrix probeset identifiers.
#' @export
#'
rToAffy <- function(x){
    output <- sub("affy.", "", x, fixed=TRUE)
    return(output)
}
