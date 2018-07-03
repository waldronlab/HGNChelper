rToAffy <-
function  #function to convert the output of affyToR back to the original Affymetrix probeset identifiers.
### This function simply strips the "affy." added by the \link{affyToR} function.
(x
 ### the character vector returned by the affyToR function.
 ){
    output <- sub("affy.", "", x, fixed=TRUE)
    return(output)
### a character vector of Affymetrix probeset identifiers.
}
