rToSymbol <-
function  #function to reverse the conversion made by symbolToR
### This function reverses the actions of the symbolToR function.
(x
### the character vector returned by the symbolToR function.
 ){
    output <- gsub("hyphen", "-", x, fixed=TRUE)
    output <- gsub("ampersand", "@", output, fixed=TRUE)
    output <- gsub("forwardslash", "/", output, fixed=TRUE)
    output <- sub("symbol.", "", output, fixed=TRUE)
    return(output)
### a character vector of HGNC gene symbols, which are not in general valid R names.
}
