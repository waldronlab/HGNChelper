symbolToR <-
function  #function to *reversibly* convert HGNC gene symbols to valid R names.
### This function reversibly converts HGNC gene symbols to valid R
### names by prepending "symbol.", and making the following
### substitutions: "-" to "hyphen", "@" to "ampersand", and "/" to
### "forwardslash".
(x
### vector of HGNC symbols
 ){
    output <- gsub("-", "hyphen", x, fixed=TRUE)
    output <- gsub("@", "ampersand", output, fixed=TRUE)
    output <- gsub("/", "forwardslash", output, fixed=TRUE)
    output <- paste("symbol.", output, sep="")
    return(output)
### a vector of valid R names, of the same length as x, which can be
### converted to the same HGNC symbols using the rToSymbol function.
}
