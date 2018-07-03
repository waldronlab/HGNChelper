findExcelGeneSymbols <-
function  #function to identify Excel-mogrified gene symbols
### This function identifies gene symbols which may have been
### mogrified by Excel or other spreadsheet programs.  If output is
### assigned to a variable, it returns a vector of the same length
### where symbols which could be mapped have been mapped.
(x,
 ### Vector of gene symbols to check for mogrified values
 mog.map=read.csv(system.file("extdata/mog_map.csv", package = "HGNChelper"), as.is=TRUE),
 ### Map of known mogrifications.  A default map is available with this
 ### package by data(mog.map), but any map may be used.  This should be
 ### a dataframe with two columns: original and mogrified, containing
 ### the correct and incorrect symbols, respectively.
 regex="[0-9]\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\.[0-9]
[0-9]E\\+[[0-9][0-9]"
 ### Regular expression, recognized by the base::grep function which
 ### is called with ignore.case=TRUE, to identify mogrified symbols.
 ### It is not necessary for all matches to have a corresponding entry
 ### in mog.map$mogrified; values in x which are matched by this regex
 ### but are not found in mog.map$mogrified simply will not be
 ### corrected.  This regex is based that provided by Zeeberg et al.,
 ### Mistaken Identifiers: Gene name errors can be introduced
 ### inadvertently when using Excel in bioinformatics.  BMC
 ### Bioinformatics 2004, 5:80.
 ){
    if(class(x) != "character"){
        x <- as.character(x)
        warning("coercing x to character.")
    }
    if(!is.null(mog.map) & class(mog.map) != "data.frame"){
        mog.map <- data.frame(mog.map)
        warning("coercing mog.map to data.frame")
    }
        
    mog.symbols <- grep(pattern=regex, x, ignore.case=TRUE, value=TRUE)
    if(!is.null(mog.map) & length(mog.symbols) > 0){
        mog.mapped <- data.frame(mogrified=mog.symbols,
                                 corrected=mog.map$original[ match(toupper(mog.symbols), toupper(mog.map$mogrified)) ],
                                 stringsAsFactors=FALSE)
        mog.mapped$abletofix <- ifelse(is.na(mog.mapped$corrected), "no", "yes")
        mog.mapped$corrected[is.na(mog.mapped$corrected)] <- mog.mapped$mogrified[is.na(mog.mapped$corrected)]
        warning.message <- paste(paste(mog.mapped$mogrified,
                                       mog.mapped$corrected,
                                       sep=" to "),
                                 collapse=", ")
        warning(paste("Transmogrified gene symbols found.  Returning the following corrections:", warning.message))
        mog.locations <- match(mog.mapped$mogrified, x)
        x[mog.locations] <- mog.mapped$corrected
    }
    invisible(x)
### if the return value of the function is assigned to a variable, the
### function will return a vector of the same length as the input,
### with corrections possible from mog.map made.
}
