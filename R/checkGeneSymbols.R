checkGeneSymbols <-
function  #function to identify outdated or Excel-mogrified gene symbols
### This function identifies gene symbols which are outdated or may have been
### mogrified by Excel or other spreadsheet programs.  If output is
### assigned to a variable, it returns a data.frame of the same number of
### rows as the input, with a second column indicating whether the symbols
### are valid and a third column with a corrected gene list.
(x,
 ### Vector of gene symbols to check for mogrified or outdated values
 unmapped.as.na=TRUE,
 ### If TRUE, unmapped symbols will appear as NA in the
 ### Suggested.Symbol column.  If FALSE, the original unmapped symbol
 ### will be kept.
 hgnc.table=NULL
 ### If hgnc.table is a data.frame with colnames(hgnc.table) identical
 ### to c("Symbol", "Approved.Symbol"), it will be used to correct
 ### gene symbols.  Otherwise, the default table data("hgnc.table",
 ### package="HGNChelper") is used.
 ){
    if(class(x) != "character"){
        x <- as.character(x)
        warning("coercing x to character.")
    }
    if(is.null(hgnc.table)){
        data("hgnc.table", package="HGNChelper", envir=environment())
    }else if(!is(hgnc.table, "data.frame") | !identical(colnames(hgnc.table), c("Symbol", "Approved.Symbol"))){
        stop("If hgnc.table is specified, it must be a dataframe with two columns named 'Symbol' and 'Approved.Symbol'")
    }
    approved <- x %in% hgnc.table$Approved.Symbol
    ##change to upper case, then change orfs back to lower case:
    x.casecorrected <- toupper(x)
    x.casecorrected <- sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", x.casecorrected)
    approvedaftercasecorrection <- x.casecorrected %in% hgnc.table$Approved.Symbol
    if (!identical(all.equal(x, x.casecorrected), TRUE))
        warning("Some lower-case letters were found and converted to upper-case.
                 HGNChelper is intended for human symbols only, which should be all
                 upper-case except for open reading frames (orf).")
    alias <- x.casecorrected %in% hgnc.table$Symbol
    df <- data.frame(x=x,
                     Approved=approved,
                     Suggested.Symbol=sapply(1:length(x), function(i)
                         ifelse(approved[i],
                                x[i],
                                ifelse(alias[i],
                                       paste(hgnc.table$Approved.Symbol[x.casecorrected[i] == hgnc.table$Symbol], collapse=" /// "),
                                       ifelse(approvedaftercasecorrection[i],
                                              x.casecorrected[i],
                                              NA)))),
                     stringsAsFactors=FALSE)
    df$Approved[is.na(df$x)] <- FALSE
    if(!unmapped.as.na){
        df[is.na(df$Suggested.Symbol), "Suggested.Symbol"] <- df[is.na(df$Suggested.Symbol), "x"]
    }
    if (sum(df$Approved) != nrow(df))
        warning("x contains non-approved gene symbols")
    return(df)
### The function will return a data.frame of the same number of rows as the input,
### with corrections possible from hgnc.table.
}
