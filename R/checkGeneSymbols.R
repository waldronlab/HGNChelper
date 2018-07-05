#' Title Identify outdated or Excel-mogrified gene symbols
#'
#' @param x Vector of gene symbols to check for mogrified or outdated values
#' @param unmapped.as.na 
#'  If TRUE (default), unmapped symbols will appear as NA in the
#'  Suggested.Symbol column.  If FALSE, the original unmapped symbol
#'  will be kept.
#' @param map 
#'  Specify if you do not want to use the default maps provided 
#'  by setting species equal to "mouse" or "human". 
#'  map can be any other data.frame with colnames(map) identical
#'  to c("Symbol", "Approved.Symbol"). The default maps can be updated 
#'  by running the interactive example below.
#' @param species
#'  A character vector of length 1, either "human" (default) or "mouse". 
#'  If NULL, or anything other than "human" or "mouse", then the map
#'  argument must be provided. 
#' @description 
#'  This function identifies gene symbols which are outdated or may have been
#'  mogrified by Excel or other spreadsheet programs.  If output is
#'  assigned to a variable, it returns a data.frame of the same number of
#'  rows as the input, with a second column indicating whether the symbols
#'  are valid and a third column with a corrected gene list.
#' @return
#' The function will return a data.frame of the same number of rows as the input,
#' with corrections possible from map.
#' @seealso \code{\link{mouse.table}} for the mouse lookup table, 
#' \code{\link{hgnc.table}} for the human lookup table
#' @importFrom methods is
#' @importFrom utils read.csv data
#' @export
#'
#' @examples
#' library(HGNChelper)
#' human = c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
#'       "OCT4-PG4", "C19ORF71", "C19orf71")
#' checkGeneSymbols(human)
#' ## mouse
#' mouse <- c("1-Feb", "Pzp", "A2m")
#' checkGeneSymbols(mouse, species="mouse")
#' if (interactive()){
#'   ##Run checkGeneSymbols with a brand-new map downloaded from HGNC:
#'   source(system.file("hgncLookup.R", package = "HGNChelper"))
#'   ## You should save this if you are going to use it multiple times,
#'   ## then load it from file rather than burdening HGNC's servers.
#'   save(hgnc.table, file="hgnc.table.rda", compress="bzip2")
#'   load("hgnc.table.rda")
#'   checkGeneSymbols(human, species=NULL, map=hgnc.table)
#'   checkGeneSymbols(human, species=NULL, map=mouse.table)
#' }
checkGeneSymbols <-function(x,
                            unmapped.as.na=TRUE,
                            map=NULL,
                            species="human"
){
  if(class(x) != "character"){
    x <- as.character(x)
    warning("coercing x to character.")
  }
  casecorrection <- FALSE
  if(identical(species, "human")){
    casecorrection <- TRUE
    if(is.null(map)){
      map <- HGNChelper::hgnc.table
    }
  }else if(identical(species, "mouse") & is.null(map)){
    map <- HGNChelper::mouse.table
  }else{
    if(is.null(map)){
      stop("If species is not 'human' or 'mouse' then map argument must be specified")
    }
  }
  if(!is(map, "data.frame") | !identical(colnames(map), c("Symbol", "Approved.Symbol")))
    stop("If map is specified, it must be a dataframe with two columns named 'Symbol' and 'Approved.Symbol'")
  approved <- x %in% map$Approved.Symbol
  if(casecorrection){
    ##change to upper case, then change orfs back to lower case:
    x.casecorrected <- toupper(x)
    x.casecorrected <- sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", x.casecorrected)
  }else{
    x.casecorrected <- x
  }
  approvedaftercasecorrection <- x.casecorrected %in% map$Approved.Symbol
  if (!identical(all.equal(x, x.casecorrected), TRUE))
    warning("Some lower-case letters were found and converted to upper-case.
                 HGNChelper is intended for human symbols only, which should be all
                 upper-case except for open reading frames (orf).")
  alias <- x.casecorrected %in% map$Symbol
  df <- data.frame(x=x,
                   Approved=approved,
                   Suggested.Symbol=sapply(1:length(x), function(i)
                     ifelse(approved[i],
                            x[i],
                            ifelse(alias[i],
                                   paste(map$Approved.Symbol[x.casecorrected[i] == map$Symbol], collapse=" /// "),
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
}
