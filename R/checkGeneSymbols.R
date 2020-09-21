#' @title Identify outdated or Excel-mogrified gene symbols
#' @description 
#'  This function identifies gene symbols which are outdated or may have been
#'  mogrified by Excel or other spreadsheet programs.  If output is
#'  assigned to a variable, it returns a data.frame of the same number of
#'  rows as the input, with a second column indicating whether the symbols
#'  are valid and a third column with a corrected gene list.
#'
#' @param x A character vector of gene symbols to check for mogrified or outdated values
#' @param unmapped.as.na 
#'  If \code{TRUE} (default), unmapped symbols will appear as NA in the
#'  \code{Suggested.Symbol} column. If \code{FALSE}, the original unmapped symbol
#'  will be kept.
#' @param map 
#'  Specify if you do not want to use the default maps provided 
#'  by setting species equal to "mouse" or "human". 
#'  Map can be any other data.frame with colnames identical
#'  to \code{c("Symbol", "Approved.Symbol")}. The default maps can be updated 
#'  by running the interactive example below.
#' @param species
#'  A character vector of length 1, either "human" (default) or "mouse". 
#'  If \code{NULL}, or anything other than "human" or "mouse", then the map
#'  argument must be provided. 
#'  
#' @return
#' The function will return a data.frame of the same number of rows as the input,
#' with corrections possible from map.
#' @seealso \code{\link{mouse.table}} for the mouse lookup table, 
#' \code{\link{hgnc.table}} for the human lookup table
#' @importFrom methods is
#' @importFrom utils read.csv data
#'
#' @examples
#' library(HGNChelper)
#' 
#' ## Human
#' human <- c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
#'       "OCT4-PG4", "C19ORF71", "C19orf71")
#' checkGeneSymbols(human)
#' 
#' ## Mouse
#' mouse <- c("1-Feb", "Pzp", "A2m")
#' checkGeneSymbols(mouse, species="mouse")
#' 
#' ## Updating the map
#' if (interactive()) {
#'     currentHumanMap <- getCurrentHumanMap()
#'     checkGeneSymbols(human, map=currentHumanMap)
#' 
#'     # You should save this if you are going to use it multiple times,   
#'     # then load it from file rather than burdening HGNC's servers.
#'     save(hgnc.table, file="hgnc.table.rda", compress="bzip2")
#'     load("hgnc.table.rda")
#'     checkGeneSymbols(human, map=hgnc.table)
#' }
#' 
#' ## Unit Tests
#' currentHumanMap <- getCurrentHumanMap()
#' 
#' ### Line 108 to 116
#' checkGeneSymbols("Sip1")  # should work as before
#' checkGeneSymbols("Sip1", chromosome=14) # should give error
#' checkGeneSymbols("Sip1", chromosome=14, map=currentHumanMap) # should give approved symbol only for chromosome 14
#' 
#' ### Line 141 to 149
#' checkGeneSymbols(rep("Sip1", 3), chromosome=c(14, 12, 3), map=currentHumanMap) # output with a warning for wrong chromosome number
#' checkGeneSymbols(rep("Sip1", 3), chromosome=c(14, 12, 2), map=currentHumanMap) # output without a warning for wrong chromosome number
#' 
#' ### Line 155 to 168
#' checkGeneSymbols("Sip1") # output as before  
#' checkGeneSymbols(rep("Sip1", 3), chromosome=c(14, 12, 2), map=currentHumanMap) # one suggested symbol for correctly specified chromosome numbers
#' checkGeneSymbols(human, c(1:5, 1, 7:11), map = currentHumanMap) # returns MTARC1 (at chromosome 1) for 1-MAR
#' checkGeneSymbols(human, c(1:5, 4, 7:11), map = currentHumanMap) # returns MARCHF1 (at chromosome 4) for 1-MAR
#' 
#'  
#'   
#' @export
checkGeneSymbols <- function(x,
                             chromosome = NULL,
                             unmapped.as.na = TRUE,
                             map = NULL,
                             species = "human") {

  lastupdate <- readLines(system.file(file.path("extdata", "date_of_last_update.txt"), 
                          package = "HGNChelper"))
  
  # check input class
  if (class(x) != "character") {
    x <- as.character(x)
    warning("coercing x to character.")
  }
  
  # load map for correct species
  if (identical(species, "human") & is.null(map)) {
    message(paste("Maps last updated on:", lastupdate, collapse = " "))
    if (!is.null(chromosome)) {
      map <- HGNChelper::hgnc.table
    } else map <- HGNChelper::hgnc.table[, 1:2]
  } else if (identical(species, "mouse") & is.null(map)) {
    message(paste("Maps last updated on:", lastupdate, collapse = " "))
    map <- HGNChelper::mouse.table
  } else {
    if (is.null(map)) {
      stop("If species is not 'human' or 'mouse' then map argument must be specified")
    }
  }
  
  if (!is.null(chromosome)) {
    if (!is(map, "data.frame") | !identical(colnames(map), c("Symbol", "Approved.Symbol", "chromosome")))
      stop("If map is specified, it must be a dataframe with three columns named 'Symbol', 'Approved.Symbol' and 'chromosome'") 
  }
  
  else {
    if (!is(map, "data.frame") | !identical(colnames(map), c("Symbol", "Approved.Symbol")))
      stop("If map is specified, it must be a dataframe with two columns named 'Symbol' and 'Approved.Symbol'")
  }
  
  approved <- x %in% map$Approved.Symbol
  
  if (identical(species, "human")) {
    # change to uppercase, then change orfs back to lowercase
    x.casecorrected <- toupper(x)
    x.casecorrected <- sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", x.casecorrected)
  } else if (identical(species, "mouse")) {
    x.casecorrected <- x
    for (i in seq_along(x)) {
      # don't correct if gene symbol starts with non-alphabet
      alphabetPrefix <- substring(x[i], 1, 1) %in% c(LETTERS, letters)
      
      # make it begin with an uppercase, followed by all lowercase
      if (isTRUE(alphabetPrefix)) {
        x.casecorrected[i] <- tolower(x[i])
        x.casecorrected[i] <- paste0(toupper(substr(x.casecorrected[i], 1, 1)), 
                                  substr(x.casecorrected[i], 2, nchar(x.casecorrected[i])))
      }
    }
  } else {
    x.casecorrected <- x
  }
  
  if (!is.null(chromosome)) {
    chromosome.check <- sapply(1:length(chromosome), function(i) 
      ifelse(x.casecorrected[i] %in% map$Symbol, 
             sum(chromosome[i] == map$chromosome[map$Symbol %in% x.casecorrected[i]])>0, 
             FALSE))
    
    if (sum(chromosome.check) != length(chromosome))
      warning("chromosome contains wrong chromosome number for some genes.")
  } else chromosome.check <- F
  
  approvedaftercasecorrection <- x.casecorrected %in% map$Approved.Symbol
  if (!identical(all.equal(x, x.casecorrected), TRUE) & identical(species, "human"))
    warning("Human gene symbols should be all upper-case except for the 'orf' in open reading frames. The case of some letters was corrected.")
  alias <- x.casecorrected %in% map$Symbol
  df <- data.frame(x=x,
                   Approved=approved,
                   Suggested.Symbol=sapply(1:length(x), function(i)
                     ifelse(approved[i],
                            x[i],
                            ifelse(alias[i],
                                   ifelse(!(is.null(chromosome)) & chromosome.check[i], 
                                          paste(map$Approved.Symbol[x.casecorrected[i] == map$Symbol
                                                                    & chromosome[i] == map$chromosome], collapse=" /// "),
                                          paste(map$Approved.Symbol[x.casecorrected[i] == map$Symbol], collapse=" /// ")), 
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
