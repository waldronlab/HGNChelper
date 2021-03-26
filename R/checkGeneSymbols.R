#' @title Identify outdated or Excel-mogrified gene symbols
#' @description This function identifies gene symbols which are outdated or may 
#' have been mogrified by Excel or other spreadsheet programs. If output is assigned 
#' to a variable, it returns a data.frame of the same number of rows as the input, 
#' with a second column indicating whether the symbols are valid and a third column 
#' with a corrected gene list.
#'
#' @importFrom methods is
#' @importFrom utils read.csv data
#' 
#' @param x A character vector of gene symbols to check for modified or outdated values
#' @param chromosome An optional integer vector containing the chromosome number of each gene
#' provided through the argument \code{x}. It should be the 
#' same length as the input for \code{x}. Currently, this argument is implemented
#' only for human gene cases.
#' @param unmapped.as.na If \code{TRUE} (default), unmapped symbols will appear as 
#' NA in the \code{Suggested.Symbol} column. If \code{FALSE}, the original unmapped 
#' symbol will be kept.
#' @param map Specify if you do not want to use the default maps provided by setting 
#' species equal to "mouse" or "human". Map can be any other data.frame with colnames 
#' identical to \code{c("Symbol", "Approved.Symbol")}. The default maps can be updated 
#' by running the interactive example below.
#' @param species A character vector of length 1, either "human" (default) or "mouse". 
#' If \code{NULL}, or anything other than "human" or "mouse", then the map argument 
#' must be provided. 
#'  
#' @return The function will return a data.frame of the same number of rows as the 
#' input, with corrections possible from map.
#' 
#' @seealso \code{\link{mouse.table}} for the mouse lookup table, \code{\link{hgnc.table}} for the human lookup table
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
#' @export
checkGeneSymbols <- function(x,
                             chromosome = NULL,
                             unmapped.as.na = TRUE,
                             map = NULL,
                             species = "human", 
                             expand.ambiguous = F) {

  lastupdate <- readLines(system.file(file.path("extdata", "date_of_last_update.txt"), 
                          package = "HGNChelper"))
  
  # check input class
  if (class(x) != "character") {
    x <- as.character(x)
    warning("coercing x to character.")
  }
  
  # chromosome argument is used only for human
  if (!identical(species, "human")) {chromosome <- NULL}
  if (!is.null(chromosome) & (length(x) != length(chromosome))) {
    stop("The length of gene symbols and chromosome lists should be same.")
  }
  
  # load map for correct species
  if (identical(species, "human") & is.null(map)) {
    message(paste("Maps last updated on:", lastupdate, collapse = " "))
    map <- HGNChelper::hgnc.table
  } else if (identical(species, "mouse") & is.null(map)) {
    message(paste("Maps last updated on:", lastupdate, collapse = " "))
    map <- HGNChelper::mouse.table
  } else {
    if (is.null(map)) {
      stop("If species is not 'human' or 'mouse' then map argument must be specified")
    }
  }
  
  # check chromosome input
  if (!is.null(chromosome)) {
    map <- map
    if (!is(map, "data.frame") | !identical(colnames(map), c("Symbol", "Approved.Symbol", "chromosome")))
      stop("If map is specified, it must be a dataframe with three columns named 'Symbol', 'Approved.Symbol' and 'chromosome'") 
  }
  
  else {
    map <- map[, 1:2]
    if (!is(map, "data.frame") | !identical(colnames(map), c("Symbol", "Approved.Symbol")))
      stop("If map is specified, it must be a dataframe with two columns named 'Symbol' and 'Approved.Symbol'")
  }
  
  approved <- x %in% map$Approved.Symbol
  
  if (!is.null(chromosome)) {
    approved.chr <- sapply(1:length(chromosome), function(i)
      ifelse(chromosome[i] %in% unique(map$chromosome[map$Approved.Symbol %in% x[i]]), 
             TRUE, 
             FALSE))
    approved <- approved & approved.chr
  }
  
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
  
  # wrong chromosome input
  if (!is.null(chromosome)) {
    chromosome.check <- sapply(1:length(chromosome), function(i) 
      ifelse(x.casecorrected[i] %in% map$Symbol, 
             sum(chromosome[i] == map$chromosome[map$Symbol %in% x.casecorrected[i]]) > 0, 
             FALSE))
    
    if (sum(chromosome.check) != length(chromosome))
      warning("Specified chromosome list contains wrong chromosome number for some genes.")
  } else chromosome.check <- FALSE
  
  approvedaftercasecorrection <- x.casecorrected %in% map$Approved.Symbol
  if (!identical(all.equal(x, x.casecorrected), TRUE) & identical(species, "human"))
    warning("Human gene symbols should be all upper-case except for the 'orf' in open reading frames. The case of some letters was corrected.")
  alias <- x.casecorrected %in% map$Symbol
  df <- data.frame(x = x,
                   Approved = approved,
                   Suggested.Symbol = sapply(1:length(x), function(i)
                     ifelse(approved[i],
                            ifelse(expand.ambiguous, 
                                   ifelse(!(is.null(chromosome)) & chromosome.check[i], 
                                          paste(unique(c(x[i], map$Approved.Symbol[x[i] == map$Symbol
                                                                                   & chromosome[i] == map$chromosome])), 
                                                collapse = " /// "),
                                          paste(unique(c(x[i], map$Approved.Symbol[x[i] == map$Symbol])), collapse = " /// ")),
                                   x[i]), 
                            ifelse(alias[i],   # format chromosome output
                                   ifelse(!(is.null(chromosome)) & chromosome.check[i], 
                                          paste(map$Approved.Symbol[x.casecorrected[i] == map$Symbol
                                                                    & chromosome[i] == map$chromosome], collapse = " /// "),
                                          paste(map$Approved.Symbol[x.casecorrected[i] == map$Symbol], collapse = " /// ")), 
                                   ifelse(approvedaftercasecorrection[i],
                                          x.casecorrected[i],
                                          NA)))),
                   stringsAsFactors=FALSE)
  df$Approved[is.na(df$x)] <- FALSE
  if (!is.null(chromosome)) {   # format chromosome output
    df$Input.chromosome = chromosome
    df$Correct.chromosome = sapply(1:length(x), function(i)
      ifelse(chromosome.check[i], 
             chromosome[i], 
             ifelse(alias[i], 
                    paste(map$chromosome[x.casecorrected[i] == map$Symbol], collapse = " /// "), 
                    NA)))
  }
  if(!unmapped.as.na) {
    df[is.na(df$Suggested.Symbol), "Suggested.Symbol"] <- df[is.na(df$Suggested.Symbol), "x"]
  }
  if (sum(df$Approved) != nrow(df))
    warning("x contains non-approved gene symbols")
  return(df)
}
