#' @name getCurrentMaps
#' @title Get the current maps for correcting gene symbols
#' @aliases getCurrentHumanMap getCurrentMouseMap
#' @description Valid human and mouse gene symbols can be updated frequently. Use
#' these functions to get the most current lists of valid symbols, which you can then use
#' as input to the "map" argument of checkGeneSymbols(). Make sure to change the default 
#' species="human" argument to checkGeneSymbols() if you are doing this for mouse.
#' getCurrentHumanMap() for HGNC human gene symbols from genenames.org
#' getCurrentMouseMap() for MGI mouse gene symbols from www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt
#' @return a `data.frame` that can be used as the checkGeneSymbols "map" argument
#' @export getCurrentMouseMap getCurrentHumanMap
#' @usage 
#' getCurrentHumanMap()
#' getCurrentMouseMap()
#' @importFrom stats complete.cases
#' @importFrom utils read.delim
#' @examples
#' \dontrun{
#' ## human
#' new.hgnc.table <- getCurrentHumanMap()
#' checkGeneSymbols(c("3-Oct", "10-3", "tp53"), map=new.hgnc.table)
#' ## mouse
#' new.mouse.table <- getCurrentMouseMap()
#' ## Set species to NULL or "mouse" 
#' ## so that human-like capitalization corrections aren't made
#' checkGeneSymbols(c("Gm46568", "1-Feb"), map=new.mouse.table, species="mouse")
#' }
#' 
getCurrentHumanMap <- function(){
  .fixttable <- function(hgnc.table) {
    ## remove withdrawn symbols with known new name
    hgnc.table <-
      hgnc.table[!(duplicated(hgnc.table$Symbol) &
                     is.na(hgnc.table$Approved.Symbol)), ]
    hgnc.table <- hgnc.table[order(hgnc.table$Symbol), ]
    
    hgnc.table$Symbol <- as.character(hgnc.table$Symbol)
    hgnc.table$Approved.Symbol <-
      as.character(hgnc.table$Approved.Symbol)
    rownames(hgnc.table) <- NULL
    
    ## In the un-approved column, convert everything but orfs to upper-case:
    hgnc.table$Symbol <- toupper(hgnc.table$Symbol)
    hgnc.table$Symbol <-
      sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", hgnc.table$Symbol)
    
    hgnc.table <- unique(hgnc.table)
    hgnc.table <- hgnc.table[complete.cases(hgnc.table), ]
    is.ascii <-
      iconv(hgnc.table[, 1], to = "ASCII", sub = ".") == hgnc.table[, 1]
    hgnc.table <- hgnc.table[is.ascii, ]
    return(hgnc.table)
  }
  url <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
  message(paste("Fetching gene symbols from", url))
  map <- read.delim(url, as.is=TRUE)
  # 2 = symbol, 3=alias_symbol, 4=prev_symbol
  # corrections
  has.corrections <- nchar(map$alias_symbol)>0 | nchar(map$prev_symbol)>0
  M  <- do.call(rbind, apply(map[has.corrections & map$status == "Approved", ], 1,
                             function(x) {
                               y <- strsplit(paste(x[c("alias_symbol", "prev_symbol")], collapse="|"), "|", fixed = TRUE)[[1]]
                               cbind(y[y!=""], x["symbol"])
                             }))
  rownames(M) <- NULL
  # valid symbols
  N <- map[, c("symbol", "symbol", "status")]
  N[N$status=="Entry Withdrawn", 2] <- NA
  N <- N[, 1:2]

  O = read.csv(system.file("extdata/mog_map.csv", package = "HGNChelper"),
               as.is=TRUE)[, 2:1]
  O$mogrified <- toupper(O$mogrified)
  
  colnames(M) <- c("Symbol", "Approved.Symbol")
  colnames(N) <- c("Symbol", "Approved.Symbol")
  colnames(O) <- c("Symbol", "Approved.Symbol")
  
  external.table <- rbind(M, N)
  ## Correct outdated symbols in extdata/mog_map.csv
  O_corrected <-
    suppressWarnings(checkGeneSymbols(O[, "Approved.Symbol"], map = external.table))
  O[, "Approved.Symbol"] <- O_corrected[, "Suggested.Symbol"]
  output <- rbind(external.table, O)
  output <- .fixttable(output)
  return(output)
}

getCurrentMouseMap <- function(){
  fname <- system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv", package = "HGNChelper")
  O = read.csv(fname, as.is=TRUE)[, 2:1]
  colnames(O) <- c("Symbol", "Approved.Symbol")
  
  map <- read.delim("http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt", as.is=TRUE, header = FALSE)
  map <- map[, 2:4]
  
  map.ok <- map[-grep("^O$", map[, 2]), c(1, 3)]
  map.ok <- data.frame(Symbol=unique(map.ok[, 1]), Approved.Symbol=unique(map.ok[, 1]), stringsAsFactors = FALSE)
  
  map.withdrawn <- map[grep("^W$", map[, 2]), c(1, 3)]
  map.withdrawn[, 2] <- sub("^withdrawn, ", "", map.withdrawn[, 2])
  map.withdrawn[, 2] <- sub("^= ", "", map.withdrawn[, 2])
  map.withdrawn <- map.withdrawn[, 2:1]
  colnames(map.withdrawn) <- c("Symbol", "Approved.Symbol")
  
  mouse.table <- rbind(O, map.withdrawn, map.ok)
  mouse.table <- unique(mouse.table)
  mouse.table <- mouse.table[complete.cases(mouse.table), ]
  rownames(mouse.table) <- NULL
  mouse.table <- mouse.table[!is.na(iconv(mouse.table[, 1], "ASCII")), ]
  return(mouse.table)
}