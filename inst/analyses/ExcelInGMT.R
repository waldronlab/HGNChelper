## Check .gmt files in the same directory for Excel-mogrified gene symbols.

library(HGNChelper)

readGMT <- function (fname, name.column = 1){
    if (!(name.column == 1 | name.column == 2))
        stop("name.column should be 1 or 2")
    x <- readLines(fname)
    x <- strsplit(x, split = "\t")
    output <- lapply(x, function(y) y[-1:-2])
    names(output) <- make.names(sapply(x, function(y) y[name.column]),
        unique = TRUE)
    return(output)
}

gmt.files <- dir(pattern="\\.gmt$")

if(length(gmt.files) == 0) stop(paste("No .gmt files found in", getwd()))

all.sigs <- lapply(gmt.files, readGMT)
all.sigs <- lapply(all.sigs, function(x) unique(unlist(x)))
##Remove LOC symbols:
all.sigs <- lapply(all.sigs, function(x) x[!grepl("^LOC", x)])
names(all.sigs) <- gmt.files
all.sigs$OVERALL <- unique(do.call(c, all.sigs))

mog <- read.csv(system.file("extdata", "mog_map.csv",package="HGNChelper"), as.is=TRUE)

simplecheck <- sapply(all.sigs, function(x){
        tmp <- checkGeneSymbols(x)
        output <- c(sum(tmp$Approved)/nrow(tmp), sum(!is.na(tmp[, 3]))/nrow(tmp), any(toupper(x) %in% toupper(mog[, 2])))
        names(output) <- c("valid.frac", "valid.after.hgnchelper.frac", "any.excel")
        return(output)
    })
simplecheck <- t(simplecheck)

simplecheck

unfixables <- sort(tmp[is.na(tmp[, 3]), 1])  ##what symbols can't be corrected?
set.seed(1)
sample(unfixables, 10)
## Genecards searches
## [1] "FLJ16126"  "FLJ37638"  "HSCRS2"    "SELO"      "DNM1DN4-1" "SCAMPER"  
## [7] "SPPM"      "MCOR"      "KTCN3"     "AUTS5"    
