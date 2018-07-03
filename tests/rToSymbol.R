library(HGNChelper)

data(hgnc.table)
hgnc.symbols <- as.character(na.omit(unique(hgnc.table[ ,2])))
if( !identical(all.equal(hgnc.symbols, rToSymbol(make.names(symbolToR(hgnc.symbols)))), TRUE))
     stop("HGNC mapping was not reversible.")
