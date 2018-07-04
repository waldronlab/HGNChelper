library(HGNChelper)

data(hgnc.table)
hgnc.symbols <- as.character(na.omit(unique(hgnc.table[ ,2])))
expect_equal(hgnc.symbols, rToSymbol(make.names(symbolToR(hgnc.symbols))))
