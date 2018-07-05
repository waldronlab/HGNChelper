data(hgnc.table)
hgnc.symbols <- as.character(na.omit(unique(hgnc.table[ ,2])))
expect_equal(hgnc.symbols, rToSymbol(make.names(symbolToR(hgnc.symbols))))

expect_equal(rToAffy(affyToR("123_at")), "123_at")

human <- read.csv(system.file("extdata/mog_map.csv", 
                              package = "HGNChelper"), as.is=TRUE)
re <- "[0-9]\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\.[0-9][0-9]E\\+[[0-9][0-9]"
expect_warning(res <- findExcelGeneSymbols(c("2-Apr", "APR2"), mog.map=human))
expect_equal(res, c("APR2", "APR2"))
expect_error(res <- findExcelGeneSymbols(c("2-Apr", "APR2"), mog.map=human[, 2:1]))
expect_warning(findExcelGeneSymbols("9-dec", mog.map=human, regex=re))

mouse <- read.csv(system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv", 
                              package = "HGNChelper"), as.is=TRUE)
expect_warning(res1 <- findExcelGeneSymbols(c("1-Feb", "1-dec"), mog.map=mouse))
expect_warning(res2 <- findExcelGeneSymbols(factor(c("1-Feb", "1-dec")), mog.map=mouse))
expect_identical(res1, res2)
