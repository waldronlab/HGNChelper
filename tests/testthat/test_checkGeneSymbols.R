## human test 1
x = c("FN1", "TP53", "UNKNOWNGENE",
      "7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct", "OCT4-PG4", 
      "C19ORF71", "C19orf71")

expect_warning(res <- checkGeneSymbols(x))

expect_equal(sum(res[,1] != x), 0)
expect_equal(res[,2], c(TRUE, TRUE, FALSE,
                        FALSE, FALSE, FALSE, FALSE, FALSE,
                        FALSE, FALSE, TRUE))

expected.correct <- c("FN1", "TP53", NA,
                      "SEPTIN7", "SEPTIN7", "MARCHF1 /// MTARC1",
                      "POU5F1", "POU5F1", "POU5F1P4", "C19orf71", "C19orf71")

expect_equal(as.character(res[,3]), expected.correct)

expect_warning(res2 <- checkGeneSymbols(x, unmapped.as.na = FALSE))
expect_identical(res[, 1:2], res2[, 1:2])
unmapped <- is.na(res[, 3])
expect_identical(res2[unmapped, 3], res2[unmapped, 1])

## human test 2
x = c("C21orf62-AS1", "c21orf62-as1",
      "MORF4L1P7", "Morf4L1P7",
      "Fn1", "FN1",
      "TP53",
      "UNKNOWNGENE",
      "7-Sep", "1-Mar", "1-MAR")

expect_warning(res <- checkGeneSymbols(x))
answers <- c(TRUE, FALSE,
             TRUE, FALSE,
             FALSE, TRUE,
             TRUE,
             FALSE,
             FALSE, FALSE, FALSE)
expect_equal(sum(res[,1] != x), 0)
expect_equal(sum(res[,2] != answers), 0)
expect_equal(as.character(res[,3]), c("C21orf62-AS1", "C21orf62-AS1",
                                      "MORF4L1P7", "MORF4L1P7",
                                      "FN1", "FN1",
                                      "TP53",
                                      NA,
                                      "SEPTIN7", "MARCHF1 /// MTARC1", "MARCHF1 /// MTARC1"))


## human test 3 - chromosome
# check chromosome input
expect_warning(res <- checkGeneSymbols("Sip1"))
answer <- "GEMIN2 /// SCAF11 /// ZEB2"
expect_equal(res[,3], answer)
expect_error(checkGeneSymbols(rep("Sip1", 2), chromosome = c(14, 12, 2)),
             "The length of gene symbols and chromosome lists should be same.")

# wrong chromosome input
expect_warning(res <- checkGeneSymbols(rep("Sip1", 3), chromosome = c(14, 12, 3)), 
               "Specified chromosome list contains wrong chromosome number for some genes.") 
expect_equal(res[,2], c(FALSE, FALSE, FALSE))
expect_equal(res[,3], c("GEMIN2", "SCAF11", "GEMIN2 /// SCAF11 /// ZEB2"))
expect_equal(res[,4], c(14, 12, 3))
expect_equal(res[,5], c("14", "12", "14 /// 12 /// 2"))

# approved symbol alias mapping on different chromosome
expect_warning(res <- checkGeneSymbols(c("TKT", "FFM","ABCA4", "TKT", "AAVS1"), chromosome = c(1, 1, 1, 3, 19)), 
               "x contains non-approved gene symbols") 
expect_equal(res[,2], c(FALSE, FALSE, TRUE, TRUE, TRUE))
expect_equal(res[,3], c("DDR2", "ABCA4", "ABCA4", "TKT", "AAVS1"))
expect_equal(res[,4], c(1, 1, 1, 3, 19))
expect_equal(res[,5], c(1, 1, 1, 3, 19))

## human test 3 - expand.ambiguous
res <- checkGeneSymbols("AAVS1", expand.ambiguous = T)
expect_identical(res$Approved, TRUE)
expect_identical(res$Suggested.Symbol, "AAVS1 /// PPP1R12C")


## mouse test 1
data(mouse.table)
orig <- c("1-Feb", "A2m", "A2mr", "lrp1", "UNKNOWNGENE", "AI893533")
correct <- c("Feb1", "A2m", "Lrp", "Lrp1", NA, NA)
expect_warning(res <- checkGeneSymbols(orig, species="mouse"))
expect_warning(res2 <- checkGeneSymbols(orig, species="mouse", chromosome = c(1, 20)))
expect_equal(res$Approved, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
expect_equal(res$Suggested.Symbol, correct)
expect_equal(res2$Approved, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
expect_equal(res2$Suggested.Symbol, correct)

## mouse test 2
orig <- c("abl2", "AbLl", "Abl2",
          "Cpamd8", "cpamd8", 
          "mug2", "Mug2")
correct <- c("Abl2", "Abl2", "Abl2",
             "Cpamd8", "Mug2", 
             "Cpamd8 /// Mug2", "Mug2")
answer <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE)
expect_warning(res <- checkGeneSymbols(orig, species="mouse"))
expect_equal(res$Approved, answer)
expect_equal(res$Suggested.Symbol, correct)

## mouse test 3 - expand.ambiguous

res <- checkGeneSymbols(c("Cpamd8", "Mug2"), species = "mouse", expand.ambiguous = T)
expect_identical(res$Approved, c(TRUE, TRUE))
expect_identical(res$Suggested.Symbol, c("Cpamd8 /// Mug2", "Mug2 /// Cpamd8"))

# check capitalization behavior
orig <- c("tp53", "TP53")
data("hgnc.table", package="HGNChelper", envir = environment())
expect_warning(res <- checkGeneSymbols(orig, species=NULL, map=hgnc.table))
expect_equal(res$Suggested.Symbol, c(NA, "TP53"))
expect_warning(res <- checkGeneSymbols(orig, species="human", map=hgnc.table))
expect_equal(res$Suggested.Symbol, c("TP53", "TP53"))
expect_warning(res <- checkGeneSymbols(orig, map=hgnc.table))
expect_equal(res$Suggested.Symbol, c("TP53", "TP53"))

# check specifying manual map
mymap <- data.frame(Symbol=c("A", "B", "BB"), Approved.Symbol=c("A", "BB", "BB"), stringsAsFactors = FALSE)

# with no human capitalization help
expect_warning(res <- checkGeneSymbols(c("a", "b", "A", "B", "BB"), map=mymap, species=NULL))
expect_equal(res$Approved, c(FALSE, FALSE, TRUE, FALSE, TRUE))
expect_equal(res$Suggested.Symbol, c(NA, NA, "A", "BB", "BB"))

# with human capitalization help
expect_warning(res <- checkGeneSymbols(c("a", "b", "A", "B", "BB"), map=mymap, species="human"))
expect_equal(res$Approved, c(FALSE, FALSE, TRUE, FALSE, TRUE))
expect_equal(res$Suggested.Symbol, c("A", "BB", "A", "BB", "BB"))

expect_error(checkGeneSymbols(c("a", "b", "A", "B", "BB"), map=as.matrix(mymap)))
expect_error(checkGeneSymbols(c("a", "b", "A", "B", "BB"), map=NULL, species=NULL))
expect_warning(res2 <- checkGeneSymbols(factor(c("a", "b", "A", "B", "BB")), map=mymap, species="human"))
expect_identical(res, res2)

# check for outdated symbols from extdata/mog_map.csv
expect_warning(res3 <- checkGeneSymbols(c("MARC1", "MARC2", "MARCH6")))
expect_equal(res3$Approved, rep(FALSE, 3))
expect_equal(res3$Suggested.Symbol, c("MTARC1", "MTARC2", "MARCHF6"))

