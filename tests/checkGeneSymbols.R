library(HGNChelper)

##test 1
x = c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
      "OCT4-PG4", "C19ORF71", "C19orf71")

res <- checkGeneSymbols(x)

stopifnot(sum(res[,1] != x) == 0)
stopifnot(sum(res[,2] !=
  c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)) == 0)

expected.correct <- c("FN1","TP53",NA,"SEPT7","SEPT7","MARC1 /// MARCH1",
                      "POU5F1","POU5F1","POU5F1P4","C19orf71","C19orf71")

if (!identical(all.equal(as.character(res[,3]), expected.correct), TRUE)){
  problem.res <- cbind(res, expected.correct)
  problem.res <- problem.res[na.omit(which(res[, 3] != expected.correct)), ]
  print(problem.res)
  stop("The above Suggested.Symbol does not match expected.correct")
}

##test 2
x = c("C21orf62-AS1", "c21orf62-as1",
      "MORF4L1P7", "Morf4L1P7",
      "Fn1", "FN1",
      "TP53",
      "UNKNOWNGENE",
      "7-Sep", "1-Mar", "1-MAR")

res <- checkGeneSymbols(x)
answers <- c(TRUE, FALSE,
             TRUE, FALSE,
             FALSE, TRUE,
             TRUE,
             FALSE,
             FALSE, FALSE, FALSE)
stopifnot(sum(res[,1] != x) == 0)
stopifnot(sum(res[,2] != answers) == 0)
stopifnot( all.equal(as.character(res[,3]), c("C21orf62-AS1", "C21orf62-AS1",
                                              "MORF4L1P7", "MORF4L1P7",
                                              "FN1", "FN1",
                                              "TP53",
                                              NA,
                                              "SEPT7", "MARC1 /// MARCH1", "MARC1 /// MARCH1"))
)
