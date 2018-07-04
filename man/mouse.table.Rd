\name{mouse.table}
\alias{mouse.table}
\docType{data}
\title{
All current and withdrawn MOUSE symbols.
}
\description{
All current and withdrawn symbols from genenames.org.
}
\format{
  A dataframe with the first column providing a gene symbol or known alias
  (including withdrawn symbols), second column
  providing the approved MOUSE symbol.
}
\source{
Extracted from http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt and 
inst/extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv
}
\examples{
data("mouse.table", package="HGNChelper", envir=environment())
head(mouse.table)
}
\keyword{datasets}
