\name{rToSymbol}
\alias{rToSymbol}
\title{function to reverse the conversion made by symbolToR}
\description{This function reverses the actions of the symbolToR function.}
\usage{rToSymbol(x)}
\arguments{
  \item{x}{the character vector returned by the symbolToR function.}
}

\value{a character vector of HGNC gene symbols, which are not in general valid R names.}

\author{Levi Waldron and Markus Riester}




\examples{library(HGNChelper)

data("hgnc.table", envir=environment())
hgnc.symbols <- as.character(na.omit(unique(hgnc.table[ ,2])))
if( !identical(all.equal(hgnc.symbols, rToSymbol(make.names(symbolToR(hgnc.symbols)))), TRUE))
     stop("HGNC mapping was not reversible.")}
