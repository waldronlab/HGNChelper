\name{affyToR}
\alias{affyToR}
\title{function to convert Affymetrix probeset identifiers to valid R names}
\description{This function simply prepends "affy." to the probeset IDs to create
valid R names.  Reverse operation is done by the \link{rToAffy}
function.}
\usage{affyToR(x)}
\arguments{
  \item{x}{vector of Affymetrix probeset identifiers, or any identifier which may with a digit.}
}

\value{a character vector that is simply x with "affy." prepended to each value.}

\author{Levi Waldron and Markus Riester}





