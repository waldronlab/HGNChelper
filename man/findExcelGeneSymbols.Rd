\name{findExcelGeneSymbols}
\alias{findExcelGeneSymbols}
\seealso{checkGeneSymbols}
\title{function to identify Excel-mogrified gene symbols}
\description{This function identifies gene symbols which may have been
mogrified by Excel or other spreadsheet programs.  If output is
assigned to a variable, it returns a vector of the same length
where symbols which could be mapped have been mapped.  Note that this
function is superceded by checkGeneSymbols, which corrects
Excel-mogrified gene symbols as well as aliases and outdated symbols.}
\usage{findExcelGeneSymbols(x, mog.map = read.csv(system.file("extdata/mog_map.csv", 
    package = "HGNChelper"), as.is = TRUE), regex = 
"[0-9]\\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\\.[0-9]
[0-9]E\\\+[[0-9][0-9]")
}
\arguments{
  \item{x}{Vector of gene symbols to check for mogrified values}
  \item{mog.map}{Map of known mogrifications.  A default map is available with this
package by data(mog.map), but any map may be used.  This should be
a dataframe with two columns: original and mogrified, containing
the correct and incorrect symbols, respectively.}
  \item{regex}{Regular expression, recognized by the base::grep function which
is called with ignore.case=TRUE, to identify mogrified symbols.
It is not necessary for all matches to have a corresponding entry
in mog.map$mogrified; values in x which are matched by this regex
but are not found in mog.map$mogrified simply will not be
corrected.  This regex is based that provided by Zeeberg et al.,
Mistaken Identifiers: Gene name errors can be introduced
inadvertently when using Excel in bioinformatics.  BMC
Bioinformatics 2004, 5:80.}
}

\value{if the return value of the function is assigned to a variable, the
function will return a vector of the same length as the input,
with corrections possible from mog.map made.}

\author{Levi Waldron and Markus Riester}





