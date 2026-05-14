# Title function to identify Excel-mogrified gene symbols

This function identifies gene symbols which may have been mogrified by
Excel or other spreadsheet programs. If output is assigned to a
variable, it returns a vector of the same length where symbols which
could be mapped have been mapped.

## Usage

``` r
findExcelGeneSymbols(
  x,
  mog.map = read.csv(system.file("extdata/mog_map.csv", package = "HGNChelper"), as.is =
    TRUE),
  regex = "impossibletomatch^"
)
```

## Arguments

- x:

  Vector of gene symbols to check for mogrified values

- mog.map:

  Map of known mogrifications. This should be a dataframe with two
  columns: original and mogrified, containing the correct and incorrect
  symbols, respectively.

- regex:

  Regular expression, recognized by the base::grep function which is
  called with ignore.case=TRUE, to identify mogrified symbols. The
  default regex will not match anything. The regex in the examples is an
  attempt to match all Excel-mogrified HGNC human gene symbols. It is
  not necessary for all matches to have a corresponding entry in
  mog.map\$mogrified; values in x which are matched by this regex but
  are not found in mog.map\$mogrified simply will not be corrected.

## Value

if the return value of the function is assigned to a variable, the
function will return a vector of the same length as the input, with
corrections possible from mog.map made.

## Examples

``` r
## Available maps from this package:
human <- read.csv(system.file("extdata/mog_map.csv", 
                              package = "HGNChelper"), as.is=TRUE)
mouse <- read.csv(system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv", 
                              package = "HGNChelper"), as.is=TRUE)
## This regex is based that provided by Zeeberg et al.,
##  Mistaken Identifiers: Gene name errors can be introduced
## inadvertently when using Excel in bioinformatics.  BMC
##  Bioinformatics 2004, 5:80.
re <- "[0-9]\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\.[0-9][0-9]E\\+[[0-9][0-9]"
findExcelGeneSymbols(c("2-Apr", "APR2"), mog.map=human, regex=re)
#> Warning: Transmogrified gene symbols found.  Returning the following corrections: 2-Apr to 2-Apr
#> [1] "2-Apr" "APR2" 
findExcelGeneSymbols(c("1-Feb", "Feb1"), mog.map=mouse)
#> Warning: Transmogrified gene symbols found.  Returning the following corrections: 1-Feb to Feb1
#> [1] "Feb1" "Feb1"
```
