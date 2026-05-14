# All current and withdrawn HGNC gene symbols and Excel-mogrified symbols

A `data.frame` with the first column providing a gene symbol or known
alias (including withdrawn symbols), second column providing the
approved HGNC human gene symbol.

- `Symbol`: All valid, Excel-mogrified, and withdrawn symbols

- `Approved.Symbol`: Approved symbols

## Usage

``` r
hgnc.table
```

## Format

An object of class `data.table` (inherits from `data.frame`) with 103939
rows and 3 columns.

## Source

Extracted from
<https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt>
and system.file("extdata/mog_map.csv", package="HGNChelper")

## Examples

``` r
data("hgnc.table", package="HGNChelper")
head(hgnc.table)
#>        Symbol Approved.Symbol chromosome
#>        <char>          <char>     <char>
#> 1:       A1BG            A1BG         19
#> 2:   A1BG-AS1        A1BG-AS1         19
#> 3:   FLJ23569        A1BG-AS1         19
#> 4:     A1BGAS        A1BG-AS1         19
#> 5:    A1BG-AS        A1BG-AS1         19
#> 6: NCRNA00181        A1BG-AS1         19
```
