# All current and withdrawn MGI mouse symbols and Excel-mogrified symbols

A `data.frame` with the first column providing a gene symbol or known
alias (including withdrawn symbols), second column providing the
approved MGI mouse gene symbol.

- `Symbol`: All valid, Excel-mogrified, and withdrawn symbols

- `Approved.Symbol`: Approved symbols

## Usage

``` r
mouse.table
```

## Format

An object of class `data.frame` with 790110 rows and 2 columns.

## Source

Extracted from
<http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt>
and system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv",
package="HGNChelper")

## Examples

``` r
data("mouse.table", package="HGNChelper")
head(mouse.table)
#>   Symbol Approved.Symbol
#> 1  1-Feb            Feb1
#> 2    2/1            Feb1
#> 3  02/01            Feb1
#> 4   2/01            Feb1
#> 5   02/1            Feb1
#> 6   02-1            Feb1
```
