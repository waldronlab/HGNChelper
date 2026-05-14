# Identify outdated or Excel-mogrified gene symbols

This function identifies gene symbols which are outdated or may have
been mogrified by Excel or other spreadsheet programs. If output is
assigned to a variable, it returns a data.frame of the same number of
rows as the input, with a second column indicating whether the symbols
are valid and a third column with a corrected gene list.

## Usage

``` r
checkGeneSymbols(
  x,
  chromosome = NULL,
  unmapped.as.na = TRUE,
  map = NULL,
  species = "human",
  expand.ambiguous = FALSE
)
```

## Arguments

- x:

  A character vector of gene symbols to check for modified or outdated
  values

- chromosome:

  An optional integer vector containing the chromosome number of each
  gene provided through the argument `x`. It should be the same length
  as the input for `x`. Currently, this argument is implemented only for
  human gene cases.

- unmapped.as.na:

  If `TRUE` (default), unmapped symbols will appear as NA in the
  `Suggested.Symbol` column. If `FALSE`, the original unmapped symbol
  will be kept.

- map:

  Specify if you do not want to use the default maps provided by setting
  species equal to "mouse" or "human". Map can be any other data.frame
  with colnames identical to `c("Symbol", "Approved.Symbol")`. The
  default maps can be updated by running the interactive example below.

- species:

  A character vector of length 1, either "human" (default) or "mouse".
  If `NULL`, or anything other than "human" or "mouse", then the map
  argument must be provided.

- expand.ambiguous:

  If `FALSE` (default), genes with multiple mapping will only map to its
  approved symbol as the correct one. If `TRUE`, genes with
  multiple/ambiguous mapping will map to all the symbols linked to it.

## Value

The function will return a data.frame of the same number of rows as the
input, with corrections possible from map.

## See also

[`mouse.table`](https://waldronlab.io/HGNChelper/reference/mouse.table.md)
for the mouse lookup table,
[`hgnc.table`](https://waldronlab.io/HGNChelper/reference/hgnc.table.md)
for the human lookup table

## Examples

``` r
library(HGNChelper)

## Human
human <- c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
      "OCT4-PG4", "C19ORF71", "C19orf71")
checkGeneSymbols(human)
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> Warning: Human gene symbols should be all upper-case except for the 'orf' in open reading frames. The case of some letters was corrected.
#> Warning: x contains non-approved gene symbols
#>              x Approved   Suggested.Symbol
#> 1          FN1     TRUE                FN1
#> 2         TP53     TRUE               TP53
#> 3  UNKNOWNGENE    FALSE               <NA>
#> 4        7-Sep    FALSE            SEPTIN7
#> 5          9/7    FALSE            SEPTIN7
#> 6        1-Mar    FALSE MARCHF1 /// MTARC1
#> 7         Oct4    FALSE             POU5F1
#> 8        4-Oct    FALSE             POU5F1
#> 9     OCT4-PG4    FALSE           POU5F1P4
#> 10    C19ORF71    FALSE            TEKTIP1
#> 11    C19orf71    FALSE            TEKTIP1

## Mouse
mouse <- c("1-Feb", "Pzp", "A2m")
checkGeneSymbols(mouse, species="mouse")
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> Warning: x contains non-approved gene symbols
#>       x Approved Suggested.Symbol
#> 1 1-Feb    FALSE             Feb1
#> 2   Pzp     TRUE              Pzp
#> 3   A2m     TRUE              A2m

## expand.ambiguous

## Human
human <- "AAVS1"
checkGeneSymbols(human, expand.ambiguous=FALSE)
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#>       x Approved Suggested.Symbol
#> 1 AAVS1     TRUE            AAVS1
checkGeneSymbols(human, expand.ambiguous=TRUE)
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#>       x Approved   Suggested.Symbol
#> 1 AAVS1     TRUE AAVS1 /// PPP1R12C

## Mouse
mouse <- c("Cpamd8", "Mug2")
checkGeneSymbols(mouse, species = "mouse", expand.ambiguous = FALSE)
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#>        x Approved Suggested.Symbol
#> 1 Cpamd8     TRUE           Cpamd8
#> 2   Mug2     TRUE             Mug2
checkGeneSymbols(mouse, species = "mouse", expand.ambiguous = TRUE)
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#>        x Approved Suggested.Symbol
#> 1 Cpamd8     TRUE  Cpamd8 /// Mug2
#> 2   Mug2     TRUE  Mug2 /// Cpamd8

## Updating the map
if (interactive()) {
    currentHumanMap <- getCurrentHumanMap()
    checkGeneSymbols(human, map=currentHumanMap)

    # You should save this if you are going to use it multiple times,   
    # then load it from file rather than burdening HGNC's servers.
    save(hgnc.table, file="hgnc.table.rda", compress="bzip2")
    load("hgnc.table.rda")
    checkGeneSymbols(human, map=hgnc.table)
}
  
```
