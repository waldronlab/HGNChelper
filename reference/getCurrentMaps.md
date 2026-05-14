# Get the current maps for correcting gene symbols

Valid human and mouse gene symbols can be updated frequently. Use these
functions to get the most current lists of valid symbols, which you can
then use as an input to the `map` argument of `checkGeneSymbols`. Make
sure to change the default `species="human"` argument to
`checkGeneSymbols` if you are doing this for mouse. Use
`getCurrentHumanMap` for HGNC human gene symbols from
<https://www.genenames.org/> and `getCurrentMouseMap` for MGI mouse gene
symbols from
<https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt>.

## Usage

``` r
getCurrentHumanMap()
getCurrentMouseMap()
```

## Value

A `data.frame` that can be used for `map` argument of `checkGeneSymbols`
function

## Examples

``` r
if (FALSE) { # \dontrun{
## human
new.hgnc.table <- getCurrentHumanMap()
checkGeneSymbols(c("3-Oct", "10-3", "tp53"), map=new.hgnc.table)

## mouse
new.mouse.table <- getCurrentMouseMap()
## Set species to NULL or "mouse" 
checkGeneSymbols(c("Gm46568", "1-Feb"), map=new.mouse.table, species="mouse")
} # }
```
