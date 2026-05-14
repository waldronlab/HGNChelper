# Title function to convert the output of affyToR back to the original Affymetrix probeset identifiers.

This function simply strips the "affy." added by the
[affyToR](https://waldronlab.io/HGNChelper/reference/affyToR.md)
function.

## Usage

``` r
rToAffy(x)
```

## Arguments

- x:

  the character vector returned by the affyToR function.

## Value

a character vector of Affymetrix probeset identifiers.
