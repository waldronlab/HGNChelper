# Title function to convert Affymetrix probeset identifiers to valid R names

This function simply prepends "affy." to the probeset IDs to create
valid R names. Reverse operation is done by the
[rToAffy](https://waldronlab.io/HGNChelper/reference/rToAffy.md)
function.

## Usage

``` r
affyToR(x)
```

## Arguments

- x:

  vector of Affymetrix probeset identifiers, or any identifier which may
  with a digit.

## Value

a character vector that is simply x with "affy." prepended to each
value.
