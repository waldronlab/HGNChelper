# Title function to \*reversibly\* convert HGNC gene symbols to valid R names.

This function reversibly converts HGNC gene symbols to valid R names by
prepending "symbol.", and making the following substitutions: "-" to
"hyphen", "@" to "ampersand", and "/" to "forwardslash".

## Usage

``` r
symbolToR(x)
```

## Arguments

- x:

  vector of HGNC symbols

## Value

a vector of valid R names, of the same length as x, which can be
converted to the same HGNC symbols using the rToSymbol function.

## See also

[`rToSymbol`](https://waldronlab.io/HGNChelper/reference/rToSymbol.md)
