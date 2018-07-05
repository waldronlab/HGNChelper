# HGNChelper: Identify and correct invalid gene symbols

[![Travis-CI Build Status](https://travis-ci.org/waldronlab/HGNChelper.svg?branch=master)](https://travis-ci.org/waldronlab/HGNChelper)
[![Coverage Status](https://codecov.io/github/waldronlab/HGNChelper/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/HGNChelper?branch=master)

To update the symbols map (data/hgnc.table.rda):

1. Download this repository, change to the "inst/" subdirectory
2. Run from the command line:
```
cd inst/
Rscript --vanilla "Rscript hgncLookup.R"
R CMD INSTALL ..
```
