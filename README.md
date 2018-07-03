# HGNChelper: Identify and correct invalid gene symbols

To update the symbols map (data/hgnc.table.rda):

1. Download this repository, change to the "inst/" subdirectory
2. Run from the command line:
```
cd inst/
Rscript --vanilla "Rscript hgncLookup.R"
R CMD INSTALL ..
```

Alternatively, Travis-CI is set up for this repository to
automatically update the map with each push to the repository. If you
are logged into GitHub.com, you can click on the DESCRIPTION file,
click on the pencil icon ("Edit this file" in the top-right-hand side
of the file panel, increment the last number in the version (e.g. add
1 to n in 0.3.8.n), and submit a pull request. Once this pull request
is merged, it Travis will update the map (I hope, this is still under
testing).

[![Travis-CI Build Status](https://travis-ci.org/waldronlab/HGNChelper.svg?branch=master)](https://travis-ci.org/waldronlab/HGNChelper)
