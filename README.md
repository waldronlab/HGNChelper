[![R build status](https://github.com/waldronlab/HGNChelper/workflows/R-CMD-check/badge.svg)](https://github.com/waldronlab/HGNChelper/actions)
[![](https://cranlogs.r-pkg.org/badges/HGNChelper)](https://cran.r-project.org/package=HGNChelper)
[![Coverage Status](https://codecov.io/github/waldronlab/HGNChelper/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/HGNChelper?branch=master)
[![Travis-CI Build Status](https://travis-ci.org/waldronlab/HGNChelper.svg?branch=master)](https://travis-ci.org/waldronlab/HGNChelper)
[![Zenodo release](https://zenodo.org/badge/139589811.svg)](https://zenodo.org/account/settings/github/repository/waldronlab/HGNChelper)


# HGNChelper: identification and correction of invalid gene symbols for human and mouse

A [pre-print is available on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.09.16.300632v2). 

## Updating

To update the symbols maps for human and mouse yourself, download this repository and run:

`./update.sh`

from its root directory. Note that this script uses the "roxygen2" 
R library to update the documentation.

Alternatively, you can use updated maps without updating the package, see `?getCurrentMaps`.

