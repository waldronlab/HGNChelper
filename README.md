[![R build status](https://github.com/waldronlab/HGNChelper/workflows/R-CMD-check/badge.svg)](https://github.com/waldronlab/HGNChelper/actions)
[![](https://cranlogs.r-pkg.org/badges/HGNChelper)](https://cran.r-project.org/package=HGNChelper)
[![Coverage Status](https://codecov.io/github/waldronlab/HGNChelper/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/HGNChelper?branch=master)
[![DOI](https://zenodo.org/badge/139589811.svg)](https://zenodo.org/badge/latestdoi/139589811)


# HGNChelper: identification and correction of invalid gene symbols for human and mouse

Please cite our software:

Oh S, Abdelnabi J, Al-Dulaimi R et al. HGNChelper: identification and
correction of invalid gene symbols for human and mouse. F1000Research 2020, 9:1493
(https://doi.org/10.12688/f1000research.28033.1)


## Updating

To update the symbols maps for human and mouse yourself, download this repository and run:

`./update.sh`

from its root directory. Note that this script uses the "roxygen2" 
R library to update the documentation.

Alternatively, you can use updated maps without updating the package, see `?getCurrentMaps`.

