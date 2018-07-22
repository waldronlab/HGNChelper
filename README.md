[![Travis-CI Build Status](https://travis-ci.org/waldronlab/HGNChelper.svg?branch=master)](https://travis-ci.org/waldronlab/HGNChelper)
[![Coverage Status](https://codecov.io/github/waldronlab/HGNChelper/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/HGNChelper?branch=master)
[![](https://cranlogs.r-pkg.org/badges/HGNChelper)](https://cran.r-project.org/package=HGNChelper)

# HGNChelper: Identify and correct invalid gene symbols

## Updating

To update the symbols maps for human and mouse yourself, download this repository and run:

`./update.sh`

from its root directory. Note that this script uses the "roxygen2" 
R library to update the documentation.

Alternatively, you can use updated maps without updating the package, see `?getCurrentMaps`.

## Updating gh-pages

Note to self - when updating the vignette, update the gh-pages website 
(https://waldronlab.github.io/HGNChelper/) like this:

```
# pip install ghp-import
Rscript -e "devtools::build_vignettes()"
ghp-import inst/doc
git push
```
