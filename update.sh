#!/bin/bash
Rscript --vanilla -e "hgnc.table <- HGNChelper::getCurrentHumanMap(); \
save(hgnc.table, file='data/hgnc.table.rda', compress='bzip2'); \
mouse.table <- HGNChelper::getCurrentMouseMap(); \
save(mouse.table, file='data/mouse.table.rda', compress='bzip2'); \
writeLines(date(), file.path('inst', 'extdata', 'date_of_last_update.txt')); \
roxygen2::roxygenize()"

# for publishing the vignette to gh-pages
# Rscript --vanilla -e "devtools::build_vignettes()"
# ghp-import inst/doc
# git push
