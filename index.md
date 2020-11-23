[![Travis-CI Build Status](https://travis-ci.org/waldronlab/HGNChelper.svg?branch=master)](https://travis-ci.org/waldronlab/HGNChelper)
[![Coverage Status](https://codecov.io/github/waldronlab/HGNChelper/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/HGNChelper?branch=master)
[![](https://cranlogs.r-pkg.org/badges/HGNChelper)](https://cran.r-project.org/package=HGNChelper)

# HGNChelper: Identification and correction of invalid gene symbols for human and mouse

A [pre-print is available on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.09.16.300632v2). 

## Why HGNChelper?

Physicians and biologists like gene symbols and bioinformaticians
hate'em. Why? For one thing, they change constantly and are given new
names or aliases. For another, some get munged into dates when
imported into spreadsheet programs - and not only Excel (Thank you
[@karawoo](https://twitter.com/kara_woo) for the
[picture](https://twitter.com/kara_woo/status/1020054225022173184)!):

![](articles/0DaysSince.png)

Myself (Levi speaking), I don't mind them. It's way easier to remember
[TP53](http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=11998)
than to remember
[7157](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=full_report&list_uids=7157)
or
[ENSG00000141510](http://www.ensembl.org/Homo_sapiens/geneview?gene=ENSG00000141510). They're
a fact of life. So Sehyun Oh, Markus Riester, and I wrote HGNChelper to make them
a little more pleasant to bioinformaticians.

* See the [Articles](articles/index.html) link above for a vignette
  describing package use
* See the [Reference](reference/index.html) for function documentation.
* See the [GitHub](https://github.com/waldronlab/HGNChelper) page for
  creating your own package update with the most current gene lists,
  rather than updating within R as described in the vignette
