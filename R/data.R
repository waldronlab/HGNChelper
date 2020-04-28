#' All current and withdrawn MGI mouse symbols and Excel-mogrified symbols
#'
#' @description A \code{data.frame} with the first column providing a gene symbol or 
#' known alias (including withdrawn symbols), second column providing the approved 
#' MGI mouse gene symbol.
#'
#' \itemize{
#'   \item \code{Symbol}: All valid, Excel-mogrified, and withdrawn symbols
#'   \item \code{Approved.Symbol}: Approved symbols
#' }
#' 
#' @examples 
#' data("mouse.table", package="HGNChelper")
#' head(mouse.table)
#' 
#' @source Extracted from \url{http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt} and 
#' system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv", package="HGNChelper")
"mouse.table"

#' All current and withdrawn HGNC gene symbols and Excel-mogrified symbols
#'
#' @description A \code{data.frame} with the first column providing a gene symbol or 
#' known alias (including withdrawn symbols), second column providing the approved
#' HGNC human gene symbol.
#'
#' \itemize{
#'   \item \code{Symbol}: All valid, Excel-mogrified, and withdrawn symbols
#'   \item \code{Approved.Symbol}: Approved symbols
#' }
#' 
#' @examples 
#' data("hgnc.table", package="HGNChelper")
#' head(hgnc.table)
#' 
#' @source Extracted from \url{ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt}
#' and system.file("extdata/mog_map.csv", package="HGNChelper")
"hgnc.table"

