#' All current and withdrawn MGI mouse symbols and Excel mogrifications.
#'
#' A dataframe with the first column providing a gene symbol or known alias
#' (including withdrawn symbols), second column providing the approved MOUSE symbol.
#'
#' \itemize{
#'   \item Symbol. All valid, Excel-mogrified, and withdrawn symbols
#'   \item Approved.Symbol. Approved symbols.
#' }
#' @examples 
#' data("mouse.table", package="HGNChelper")
#' head(mouse.table)
#' @source Extracted from \url{http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt} and 
#' system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv", package="HGNChelper")
"mouse.table"

#' All current and withdrawn HGNC gene symbols and Excel mogrifications.
#'
#' A dataframe with the first column providing a gene symbol or known alias
#' (including withdrawn symbols), second column providing the approved HGNC 
#' human gene symbol.
#'
#' \itemize{
#'   \item Symbol. All valid, Excel-mogrified, and withdrawn symbols
#'   \item Approved.Symbol. Approved symbols.
#' }
#' @examples 
#' data("hgnc.table", package="HGNChelper")
#' head(hgnc.table)
#' @source Extracted from 
#' \url{http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&hgnc_dbtag=on&submit=submit}
#'  and 
#' system.file("extdata/mog_map.csv", package="HGNChelper")
"hgnc.table"

