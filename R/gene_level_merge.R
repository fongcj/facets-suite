
##########################################################################################
##########################################################################################
# MSKCC CMO
# Merging gene-level FACETS calls from purity, hisens runs
# TODO: Extend to handle 3+ Cvals
##########################################################################################
##########################################################################################


#' @name merge_calls
#' @title don't know
#' @description
#'
#' don't know
#'
#' @param hisens numeric first argument
#' @param purity numeric first argument
#' @return don't know
merge_calls = function(hisens, purity) {
  
  hisens[, id := basename(gsub("_hisens.cncf.txt", "", Tumor_Sample_Barcode))]
  purity[, id := basename(gsub("_purity.cncf.txt", "", Tumor_Sample_Barcode))]
  
  hisens_fails <- setdiff(unique(purity$id), unique(hisens$id))
  purity <- purity[id %in% hisens_fails]
  merged <- rbind(hisens, purity)
  merged[, id := NULL]
  
  return(merged)
  
}

#' @name gene_level_merge_function
#' @title don't know
#' @description
#'
#' don't know
#'
#' @param hisens character FACETS hisens gene-level calls
#' @param purity character FACETS purity gene-level calls
#' @param outfile character Output filename
#' @return don't know
gene_level_merge_function = function(hisens, purity, outfile){

  ## if( ! interactive() ) {
  ##  
  ##  pkgs = c('data.table', 'argparse')
  ##  lapply(pkgs, require, character.only = T)
  ##  
  ##  parser=ArgumentParser()
  ##  parser$add_argument('-s', '--hisens', type='character', help='FACETS hisens gene-level calls')
  ##  parser$add_argument('-p', '--purity', type='character', help='FACETS purity gene-level calls')  
  ##  parser$add_argument('-o', '--outfile', type='character', help='Output filename.')
  ##  ## args=parser$parse_args()
  
  ## hisens <- fread(args$hisens)
  ## purity <- fread(args$purity)
  ## outfile <- args$outfile
  ## }
  
  merged_calls <- merge_calls(hisens, purity)
  write.table(merged_calls, outfile, quote = FALSE, col.names = TRUE, row.names = FALSE,sep = "\t")
   
}




