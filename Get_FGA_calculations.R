###########################################################################
# Get_FGA_calculations
#
#
# By Chris Fong, MSKCC 2018
#
###########################################################################
# if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
# if("facetsSuite" %in% rownames(installed.packages()) == FALSE) {install.packages("facetsSuite")}
# if("usethis" %in% rownames(installed.packages()) == FALSE) {install.packages("usethis", url='http://ftp.osuosl.org/pub/cran/')}
# if("readr" %in% rownames(installed.packages()) == FALSE) {install.packages("readr")}

library("data.table")
library("dplyr")
library("methods")
library("readr")
library("usethis")
# library("facetsSuite")

# Source facets functions
print('Loading facets suite')
source("R/copy-number-scores.R")
source("data-raw/sysdata.R")
source("FGA.R")
source("segment_stdev.R")

printf <- function(...)print(sprintf(...))

print('Loading facets locations')
fname_out_prefix <- 'facets_0.5.14_fga_all_20190422_with_multi_feature'
extension <- '.csv'
path_facets_loc <- '/Users/fongc2/Desktop/luna_transfer/facets_0.5.14/'
path_facets_orig_path <- '/ifs/res/taylorlab/impact_facets/facets_0.5.14/'
# path_facets_loc <- '/ifs/res/taylorlab/impact_facets/facets_0.5.14/all/P-00342/P-0034223-T01-IM6_P-0034223-N01-IM6/default/'
path_facets_loc_all <- 'all'
path_facets_loc_manifests <- 'manifests'
fname_loc = 'impact_facets_manifest_2019_04_21.txt'

# Load manifest file
filenamePath <- file.path(path_facets_loc, path_facets_loc_manifests, fname_loc)
df_facets_loc <- read.csv(filenamePath, header = TRUE, sep = '\t')
head(df_facets_loc)

cols_facets <- c('patient', 'tumor_sample', 'run_output_dir', 'tag', 'has_hisens_run', 'has_purity_run')
df_facets_output <- df_facets_loc[ , cols_facets]

# p <- '/Users/fongc2/Desktop/luna_transfer/all/P-00342/P-0034223-T01-IM6_P-0034223-N01-IM6/default/P-0034223-T01-IM6_P-0034223-N01-IM6_purity.Rdata'
# load(p)
cols_cb <- c('genome_doubled', 'fraction_cna', 'genome_doubled_bn_algo', 'fraction_cna_bn_algo', 'loglik', 'purity', 'ploidy', 'dipLogR', 'timestamp')
cols_bn_chr <- c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11',  'Chr12', 
                 'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21', 'Chr22', 'Chr23')
cols_bn_fga <- c('FGA', 'GAIN', 'LOSS', 'LOH', 'IS_WGD', 'EM_FLAGS')
cols_output <- c(cols_cb, cols_bn_fga, cols_bn_chr)

df_facets_output[cols_output] <- NA

# Compute Fraction CNA
print('Calculating Fraction of Genome Altered')
for (k in 1:1) {
  if (k == 1) {
    logic_run_type <- df_facets_output$has_hisens_run == 'TRUE'
    fname_suffix <- '_hisens'
    rdata_suffix <- '_hisens.Rdata'
  } else {
    logic_run_type <- df_facets_output$has_purity_run == 'TRUE'
    fname_suffix <- '_purity'
    rdata_suffix <- '_purity.Rdata'
  }
  df_facets_output_current <- df_facets_output[logic_run_type, ]
  fname_out <- paste(fname_out_prefix, fname_suffix, extension, sep = '')
  
  # for (i in 1:2) {
  for (i in 20029:nrow(df_facets_output_current)) {
    # Build Pathname
    folder_facets <- as.character(df_facets_output_current$run_output_dir[i])
    folder_facets_append <- gsub(path_facets_orig_path, path_facets_loc, folder_facets)
    # Build filename
    filename_facets_pre <- as.character(df_facets_output_current$tag[i])
    filename_facets <- paste(filename_facets_pre, '_hisens.Rdata', sep = '')
      
    pathfilename <- file.path(folder_facets_append, filename_facets)
    # Load R data
    if (file.exists(pathfilename)) {
      # Load Rdata
      printf("Loading: %s", filename_facets)
      load(pathfilename)
      
      # Check if columns exist and add other facet info
      if("loglik" %in% names(fit)) {
        loglik <- fit$loglik
        df_facets_output_current[i, 'loglik'] <- loglik
      }
      if("purity" %in% names(fit)) {
        purity <- fit$purity
        df_facets_output_current[i, 'purity'] <- purity
      }
      if("ploidy" %in% names(fit)) {
        ploidy <- fit$ploidy
        df_facets_output_current[i, 'ploidy'] <- ploidy
      }
      if("dipLogR" %in% names(fit)) {
        dipLogR <- fit$dipLogR
        df_facets_output_current[i, 'dipLogR'] <- dipLogR
      }
      
      # Get time stamp of the .Rdata file
      time_stamp <- file.info(pathfilename)$ctime
      df_facets_output_current[i, 'timestamp'] <- time_stamp
      
      # Compute CNA via PJ/CG's FGA script
      facets_cna <- calculate_fraction_cna(segs=fit$cncf, ploidy=fit$ploidy, "hg19", "em")
      # Place into df
      gd <- facets_cna$genome_doubled
      fga <- facets_cna$fraction_cna
      df_facets_output_current[i, 'genome_doubled'] <- gd
      df_facets_output_current[i, 'fraction_cna'] <- fga
      
      # Compute CNA via PJ/CG's FGA script
      fga_output <- get_FGA_BN(fit = fit, out = out)
      chr_output <- get_FChrA(fit = fit, out = out)
  
      # Place into df
      cols_bn_fga_current <- intersect(colnames(fga_output), colnames(df_facets_output_current))
      cols_bn_chr_current <- intersect(colnames(chr_output), colnames(df_facets_output_current))
      df_facets_output_current[i, cols_bn_fga_current] <- fga_output[ , cols_bn_fga_current]
      df_facets_output_current[i, cols_bn_chr_current] <- chr_output[ , cols_bn_chr_current]
      
    } else {
      print('File does not exist')
    }
  }
  # colnames(df_facets_output)[colnames(df_facets_output)=="FGA"] <- "CNA.FACETS.BN"
  write.csv(df_facets_output_current, file = fname_out, row.names = FALSE)
}



