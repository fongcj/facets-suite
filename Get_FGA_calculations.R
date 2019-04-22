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

print('Loading facets locations')
fname_out <- 'facets_fga_all_20190417_with_multi_feature.csv'
path_facets_loc <- '/Users/fongc2/Desktop/luna_transfer'
# path_facets_loc <- '/ifs/res/taylorlab/impact_facets/facets_0.5.14/all/P-00342/P-0034223-T01-IM6_P-0034223-N01-IM6/default/'
path_facets_loc_all <- 'all'
path_facets_loc_manifests <- 'manifests'
fname_loc = 'impact_facets_manifest_2019_04_21.txt'
filenamePath <- file.path(path_facets_loc, path_facets_loc_manifests, fname_loc)

df_facets_loc <- read.csv(filenamePath, header = TRUE, sep = '\t')
head(df_facets_loc)
path_facets_loc <- '/Users/fongc2/Desktop/luna_transfer'

# p <- '/Users/fongc2/Desktop/luna_transfer/all/P-00342/P-0034223-T01-IM6_P-0034223-N01-IM6/default/P-0034223-T01-IM6_P-0034223-N01-IM6_purity.Rdata'
# load(p)
cols_bn_chr <- c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11',  'Chr12', 
                 'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21', 'Chr22', 'Chr23')
cols_bn_fga <- c('FGA', 'GAIN', 'LOSS', 'LOH', 'IS_WGD', 'EM_FLAGS')

df_facets_loc$genome_doubled <- NaN
df_facets_loc$fraction_cna <- NaN
df_facets_loc$genome_doubled_bn_algo <- NaN
df_facets_loc$fraction_cna_bn_algo <- NaN
df_facets_loc$loglik <- NaN
df_facets_loc$purity <- NaN
df_facets_loc$ploidy <- NaN
df_facets_loc$dipLogR <- NaN
df_facets_loc$timestamp <- NaN
df_facets_loc$TCN_mean <- NaN
df_facets_loc$TCN_stdev <- NaN

df_facets_loc[cols_bn_fga] <- NA
df_facets_loc[cols_bn_chr] <- NA

# Compute Fraction CNA
print('Calculating Fraction of Genome Altered')
for (i in 1:nrow(df_facets_loc)) {
# for (i in 1:2) {
  # Build Pathname
  folder_sample <- as.character(df_facets_loc$FOLDER_NAME_TOP[i])
  folder_facets <- as.character(df_facets_loc$FOLDER_NAME_FACETS[i])
  filename_facets <- as.character(df_facets_loc$FILENAME[i])
  pathfilename <- file.path(path_facets_loc, folder_sample, folder_facets, filename_facets)
  # Load R data
  if (file.exists(pathfilename)) {
    # Load Rdata
    load(pathfilename)
    
    # Check if columns exist and add other facet info
    if("loglik" %in% names(fit)) {
      loglik <- fit$loglik
      df_facets_loc[i, 'loglik'] <- loglik
    }
    if("purity" %in% names(fit)) {
      purity <- fit$purity
      df_facets_loc[i, 'purity'] <- purity
    }
    if("ploidy" %in% names(fit)) {
      ploidy <- fit$ploidy
      df_facets_loc[i, 'ploidy'] <- ploidy
    }
    if("dipLogR" %in% names(fit)) {
      dipLogR <- fit$dipLogR
      df_facets_loc[i, 'dipLogR'] <- dipLogR
    }
    
    # Get time stamp of the .Rdata file
    time_stamp <- file.info(pathfilename)$ctime
    df_facets_loc[i, 'timestamp'] <- time_stamp
    
    # Compute CNA via PJ/CG's FGA script
    facets_cna <- calculate_fraction_cna(segs=fit$cncf, ploidy=fit$ploidy, "hg19", "em")
    # Place into df
    gd <- facets_cna$genome_doubled
    fga <- facets_cna$fraction_cna
    df_facets_loc[i, 'genome_doubled'] <- gd
    df_facets_loc[i, 'fraction_cna'] <- fga
    
    # Compute CNA via PJ/CG's FGA script
    fga_output <- get_FGA_BN(fit = fit, out = out)
    chr_output <- get_FChrA(fit = fit, out = out)

    # Place into df
    df_facets_loc[i, cols_bn_fga] <- fga_output[ , cols_bn_fga]
    df_facets_loc[i, cols_bn_chr] <- chr_output[ , cols_bn_chr]
    
    # Compute mean and stdev of TCN data
    tcn_mu_sig <- facet_seg_stdev(cncf=fit)
    # Place into df
    df_facets_loc[i, 'TCN_mean'] <- tcn_mu_sig$tcn_mean
    df_facets_loc[i, 'TCN_stdev'] <- tcn_mu_sig$tcn_stdev
    
  }
}
colnames(df_facets_loc)[colnames(df_facets_loc)=="FGA"] <- "CNA.FACETS.BN"

write.csv(df_facets_loc, file = fname_out)

