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

print('Loading facets locations')
pathname = "data"
fname = 'facets_file_locations_20190309.csv'
fname_out <- "facets_fga_all_20190311.csv"
filenamePath <- file.path(pathname, fname)
df_facets_loc <- read.csv(filenamePath, header = TRUE, sep = ',')
head(df_facets_loc)
# path_facets_loc <- '/Users/fongc2/Desktop/luna_transfer/all'
path_facets_loc <- '/ifs/res/taylorlab/impact_facets/all'

# load('/Users/fongc2/Desktop/luna_transfer/all/P-0000004-T01-IM3_P-0000004-N01-IM3/facets_R0.5.6c50p150m15p15/P-0000004-T01-IM3_P-0000004-N01-IM3_purity.Rdata')
df_facets_loc$genome_doubled <- NaN
df_facets_loc$fraction_cna <- NaN
df_facets_loc$loglik <- NaN
df_facets_loc$purity <- NaN
df_facets_loc$ploidy <- NaN
df_facets_loc$dipLogR <- NaN

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
    load(pathfilename)
    loglik <- fit$loglik
    purity <- fit$purity
    ploidy <- fit$ploidy
    dipLogR <- fit$dipLogR
    
    # Compute CNA
    facets_cna <- calculate_fraction_cna(segs=fit$cncf, ploidy=fit$ploidy, "hg19", "em")
    
    # Place into df
    gd <- facets_cna$genome_doubled
    fga <- facets_cna$fraction_cna
    df_facets_loc[i, 'genome_doubled'] <- gd
    df_facets_loc[i, 'fraction_cna'] <- fga 
    
    # Other facet info
    df_facets_loc[i, 'loglik'] <- loglik
    df_facets_loc[i, 'purity'] <- purity
    df_facets_loc[i, 'ploidy'] <- ploidy
    df_facets_loc[i, 'dipLogR'] <- dipLogR

  }
}

write.csv(df_facets_loc, file = fname_out)

