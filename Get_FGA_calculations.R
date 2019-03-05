###########################################################################
# Get_FGA_calculations
#
#
# By Chris Fong, MSKCC 2018
#
###########################################################################
# if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
# if("facetsSuite" %in% rownames(installed.packages()) == FALSE) {install.packages("facetsSuite")}
# if("usethis" %in% rownames(installed.packages()) == FALSE) {install.packages("usethis")}
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
fname = 'facets_file_locations_AACR_2019.csv'
filenamePath <- file.path(pathname, fname)
df_facets_loc <- read.csv(filenamePath, header = TRUE, sep = '\t')

path_facets_loc <- '/ifs/res/taylorlab/impact_facets/'

# Compute Fraction CNA
print('Calculating Fraction of Genome Altered')
for (i in 1:nrow(df_facets_loc)) {
  # Build Pathname
  folder_sample <- df$FOLDER_NAME_TOP[i]
  folder_facets <- df$FOLDER_NAME_FACETS[i]
  filename_facets <- df$FILENAME[i]
  pathfilename <- file.path(path_facets_loc, folder_sample, folder_facets, filename_facets)
  # Load R data
  load(pathfilename)
  
  # Compute CNA
  facets_cna <- calculate_fraction_cna(segs=fit$cncf, ploidy=fit$ploidy, "hg19", "em")
  
  # Place into df
  gd <- fraction_cna$genome_doubled
  fga <- fraction_cna$fraction_cna
  df_facets_loc$genome_doubled[i] <- gd
  df_facets_loc$fraction_cna[i] <- fga
}


write.csv(df_facets_loc, file = "facets_fga_AACR_2019.csv")

