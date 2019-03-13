###########################################################################
# view_facets.R
#
#
# By Chris Fong, MSKCC 2018
#
###########################################################################
if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyr")}
if("gplots" %in% rownames(installed.packages()) == FALSE) {install.packages("gplots")}

library("ggplot2")


# Load data raw
pathname = "/Users/fongc2/Documents/github/MSK/metastatic_tropisms"
fname = "facets_fga_all_20190311_with_fga.csv"
filenamePath <- file.path(pathname, fname)
df1 <- read.csv(filenamePath, header = TRUE)
df1$genome_doubled <- as.logical(df1$genome_doubled)
head(df1)

# Filter out cases unknown to be hisens or purity
df <- df1[df1$FACETS_RUN_TYPE != 'Other', ]

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TMB  
p1 <- ggplot(df, aes(x=CNA.Portal, y=CNA.FACETS)) +
  geom_point(aes(color=PURITY_LEVEL), alpha=0.6) +
  facet_grid(vars(FACETS_RUN_TYPE), vars(genome_doubled), labeller=label_both)

p1
