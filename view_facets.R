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
pathname = "/Users/fongc2/Documents/github/MSK/facets-suite/data"
fname = "facets_fga_all_20190314_with_multi_feature.csv"
filenamePath <- file.path(pathname, fname)
df1 <- read.csv(filenamePath, header = TRUE)
df1$genome_doubled <- as.logical(df1$genome_doubled)
df1$genome_doubled_bn_algo <- as.logical(df1$genome_doubled_bn_algo)
head(df1)

# Filter out cases unknown to be hisens or purity
df <- df1[df1$FACETS_RUN_TYPE != 'Other', ]

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Correlation of facets FGA and Portal
p1 <- ggplot(df, aes(x=CNA.Portal, y=CNA.FACETS)) +
  geom_point(aes(color=PURITY_LEVEL), alpha=0.6) +
  facet_grid(vars(FACETS_RUN_TYPE), vars(genome_doubled), labeller=label_both)

p1

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Correlation of facets FGA and Portal
p_b <- ggplot(df, aes(x=CNA.Portal, y=fraction_cna_bn_algo)) +
  geom_point(aes(color=PURITY_LEVEL), alpha=0.6) +
  facet_grid(vars(FACETS_RUN_TYPE), vars(genome_doubled_bn_algo), labeller=label_both)

p_b

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Correlation of facets FGA and Bastien's
p_b2 <- ggplot(df[df$FACETS_RUN_TYPE == 'hisens', ], aes(x=fraction_cna_bn_algo, y=CNA.FACETS)) +
  geom_point(aes(color=PURITY_LEVEL), alpha=0.6) +
  facet_grid(vars(genome_doubled_bn_algo), vars(genome_doubled), labeller=label_both) +
  theme(axis.line = element_line(size = 3, colour = "grey80"))

p_b2

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Box plot of mean and std 
p_stdev <- ggplot(df[df$FACETS_RUN_TYPE == 'hisens', ], aes(x=PURITY_LEVEL, y=TCN_stdev)) +
  geom_boxplot() +
  ylim(0, 5)


p_stdev


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load summary data
pathname = "/Users/fongc2/Documents/github/MSK/facets-suite/data"
fname = "facets_p_vs_h_summary.csv"
filenamePath <- file.path(pathname, fname)
df2 <- read.csv(filenamePath, header = TRUE)
df2$genome_doubled <- as.logical(df1$genome_doubled)
head(df2)

# Correlation of facets FGA and Portal
p2 <- ggplot(df2, aes(x=CNA.FACETS.hisens, y=CNA.FACETS.purity)) +
  geom_point(aes(color=FOLDER_NAME_FACETS), alpha=0.6) +
  geom_abline(alpha=0.4)

p2

p3 <- ggplot(df2, aes(x=purity.hisens, y=purity.purity)) +
  geom_point(aes(color=FOLDER_NAME_FACETS), alpha=0.6) +
  geom_abline(alpha=0.4)

p3

p4 <- ggplot(df2, aes(x=ploidy.hisens, y=ploidy.purity)) +
  geom_point(aes(color=FOLDER_NAME_FACETS), alpha=0.6) +
  geom_abline(alpha=0.4) +
  xlim(0,15) +
  ylim(0,15) 

p4
