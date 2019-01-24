#!/usr/bin/Rscript
suppressPackageStartupMessages({
    library(argparse)
    })

parser = ArgumentParser(description = 'Wrapper to execute Facets and generate various output from the output from snp-pileup')

parser$add_argument('-v', '--verbose', action="store_true", default = TRUE,
                    help = 'Print run info')
parser$add_argument('-f', '--counts-file', required = T,
                    help = 'Merged, gzipped tumor-normal output from snp-pileup')
parser$add_argument('-s', '--sample-id', required = F,
                    help = 'Sample ID, preferrable Tumor_Normal to keep track of the normal used')
parser$add_argument('-D', '--directory', required = F,
                    default = getwd(), help = 'Output directory [default %(default)s]')
parser$add_argument('-c', '--cval', required = F,
                    default = 100, help = 'Segmentation parameter (cval) [default %(default)s]')
parser$add_argument('-pc', '--purity-cval', required = F,
                    default = NULL, help = 'If two-pass, purity segmentation parameter (cval)')
parser$add_argument('-d', '--diplogr', required = F,
                    default = NULL, help = 'Manual diplogr')
parser$add_argument('-g', '--genome', required = F,
                    choices = c('hg19', 'hg38', 'hg18', 'mm9', 'mm10'),
                    default = 'hg19', help = 'Reference genome [default %(default)s]')
parser$add_argument('-s', '--seed', required = F,
                    default = NULL, help = 'Manual seed value')


# 
# parser$add_argument('-q', '--quietly', action='store_false', 
#                     dest='verbose', help='Print little output')
# parser$add_argument('-c', '--count', type='integer', default=5, 
#                     help='Number of random normals to generate [default %(default)s]',
#                     metavar='number')
# parser$add_argument('--generator', default='rnorm', 
#                     help = 'Function to generate random deviates [default \'%(default)s\']')
# parser$add_argument('--mean', default=0, type='double',
#                     help='Mean if generator == \'rnorm\' [default %(default)s]')
# parser$add_argument('--sd', default=1, type='double',
#                     metavar='standard deviation',
#                     help='Standard deviation if generator == \'rnorm\' [default %(default)s]')


args = parser$parse_args()
if (is.null(args$sample_id)) args$sample_id = basename(args$counts_file)

# Determine if running two-pass
two_pass = FALSE
if (is.null(args$purity_cval)) {
    two_pass = TRUE
    purity_output = run_facets(read_counts = args$counts_file,
                               sample_id = args$sample_id,
                               directory = args$directory,
                               cval = args$purity_cval,
                               diplogr = args$diplogr,
                               ndepth = args$ndepth,
                               snp_nbhd = args$snp_nbhd,
                               min_nhet = args$purity_min_nhet,
                               genome = args$genome,
                               seed = args$seed,
                               rlib_path = NULL)
    diplogr = purity_output$diplogr
} else {
    diplogr = args$diplogr
}

# Run Facets
facets_output = run_facets(read_counts = args$counts_file,
                           sample_id = args$sample_id,
                           directory = args$directory,
                           cval = args$cval,
                           diplogr = diplogr,
                           ndepth = args$ndepth,
                           snp_nbhd = args$snp_nbhd,
                           min_nhet = args$min_nhet,
                           genome = args$genome,
                           seed = args$seed,
                           rlib_path = NULL)

# Generate output files
plot_facets()

# Write output

# Exit


