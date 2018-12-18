#' Example SNP counts
#'
#' Output generated with snp-pileup for testing, processed with \code{read_snp_matrix}.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{Chromosome}
#'   \item{Position}
#'   \item{NOR.DP}{Depth in normal sample}
#'   \item{TUM.DP}{Depth in tumor sample}
#'   \item{NOR.RD}{Reads supporting reference allele in normal sample}
#'   \item{TUM.RD}{Reads supporting reference allele in tumor sample}
#' }
#' @source \url{portal.gdc.cancer.gov}
'example_read_counts'