###########################################################################
# segment_stdev
#
#
# By Chris Fong and Francisco Sanchez-Vega, MSKCC 2019
#
###########################################################################

# Compute the standard deviation of the facets segmentation 
# Input -- For a given Rdata file from facets, cncf = fit
facet_seg_stdev = function(cncf) {
  cncf <- fit$cncf

  # Compute the mean (mu) to be the sum of the segment length * total copy number, all divided by the total length of segments
  seg_len <- cncf$end - cncf$start
  G <- sum(seg_len)
  G_inv <- 1/G
  mu <- sum(cncf$tcn.em * seg_len) * G_inv
  
  # Compute the standard deviation to be the sum of squares of total copy number - mean, times the segment length, all divided by the total length
  stdev <- sqrt(G_inv * sum((cncf$tcn.em - mu)^2 * seg_len))

  list(
    tcn_mean = mu,
    tcn_stdev = stdev
  )
}


