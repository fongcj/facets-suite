#' Run FACETS
#'
#' Runs FACETS on an input count file, with specified parameter settings.
#'
#' @param counts Counts object.
#'
#' @return A list object with \code{out} and \code{fit} objects from Facets run.
#'
#' @examples
#' facets_closeup('path_to_sample.Rdata', 'BRCA1', 17)
#' 
#' @import facets
#' 
#' @name run_facets
NULL

#' @export run_facets
run_facets = function(read_counts,
                      directory = getwd(),
                      cval = 100,
                      diplogr = NULL,
                      ndepth = 35,
                      snp_nbhd = 250,
                      min_nhet = 15,
                      genome = c('hg19', 'hg38', 'hg18', 'mm9', 'mm10'),
                      seed = 100,
                      rlib_path = NULL
) {
    set.seed(seed)
    genome = match.arg(genome)
    
    # run FACETS
    dat = facets::preProcSample(read_counts, ndepth = ndepth, het.thresh = 0.25, snp.nbhd = snp_nbhd, cval = cval,
                                gbuild = genome, hetscale = TRUE, unmatched = FALSE, ndepthmax = 1000)
    out = facets::procSample(dat, cval = cval, min.nhet = min_nhet, dipLogR = diplogr)
    fit = facets::emcncf(out)
    
    # generate output
    fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
    
    list(
        snps = out$jointseg,
        segs = fit$cncf,
        purity = fit$purity,
        ploidy = fit$ploidy,
        diplogr = out$dipLogR,
        alballogr = out$alBalLogR,
        flags = out$flags,
        em_flags = fit$emflags,
        loglik = fit$loglik,
        mafr_thresh = out$mafR.thresh
    )
}

