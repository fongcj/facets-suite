#' Format .seg file
#'
#' Format Facets output for viewing in IGV.
#'
#' @param out Facets out objet.
#' @param sample_id Sample name.
#' 
#' @import dplyr
#'
#' @return Segmentation output formatted for IGV.

format_seg = function(out, sample_id) {
    group_by(out$jointseg, chrom, seg) %>% 
        summarize(loc.start = min(maploc),
                  loc.end = max(maploc)) %>% 
        left_join(., select(out$out, chrom, seg, num.mark, seg.mean = cnlr.median)) %>% 
        mutate(ID = sample_id) %>% 
        select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
}
