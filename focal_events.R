#!/opt/common/CentOS_6-dev/R/R-3.4.1/bin/Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
    library(gridExtra)
    library(argparse)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(tidyr)
})


# Define functions ------------------------------------------------------------------------------------------------
load_genes = function(ivs, genes = NULL) {
    first_line = readLines(ivs, n = 1) %>% 
        strsplit(., '\t') %>% 
        unlist
    
    strand_ = which(first_line %in% c('+', '-'))[1]
    info_ = which(first_line %like% 'NM|ENS|[A-Z]+')
    
    if (strand_ == 4 & info_ == 5) {
        cols = c('chrom', 'start', 'end', 'strand', 'info')
    } else if (strand_ == 6 & info_ == 4) {
        cols = c('chrom', 'start', 'end', 'info', 'empty', 'strand')
    }
    
    ivs = fread(ivs, header = F, col.names = cols) %>% 
        mutate(gene = str_extract(info, '^[^\\:]+'),
               gene = ifelse(gene %like% 'CDKN2A|CDKN2B', str_extract(gene, 'CDKN2A|CDKN2B'), gene)) %>% 
        group_by(gene) %>% 
        filter(chrom %in% c(1:22, 'X')) %>% 
        mutate(exon = ifelse(strand == '+',
                             row_number(),
                             max(row_number()) - (row_number() - 1)),
               chrom = ifelse(chrom == 'X', as.numeric(23), as.numeric(chrom))) %>% 
        ungroup() 
    if (is.null(genes)) {
        ivs
    } else {
        filter(ivs, gene %in% genes)
    }
}

chrom_info = function() { # hg19
    data.frame(
        chr = 1:24,
        size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                 59128983, 63025520, 48129895, 51304566, 155270560, 59373566),
        centromere = c(123035434, 93826171, 92004854, 51160117, 47905641, 60330166, 59554331, 45338887, 48867679,
                       40754935, 53144205, 36356694, 17500000, 17500000, 18500000, 36835801, 23763006, 16960898,
                       26181782, 27869569, 12788129, 14500000, 60132012, 11604553)
    )
}

maploc_cum = function(mat) {
    
    chroms = chrom_info()
    chroms = chroms[chroms$chr != 24,]
    cum.chrom.lengths = cumsum(as.numeric(chroms$size))
    mid = cum.chrom.lengths - chroms$size/2
    names(mid) = c(1:22, 'X')
    centromeres = chroms$centromere + c(0, cum.chrom.lengths[-length(cum.chrom.lengths)])
    
    chr.maploc.to.gen.maploc = function(x) {
        mat[mat$chrom == x, ]$maploc + cum.chrom.lengths[x - 1]
    }
    chr.maploc = sapply(2:23, chr.maploc.to.gen.maploc)
    chr.maploc = unlist(chr.maploc)
    chr.maploc = c(mat[mat$chrom == 1,]$maploc, chr.maploc)
    mat = cbind(mat, chr.maploc)
    list(mat = mat, mid = mid, centromeres = centromeres)
}

plot_events = function(snps, segs, gmat, diplogr) {
    
    snps = maploc_cum(snps)
    snps$mat$chrom[which(snps$mat$chrom == 23)] = 'X'
    mid = snps$mid[names(snps$mid) %in% snps$mat$chrom]
    centromeres = snps$centromeres
    snps = snps$mat
    
    cnlr.median = rep(segs$cnlr.median, segs$num.mark)
    mat = cbind(snps, cnlr.median)
    
    starts = cumsum(c(1, segs$num.mark))[1:length(segs$num.mark)]
    ends = cumsum(c(segs$num.mark))
    my.starts = mat[starts, c('chr.maploc', 'cnlr.median')]
    my.ends = mat[ends, c('chr.maploc', 'cnlr.median')]
    
    if (any(gmat$signif) == T) {
        pts = filter(gmat, signif == T) %>%
            separate_rows(pos, sep = ',') %>% 
            mutate(up = fold_change_gene > 0,
                   down = fold_change_gene < 0,
                   maploc = as.numeric(pos),
                   chrom = ifelse(chrom %like% '23', 'X', gsub('(p|q)', '', chrom)),
                   chrom = as.character(chrom)) %>% 
            select(chrom, maploc, up, down) 
        
        mat = left_join(mat, pts, by = c('chrom', 'maploc')) %>% 
            mutate(
                type = case_when(
                    up == T ~ 'up',
                    down == T ~ 'down',
                    chrom %in% seq(2, 22, 2) ~ 'even',
                    chrom %in% c(seq(1, 21, 2), 'X') ~ 'odd')
            )
    } else {
        mat = mutate(mat, type = ifelse(chrom %in% seq(2, 22, 2), 'even', 'odd'))
    }
    
    ymin = floor(min(range(segs$cnlr.median), na.rm = T))
    if (ymin > -3) ymin = -3
    ymax = max(ceiling(max(range(segs$cnlr.median), na.rm = T)), mat$cnlr[mat$type == 'up'])
    if (ymax < 3) ymax = 3
    
    cn0_snps = snps[cn0 == T]
    # lower = mean(cn0_snps$cnlr)-log(2, base = 2)
    # upper = mean(cn0_snps$cnlr)+log(3, base = 2)
    
    # adjust for diplogr
    mat$cnlr = mat$cnlr - diplogr
    my.starts$cnlr.median = my.starts$cnlr.median - diplogr
    my.ends$cnlr.median = my.ends$cnlr.median - diplogr
    
    ggplot(mat, aes(x = chr.maploc, y = cnlr, col = type)) +
        geom_vline(xintercept = centromeres, linetype = 'dashed', color = 'darkgrey', size = .25) +
        geom_point(data = filter(mat, type %in% c('odd', 'even')), pch = 19, size = .4) +
        geom_point(data = filter(mat, type %in% c('up', 'down')), pch = 19, size = .6) +
        scale_color_manual(values = c('even' = '#bcbcbc', 'odd' = '#dbdbdb', 'up' = '#d7301f', 'down' = '#d7301f'),
                           guide = F) +
        scale_x_continuous(breaks = mid, labels = names(mid), expand = c(0.01, 0)) +
        # geom_vline(xintercept = as.numeric(mid[seq(1, 23)])) +
        scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(ymin, ymax)) +
        geom_hline(yintercept = c(-log2(2), log2(3)), linetype = c('dashed', 'dashed'),
                   color = c('grey35', 'grey35'), size = .5) +
        geom_segment(data = segs, col = 'dodgerblue2', size = .5,
                     aes(
                         x = my.starts$chr.maploc,
                         xend = my.ends$chr.maploc,
                         y = my.starts$cnlr.median,
                         yend = my.ends$cnlr.median
                     )) +
        theme_bw() +
        theme(text = element_text(size = 10),
              axis.text.x = element_text(angle = 0, size = 8, color = 'black'),
              axis.text.y = element_text(angle = 0, size = 8, color = 'black'),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              #panel.grid.minor.x = element_line(colour = 'grey', size = .2),
              #panel.grid.major.x = element_line(colour = 'grey', size = 0),
              plot.margin = unit(c(0, 1, 0, 0), 'lines')) +
        labs(x = NULL, y = 'Copy number log ratio')
}

plot_table = function(gmat, n) {
    up = filter(gmat, fold_change_gene > 0) %>% 
        arrange(desc(signif), desc(fold_change_gene)) %>% 
        head(n) %>% 
        select(gene, chrom, fc = fold_change_gene, p = pval_gene, p_adj = padj_gene, signif) %>% 
        mutate(signif = ifelse(signif == T, 'yes', ''))
    down = filter(gmat, fold_change_gene < 0) %>% 
        arrange(desc(signif), fold_change_gene) %>% 
        head(n) %>% 
        select(gene, chrom, fc = fold_change_gene, p = pval_gene, p_adj = padj_gene, signif) %>% 
        mutate(signif = ifelse(signif == T, 'yes', ''))
    
    table_theme = ttheme_minimal(base_size = 6)
    
    up = tableGrob(up, theme = table_theme, rows = NULL)
    down = tableGrob(down, theme = table_theme, rows = NULL)
    
    list(up = up, down = down)
}

find_events = function(rdata,
                       gene_intervals = NULL,
                       gene_list = NULL)
{
    
    # load gene annotation
    if (is.null(gene_intervals)) {
        gene_intervals = '/ifs/depot/resources/dmp/data_20150226/mskdata/interval-lists/VERSIONS/cv6/genelist.with_aa.interval_list'
    }
    genes = load_genes(gene_intervals, gene_list)
    genes = as.data.table(genes)
    setkey(genes, 'chrom', 'start', 'end')
    
    # load Facets output
    load(rdata)
    diplogr = out$dipLogR
    ploidy = fit$ploidy
    segs = as.data.table(fit$cncf)
    snps = as.data.table(out$jointseg)
    if (nrow(snps) < 1000) stop('Not enough SNPs, is something wrong with FACETS output?')
    snps$maploc2 = snps$maploc
    snps[, chrom := as.numeric(chrom)]
    setcolorder(snps, c('chrom', 'maploc', 'maploc2'))
    setkey(snps, 'chrom', 'maploc', 'maploc2')
    snps = merge(snps,
                 segs[, .(seg, tcn.em)],
                 by = 'seg')
    
    # set baseline for copy-number unaltered
    if (!any(segs$cnlr.median.clust == diplogr)) {
        new_dlr = unique(segs$cnlr.median.clust)[which.min(abs(unique(segs$cnlr.median.clust)-diplogr))]
        diplogr = new_dlr
    }
    cn0_segs = segs[cnlr.median.clust == diplogr]
    cn0_snps_all = snps[segclust %in% cn0_segs$segclust | tcn.em == unique(cn0_segs$tcn.em)] # baseline can be any segment on diplogr or segment with same tcn // this might need to be adjusted for GD tumors
    cn0_snps = cn0_snps_all[between(cnlr, quantile(cn0_snps_all$cnlr, .25), quantile(cn0_snps_all$cnlr, .75))] # remove noise
    snps[, cn0 := paste0(chrom, ':', maploc) %in% paste0(cn0_snps$chrom, ':', cn0_snps$maploc)]
    
    chroms = chrom_info()
    
    gmat = mclapply(unique(genes$gene), function(g) {

        g_coords = genes[gene == g]
        g_chrom = unique(g_coords$chrom)
        g_start = min(g_coords$start)
        g_end = max(g_coords$end)
        
        # g_snps = data.table::foverlaps(
        #     y = g_coords,
        #     x = snps,
        #     by.y = c('chrom', 'start', 'end'),
        #     by.x = c('chrom', 'maploc', 'maploc2'),
        #     type = 'within',
        #     mult = 'all',
        #     nomatch = 0L
        # )
        g_snps = snps[chrom == g_chrom & between(maploc, g_start, g_end)]
        
        if (nrow(g_snps) >= 5) { # require at least five SNPs
            
            # test whole gene
            if (mean(g_snps$cnlr) > mean(cn0_snps$cnlr)) {
                p_diplogr = pnorm(
                    mean(g_snps$cnlr),
                    mean = mean(cn0_snps$cnlr),
                    sd = sd(cn0_snps$cnlr),
                    lower.tail = F)
            } else {
                p_diplogr = pnorm(
                    mean(g_snps$cnlr),
                    mean(cn0_snps$cnlr),
                    sd(cn0_snps$cnlr),
                    lower.tail = T)
            }
            logr_diplogr = mean(g_snps$cnlr) - mean(cn0_snps$cnlr)
            fc_diplogr = ifelse(logr_diplogr < 0, -2^(-logr_diplogr), 2^(logr_diplogr))
            
            # test chromosome arm against baseline
            arm = ifelse(max(g_snps$maploc) > chroms$centromere[chroms$chr==unique(g_snps$chrom)], 'q', 'p')
            chr = paste0(unique(g_coords$chrom), arm)
            arm_snps = snps[if (arm == 'p') {
                chrom %in% unique(g_snps$chrom) &
                    !maploc %in% g_snps$maploc &
                    maploc < chroms$centromere[chroms$chr == unique(g_snps$chrom)]
            } else {
                {
                    chrom %in% unique(g_snps$chrom) &
                        !maploc %in% g_snps$maploc &
                        maploc > chroms$centromere[chroms$chr == unique(g_snps$chrom)]
                }
            }]
            arm_snps = arm_snps[between(cnlr, quantile(arm_snps$cnlr, .1), quantile(arm_snps$cnlr, .9))] # remove noise
            
            if (mean(arm_snps$cnlr) > mean(cn0_snps$cnlr)) {
                p_arm = pnorm(
                    mean(arm_snps$cnlr),
                    mean(cn0_snps$cnlr),
                    sd(cn0_snps$cnlr),
                    lower.tail = F)
            } else {
                p_arm = pnorm(
                    mean(arm_snps$cnlr),
                    mean(cn0_snps$cnlr),
                    sd(cn0_snps$cnlr),
                    lower.tail = T)
            }
            logr_arm = mean(arm_snps$cnlr) - mean(cn0_snps$cnlr)
            fc_arm = ifelse(logr_arm < 0, -2^(logr_arm), 2^(logr_arm))
            
            # if arm is different from baseline, test gene against arm
            if (mean(g_snps$cnlr) > mean(arm_snps$cnlr)) {
                p_g_arm = pnorm(
                    mean(g_snps$cnlr),
                    mean(arm_snps$cnlr),
                    sd(arm_snps$cnlr),
                    lower.tail = F)
            } else {
                p_g_arm = pnorm(
                    mean(g_snps$cnlr),
                    mean(arm_snps$cnlr),
                    sd(arm_snps$cnlr),
                    lower.tail = T)
            }
            fc_g_arm = mean(g_snps$cnlr)/mean(arm_snps$cnlr)
                
            # also, test intragenic
            intragenic_snps = g_snps[between(cnlr, quantile(g_snps$cnlr, .1), quantile(g_snps$cnlr, .9))]
            if (nrow(intragenic_snps) >= 5) {
                if (mean(intragenic_snps$cnlr) > mean(cn0_snps$cnlr)) {
                    p_intragenic = 1 - pnorm(
                        mean(intragenic_snps$cnlr),
                        mean = mean(cn0_snps$cnlr),
                        sd = sd(cn0_snps$cnlr),
                        lower.tail = T)
                } else {
                    p_intragenic = pnorm(
                        mean(intragenic_snps$cnlr),
                        mean = mean(cn0_snps$cnlr),
                        sd = sd(cn0_snps$cnlr),
                        lower.tail = T)
                }
                
                logr_diplogr_intragenic = mean(intragenic_snps$cnlr) - mean(cn0_snps$cnlr)
                fc_diplogr_intragenic = ifelse(logr_diplogr_intragenic < 0,
                                               -2^(-logr_diplogr_intragenic),
                                               2^(logr_diplogr_intragenic))
        
                tibble(gene = g,
                       chrom = chr,
                       fold_change_gene = fc_diplogr,
                       pval_gene = p_diplogr,
                       fold_change_arm = fc_arm,
                       pval_arm = p_arm,
                       fold_change_gene_arm = fc_g_arm,
                       pval_gene_arm = p_g_arm,
                       fold_change_intragenic = fc_diplogr_intragenic,
                       pval_intragenic = p_intragenic,
                       pos = paste0(g_snps$maploc, collapse = ','),
                       pos_intragenic = paste0(intragenic_snps$maploc, collapse = ','))
            } else {
                tibble(gene = g,
                       chrom = chr,
                       fold_change_gene = fc_diplogr,
                       pval_gene = p_diplogr,
                       fold_change_arm = fc_arm,
                       pval_arm = p_arm,
                       fold_change_gene_arm = p_g_arm,
                       pval_gene_arm = fc_g_arm,
                       pos = paste0(g_snps$maploc, collapse = ','))
            }
        }
    }, mc.cores = detectCores())
    gmat = rbindlist(gmat, fill = T)
    
    # multiple hypothesis testing correction and mark significant changes
    gmat = mutate(gmat,
                  pval_intragenic = ifelse(is.na(pval_intragenic), 1, pval_intragenic),
                  pval = ifelse(pval_gene < pval_intragenic, pval_gene, pval_intragenic),
                  fold_change = ifelse(pval_gene < pval_intragenic, fold_change_gene, fold_change_intragenic))
    gmat$padj = p.adjust(gmat$pval, method = 'BH') 
    gmat$padj_gene = p.adjust(gmat$pval_gene, method = 'BH') 
    gmat$padj_intragenic = p.adjust(gmat$pval_intragenic, method = 'BH')
    gmat = mutate(gmat,
                  signif = (fold_change >= 6 | fold_change < -2) & padj < .05,
                  signif = ifelse(pval_arm < .01 & abs(fold_change_gene_arm) < 2 & sign(fold_change) == sign(fold_change_arm),
                                  FALSE, signif)) %>% 
        select(gene, chrom, padj, fold_change, signif, everything())
    
    # plot, write
    sample_name = str_extract(rdata, '(P-[0-9]{7}-T[0-9]{2}-IM[0-9]{1}|TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})')
    if (is.na(sample_name)) sample_name = str_replace(basename(rdata), '(_purity.Rdata|_purity.Rdata)', '')
    
    p = plot_events(snps, segs, gmat, diplogr) + ggtitle(sample_name)
    gmat = select(gmat, -pos, -pos_intragenic) %>% 
        mutate_at(vars(matches('fold_change|pval|padj')), funs(signif(., 2)))
    pt = plot_table(gmat, 10)
    
    write.table(gmat, paste0(sample_name, '_focal_events.txt'), sep = '\t', quote = F, col.names = T, row.names = F)
    ggplot() +
        coord_equal(xlim = c(0, 40), ylim = c(0, 30), expand = 0) +
        annotation_custom(ggplotGrob(p), xmin = 0, xmax = 40, ymin = 30, ymax = 15) +
        annotation_custom(pt$up, xmin = 1, xmax = 20, ymin = 0, ymax = 14) +
        annotation_custom(pt$down, xmin = 21, xmax = 40, ymin = 0, ymax = 14) +
        theme(panel.background = element_rect(fill = 'white'))
    ggsave(paste0(sample_name, '_focal_events.pdf'), w = 8.5, h = 5)
}

# Run as a script -------------------------------------------------------------------------------------------------
if (!interactive()) {
    args = commandArgs(F)
    if (is.null(args) | length(args)<1) {
        message('No input argument, run focal_events.R -h for help.')
        quit()
    }
    
    parser = ArgumentParser()
    parser$add_argument('-r', '--rdata', type = 'character', help = 'FACETS Rdata file')
    parser$add_argument('-g', '--gene_bed', type = 'character', help = 'Gene or gene-exon BED file')
    args = parser$parse_args()

    find_events(args$rdata, args$gene_bed)
    
}

