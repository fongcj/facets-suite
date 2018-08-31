
## library(argparse)
## library(data.table)

## parser = ArgumentParser()


## parser$add_argument('-t', '--tumor_counts', type='character', help='Tumor counts file name')
## parser$add_argument('-n', '--normal_counts', type='character', help='Normal counts file name')
## parser$add_argument('-o', '--outfile', type='character', help='output file (gzipped')
## args=parser$parse_args()


#' @import data.table


#' @name read_counts
#' @title don't know
#' @description
#'
#' don't know
#'
#' @param file numeric first argument
#' @return don't know
read_counts = function(file){
  dt = fread(paste0("gunzip --stdout ", file), showProgress=FALSE)
  setkeyv(dt, c("Chrom", "Pos", "Ref", "Alt"))
  dt[,Refidx := NULL]
  dt[,TOTAL_depth := NULL]
  dt[,MAPQ_depth := NULL]
  dt[,ID := NULL]
  dt[,INS := NULL]
  dt[,DEL := NULL]
  dt[Ref != "N" & nchar(Ref) == 1 & nchar(Alt) == 1]
}


#' @name annotate_integer_copy_number
#' @title don't know
#' @description
#'
#' don't know
#'
#' @param tumor_counts numeric first argument
#' @param normal_counts numeric first argument
#' @param outfile numeric first argument
#' @param MINCOV_NORMAL numeric first argument
#' @return don't know
mergeTN_function = function(tumor_counts, normal_counts, outfile, MINCOV_NORMAL=25){


TUMOR = tumor_counts
NORMAL = normal_counts
GZOUT = gzfile(outfile)
## MINCOV_NORMAL=25

write("Reading normal ...", stdout())
NORMAL_dt = read_counts(NORMAL)
NORMAL_dt = NORMAL_dt[BASEQ_depth >= MINCOV_NORMAL]
write("done ...", stderr())

write("Reading tumor ...", stderr())
TUMOR_dt = read_counts(TUMOR)
write("done ...", stderr())

mergeTN = merge(TUMOR_dt, NORMAL_dt, suffixes = c(".TUM", ".NOR"))

setnames(mergeTN,
         c("Chrom", "Pos", "Ref", "Alt", "TUM.DP", "TUM.Ap", "TUM.Cp", 
           "TUM.Gp", "TUM.Tp", "TUM.An", "TUM.Cn", "TUM.Gn", "TUM.Tn", "NOR.DP", 
           "NOR.Ap", "NOR.Cp", "NOR.Gp", "NOR.Tp", "NOR.An", "NOR.Cn", "NOR.Gn", 
           "NOR.Tn"))

mergeTN[, Chrom := factor(Chrom, levels=c(c(1:22, "X", "Y", "MT"), paste0("chr", c(1:22, "X", "Y", "M"))))]
mergeTN = mergeTN[order(Chrom, Pos)]

write.table(mergeTN, file=GZOUT,
            quote = FALSE, 
            col.names = TRUE, 
            row.names = FALSE, 
            sep = "\t")


}
