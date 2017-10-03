### write tab delimited
write.tab <- function(...){
  write.table(..., quote = F,
              col.names=T,
              row.names=F,
              sep='\t')
}

### fread a bunch of files and rbind them
fread_rbind <- function(filenames, fn = fread){
  rbindlist(idcol = "filename",
            sapply(USE.NAMES = T,
                   simplify = F,
                   filenames,
                   fn))
}

#' Title
#'
#' @param arg_line
#' @import data.table
#' @export
facets.somatic <- function(arg_line = NA){

  ### process args
  if(!is.na(arg_line) | interactive()) {
    # print("reading from arg_line")
    raw_args <- unlist(stringr::str_split(arg_line, " "))
  } else {
    # print("batch")
    raw_args <- commandArgs(TRUE)
  }
  # print(raw_args)

  option_list <- list(
    optparse::make_option(c("-s", "--samplefile"),
                          type="character", default="samples.txt",
                          help="sample file"),
    optparse::make_option(c("-o", "--outdir"),
                          type="character", default=NA,
                          help="output directory"),
    optparse::make_option(c('-t', '--targetFile'),
                          type='character', default='IMPACT468',
                          help="IMPACT341/410/468, or a Picard interval list file of gene target coordinates [default IMPACT468]")
  )
  if(any(sapply(
    option_list,
    function(option){
      option@short_flag == "-g"
    }))){
    stop("cannot use short option '-g', conflicts with Rscript --gui")
  }

  opts <- optparse::parse_args(
    optparse::OptionParser(option_list=option_list),
    args = raw_args,
    positional_arguments = TRUE
  )

  samplefile <- opts$options$samplefile
  outdir <- opts$options$outdir
  targetFile <- opts$options$targetFile

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  s2c <- suppressWarnings(fread(samplefile))

  original_directory <- getwd()
  setwd(dirname(samplefile))
  if(!all(file.exists(s2c[, cncf]))) {
    warning("Cannot find all paths in sample table\n")
  } else {
    s2c[, cncf := tools::file_path_as_absolute(cncf), cncf]
  }
  setwd(original_directory)

  if (targetFile == "IMPACT341") {
    gene_targets = facets.somatic::IMPACT341_targets
  } else if (targetFile == "IMPACT410") {
    gene_targets = facets.somatic::IMPACT410_targets
  } else if (targetFile == "IMPACT468") {
    gene_targets = facets.somatic::IMPACT468_targets
  } else {
    # Note the target file needs to not only be in the PICARD interval list format
    # But the names must match the regex: /GENESYMBOL_.*/ (e.g. TP53_target_02)
    gene_targets <- suppressWarnings(fread(paste0('grep -v "^@" ', targetFile)))
    setnames(gene_targets, c("chr", "start", "end", "strand", "name"))
    setkey(gene_targets, chr, start, end)
  }

  cncf <- generate_cncf_file(s2c)
  write.tab(cncf, file.path(outdir, "cncf.txt"))

  seg <- generate_seg(cncf)
  write.tab(seg, file.path(outdir, "cncf.seg"))

  arm_level_calls = get_arm_level_calls(cncf)
  write.tab(arm_level_calls, file.path(outdir, "armLevel.txt"))

  gene_level_calls = get_gene_level_calls(cncf, gene_targets = gene_targets)
  write.tab(gene_level_calls, file.path(outdir, "geneLevel.txt"))
}
