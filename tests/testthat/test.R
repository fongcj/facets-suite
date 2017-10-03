test_that("run facets.somatic", {
  samplefile <- system.file(package = "facets.somatic", "example_data/samples.txt")
  maffile <- system.file(package = "facets.somatic", "example_data/tcga.maf")
  outdir <- system.file(package = "facets.somatic", "example_data/example_output")
  arg_line <- paste("-s", samplefile, "-m", maffile, "-o", outdir)
  cat("Rscript -e 'facets.somatic::facets.somatic()'", arg_line, "\n")
  facets.somatic(arg_line)
})
