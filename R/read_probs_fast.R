# Read genotype probability object from file
read_probs_fast <- function(chr, datapath, allele = TRUE,
                            probdir = "genoprob") {

  genoprob_dir <- file.path(datapath, probdir)
  if(!dir.exists(genoprob_dir))
    return(NULL)
  
  allele_rds <- paste0("fst", ifelse(allele, "_aprobs.rds", "_probs.rds"))
  if(!file.exists(probfile <- file.path(genoprob_dir, allele_rds)))
    return(NULL)

  ## Read in fst_genotype object (small).
  probs <- readRDS(probfile)

  ## Modify fst directory to match current datapath.
  pr <- unclass(probs)
  pr$fst <- file.path(genoprob_dir, basename(pr$fst))
  probs <- modify_object(probs, pr)
  
  # subset to desired chromosome(s)
  subset_probs_fast(probs, chr = chr)
}
subset_probs_fast <- function(probs, chr=NULL, mar=NULL) {
  qtl2fst::subset_fst_genoprob(probs, chr = chr, mar = mar)
}
