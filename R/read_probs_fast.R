# Read genotype probability object from file
read_probs_fast <- function(chr, datapath, allele = TRUE,
                            fast = c("fst","feather"),
                            dirname = "genoprob") {

  ## Temporary kludge while directory names are changed.
  genoprob_dir <- file.path(datapath, dirname)
  if(!dir.exists(genoprob_dir)) {
    dirname <- "genoprob"
    genoprob_dir <- file.path(datapath, dirname)
  }
  
  fast <- match.arg(fast)
  
  ## Read in feather_genotype object (small).
  probs <- readRDS(file.path(genoprob_dir,
                             paste0(fast, ifelse(allele, "_aprobs.rds", "_probs.rds"))))

  ## Modify feather directory to match current datapath.
  pr <- unclass(probs)
  pr[[fast]] <- file.path(genoprob_dir, basename(pr[[fast]]))
  probs <- modify_object(probs, pr)
  
  # subset to desired chromosome(s)
  subsetfn <- switch(fast,
                     feather = qtl2feather::subset_feather_genoprob,
                     fst     =     qtl2fst::subset_fst_genoprob)
  subsetfn(probs, chr = chr)
}
