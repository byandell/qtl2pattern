# Read genotype probability object from file
read_probs_fast <- function(chr, datapath, allele = TRUE,
                            probdir = "genoprob") {

  ## Temporary kludge while directory names are changed.
  genoprob_dir <- file.path(datapath, probdir)
  if(!dir.exists(genoprob_dir)) {
    probdir <- "genoprob"
    genoprob_dir <- file.path(datapath, probdir)
    if(!dir.exists(genoprob_dir))
      return(NULL)
  }
  
  allele_rds <- ifelse(allele, "_aprobs.rds", "_probs.rds")
  if(file.exists(probfile <- 
    file.path(genoprob_dir, paste0("fst", allele_rds)))) {
    fast <- "fst"
  } else {
    if(!file.exists(probfile <- 
                    file.path(genoprob_dir, paste0("feather", allele_rds))))
      return(NULL)
    fast <- "feather"
  }

  ## Read in feather_genotype object (small).
  probs <- readRDS(probfile)

  ## Modify feather directory to match current datapath.
  pr <- unclass(probs)
  pr[[fast]] <- file.path(genoprob_dir, basename(pr[[fast]]))
  probs <- modify_object(probs, pr)
  
  # subset to desired chromosome(s)
  subset_probs_fast(probs, chr = chr)
}
subset_probs_fast <- function(probs, chr=NULL, mar=NULL) {
  subsetfn <- switch(class(probs)[1],
                     feather_genoprob = qtl2feather::subset_feather_genoprob,
                     fst_genoprob     =     qtl2fst::subset_fst_genoprob)
  subsetfn(probs, chr = chr, mar = mar)
}
