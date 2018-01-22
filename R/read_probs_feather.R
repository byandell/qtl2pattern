# Read genotype probability object from file
read_probs_feather <- function(chr, datapath, allele = TRUE) {

  ## Temporary kludge while directory names are changed.
  base_dir <- "genoprob"
  feather_dir <- file.path(datapath, base_dir)
  if(!dir.exists(feather_dir)) {
    base_dir <- "feather"
    feather_dir <- file.path(datapath, base_dir)
  }

  ## Read in feather_genotype object (small).
  probs <- readRDS(file.path(feather_dir,
                             ifelse(allele, "faprobs.rds", "fprobs.rds")))

  ## Modify feather directory to match current datapath.
  pr <- unclass(probs)
  oldpath <- pr$feather
  pr$feather <-
    stringr::str_replace(oldpath,
                         file.path(paste0(".*", base_dir), basename(oldpath)),
                         file.path(feather_dir, basename(oldpath)))
  probs <- modify_object(probs, pr)
  qtl2feather::subset_feather_genoprob(probs, chr = chr)
}
