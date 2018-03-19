# Read genotype probability object from file
read_probs_calc <- function(chr, datapath, allele = TRUE, probdir = "genoprob") {

  if(!allele & length(chr) > 1)
    stop("must supply at most one chr")

  ## Read in probs for selected chromosomes and cbind.
  prefix <- ifelse(allele, "aprobs_", "probs_")
  probs <- convert_probs(
    readRDS(file.path(datapath, probdir,
                      paste0(prefix, chr[1], ".rds"))))

  if(!allele) {

    ## Fix rownames of probs. Begin with "DO-".
    pr <- probs
    tmp <- substring(rownames(pr[[chr]]), 4)
    rownames(pr[[chr]]) <- tmp
    ## Sort them in increasing number order.
    pr[[chr]] <- pr[[chr]][order(as.numeric(tmp)),,, drop = FALSE]
    probs <- modify_object(probs, pr)
  }

  if(length(chr) > 1) {
    for(chri in chr[-1]) {
      probs <- cbind(probs,
                     convert_probs(
                       readRDS(file.path(datapath, probdir,
                                         paste0(prefix, chri, ".rds")))))
    }
  }
  probs
}
dimnames.calc_genoprob <- function (x)
{
  dnames <- lapply(x, dimnames)
  list(ind = dnames[[1]][[1]],
       gen = lapply(dnames, "[[", 2),
       mar = lapply(dnames, "[[", 3))
}
