rm(list = ls())
library(tidyverse)
library(runjags)
library(swfscMisc)

load("results/221229_0535.posterior.rdata")

p <- myFuncs::runjags2list(post)

dimnames(p$group)[[1]] <- list(sites.to.keep)
dimnames(p$b.site)[[1]] <- list(sites.to.keep)
dimnames(p$pred.age)[[1]] <- list(ages$swfsc.id)
dimnames(p$p.meth)[1:2] <- list(ages$swfsc.id, sites.to.keep)
dimnames(p$obs.p.meth)[1:2] <- list(ages$swfsc.id, sites.to.keep)


pctSameGroup <- function(mat, pct = TRUE) {
  nr <- nrow(mat)
  num.same <- matrix(ncol(mat), nrow = nr, ncol = nr)
  for(i in combn(1:nr, 2, simplify = FALSE)) {
    x <- apply(mat[i, ], 2, function(mat.pairs) mat.pairs[1] == mat.pairs[2])
    num.same[i[1], i[2]] <- num.same[i[2], i[1]] <- sum(x)
  }
  if(pct) num.same / ncol(mat) else num.same
}


by.locus <- mclapply(
  split(locus.map$loc.site, locus.map$locus), 
  function(loc.site) {
    pctSameGroup(p$group[rownames(p$group) %in% loc.site, ])
  },
  mc.cores = 14
)
