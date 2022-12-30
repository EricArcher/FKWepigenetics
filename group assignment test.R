rm(list = ls())
library(tidyverse)
library(runjags)
library(swfscMisc)

load("results/221229_0535.posterior.rdata")

p <- myFuncs::runjags2list(post)

dimnames(p$group)[[1]] <- sites.to.keep
dimnames(p$b.site)[[1]] <- sites.to.keep
dimnames(p$pred.age)[[1]] <- ages$swfsc.id
dimnames(p$p.meth)[1:2] <- list(ages$swfsc.id, sites.to.keep)
dimnames(p$obs.p.meth)[1:2] <- list(ages$swfsc.id, sites.to.keep)

by.locus <- lapply(
  split(locus.map$loc.site, locus.map$locus), 
  function(loc.site) t(p$group[rownames(p$group) %in% loc.site, ])
)


numSame <- function(mat) {
  nc <- ncol(mat)
  num.same <- matrix(nrow(mat), nrow = nc, ncol = nc)
  for(i in combn(1:nc, 2, simplify = FALSE)) {
    x <- sum(apply(mat[, i], 1, function(x) x[1] == x[2]))
    num.same[i[1], i[2]] <- num.same[i[2], i[1]] <- x
  }
  num.same / nrow(mat)
}


x <- replicate(4, {
  matrix(sample(1:6, 500, T), ncol = 10)
}, simplify = FALSE)

x2 <- lapply(x, numSame)

x3 <- do.call(abind::abind, c(x2, list(along = 3)))

x4 <- apply(x3, 1:2, mean)

x4