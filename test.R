rm(list = ls())
library(tidyverse)
library(runjags)
library(swfscMisc)

load("results/221229_0535.posterior.rdata")

p <- myFuncs::runjags2list(post, FALSE)

dimnames(p$group)[[1]] <- list(sites.to.keep)
dimnames(p$b.site)[[1]] <- list(sites.to.keep)
dimnames(p$pred.age)[[1]] <- list(ages$swfsc.id)
dimnames(p$p.meth)[1:2] <- list(ages$swfsc.id, sites.to.keep)
dimnames(p$obs.p.meth)[1:2] <- list(ages$swfsc.id, sites.to.keep)

obs.p.meth <- model.data$freq.meth / model.data$meth.cov
p.conv <- model.data$conversion[, 1] / model.data$conversion[, 2]
p.meth <- 1 - ((1 - obs.p.meth) / p.conv)
lo.meth <- ifelse(p.meth <= 0, logOdds(median.p), logOdds(p.meth))

cbind(lo.meth[1, ], p$lo.meth[1, , 1], lo.meth[1, ] - p$lo.meth[1, , 1])

pred.age <- t(apply(lo.meth, 1, function(x) {
  age.post <- p$intercept + colSums(p$b.site * x)
  c(
    median = median(age.post), 
    HDInterval::hdi(age.post), 
    setNames(range(age.post), c("min", "max"))
  )
})) %>% 
  as.data.frame() %>% 
  rownames_to_column("swfsc.id") %>% 
  left_join(
    select(ages, swfsc.id, age.best, age.min, age.max, age.confidence),
    by = "swfsc.id"
  ) %>% 
  mutate(age.diff = median - age.best)