rm(list = ls())
library(tidyverse)
source("../00_misc_funcs.R")
load("../age_and_methylation_data.rdata")

set.seed(1234)
n.ind <- 5

# get ids present in both age and methylation data
ids <- sort(intersect(cpg$id, age.df$swfsc.id))
#ids <- sample(ids, 5, replace = FALSE)

# filter ages for individuals with methylation data
ages <- age.df %>% 
  filter(swfsc.id %in% ids) %>% 
  arrange(swfsc.id)

# draw random samples of methylation across individuals and sites
lo.pct.meth <- pctMethSamples(cpg, n.ind) %>% 
  mutate(across(-c("id", "rep"), swfscMisc::logOdds)) %>% 
  filter(id %in% ids) %>% 
  select(-rep) %>% 
  arrange(id)

rep.ids <- match(lo.pct.meth$id, ids)

# match locus.site combinations to loci
locus <- setNames(locus.map$locus.num, locus.map$loc.site)

# create map of site groups (columns) by locus (rows) for coefficient clustering
site.group.map <- locus.map %>%
  select(locus.num, locus.site.num) %>%
  mutate(present = 1) %>%
  arrange(locus.site.num, locus.num) %>%
  pivot_wider(
    locus.num,
    names_from = "locus.site.num",
    values_from = "present",
    values_fill = 0
  ) %>%
  select(-locus.num) %>%
  as.matrix()

model.data <- list(
  num.ind = nrow(lo.pct.meth),
  num.sites = ncol(lo.pct.meth) - 1,
  num.loci = length(unique(locus.map$locus.num)),
  locus = locus[colnames(lo.pct.meth)[-1]],
  lo.pct.meth = as.matrix(lo.pct.meth[, -1]),
  num.groups = ncol(site.group.map),
  group.prior = site.group.map,
  age = ages$age.best[rep.ids],
  age.min = ages$age.min[rep.ids],
  age.max = ages$age.max[rep.ids],
  scale = ages$sn.scale[rep.ids],
  shape = ages$sn.shape[rep.ids],
  scale.m0 = ages$sn.scale[rep.ids] * swfscMisc::sn.m0(ages$sn.shape[rep.ids]),
  const = 1.1 * max(sapply(1:nrow(ages), function(i) {
    sn.pars <- unlist(ages[i, c("sn.location", "sn.scale", "sn.shape")])
    sn::dsn(swfscMisc::sn.mode(sn.pars), dp = sn.pars)
  })),
  ones = rep(1, nrow(lo.pct.meth))
)

saveRDS(model.data, "full_absolute_data.rds")