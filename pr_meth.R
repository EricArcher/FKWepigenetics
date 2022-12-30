rm(list = ls())
library(runjags)
library(tidyverse)
load("age_and_methylation_data.rdata")

conversion <- non.cpg %>% 
  select(id, converted, coverage) %>% 
  arrange(id) %>% 
  column_to_rownames("id") %>% 
  as.matrix()

freq.meth <- cpg %>% 
  select(id, loc.site, freq.meth) %>% 
  pivot_wider(names_from = loc.site, values_from = freq.meth) %>% 
  column_to_rownames("id") %>% 
  as.matrix()

meth.cov <- cpg %>% 
  select(id, loc.site, corrected.cov) %>% 
  pivot_wider(names_from = loc.site, values_from = corrected.cov) %>% 
  column_to_rownames("id") %>% 
  as.matrix()

sites.to.keep <- apply(freq.meth, 2, function(x) !any(is.na(x)))
freq.meth <- freq.meth[, sites.to.keep]
meth.cov <- meth.cov[rownames(freq.meth), sites.to.keep]

post <- run.jags(
  model = "model {
    for(i in 1:num.ind) {
      # Probability of conversion
      p.conv[i] ~ dunif(0, 1)
      conversion[i, 1] ~ dbinom(p.conv[i], conversion[i, 2])
      
      for(s in 1:num.sites) {
        # Probability of observed methylation
        obs.p.meth[i, s] ~ dunif(0, 1)
        freq.meth[i, s] ~ dbinom(obs.p.meth[i, s], meth.cov[i, s])
        
        # Probability of actual methylation
        p.meth[i, s] <- 1 - ((1 - obs.p.meth[i, s]) / p.conv[i])
        # lo.meth[i, s] <- ifelse(p.meth[i, s] <= 0, -1e308, logit(p.meth[i, s]))
      }
    }
  }",
  monitor = c("deviance", "p.meth"), 
  data = list(
    num.ind = nrow(freq.meth),
    num.sites = ncol(freq.meth),
    freq.meth = freq.meth,
    meth.cov = meth.cov,
    conversion = conversion[rownames(freq.meth), ]
  ),
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  jags.refresh = 10,
  method = "parallel",
  n.chains = 10,
  adapt = 100,
  burnin = 10000,
  sample = 1000,
  thin = 10
)
post$timetaken <- swfscMisc::autoUnits(post$timetaken)

p <- myFuncs::runjags2list(post)
dimnames(p$p.meth)[1:2] <- list(rownames(meth.cov), colnames(meth.cov))

save.image("results/pr.meth.posterior.rdata")

post$timetaken