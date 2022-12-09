rm(list = ls())
library(runjags)
library(tidyverse)
source("../00_misc_funcs.R")
load("../age_and_methylation_data.rdata")

#set.seed(1234)

chains <- 14
adapt <- 100
burnin <- 5000000
total.samples <- 10000
thin <- 100
method <- "parallel"


# Extract conversion,  methylation,  and coverage matrices ----------------

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


# IDs with no missing sites -----------------------------------------------

ids.to.keep <- cpg %>% 
  filter(id %in% age.df$swfsc.id & type == "CpG") %>% 
  group_by(id) %>% 
  summarize(
    is.complete = !any(is.na(freq.meth)),
    median.cov = median(coverage), 
    .groups = "drop"
  ) %>% 
  filter(median.cov >= 1000 & is.complete) %>% 
  pull("id")


# Summarize and select sites ----------------------------------------------

site.smry <- cpg %>% 
  filter(id %in% ids.to.keep & type == "CpG") %>% 
  mutate(
    lo.meth = swfscMisc::logOdds(freq.meth / corrected.cov),
    lo.meth = ifelse(is.infinite(lo.meth), NA, lo.meth)
  ) %>% 
  left_join(age.df, by = c(id = "swfsc.id")) %>% 
  group_by(loc.site) %>% 
  summarize(
    age.cor = cor(lo.meth, age.best, use = "complete.obs"),
    lo.var = var(lo.meth, na.rm = TRUE),
    cov.var = var(corrected.cov, na.rm = TRUE),
    cov.median = median(corrected.cov, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  arrange("loc.site")

plot(sort(site.smry$lo.var), type = "l")

plot(sort(abs(site.smry$age.cor)), type = "l")

plot(sort(site.smry$cov.median), type = "l")

plot(site.smry$cov.median, site.smry$lo.var)

sites.to.keep <- site.smry %>% 
  filter(cov.median >= 1000) %>% 
  #filter(abs(age.cor) >= 0.5 & cov.median >= 1000) %>% 
  pull("loc.site")


# Select age data ---------------------------------------------------------

#ids.to.keep <- sample(ids.to.keep, 50)

ages <- age.df %>% 
  filter(swfsc.id %in% ids.to.keep & age.confidence == 5) %>% 
  arrange(swfsc.id)


# Create model data list --------------------------------------------------

model.data <- list(
  num.ind = nrow(ages),
  num.sites = length(sites.to.keep),
  freq.meth = freq.meth[ages$swfsc.id, sites.to.keep],
  meth.cov = meth.cov[ages$swfsc.id, sites.to.keep],
  conversion = conversion[ages$swfsc.id, ],
  age = ages$age.best,
  age.min = ages$age.min,
  age.max = ages$age.max,
  scale = ages$sn.scale,
  shape = ages$sn.shape,
  scale.m0 = ages$sn.scale * swfscMisc::sn.m0(ages$sn.shape),
  const = 1.1 * max(sapply(1:nrow(ages), function(i) {
    sn.pars <- unlist(ages[i, c("sn.location", "sn.scale", "sn.shape")])
    sn::dsn(swfscMisc::sn.mode(sn.pars), dp = sn.pars)
  })),
  ones = rep(1, nrow(ages))
)


# Run model ---------------------------------------------------------------

post <- run.jags(
  model = "model {
    intercept ~ dunif(-100, 500) #dnorm(0, 1 / sqrt(40))
    
    # Model site coefficient
    for(s in 1:num.sites) {
      b[s] ~ dnorm(0, 0.02) #1 / sqrt(50))
    }
    
    for(i in 1:num.ind) {
      p.conv[i] ~ dunif(0, 1)
      conversion[i, 1] ~ dbinom(p.conv[i], conversion[i, 2])
      
      for(s in 1:num.sites) {
        p.meth[i, s] ~ dunif(0, 1)
        freq.meth[i, s] ~ dbinom(p.meth[i, s], meth.cov[i, s])
        pr.meth.hat[i, s] <- 1 - ((1 - p.meth[i, s]) / p.conv[i])
        lo.meth[i, s] <- ifelse(pr.meth.hat[i, s] < 0, -1000, logit(pr.meth.hat[i, s]))
      }
      
      # Linear model predicting age
      pred.age[i] <- intercept + inprod(b[], lo.meth[i, ])
      
      # Normal age using predicted age as mode of skew normal
      norm.age[i] <- (pred.age[i] - (age[i] - scale.m0[i])) / scale[i]
        
      # Skew normal likelihood, restricted to min and max ages
      dsn[i] <- ifelse(
        pred.age[i] >= age.min[i] && pred.age[i] <= age.max[i],
        (2 / scale[i]) * dnorm(norm.age[i], 0, 1) * pnorm(norm.age[i] * shape[i], 0, 1),
        1e-308
      )
      
      # Ones likelihood (Kruschke et al 2010, DBDA2E Section 8.6.1 pp. 214-215)
      ones[i] ~ dbern(dsn[i] / const)
    }
  }",
  monitor = c("deviance", "intercept", "b", "pred.age", "pr.meth.hat", "p.meth", "dsn"), 
  data = model.data,
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  jags.refresh = 10,
  method = method,
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.samples / chains),
  thin = thin
)
units(post$timetaken) <- "hours"

postSummary(post, c("deviance", "intercept", "b", "pred.age", "dsn"), "sn.age_binom.meth")
