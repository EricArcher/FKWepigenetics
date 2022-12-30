rm(list = ls())
library(runjags)
library(tidyverse)
load("age_and_methylation_data.rdata")

chains <- 14
adapt <- 100
burnin <- 100000
total.samples <- 1000
thin <- 100


# Extract conversion, methylation, and coverage matrices ----------------

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

sites.to.keep <- site.smry %>% 
  filter(cov.median >= 1000) %>% 
  pull("loc.site")

locus <- locus.map %>% 
  select(loc.site, locus.num) %>% 
  deframe()


# Prior for groups x loci for coefficient clustering ----------------------

group.prior <- locus.map %>%
  filter(loc.site %in% sites.to.keep) %>% 
  group_by(locus) %>% 
  mutate(locus.site.num = as.numeric(factor(site))) %>% 
  ungroup() %>% 
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


# Select age data ---------------------------------------------------------

ages <- age.df %>% 
  filter(swfsc.id %in% ids.to.keep & age.confidence >= 1) %>% 
  #slice_sample(n = 1) %>% 
  arrange(swfsc.id)

# Create model data list --------------------------------------------------

model.data <- list(
  num.ind = nrow(ages),
  num.sites = length(sites.to.keep),
  num.loci = length(unique(locus.map$locus.num)),
  freq.meth = freq.meth[ages$swfsc.id, sites.to.keep, drop = FALSE],
  meth.cov = meth.cov[ages$swfsc.id, sites.to.keep, drop = FALSE],
  conversion = conversion[ages$swfsc.id, , drop = FALSE],
  locus = locus[sites.to.keep],
  num.groups = ncol(group.prior),
  group.prior = group.prior,
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
  ones = rep(1, nrow(ages)),
  zero.lo.meth = median.p[ages$swfsc.id, sites.to.keep]
)


# Run model ---------------------------------------------------------------

post <- run.jags(
  model = "model {
    intercept ~ dunif(-100, 500)
    
    for(l in 1:num.loci) {
      # Prior for probability of site group membership (Curtis and Ghosh 2011)
      pr.group[l, 1:num.groups] ~ ddirich(group.prior[l, 1:num.groups])
      
      # First group coefficient is always 0 (toggles site off)
      b.group[l, 1] <- 0
      
      # Prior for site group coefficients
      for(g in 2:num.groups) {
        b.group[l, g] ~ dnorm(0, 3e-3)
      }
    }
    
    # Model site coefficient
    for(s in 1:num.sites) {      
      # Prior for site group membership
      group[s] ~ dcat(pr.group[locus[s], ])
    
      # Model site coefficient (column 2)
      b.site[s] <- b.group[locus[s], group[s]]
    }
    
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
        lo.meth[i, s] <- ifelse(p.meth[i, s] <= 0, zero.lo.meth[i, s], logit(p.meth[i, s]))
      }
      
      # Linear model predicting age
      pred.age[i] <- intercept + inprod(b.site[], lo.meth[i, ])
      
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
  monitor = c(
    "deviance", "intercept", "pred.age", "group", 
    "b.site", "b.group", "lo.meth", "obs.p.meth"
  ), 
  data = model.data,
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  jags.refresh = 10,
  method = "parallel",
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.samples / chains),
  thin = thin
)
end.time <- Sys.time()
post$timetaken <- swfscMisc::autoUnits(post$timetaken)

save.image(format(end.time, "%y%m%d_%H%M.posterior.rdata"))
  
plot(
  post, 
  vars = c("deviance", "intercept", "pred.age", "group", "b.site"), 
  plot.type = c("trace", "histogram"),
  file = format(end.time, "%y%m%d_%H%M.plots.pdf")
)

post$timetaken