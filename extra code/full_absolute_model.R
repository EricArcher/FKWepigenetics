source("full_absolute_format.R")

rm(list = ls())
library(runjags)

chains <- 5
adapt <- 100
burnin <- 10
total.samples <- 100
thin <- 1
method <- "parallel"

nrep <- 50
params <- data.frame(
  adapt = sample(100:300, nrep, replace = TRUE),
  chains = sample(3:7, nrep, replace = TRUE),
  burnin = sample(10:50, nrep, replace = TRUE),
  total.samples = sample(10:50, nrep, replace = TRUE),
  thin = sample(1:20, nrep, replace = TRUE)
)

for(i in 1:nrow(params)) {
  
adapt <- params$adapt[i]
chains <- params$chains[i]
burnin <- params$burnin[i]
total.samples <- params$total.samples[i]
thin <- params$thin[i]

post <- run.jags(
  model = "model {
    # Prior for precision of coefficients
    # b.tau ~ dunif(0, 1e-4)
    
    for(l in 1:num.loci) {
      # Prior for locus inclusion
      pr.w.locus[l] ~ dunif(0, 1)
      w.locus[l] ~ dbern(pr.w.locus[l])
      
      # Prior for probability of site group membership (Curtis and Ghosh 2011)
      pr.group[l, 1:num.groups] ~ ddirich(group.prior[l, 1:num.groups])
    
      # Prior for site group coefficients
      for(g in 1:num.groups) {
        b[l, g] ~ dnorm(0, 1e-4)
      }
    }
    
    for(s in 1:num.sites) {
      # Prior for site group membership 
      b.group[s] ~ dcat(pr.group[locus[s], ])
    
      # Prior for site inclusion
      pr.w.site[s] ~ dunif(0, 1)
      w.site[s] ~ dbern(pr.w.site[s])
    
      # Model site coefficient
      b.prime[s] <- b[locus[s], b.group[s]] * w.locus[locus[s]] * w.site[s]
    }
    
    for(i in 1:num.ind) {
      # Linear model predicting age
      pred.age[i] <- 32.5 + inprod(b.prime[], lo.pct.meth[i, ])
      
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
  monitor = c("w.locus", "w.site", "b.group", "pr.group", "b.prime", "pred.age", "dsn"), 
  data = readRDS("full_absolute_data.rds"),
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
end.time <- Sys.time()
fname <- paste("test", i, "posterior_full_absolute_%y%m%d_%H%M.rdata")
fname <- format(end.time, fname)
save(adapt, chains, burnin, total.samples, thin, post, fname, end.time, file = fname)

}