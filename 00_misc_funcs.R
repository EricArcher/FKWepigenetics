# estimate skew-normal parameters from best, min, and max
sn.params <- function(
    age.mode, age.min, age.max, p, shape.const = 5, scale.const = 0.5, 
    maxit = 5000
) {
  age.mean <- mean(c(age.min, age.max))
  shape <- shape.const * (age.mean - age.mode) / (age.mean - age.min)
  scale <- scale.const * (age.max - age.min)
  optim(
    par = c(age.mode, scale, shape), 
    fn = function(dp, obs, p) {
      if(dp[2] <= 0) dp[2] <- .Machine$double.eps 
      abs(obs[1] - swfscMisc::sn.mode(dp)) +
        abs(obs[2] - sn::qsn(1 - p, dp = dp, solver = "RFB")) +
        abs(obs[3] - sn::qsn(p, dp = dp, solver = "RFB"))
    }, 
    obs = c(age.mode, age.min, age.max), 
    p = p,
    control = list(
      maxit = maxit, 
      abstol = .Machine$double.eps,
      reltol = 1e-12
    )
  )
}

confidence.dens <- function(dens.sn, p, dens.unif) {
  dens.unif - ((dens.unif - dens.sn) * p)
}

pctMethSamples <- function(df, n) { 
  p <- do.call(
    rbind,
    parallel::mclapply(1:nrow(df), function(i) {
      with(df, {
        p.meth <- rbeta(n, meth.shape1[i], meth.shape2[i])
        p.conv <- rbeta(n, conv.shape1[i], conv.shape2[i])
        1 - ((1 - p.meth) / p.conv)
      })
    }, mc.cores = 8)
  )
  lt.0 <- p <= 0
  p[lt.0] <- runif(sum(lt.0), 1e-300, min(p[!lt.0]) / 10)
  
  p %>% 
    as.data.frame() %>% 
    setNames(1:ncol(p)) %>%
    mutate(id.site = df$id.site) %>% 
    pivot_longer(-id.site, names_to = "rep", values_to = "pct.meth") %>% 
    mutate(rep = as.numeric(rep)) %>% 
    left_join(select(df, id.site, id, loc.site), by = "id.site") %>% 
    select(-id.site) %>% 
    pivot_wider(
      all_of(c("id", "rep")), 
      names_from = "loc.site", 
      values_from = "pct.meth"
    ) %>% 
    arrange(id, rep)
}

postSummary <- function(post, vars, label) {
  units(post$timetaken) <- "hours"
  end.time <- Sys.time()
  time.stamp <- format(end.time, "%y%m%d_%H%M")
  fname <- paste0("posterior.summary_", label, "_", time.stamp, ".rdata")
  
  post.smry <- summary(post, vars = vars)
  save(end.time, post, post.smry, fname, file = fname)
  
  graphics.off()
  pdf(paste0("plots_", label, "_", time.stamp, ".pdf"))
  plot(post, vars = vars, plot.type = c("trace", "histogram"))
  dev.off()
  
  View(post.smry)
}