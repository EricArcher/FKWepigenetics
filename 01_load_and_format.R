rm(list = ls())
library(tidyverse)
library(sn)
library(swfscMisc)
library(runjags)
library(parallel)
library(abind)
library(gridExtra)
source("00_misc_funcs.R")
set.seed(1)


# Confidence rating mapping -----------------------------------------------

cr.mapping <- rbind(
  mk = c(0.05, 0.2, 0.55, 0.75, 1),
  sm = c(0.05, 0.2, 0.55, 0.75, 1),
  rb = c(0.2, 0.4, 0.6, 0.8, 1)
)


# Load and format age data ------------------------------------------------

crc.data <- readxl::read_xlsx(
  "data/Pseudorca_AgeEstimates_ChosenSamples_lasertagdates_2022JUNv5.xlsx",
  na = c(NA, "NA", "")
) 

# rename necessary columns
cols <- c(1, 4, 5, 6, 7, 21, 22, 25, 27, 29, 33)
colnames(crc.data)[cols] <- c(
  "crc.id", "sex", "swfsc.id", "biopsy.id", "date.biopsy",
  "age.best", "age.confidence", "age.min", "age.max.CRR", "ADULT.max.best", 
  "pair.id"
)

# split out and replicate ids
crc.data <- do.call(
  rbind,
  lapply(1:nrow(crc.data), function(i) {
    id <- strsplit(crc.data$swfsc.id[i], ",")[[1]] %>% 
      stringi::stri_trim()
    df <- crc.data[rep(i, length(id)), ]
    df$swfsc.id <- id
    df
  })
)

# format labid
crc.data <- crc.data %>% 
  mutate(
    swfsc.id = zero.pad(as.numeric(swfsc.id)),
    swfsc.id = ifelse(is.na(swfsc.id), NA, paste0("z0", swfsc.id))
  )

# choose best maximum age
crc.data$age.max <- sapply(1:nrow(crc.data), function(i) {
  min(crc.data$ADULT.max.best[i], crc.data$age.max.CRR[i], na.rm = TRUE)
})


# format age class
age.class.map <- data.frame(
  class = c("calf", "juvenile", "sub.adult", "adult.f.sm", "adult.pm"),
  f.ages = c(2, 5, 9, max(crc.data$age.best), NA),
  m.ages = c(2, 8, 14, 24, max(crc.data$age.best))
)
crc.data <- crc.data %>% 
  mutate(
    age.class = ifelse(
      sex == "Female", 
      as.numeric(cut(age.best, c(0, age.class.map$f), include = T, right = T)),
      as.numeric(cut(age.best, c(0, age.class.map$m), include = T, right = T))
    ),
    age.class = factor(
      age.class.map$class[age.class], 
      levels = age.class.map$class
    )
  )

crc.cols <- c(
  "crc.id", "swfsc.id", "biopsy.id", "date.biopsy", "sex",
  "age.best", "age.confidence", "age.min", "age.max", "age.class", "pair.id"
)
age.df <- crc.data[, crc.cols] %>% 
  mutate( 
    age.range = age.max - age.min,
    confidence.p = cr.mapping["mk", age.confidence],
    dens.unif = 1 / age.range,
    date.biopsy = as.POSIXct(date.biopsy, format = "%m/%d/%y")
  ) %>%
  filter(!is.na(age.min)) %>% 
  arrange(swfsc.id)
attributes(age.df$date.biopsy)$tzone <- "HST"

# read glgs and compute minimum glg
glgs <- readxl::read_xlsx("data/FKW GLG Age Results_KMR.xlsx") %>% 
  setNames(c("biopsy.id", "kmr", "min.kmr", "sjc", "crc", "notes"))
glgs$min.glgs <- sapply(1:nrow(glgs), function(i) {
  x <- min(glgs$kmr[i], glgs$min.kmr[i], glgs$sjc[i], glgs$crc[i], na.rm = T)
  if(is.finite(x)) x else NA
})

# make minimum age = minimum glg when present
age.df <- age.df %>% 
  left_join(select(glgs, biopsy.id, min.glgs), by = "biopsy.id") %>% 
  mutate(age.min = ifelse(!is.na(min.glgs), min.glgs, age.min))


age.sn.optim <- mclapply(1:nrow(age.df), function(i) {
  sn.params(
    age.df$age.best[i],
    age.df$age.min[i],
    age.df$age.max[i],
    p = 0.975,
    shape.const = 4
  )
}, mc.cores = 14)
age.sn.params <- t(sapply(age.sn.optim, function(x) x$par))
colnames(age.sn.params) <- c("sn.location", "sn.scale", "sn.shape")
age.df <- cbind(age.df, age.sn.params)


# Load and format methylation data ----------------------------------------

load("data/Pcra.epi.data.for.Eric.Rdata")
load("data/res.sum.Rdata")

epi.df <- do.call(
  rbind,
  lapply(names(res.sum), function(locus) {
    do.call(
      rbind,
      lapply(names(res.sum[[locus]]), function(type) {
        x <- res.sum[[locus]][[type]]
        if(is.null(x)) return(NULL)
        x$coverage %>% 
          select(-avg) %>% 
          pivot_longer(-id, names_to = "site", values_to = "coverage") %>% 
          left_join(
            pivot_longer(x$freq.meth, -id, names_to = "site", values_to = "freq.meth"),
            by = c("id", "site")
          ) %>% 
          left_join(
            pivot_longer(x$errors, -id, names_to = "site", values_to = "errors"),
            by = c("id", "site")
          ) %>% 
          mutate(locus = locus, type = type)
      })
    )
  })
) %>%   
  mutate(
    site = as.numeric(site),
    loc.site = paste0(locus, "_", zero.pad(site)),
    id.site = paste0(id, "_", loc.site),
    corrected.cov = coverage - errors,
    pct.meth = freq.meth / corrected.cov,
    type = gsub(".sum", "", type)
  ) %>% 
  arrange(type, id, locus, site) 

locus.map <- epi.df %>% 
  select(loc.site, locus, site) %>% 
  unique() %>% 
  mutate(locus.num = as.numeric(factor(locus))) %>% 
  arrange(locus, site) %>% 
  group_by(locus) %>% 
  mutate(locus.site.num = as.numeric(factor(site))) %>% 
  ungroup()


# Calculate conversion rate -----------------------------------------------

non.cpg <- epi.df %>%
  filter(type == "non.CpG") %>%
  group_by(id) %>%
  summarize(
    coverage = sum(corrected.cov, na.rm = TRUE),
    freq.meth = sum(freq.meth, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    converted = coverage - freq.meth,
    pct.conv = converted / coverage,
    conv.shape1 = converted + 1,
    conv.shape2 = (coverage - converted) + 1
  )

cpg <- epi.df %>% 
  filter(type == "CpG") %>% 
  mutate(
    meth.shape1 = ifelse(
      is.na(freq.meth) | is.na(corrected.cov), 
      1, freq.meth + 1
    ),
    meth.shape2 = ifelse(
      is.na(freq.meth) | is.na(corrected.cov), 
      1, (corrected.cov - freq.meth) + 1
    )
  ) %>% 
  left_join(
    select(non.cpg, id, pct.conv, conv.shape1, conv.shape2),
    by = "id"
  )


# Minimum Pr(methylation) Bayesian model ----------------------------------

meth.cov <- cpg %>% 
  select(id, loc.site, corrected.cov) %>% 
  pivot_wider(names_from = loc.site, values_from = corrected.cov) %>% 
  column_to_rownames("id") %>% 
  as.matrix()

sites.to.keep <- apply(meth.cov, 2, function(x) !any(is.na(x)))
meth.cov <- meth.cov[, sites.to.keep]

post <- run.jags(
  model = "model {
    for(i in 1:num.ind) {
      for(s in 1:num.sites) {
        p.meth[i, s] ~ dunif(0, 1)
        zeroes[i, s] ~ dbinom(p.meth[i, s], meth.cov[i, s])
      }
    }
  }",
  monitor = c("deviance", "p.meth"), 
  data = list(
    num.ind = nrow(meth.cov),
    num.sites = ncol(meth.cov),
    meth.cov = meth.cov,
    zeroes = matrix(0, nrow = nrow(meth.cov), ncol = ncol(meth.cov))
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
save(post, file = "results/min.pr.meth.posterior.rdata")

p <- myFuncs::runjags2list(post)
dimnames(p$p.meth)[1:2] <- list(rownames(meth.cov), colnames(meth.cov))
median.p <- apply(p$p.meth, 1:2, median)


# Save data ---------------------------------------------------------------

save(
  age.df, age.class.map, epi.df, locus.map, cpg, non.cpg, median.p,
  file = "age_and_methylation_data.rdata"
)


# Plot age distributions --------------------------------------------------

graphics.off()
pdf(paste0("results/plots - estimated age priors.pdf"), height = 12, width = 12)

p <- age.df %>% 
  mutate(
    id = factor(swfsc.id, levels = .$swfsc.id[order(-age.best, -age.min, -age.max)]),
    Confidence = factor(age.confidence)
  ) %>% 
  ggplot(aes(y = id)) +
  geom_segment(
    aes(x = age.min, xend = age.max, yend = id, color = Confidence), 
    size = 2, alpha = 0.8
  ) +
  geom_point(
    aes(x = age.best, color = Confidence), 
    shape = 21, size = 3.5, fill = "white"
  ) + 
  labs(x = "Age", y = NULL) + 
  scale_x_continuous(breaks = seq(0, max(age.df$age.max), 5)) +
  theme(
    text = element_text(size = 10),
    legend.position = "top"
  )
print(p)

p <- age.df %>% 
  mutate(Confidence = factor(age.confidence)) %>% 
  ggplot(aes(age.best, age.range)) +
  geom_point(aes(fill = Confidence), color = "white", shape = 21, size = 4) +
  labs(x = "CRC best age", y = "Maximum - minimum age") +
  theme(legend.position = "top")
print(p)

age.df <- arrange(age.df, swfsc.id)

for(age.cr in split(age.df, age.df$age.confidence)) {
  x.age <- seq(min(age.cr$age.min) * 0.8, max(age.cr$age.max) * 1.2, length.out = 1000)
  conf <- unique(age.cr$age.confidence)
  pct <- c(0.05, 0.1, 0.5, 0.75, 1)[conf]
  color <- scales::hue_pal()(5)[conf]
  
  p <- lapply(1:nrow(age.cr), function(i) {
    adj.dens <- confidence.dens(
      dsn(x.age, dp = unlist(age.cr[i, c("sn.location", "sn.scale", "sn.shape")])), 
      pct, age.cr$dens.unif[i]
    )
    sn <- data.frame(age = x.age, dens = adj.dens) %>% 
      mutate(dens = ifelse(age < age.cr$age.min[i] | age > age.cr$age.max[i], 0, dens))
    
    ggplot() +
      geom_vline(xintercept = age.cr$age.best[i]) +
      geom_line(aes(x = age, y = dens), data = sn, color = color, size = 1.5) + 
      xlim(range(x.age)) +
      ggtitle(age.cr$swfsc.id[i]) +
      theme(axis.title = element_blank())
  })
  p$top <- paste0(
    "Confidence rating = ", unique(age.cr$age.confidence), 
    ", % from uniform = ", pct
  )
  p$bottom <- "Age"
  p$left <- "Density"
  do.call(grid.arrange, p)
}

dev.off()


# Plot percent methylation posteriors -------------------------------------

smry <- pctMethSamples(cpg, 1000) %>% 
  pivot_longer(-c(id, rep), names_to = "loc.site", values_to = "pct.meth") %>% 
  group_by(id, loc.site) %>% 
  summarize(
    mode = modeest::mlv(pct.meth, method = "Venter"),
    lci = HDInterval::hdi(pct.meth, credMass = 0.95)[1],
    uci = HDInterval::hdi(pct.meth, credMass = 0.95)[2],
    .groups = "drop"
  ) %>% 
  left_join(locus.map, by = "loc.site") %>% 
  mutate(id = factor(id, levels = sort(unique(.$id), d = TRUE)))

ids.by.age <- age.df %>%
  arrange(age.best, age.min, age.max, age.confidence, swfsc.id) %>%
  pull("swfsc.id") %>% 
  rev()

graphics.off()
pdf("results/plots - locus posteriors.pdf", width = 30, height = 20)
for(loc in split(smry, smry$locus)) {
  for(i in 1:2) {
    if(i == 2) loc <- mutate(loc, id = factor(id, levels = ids.by.age))
    p <- loc %>% 
      ggplot(aes(y = id)) +
      geom_segment(aes(x = lci, xend = uci, yend = id)) +
      geom_point(aes(x = mode)) +
      scale_x_continuous(n.breaks = 4) +
      facet_grid(~ site, scale = "free_x") +
      labs(
        x = "% methylation", 
        title = paste(unique(loc$locus), ifelse(i == 1, "by ID", "by age"))
      ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  }
}
dev.off()