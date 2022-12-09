rm(list = ls())
library(runjags)
library(tidyverse)
load("posterior_sn.age_binom.meth_221122_1057.rdata")

p <- myFuncs::runjags2list(post)

pr.meth <- t(as.matrix(p$pr.meth))

pr.meth.arr <- array(
  as.vector(pr.meth),
  dim = c(
    nrow(model.data$freq.meth), 
    ncol(model.data$freq.meth), 
    ncol(pr.meth)
  ),
  dimnames = list(
    rownames(model.data$freq.meth),
    colnames(model.data$freq.meth),
    1:ncol(pr.meth)
  )
)

post.compare <- do.call(rbind, apply(
  expand.grid(id = rownames(model.data$freq.meth), site = colnames(model.data$freq.meth)), 1,
  function(x) { 
    id <- x[1]
    site <- x[2]
    freq.meth <- model.data$freq.meth[id, site]
    coverage <- model.data$coverage[id, site]
    obs.p <- freq.meth / coverage
    post.x <- pr.meth.arr[id, site, ]
    p.gte <- mean(pbinom(freq.meth, coverage, post.x) >= obs.p)
    data.frame(
      id = id, loc.site = site, freq.meth = freq.meth, coverage = coverage, obs.p = obs.p, 
      obs.p.lci = qbinom(0.025, coverage, obs.p) / coverage,
      obs.p.hci = qbinom(0.975, coverage, obs.p) / coverage,
      post.median = median(post.x),
      post.lci = quantile(post.x, 0.025),
      post.hci = quantile(post.x, 0.975),
      p.gte = p.gte
    )
  },
  simplify = FALSE
)) %>% 
  left_join(locus.map, by = "loc.site")

save.image("meth compare ws.rdata")

graphics.off()
pdf("test.pdf", height = 14, width = 14)
for(df in split(post.compare, post.compare$locus)) {
  p <- ggplot(df) +
    geom_segment(
      aes(x = obs.p.lci, xend = obs.p.hci, y = id, yend = id), 
      color = "red", linewidth = 2
    ) +
    geom_point(aes(x = obs.p, y = id), fill = "white", color = "red", shape = 21, size = 2.5) +
    geom_segment(
      aes(x = post.lci, xend = post.hci, y = id, yend = id)
    ) +
    geom_point(aes(x = post.median, y = id), shape = 20) +
    scale_x_continuous(n.breaks = 3) +
    labs(x = "% methylated", title = unique(df$locus)) +
    facet_wrap(~ site, nrow = 1, scale = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank()
    )
  print(p)
}
dev.off()


min.cov.by.id <- post.compare %>% 
  group_by(id) %>% 
  summarize(min.cov = min(coverage), .groups = "drop") %>% 
  arrange(min.cov) 

cov.by.id <- post.compare %>% 
  group_by(id) %>% 
  summarize(median.cov = median(coverage), .groups = "drop") %>% 
  arrange(median.cov) 

ids.to.delete <- cov.by.id %>% 
  filter(median.cov < 1000) %>% 
  pull("id")

cov.by.site <- post.compare %>% 
  filter(!id %in% ids.to.delete) %>% 
  group_by(locus, site) %>% 
  summarize(median.cov = median(coverage), .groups = "drop") %>% 
  arrange(median.cov) 

sites.to.delete <- cov.by.site %>% 
  filter(median.cov < 500)
