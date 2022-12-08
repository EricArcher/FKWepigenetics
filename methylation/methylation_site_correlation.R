rm(list = ls())
library(tidyverse)
library(swfscMisc)
load("meth posterior.rdata")

df <- cbind(df, post.smry)%>% 
  mutate(lo.meth = logOdds(median)) %>% 
  group_by(locus, site) %>% 
  mutate(meth.scale = scale(lo.meth)[, 1]) %>% 
  ungroup()

df.split <- split(df, list(df$locus, df$id))
site.meth.diff <- do.call(
  rbind,
  lapply(df.split, function(x) {
    if(length(unique(x$site)) < 2) return(NULL)
    site.pair <- t(combn(x$site, 2))
    colnames(site.pair) <- c("site.1", "site.2")
    cbind(
      data.frame(locus = unique(x$locus), id = unique(x$id)),
      site.pair,
      t(apply(site.pair, 1, function(i) { 
        pair.x <- filter(x, site %in% i)
        c(
          site.diff = abs(diff(pair.x$site)),
          meth.diff = abs(diff(pair.x$lo.meth))
        )
      }))
    )
  })
)

save.image("correlation ws.rdata")

pdf("site diff correlation.pdf", width = 15, height = 10)
diff.split <- split(site.meth.diff, site.meth.diff$locus)
for(x in diff.split) {
  p <- x %>% 
    ggplot(aes(site.diff, meth.diff)) +
    geom_point(size = 0.3, shape = 21, color = "black", fill = "white") +
    stat_summary(fun.data = mean_cl_normal, color = "red") + 
    geom_smooth(method = 'lm', formula = y~x) +
    labs(
      x = "Absolute distance between sites", 
      y = "Absolute difference between LO(methylation)",
      title = unique(x$locus)
    )
  print(p)
}
dev.off()
