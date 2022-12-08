rm(list = ls())
library(mgcv)
library(tidyverse)
library(swfscMisc)

load("methylation data.rdata")
age.df <- readRDS("age df.rds")

df <- complete.epi.df %>% 
  mutate(
    loc.site = paste0(locus, "_", site),
    pct.meth = freq.meth / (coverage - errors),
    lo.meth = logOdds(pct.meth)
  ) %>% 
  filter(pct.meth > 0 & type == "CpG") %>% 
  select(id, loc.site, lo.meth) %>% 
  pivot_wider(
    names_from = "loc.site",
    values_from = "lo.meth"
  ) %>% 
  left_join(
    select(age.df, id, age.point),
    by = "id"
  ) %>% 
  select(age.point, everything()) %>% 
  column_to_rownames("id") %>% 
  as.data.frame()

epi.lm.gam <- sapply(setdiff(colnames(df), "age.point"), function(x) {
  if(length(unique(df[, x])) < 6) return(NULL)
  site.term <- paste0("s(", x, ", k = 4, bs = 'ts')", collapse = " + ")
  list(
    lm = lm(
      age.point ~ ., 
      data = df[, c("age.point", x)], 
      na.action = "na.omit"
    ),
    gam = gam(
      as.formula(paste0("age.point ~ ", site.term)),
      data = df,
      na.action = "na.omit"
    )
  )
}, simplify = FALSE)

lm.coeffs <- sapply(
  epi.lm.gam, 
  function(x) if(!is.null(x)) coef(x$lm) else NULL,
  simplify = FALSE
)
lm.coeffs <- do.call(
  rbind,
  lm.coeffs[!sapply(lm.coeffs, is.null)]
)

pdf("gam site plots.pdf")
for(x in epi.gams) if(!is.null(x)) plot(x$gam)
dev.off()
