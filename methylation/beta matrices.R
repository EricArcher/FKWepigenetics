# meth.shape <- abind::abind(
#   meth.shape1 = cpg %>% 
#     filter(id %in% ids & loc.site %in% sites) %>% 
#     select(id, loc.site, meth.shape1) %>% 
#     pivot_wider(names_from = loc.site, values_from = meth.shape1) %>% 
#     column_to_rownames("id") %>% 
#     as.matrix(),
#   meth.shape2 = cpg %>% 
#     filter(id %in% ids & loc.site %in% sites) %>% 
#     select(id, loc.site, meth.shape2) %>% 
#     pivot_wider(names_from = loc.site, values_from = meth.shape2) %>% 
#     column_to_rownames("id") %>% 
#     as.matrix(),
#   along = 3
# )

# conv.shape <- non.cpg %>% 
#   filter(id %in% ids) %>% 
#   select(id, conv.shape1, conv.shape2) %>% 
#   arrange(id) %>% 
#   column_to_rownames("id") %>% 
#   as.matrix()