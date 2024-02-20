# Figuring out the deal with Canadian Eastern Arctic - West Greenland 

sub <- df %>% filter(LME_NAME == "Canadian Eastern Arctic - West Greenland")
month_3_sub <- sub %>% filter(month == 3)
View(month_3_sub) # Thus, just 1 observation. That solves the matter

