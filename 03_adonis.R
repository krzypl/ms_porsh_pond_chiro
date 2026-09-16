library(tidyverse)
library(vegan)
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 9), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")),
         counts = rowSums(select(., Ablabesmyia:Zavreliella))) %>%
  filter(!counts < 12)

meta_cnts <- data.frame(
  sample_type = cnts$sample_type,
  core_id = cnts$core_id)

spe <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

spe_rowsums <- rowSums(spe)

spe_perc <- (spe/spe_rowsums)*100

#Test whether there is a lake-scale difference in chironomid assamblages between the sample classes
set.seed(12)
pw_blocked <- pairwise.adonis2(
  spe_perc ~ sample_type,
  data = meta_cnts,
  strata = "core_id",
  nperm = 9999
)
