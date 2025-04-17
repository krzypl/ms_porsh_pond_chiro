library(tidyverse)
library(vegan)

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(c(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 5), "post-tsunami", "pre-tsunami top", "pre-tsunami bot", rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 3)), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")))

spe <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

spe_rowsums <- rowSums(spe)

spe_perc <- (spe/spe_rowsums)*100

sample_type <- factor(c(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 5), "post-tsunami", "pre-tsunami top", "pre-tsunami bot", rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 3)), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern"))

spe_pre_bot_vs_post <- spe_perc %>% 
  mutate(sample_type = sample_type) %>% 
  filter(!sample_type %in% c("modern", "pre-tsunami top"))

distance_pre_bot_vs_post <- vegdist(spe_pre_bot_vs_post[,-length(spe_pre_bot_vs_post)], method = "bray")

set.seed(12)
(adonis_pre_bot_vs_post <- adonis2(distance_pre_bot_vs_post ~ spe_pre_bot_vs_post$sample_type))

spe_pre_top_vs_post <- spe_perc %>% 
  mutate(sample_type = sample_type) %>% 
  filter(!sample_type %in% c("modern", "pre-tsunami bot"))

distance_pre_top_vs_post <- vegdist(spe_pre_top_vs_post[,-length(spe_pre_bot_vs_post)], method = "bray")

set.seed(12)
(adonis_pre_bot_vs_post <- adonis2(distance_pre_top_vs_post ~ spe_pre_top_vs_post$sample_type))
