library(tidyverse)
library(vegan)
library(coin)

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

spe_post <- spe_perc %>% 
  mutate(sample_type = sample_type) %>% 
  filter(sample_type == "post-tsunami") %>% 
  select(!sample_type)

spe_pre_top <- spe_perc %>% 
  mutate(sample_type = sample_type) %>% 
  filter(sample_type == "pre-tsunami top") %>% 
  select(!sample_type)

dm_post_vs_pre_top <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_post_vs_pre_top[i] <- vegdist(rbind(spe_post[i, ], spe_pre_top[i, ]), method = "bray")[[1]]
}

spe_pre_bot <- spe_perc %>% 
  mutate(sample_type = sample_type) %>% 
  filter(sample_type == "pre-tsunami bot") %>% 
  select(!sample_type)

dm_post_vs_pre_bot <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_post_vs_pre_bot[i] <- vegdist(rbind(spe_post[i, ], spe_pre_bot[i, ]), method = "bray")[[1]]
}

dm_pre_bot_vs_pre_top <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_pre_bot_vs_pre_top[i] <- vegdist(rbind(spe_pre_top[i, ], spe_pre_bot[i, ]), method = "bray")[[1]]
}

spe_modern <- spe_perc %>% 
  mutate(sample_type = sample_type) %>% 
  filter(sample_type == "modern") %>% 
  select(!sample_type)

dm_post_vs_modern <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_post_vs_modern[i] <- vegdist(rbind(spe_post[i, ], spe_modern[i, ]), method = "bray")[[1]]
}

dm_pre_top_vs_modern <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_pre_top_vs_modern[i] <- vegdist(rbind(spe_pre_top[i, ], spe_modern[i, ]), method = "bray")[[1]]
}

dm_pre_bot_vs_modern <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_pre_bot_vs_modern[i] <- vegdist(rbind(spe_pre_bot[i, ], spe_modern[i, ]), method = "bray")[[1]]
}

corewise_diss <- tibble(
  core_id = unique(cnts$core_id),
  post_vs_pre_top = dm_post_vs_pre_top,
  post_vs_pre_bot = dm_post_vs_pre_bot,
  pre_bot_vs_pre_top = dm_pre_bot_vs_pre_top,
#  dm_post_vs_modern,
#  dm_pre_top_vs_modern,
#  dm_pre_bot_vs_modern
) %>% 
  pivot_longer(!core_id, names_to = "samples_compared", values_to = "dissimilarity") %>% 
  mutate(samples_compared = factor(samples_compared))

diss_plot <- ggplot(corewise_diss) +
  geom_boxplot(aes(y = dissimilarity)) +
  facet_wrap(.~samples_compared, scales = "fixed", nrow = 1)

kruskal_post_vs_pre_top_prep <- corewise_diss %>% 
  filter(!samples_compared == "pre-tsunami bot")
  
kruskal_post_vs_pre_top <- kruskal_test(dissimilarity ~ samples_compared, data = kruskal_post_vs_pre_top_prep,
                              distribution = approximate(nresample = 9999))

kruskal_post_vs_pre_bot_prep <- corewise_diss %>% 
  filter(!samples_compared == "pre-tsunami top")

kruskal_post_vs_pre_bot <- kruskal_test(dissimilarity ~ samples_compared, data = kruskal_post_vs_pre_bot_prep,
                                        distribution = approximate(nresample = 9999))


