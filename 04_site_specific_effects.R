library(tidyverse)
library(tidypaleo)
library(vegan)
library(coin)
theme_set(theme_paleo(12))

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 9), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")),
         counts = rowSums(select(., Ablabesmyia:Zavreliella))) %>%
  filter(!counts < 12)

spe <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

spe_rowsums <- rowSums(spe)

spe_perc <- (spe/spe_rowsums)*100

sample_type <- factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 9), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern"))

spe_post <- spe_perc %>% 
  mutate(sample_type = cnts$sample_type) %>% 
  filter(sample_type == "post-tsunami") %>% 
  select(!sample_type)

spe_pre_top <- spe_perc %>% 
  mutate(sample_type = cnts$sample_type) %>% 
  filter(sample_type == "pre-tsunami top") %>% 
  select(!sample_type)

dm_post_vs_pre_top <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_post_vs_pre_top[i] <- vegdist(rbind(spe_post[i, ], spe_pre_top[i, ]), method = "bray")[1]
}

spe_pre_bot <- spe_perc %>% 
  mutate(sample_type = cnts$sample_type) %>% 
  filter(sample_type == "pre-tsunami bot") %>% 
  select(!sample_type)

dm_post_vs_pre_bot <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_post_vs_pre_bot[i] <- vegdist(rbind(spe_post[i, ], spe_pre_bot[i, ]), method = "bray")[1]
}

dm_pre_bot_vs_pre_top <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_pre_bot_vs_pre_top[i] <- vegdist(rbind(spe_pre_top[i, ], spe_pre_bot[i, ]), method = "bray")[1]
}

spe_modern <- spe_perc %>% 
  mutate(sample_type = cnts$sample_type) %>% 
  filter(sample_type == "modern") %>% 
  select(!sample_type)

dm_post_vs_modern <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_post_vs_modern[i] <- vegdist(rbind(spe_post[i, ], spe_modern[i, ]), method = "bray")[1]
}

dm_pre_top_vs_modern <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_pre_top_vs_modern[i] <- vegdist(rbind(spe_pre_top[i, ], spe_modern[i, ]), method = "bray")[1]
}

dm_pre_bot_vs_modern <- vector("numeric", length = nrow(spe_post))

for (i in 1:nrow(spe_post)) {
  dm_pre_bot_vs_modern[i] <- vegdist(rbind(spe_pre_bot[i, ], spe_modern[i, ]), method = "bray")[1]
}

corewise_diss <- tibble(
  no = letters[1:8],
  post_vs_pre_top = dm_post_vs_pre_top,
  post_vs_pre_bot = dm_post_vs_pre_bot,
  pre_bot_vs_pre_top = dm_pre_bot_vs_pre_top,
  post_vs_modern = dm_post_vs_modern,
  pre_top_vs_modern = dm_pre_top_vs_modern,
  pre_bot_vs_modern = dm_pre_bot_vs_modern
) %>% 
  pivot_longer(!no, names_to = "samples_compared", values_to = "dissimilarity") %>% 
  mutate(samples_compared = factor(samples_compared))

set.seed(12)
kruskal_diss <- kruskal_test(dissimilarity ~ samples_compared, data = corewise_diss,                                               distribution = approximate(nresample = 9999))

diss_plot <- ggplot(corewise_diss) +
  geom_boxplot(aes(y = dissimilarity, x = samples_compared)) +
#  facet_wrap(.~samples_compared, scales = "fixed", nrow = 1) +
#  scale_x_discrete(labels = NULL) +
  labs(y = "Dissimilarity") +
  geom_text(
    aes(x = Inf, y = Inf, 
        label = paste("p =", round(pvalue(kruskal_diss)[1], digits = 3))),
    hjust = 1.5, vjust = 2,
    inherit.aes = FALSE) +
  scale_x_discrete(labels = c(
    "post_vs_modern"     = "post-tsunami vs.\n modern",
    "post_vs_pre_bot"    = "post-tsunami vs.\n pre-tsunami bot",
    "post_vs_pre_top"    = "post-tsunami vs.\n pre-tsunami top",
    "pre_bot_vs_modern"  = "pre-tsunami top vs.\n modern",
    "pre_bot_vs_pre_top" = "pre-tsunami bot vs.\n pre-tsunami top",
    "pre_top_vs_modern"  = "pre-tsunami top vs.\n modern"
  )) +
  labs(x = NULL, title = "(D)")

saveRDS(diss_plot, "figures/dissimilarity_boxplot.rds")




# ta czesc ponizej byla dla pairwise tests i z tego zrezygnowalem, zeby nie mnozyc tych testow
# kruskal_post_vs_pre_top_and_pre_bot_vs_pre_top_prep <- corewise_diss %>% 
#   filter(samples_compared == "post_vs_pre_top" | samples_compared == "pre_bot_vs_pre_top")
# 
# set.seed(12)
# kruskal_post_vs_pre_top_and_pre_bot_vs_pre_top <- kruskal_test(dissimilarity ~ samples_compared, data = kruskal_post_vs_pre_top_and_pre_bot_vs_pre_top_prep,
#                                         distribution = approximate(nresample = 9999))
# 
# kruskal_post_vs_pre_bot_and_pre_bot_vs_pre_top_prep <- corewise_diss %>% 
#   filter(samples_compared %in% c("post_vs_pre_bot", "pre_bot_vs_pre_top"))
# 
# set.seed(12)
# kruskal_post_vs_pre_bot_and_pre_bot_vs_pre_top <- kruskal_test(dissimilarity ~ samples_compared, data = kruskal_post_vs_pre_bot_and_pre_bot_vs_pre_top_prep,
#                                         distribution = approximate(nresample = 9999))
# 
# kruskal_post_vs_pre_bot_and_post_vs_pre_top_prep <- corewise_diss %>% 
#   filter(samples_compared %in% c("post_vs_pre_bot", "post_vs_pre_top"))
# 
# set.seed(12)
# kruskal_post_vs_pre_bot_and_post_vs_pre_top <- kruskal_test(dissimilarity ~ samples_compared, data = kruskal_post_vs_pre_bot_and_post_vs_pre_top_prep,
#                                         distribution = approximate(nresample = 9999))
# 
# 
# # kruskal_post_vs_modern_and_post_vs_pre_bot_prep <- corewise_diss %>% 
# #   filter(samples_compared %in% c("post_vs_modern", "post_vs_pre_bot"))
# # 
# # set.seed(12)
# # kruskal_post_vs_modern_and_post_vs_pre_bot <- kruskal_test(dissimilarity ~ samples_compared, data = kruskal_post_vs_modern_and_post_vs_pre_bot_prep,
# #                                                             distribution = approximate(nresample = 9999))
# 
