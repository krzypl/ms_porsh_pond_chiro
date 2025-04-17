library(tidyverse)
library(tidypaleo)
library(ggrepel)
library(vegan)
theme_set(theme_paleo(12))


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

set.seed(12)
spe_perc_nmds <- metaMDS(spe_perc, trymax = 20, trace = TRUE, autotransform = FALSE, plot = TRUE) #with the default wisconsing transformation stress is larger and shepard plot reveal weaker fit

stressplot(spe_perc_nmds, main = "Shepard plot")

nmds_4ploting <- cnts %>% 
  select(core_id, depth_mid) %>% 
  mutate(NMDS1 = spe_perc_nmds$points[,1],
         NMDS2 = spe_perc_nmds$points[,2],
         sample_type = factor(c(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 5), "post-tsunami", "pre-tsunami top", "pre-tsunami bot", rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 3)), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")))

nmds_plot <- ggplot(nmds_4ploting) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = sample_type), size = 4) +
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 0, 8, 15)) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = core_id), size = 3) +
  coord_equal(ratio = 1) +
  labs(color = "Sample type") +
  scale_color_manual(values = c(
    "modern" = "magenta",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ))

gof <- goodness(spe_perc_nmds)
plot(spe_perc_nmds, type = "t", main = "Goodness of fit")
points(spe_perc_nmds, display = "sites", cex = gof * 300)