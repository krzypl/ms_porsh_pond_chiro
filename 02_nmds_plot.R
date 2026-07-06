library(tidyverse)
library(tidypaleo)
library(ggrepel)
library(vegan)
library(pairwiseAdonis)
theme_set(theme_paleo(12))


cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 9), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")),
         counts = rowSums(select(., Ablabesmyia:Zavreliella))) %>% 
  filter(!counts < 12) #remove two samples with very low counts (8 and 11)

spe <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

spe_rowsums <- rowSums(spe)

spe_perc <- (spe/spe_rowsums)*100

set.seed(12)
spe_perc_nmds <- metaMDS(spe_perc, trymax = 20, trace = TRUE, autotransform = TRUE, plot = TRUE)

stressplot(spe_perc_nmds, main = "Shepard plot")

nmds_4ploting <- cnts %>% 
  select(core_id, depth_mid, sample_type) %>% 
  mutate(NMDS1 = spe_perc_nmds$points[,1],
         NMDS2 = spe_perc_nmds$points[,2]) %>% 
  mutate(sample_type = factor(sample_type,
                       levels = c("modern",
                                  "post-tsunami",
                                  "pre-tsunami top",
                                  "pre-tsunami bot")))

meta_cnts <- data.frame(
  sample_type = cnts$sample_type,
  core_id = cnts$core_id)

set.seed(12)
pw_adonis <- pairwise.adonis2(
  spe_perc ~ sample_type,
  data = meta_cnts,
  nperm = 9999
)

pw_adonis_pvals <- c(
  paste("p =", pw_adonis$`modern_vs_post-tsunami`$`Pr(>F)`[1]),
  paste("p =", pw_adonis$`modern_vs_pre-tsunami top`$`Pr(>F)`[1]),
  paste("p =", pw_adonis$`modern_vs_pre-tsunami bot`$`Pr(>F)`[1]),
  paste("p =", pw_adonis$`post-tsunami_vs_pre-tsunami top`$`Pr(>F)`[1]),
  paste("p =", pw_adonis$`post-tsunami_vs_pre-tsunami bot`$`Pr(>F)`[1]),
  paste("p =", pw_adonis$`pre-tsunami top_vs_pre-tsunami bot`$`Pr(>F)`[1])
)

nmds_plot <- ggplot(nmds_4ploting) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = sample_type), size = 4) +
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 0, 8, 15)) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = core_id), size = 3) +
  coord_equal(ratio = 1) +
  labs(color = "Sample class") +
  scale_color_manual(values = c(
    "modern" = "magenta",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  )) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste(pw_adonis_pvals, collapse = "\n"),
    hjust = 1.5, vjust = 2
  )

ggsave(filename="figures/fig6_nmds_plot.svg",
       plot = nmds_plot,
       device = svg,
       width = 10,
       height = 7,
       units = "in")

ggsave(filename="figures/fig6_nmds_plot.pdf",
       plot = nmds_plot,
       device = pdf,
       width = 10,
       height = 7,
       units = "in")

ggsave(filename="figures/fig6_nmds_plot.jpeg",
       plot = nmds_plot,
       device = jpeg,
       width = 10,
       height = 7,
       units = "in")

gof <- goodness(spe_perc_nmds)
plot(spe_perc_nmds, type = "t", main = "Goodness of fit")
points(spe_perc_nmds, display = "sites", cex = gof * 300)
