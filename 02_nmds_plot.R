library(tidyverse)
library(vegan)


cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from)

spe <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

set.seed(12)
spe_nmds <- metaMDS(spe, trace = FALSE)
set.seed(12)
spe_nmds <- metaMDS(spe, previous.best = spe_nmds, trace = FALSE)

layout(matrix(1:2, ncol = 2))
plot(spe_nmds, main = "Chironomid NMDS plot"); stressplot(sol, main = "Shepard plot")

gof <- goodness(spe_nmds)
plot(spe_nmds, type = "t", main = "Goodness of fit")
points(spe_nmds, display = "sites", cex = gof * 300)
