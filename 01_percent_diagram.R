library(tidyverse)
library(tidypaleo)
theme_set(theme_paleo(10))

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 8), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")))

chiro_long <- cnts %>% 
  select(!depth_from & !depth_to & !sample_thick) %>% 
  pivot_longer(!core_id & !depth_mid & !volume & !sample_type, names_to = "taxon", values_to = "count")

chiro_sum <- chiro_long %>% 
  group_by(core_id, depth_mid) %>% 
  summarise(count_sum = sum(count)) %>%
  ungroup() %>% 
  mutate(volume = cnts$volume)

chiro_conc <- chiro_sum %>% 
  mutate(concentration = count_sum/volume) %>% 
  pivot_longer(!core_id & !depth_mid, names_to = "param", values_to = "value") %>% 
  filter(param != "volume")

chiro_perc <- chiro_long %>%
  left_join(chiro_sum) %>% 
  left_join(chiro_conc[which(chiro_conc$param == "concentration"),]) %>% 
  mutate(rel_abund = count/count_sum*100) %>% 
  select(!param) %>% 
  rename(concentration = value)

chiro_zero <- chiro_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(cnts$depth_mid)) #filter out missing taxa

chiro_red <- chiro_perc %>% 
  filter(!taxon %in% chiro_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 4) %>% 
  ungroup()

chiro_order <- unique(chiro_red$taxon)

chiro_plot <- chiro_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, chiro_order)) %>% 
  ggplot(aes(x = rel_abund, y = core_id, fill = sample_type)) +
  geom_colh(width = 0.5, position = "dodgev") + 
  facet_abundanceh(vars(taxon)) +
  facet_abundanceh(vars(taxon), rotate_facet_labels = 90) +
  labs(x = "Relative abundance (%)", y = NULL, fill = "Sample Type") +
  scale_fill_manual(values = c(
    "modern" = "green",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
    breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"))

chiro_conc_plot <- ggplot(chiro_red, aes(x = concentration, y = core_id, fill = sample_type)) +
  geom_lineh() +
  geom_point() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(component), grouping = vars(location)) +
  labs(x = NULL)

chiro_red_rev <- chiro_perc %>% 
  filter(!taxon %in% chiro_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) <= 4) %>% 
  ungroup()

chiro_plot_rare <- chiro_red_rev %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, chiro_order)) %>% 
  ggplot(aes(x = rel_abund, y = core_id, fill = sample_type)) +
  geom_colh(width = 0.5, position = "dodgev") + 
  facet_abundanceh(vars(taxon)) +
  facet_abundanceh(vars(taxon), rotate_facet_labels = 90) +
  labs(x = "Relative abundance (%)", y = NULL, fill = "Sample Type") +
  scale_fill_manual(values = c(
    "modern" = "green",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
  breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"))
