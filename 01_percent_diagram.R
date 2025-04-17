library(tidyverse)
library(tidypaleo)
theme_set(theme_paleo(10))

#read data----

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(c(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 5), "post-tsunami", "pre-tsunami top", "pre-tsunami bot", rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 3)), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")))

#chironomid percent diagram----
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
  mutate(rel_abund = count/count_sum*100)

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
  labs(x = "Relative abundance (%)", y = NULL, fill = "Sample type") +
  scale_fill_manual(values = c(
    "modern" = "magenta",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
    breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot")) +
  theme(legend.position = "bottom")

#additional data----
counts_raw_prep <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

counts_raw <- rowSums(counts_raw_prep)

loi_and_sand <- read_csv("data/loi_and_sand.csv")

chiro_adds <- cnts %>% 
  select(core_id, volume, sample_type, depth_from) %>% 
  mutate(counts = counts_raw,
         concentration = counts/volume) %>% 
  left_join(loi_and_sand, by = c("core_id", "depth_from"))

chiro_adds_long <- chiro_adds %>% 
  select(!volume & !depth_from & !depth_to) %>% 
  pivot_longer(!core_id & !sample_type, names_to = "param", values_to = "value")

chiro_aads_plot <- ggplot(chiro_adds_long, aes(x = value, y = core_id, fill = sample_type)) +
  geom_colh(width = 0.5, position = "dodgev") +
  facet_wrap(.~param, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c(
    "modern" = "green",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
  breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot")) +
  labs(x = "Concentration")

#identifing taxa missing in post-tusnami samples

chiro_pre_tsunami_bot <- chiro_perc %>% 
  filter(sample_type == "pre-tsunami bot") %>% 
  group_by(taxon) %>% 
  summarise(tot_sum = sum(count)) %>% 
  filter(!tot_sum == 0) %>% 
  pull(taxon)

chiro_pre_tsunami_top <- chiro_perc %>% 
  filter(sample_type == "pre-tsunami top") %>% 
  group_by(taxon) %>% 
  summarise(tot_sum = sum(count)) %>% 
  filter(!tot_sum == 0) %>% 
  pull(taxon)

chiro_post_tsunami <- chiro_perc %>% 
  filter(sample_type == "post-tsunami") %>% 
  group_by(taxon) %>% 
  summarise(tot_sum = sum(count)) %>% 
  filter(!tot_sum == 0) %>% 
  pull(taxon)

chiro_modern <- chiro_perc %>%
  filter(sample_type == "modern") %>% 
  group_by(taxon) %>% 
  summarise(tot_sum = sum(count)) %>% 
  filter(!tot_sum == 0) %>% 
  pull(taxon)

setdiff(chiro_pre_tsunami_bot, chiro_pre_tsunami_top) # "Cricotopus sylvestris-type", "Eukiefferiella claripennis-type", "Lasiodiamesa",  "Mnodiamesa", "Parachironomus varus-type", "Tanypus" 

setdiff(chiro_pre_tsunami_top, chiro_pre_tsunami_bot) #"Paratendipes nudisquama-type", "Xenoxchironomus"

setdiff(chiro_post_tsunami, chiro_pre_tsunami_bot) #no missing in pre-tsunami bot compared to post-tsunami

setdiff(chiro_post_tsunami, chiro_pre_tsunami_top) #no missing in pre-tsunami bot compared to post-tsunami

setdiff(chiro_pre_tsunami_bot, chiro_post_tsunami) #Althoghether 10 missing in post-tsunami

setdiff(chiro_pre_tsunami_top, chiro_post_tsunami) #

