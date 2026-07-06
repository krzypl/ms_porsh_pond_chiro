library(tidyverse)
library(tidypaleo)
library(vegan)
library(patchwork)
theme_set(theme_paleo(10))

#read data----

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 9), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")))

#tsunami layer thickness ----

ts_thick_pre_top <- cnts %>% 
  select(core_id, sample_type, depth_to) %>% 
  filter(sample_type == "pre-tsunami top")

ts_thick_post_bot <- cnts %>% 
  select(core_id, sample_type, depth_from) %>% 
  filter(sample_type == "post-tsunami")

ts_thick <- tibble(
  core_id = unique(cnts$core_id),
  thick = ts_thick_pre_top$depth_to - ts_thick_post_bot$depth_from
)

ts_thick_plot <- ggplot(ts_thick) + 
  geom_col(aes(x = thick, y = core_id))

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
  filter(max(rel_abund) > 8) %>% 
  ungroup()

chiro_order <- unique(chiro_red$taxon)

chiro_plot <- chiro_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, chiro_order)) %>% 
  ggplot(aes(x = rel_abund, y = core_id, fill = sample_type)) +
  geom_colh(width = 0.5, position = "dodgev") +
  facet_abundanceh(vars(taxon), rotate_facet_labels = 45,
                   dont_italicize = c("\\btype\\b", "\\bgroup\\b")) +
  labs(x = "Relative abundance (%)", y = "Core ID", fill = "Sample class") +
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

spe4tr_prep <- counts_raw_prep %>% 
  mutate(core_id = cnts$core_id, sample_type = cnts$sample_type)

spe4tr <- spe4tr_prep %>% 
  select(!core_id & !sample_type) %>% 
  filter(!rowSums(.) < 12)

spe4tr_add <- spe4tr_prep %>% 
  select(core_id, sample_type) %>% 
  filter(!rowSums(counts_raw_prep) < 12)

tr_spe <- rarefy(spe4tr, min(rowSums(spe4tr)))

rarecurve(spe4tr, min(rowSums(spe4tr)))

tr_final <- spe4tr_add %>% 
  mutate(value = tr_spe,
         param = "taxon richness") %>% 
  add_row(tibble("core_id" = c("TL09_C_30", "TL09_B_60"), 
                 "sample_type" = c("post-tsunami", "modern"),
                 "value" = c(0, 0),
                 "param" = c("taxon richness", "taxon richness")))

chiro_adds_long <- chiro_adds %>% 
  select(!volume & !depth_from & !depth_to) %>% 
  pivot_longer(!core_id & !sample_type, names_to = "param", values_to = "value") %>% 
  filter(!param == "loi550") %>% # nie ufam tym danym, do usuniecia
  add_row(tr_final) %>% 
  mutate(param = factor(param, levels = c("taxon richness", "counts", "concentration", "no_sand")),
         sample_type = factor(sample_type, levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")))

chiro_adds_plot <- ggplot(chiro_adds_long, aes(x = value, y = core_id, fill = sample_type)) +
  geom_colh(width = 0.5, position = "dodgev") +
  facet_wrap(
    . ~ param,
    scales = "free_x",
    nrow = 1,
    labeller = as_labeller(c(
      "counts" = "Counts~(n~hc)",
      "concentration" = "Concentration~(hc~cm^{-3})",
      "no_sand" = "Sand~'>0.25'~mm~(n~cm^{-3})",
      "taxon richness" = "Taxon~richness"
    ), label_parsed)
  ) +
  scale_fill_manual(
    values = c(
      "modern" = "magenta",
      "post-tsunami" = "orange",
      "pre-tsunami top" = "darkblue",
      "pre-tsunami bot" = "blue"
    ),
    breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot")
  ) +
  labs(x = NULL)

plots_wraped <- wrap_plots(
    chiro_plot,
    chiro_adds_plot +
      theme(axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(angle = 90, hjust = 0)) +
      labs(y = NULL),
    nrow = 1,
    widths = c(7, 2)
) 


ggsave(filename="figures/fig4_chiro_plot.svg",
       plot = plots_wraped,
       device = svg,
       width = 12.5,
       height = 10,
       units = "in")

ggsave(filename="figures/fig4_chiro_plot.pdf",
       plot = plots_wraped,
       device = pdf,
       width = 12.5,
       height = 10,
       units = "in")

ggsave(filename="figures/fig4_chiro_plot.jpeg",
       plot = plots_wraped,
       device = jpeg,
       width = 12.5,
       height = 10,
       units = "in")

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

diff_post_pre <- setdiff(c(chiro_pre_tsunami_top, chiro_pre_tsunami_bot), chiro_post_tsunami) #Althoghether 14 taxa missing in post-tsunami

diff_pre_post <- setdiff(chiro_post_tsunami, c(chiro_pre_tsunami_top, chiro_pre_tsunami_bot)) # no missing in pre-tsunami bot compared to post-tsunami

diff_pre_modern <- setdiff(c(chiro_pre_tsunami_top, chiro_pre_tsunami_bot), chiro_modern) # "Paratendipes nudisquama-type", "Psectrocladius calcaratus-type", "Tanytarsus chinyensis-type", "Zavreliella", "Cricotopus sylvestris-type", "Eukiefferiella claripennis-type", "Mnodiamesa",     "Tanypus"

diff_modern_pre <- setdiff(chiro_modern, c(chiro_pre_tsunami_top, chiro_pre_tsunami_bot)) #Paracladius

