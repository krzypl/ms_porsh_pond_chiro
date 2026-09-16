library(tidyverse)
library(tidypaleo)
library(analogue)
library(patchwork)
library(vegan)
library(readxl)
theme_set(theme_paleo(12))


chiro_core_raw <-
  read_excel("data/chiro_tl09_c_110_counts_raw.xlsx",
             skip = 2) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))
  )

chiro_core_spe <- chiro_core_raw %>%
  mutate(depth = (depth_from + depth_to)/2) %>%
  select(depth, volume, Ablabesmyia:Zavrelia) %>%
  arrange(depth) %>%
  pivot_longer(Ablabesmyia:Zavrelia, names_to = "taxon", values_to = "count") %>%
  group_by(depth) %>%
  mutate(count_sum = sum(count),
         Concentration = count_sum/volume,
         perc_abund = 100*count/count_sum) %>%
  ungroup()

chiro_to_retain <- chiro_core_spe %>%
  group_by(taxon) %>%
  summarise(n_zeros = sum(perc_abund == 0),
            max_abund = max(perc_abund)) %>%
  ungroup() %>%
  filter(n_zeros < length(unique(chiro_core_spe$depth)) - 2 & max_abund >= 2) %>%
  pull(taxon)

spe_red <- chiro_core_spe %>%
  filter(taxon %in% chiro_to_retain)

chiro_to_retain_plot <- chiro_core_spe %>%
  group_by(taxon) %>%
  summarise(n_zeros = sum(perc_abund == 0),
            max_abund = max(perc_abund)) %>%
  ungroup() %>%
  filter(n_zeros < length(unique(chiro_core_spe$depth)) - 4 & max_abund >= 4) %>%
  pull(taxon)

spe_red_plot <- chiro_core_spe %>%
  filter(taxon %in% chiro_to_retain_plot)

length(unique(spe_red_plot$taxon))

spe_red %>%
  group_by(taxon) %>%
  summarise(x = sum(count))


spe_red %>%
  group_by(taxon) %>%
  summarise(x = sum(count))

tl09_coniss <- spe_red %>%
  mutate(perc_abund = perc_abund/100) %>%
  nested_data(qualifiers = depth, key = taxon, value = perc_abund) %>%
  nested_chclust_coniss()

y <- tibble(n_gr = tl09_coniss[[14]][[1]]$n_groups,
            dispersion = tl09_coniss[[14]][[1]]$dispersion,
            bs_dispersion = tl09_coniss[[14]][[1]]$broken_stick_dispersion)

ggplot(y) +
  geom_line(aes(x = n_gr, y = dispersion), color = "black") +
  geom_line(aes(x = n_gr, y = bs_dispersion), color = "red")

spe_plot <- ggplot(spe_red_plot, aes(x = perc_abund, y = depth)) +
  geom_colh(width = 1) +
  geom_colh(data = filter(spe_red_plot, depth %in% c(0.25)),
            aes(x = perc_abund, y = depth),
            fill = "magenta") +
  geom_colh(data = filter(spe_red_plot, depth %in% c(36.75)),
            aes(x = perc_abund, y = depth),
            fill = "orange") +
  geom_colh(data = filter(spe_red_plot, depth %in% c(43.25)),
            aes(x = perc_abund, y = depth),
            fill = "darkblue") +
  geom_colh(data = filter(spe_red_plot, depth %in% c(43.75)),
            aes(x = perc_abund, y = depth),
            fill = "blue") +
  scale_y_reverse() +
  facet_abundanceh(vars(taxon), rotate_facet_labels = 90,
                   dont_italicize = c("\\btype\\b", "\\bgroup\\b")) +
  geom_rect(
    mapping = aes(ymin = 37.5, ymax = 42.5, xmin = -Inf, xmax = Inf),
    alpha = 0.2,
    fill = "gold",
    inherit.aes = FALSE
  ) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)")

chiro_prc_prep <- spe_red %>%
  select(depth, taxon, perc_abund) %>%
  pivot_wider(id_cols = depth, names_from = taxon, values_from = perc_abund) %>%
  select(!depth)

set.seed(12)
chiro_prc <- prcurve(
  sqrt(chiro_prc_prep),
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
)

chiro_prc# variation explained by prc = 24%

chiro_prc_plot_prep <- tibble(
  depth = unique(spe_red$depth),
  param = "PrC score",
  value = chiro_prc$lambda)

count_plot_prep <- spe_red %>%
  select(depth, count_sum, Concentration) %>%
  pivot_longer(count_sum:Concentration, names_to = "param",
               values_to = "value") %>%
  add_row(chiro_prc_plot_prep) %>%
  mutate(param = gsub("count_sum", "Counts", param)) %>%
  group_by(param) %>%
  distinct(depth, .keep_all = TRUE) %>%
  ungroup()

# spe4tr <- chiro_core_raw %>%
#   arrange(depth_from) %>%
#   select(Ablabesmyia:Zavrelia)
# 
# tr_spe <- rarefy(spe4tr, min(rowSums(spe4tr)))
# 
# rarecurve(spe4tr, min(rowSums(spe4tr)))
# 
# tr_4plot <- tibble(depth = unique(count_plot_prep$depth),
#                    param = "Taxon richness",
#                    value = tr_spe)

count_plot_prep <- count_plot_prep %>%
#  add_row(tr_4plot) %>%
  mutate(param =
           factor(param, levels = c(#"Taxon richness",
                                    "Counts",
                                    "Concentration",
                                    "PrC score")))

count_plot <- ggplot(data = count_plot_prep, aes(x = value, y = depth)) +
  geom_lineh() +
  geom_point(size = 2) +
  geom_point(data = filter(count_plot_prep, depth %in% c(0.25)),
             aes(x = value, y = depth),
             color = "magenta", size = 2) +
  geom_point(data = filter(count_plot_prep, depth %in% c(36.75)),
             aes(x = value, y = depth),
             color = "orange", size = 2) +
  geom_point(data = filter(count_plot_prep, depth %in% c(43.25)),
             aes(x = value, y = depth),
             color = "darkblue", size = 2) +
  geom_point(data = filter(count_plot_prep, depth %in% c(43.75)),
             aes(x = value, y = depth),
             color = "blue", size = 2) +
  facet_geochem_gridh(vars(param), scales = "free",
                      units = c(
                     #   "Taxon richness" = NA,
                        "Counts" = "n hc",
                        "Concentration" = "hc cm⁻³",
                        "PrC score" = "24 %",
                        "CONISS" = "Total sum \n of squares")
  ) +
  scale_y_reverse() +
  layer_dendrogram(tl09_coniss, aes(y = depth),
                   param = "CONISS") +
  layer_zone_boundaries(tl09_coniss, aes(y = depth)) +
  geom_rect(
    mapping = aes(ymin = 37.5, ymax = 42.5, xmin = -Inf, xmax = Inf),
    alpha = 0.2,
    fill = "gold",
    inherit.aes = FALSE
  ) +
  labs(y = NULL, x = NULL) +
  rotated_facet_labels(
    angle = 90,
    direction = "x",
    remove_label_background = TRUE
  )


plots_wraped <- wrap_plots(
  spe_plot +
    theme(strip.background = element_blank()),
  count_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()),
  nrow = 1,
  widths = c(10, 4)
)

ggsave(filename="figures/fig5_chiro_tl09_c_110.svg",
       plot = plots_wraped,
       device = svg,
       width = 12.5,
       height = 7,
       units = "in")

ggsave(filename="figures/fig5_chiro_tl09_c_110.pdf",
       plot = plots_wraped,
       device = pdf,
       width = 12.5,
       height = 7,
       units = "in")

ggsave(filename="figures/fig5_chiro_tl09_c_110.jpg",
       plot = plots_wraped,
       device = jpeg,
       width = 12.5,
       height = 7,
       units = "in")