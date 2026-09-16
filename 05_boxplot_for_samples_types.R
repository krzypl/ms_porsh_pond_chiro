library(tidyverse)
library(tidypaleo)
library(vegan)
library(patchwork)

theme_set(theme_paleo(12))

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(
    depth_mid = (depth_from + depth_to) / 2,
    sample_thick = depth_to - depth_from,
    core_id = factor(core_id)
  ) %>%
  arrange(core_id, depth_from) %>% 
  mutate(
    sample_type = factor(
      rep(c("modern", "post-tsunami",
            "pre-tsunami top", "pre-tsunami bot"), 9),
      levels = c(
        "pre-tsunami bot",
        "pre-tsunami top",
        "post-tsunami",
        "modern"
      )
    ),
    counts = rowSums(select(., Ablabesmyia:Zavreliella))
  )

loi_and_sand <- read_csv("data/loi_and_sand.csv")

counts_raw_prep <- cnts %>% 
  select(Ablabesmyia:Zavreliella)

counts_raw <- rowSums(counts_raw_prep)

chiro_adds <- cnts %>% 
  select(core_id, volume, sample_type, depth_from) %>% 
  mutate(
    counts = counts_raw,
    concentration = counts / volume
  ) %>% 
  left_join(
    loi_and_sand,
    by = c("core_id", "depth_from")
  )

blocked_kw <- function(data, response, block = "core_id",
                       group = "sample_type", nperm = 9999) {
  
  # observed KW statistic
  observed <- kruskal.test(
    data[[response]] ~ data[[group]]
  )$statistic
  
  # split data into blocks
  blocks <- split(data, data[[block]])
  
  one_permutation <- function() {
    
    permuted <- map_dfr(blocks, function(dat) {
      
      dat[[group]] <- sample(dat[[group]])
      
      dat
    })
    
    kruskal.test(
      permuted[[response]] ~ permuted[[group]]
    )$statistic
  }
  
  set.seed(12)
  
  perm_stats <- replicate(
    nperm,
    one_permutation()
  )
  
  p_value <- (sum(perm_stats >= observed) + 1) /
    (nperm + 1)
  
  list(
    statistic = unname(observed),
    p.value = p_value,
    permutations = perm_stats
  )
}

kruskal_concentration <- blocked_kw(
  data = chiro_adds,
  response = "concentration",
  block = "core_id",
  group = "sample_type",
  nperm = 9999
)

kruskal_concentration$statistic
kruskal_concentration$p.value

kruskal_sand <- blocked_kw(
  data = chiro_adds,
  response = "no_sand",
  block = "core_id",
  group = "sample_type",
  nperm = 9999
)

kruskal_sand$statistic
kruskal_sand$p.value

spe4tr_prep <- counts_raw_prep %>% 
  mutate(
    core_id = cnts$core_id,
    sample_type = cnts$sample_type
  )

spe4tr <- spe4tr_prep %>% 
  select(!core_id & !sample_type) %>% 
  filter(!rowSums(.) < 12)

spe4tr_add <- spe4tr_prep %>% 
  select(core_id, sample_type) %>% 
  filter(!rowSums(counts_raw_prep) < 12)

tr_spe <- rarefy(
  spe4tr,
  min(rowSums(spe4tr))
)

tr_df <- tibble(
  core_id = cnts$core_id[which(counts_raw >= 12)],
  sample_type = cnts$sample_type[which(counts_raw >= 12)],
  tr = tr_spe
)

tr_test <- blocked_kw(
  data = tr_df,
  response = "tr",
  block = "core_id",
  group = "sample_type",
  nperm = 9999
)

tr_test$statistic
tr_test$p.value

kw_results <- tibble(
  variable = c(
    "Rarefied taxon richness",
    "Chironomid concentration",
    "Sand-grain counts",
    "Bray–Curtis dissimilarity"
  ),
  H = c(
    tr_test$statistic,
    kruskal_concentration$statistic,
    kruskal_sand$statistic,
    H_observed
  ),
  p = c(
    tr_test$p.value,
    kruskal_concentration$p.value,
    kruskal_sand$p.value,
    p_perm
  )
)

kw_results

# plots ----

tr_plot <- tr_df %>% 
  ggplot(aes(y = tr, x = sample_type, color = sample_type)) +
  geom_boxplot() +
  scale_color_manual(values = c(
    "modern" = "magenta",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
  breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot")) +
  geom_text(
    aes(x = Inf, y = Inf, 
        label = paste("p =", round(kw_results$p[kw_results$variable ==
                                                  "Rarefied taxon richness"], digits = 3))),
    hjust = 1.5, vjust = 2,
    inherit.aes = FALSE) +
  labs(x = NULL, y = "Taxon richness", title = "(A)")

conc_disb_plot <- chiro_adds %>% 
  ggplot(aes(y = concentration, x = sample_type, color = sample_type)) +
  geom_boxplot() +
  scale_color_manual(values = c(
    "modern" = "magenta",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
  breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot")) +
  geom_text(
    aes(x = Inf, y = Inf, 
        label = paste("p =",
                      round(kw_results$p[kw_results$variable ==
                                           "Chironomid concentration"],
                            digits = 3))),
    hjust = 1.5, vjust = 2,
    inherit.aes = FALSE) +
  labs(x = NULL, y= expression(Concentration~(hc~cm^{-3})), title = "(B)")

sand_disb_plot <- chiro_adds %>% 
  ggplot(aes(y = no_sand, x = sample_type, color = sample_type)) +
  geom_boxplot() +
  scale_color_manual(values = c(
    "modern" = "magenta",
    "post-tsunami" = "orange",
    "pre-tsunami top" = "darkblue",
    "pre-tsunami bot" = "blue"
  ),
  breaks = c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot")) +
  geom_text(
    aes(x = Inf, y = Inf, 
        label = paste("p ", format.pval(kw_results$p[kw_results$variable ==
                                                       "Sand-grain counts"],
                                        digits = 3, eps = 0.001))),
    hjust = 1.5, vjust = 2,
    inherit.aes = FALSE) +
  labs(x = NULL, y = expression(Sand~'>0.25'~mm~(n~cm^{-3})), color = "Sample class", title = "(C)")

diss_plot <- readRDS("figures/dissimilarity_boxplot.rds")

wrapped_plots <- wrap_plots(
  tr_plot +
    theme(legend.position = "none"),
  conc_disb_plot +
    theme(legend.position = "none"),
  sand_disb_plot +
    theme(legend.position = "bottom"),
  diss_plot +
    theme(legend.position = "botoom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)),
  nrow = 2
)


ggsave(filename="figures/fig5_boxplots.svg",
       plot = wrapped_plots,
       device = svg,
       width = 12.5,
       height = 9,
       units = "in")

ggsave(filename="figures/fig5_boxplots.pdf",
       plot = wrapped_plots,
       device = pdf,
       width = 12.5,
       height = 9,
       units = "in")

ggsave(filename="figures/fig5_boxplots.jpeg",
       plot = wrapped_plots,
       device = jpeg,
       width = 12.5,
       height = 9,
       units = "in")