library(tidyverse)
library(tidypaleo)
library(vegan)

theme_set(theme_paleo(12))

cnts <- read_csv("data/chiro_counts_raw.csv") %>% 
  mutate(depth_mid = (depth_from + depth_to)/2,
         sample_thick = depth_to - depth_from,
         core_id = factor(core_id)) %>%
  arrange(core_id, depth_from) %>% 
  mutate(sample_type = factor(rep(c("modern", "post-tsunami", "pre-tsunami top", "pre-tsunami bot"), 9), levels = c("pre-tsunami bot", "pre-tsunami top", "post-tsunami", "modern")),
         counts = rowSums(select(., Ablabesmyia:Zavreliella))) %>%
  filter(!counts < 12)

chiro <- cnts %>%
  select(core_id, sample_type, Ablabesmyia:Zavreliella) %>%
  mutate(
    total = rowSums(across(Ablabesmyia:Zavreliella))
  ) %>%
  mutate(
    across(Ablabesmyia:Zavreliella, ~ .x / total * 100)
  ) %>%
  select(-total)

pairwise_diss <- chiro %>%
  group_by(core_id) %>%
  group_modify(~ {
    
    dat <- .x
    
    comm <- dat %>%
      select(Ablabesmyia:Zavreliella) %>%
      as.matrix()
    
    d <- as.matrix(vegdist(comm, method = "bray"))
    
    pairs <- combn(seq_len(nrow(dat)), 2)
    
    tibble(
      i = pairs[1, ],
      j = pairs[2, ],
      sample1 = as.character(dat$sample_type[pairs[1, ]]),
      sample2 = as.character(dat$sample_type[pairs[2, ]]),
      dissimilarity = d[cbind(pairs[1, ], pairs[2, ])]
    )
    
  }) %>%
  ungroup()

get_comparison <- function(x, y) {
  
  paste(sort(c(x, y)), collapse = " | ")
}

comparison_labels <- c(
  "modern | post-tsunami" =
    "post_vs_modern",
  "modern | pre-tsunami bot" =
    "pre_bot_vs_modern",
  "modern | pre-tsunami top" =
    "pre_top_vs_modern",
  "post-tsunami | pre-tsunami bot" =
    "post_vs_pre_bot",
  "post-tsunami | pre-tsunami top" =
    "post_vs_pre_top",
  "pre-tsunami bot | pre-tsunami top" =
    "pre_bot_vs_pre_top"
)

pairwise_diss <- pairwise_diss %>%
  mutate(
    pair_name = map2_chr(sample1, sample2, get_comparison),
    samples_compared = unname(comparison_labels[pair_name]),
    samples_compared = factor(
      samples_compared,
      levels = c(
        "post_vs_modern",
        "post_vs_pre_bot",
        "post_vs_pre_top",
        "pre_bot_vs_modern",
        "pre_bot_vs_pre_top",
        "pre_top_vs_modern"
      )
    )
  )

kw_observed <- kruskal.test(
  dissimilarity ~ samples_compared,
  data = pairwise_diss
)

kw_observed

H_observed <- unname(kw_observed$statistic)

core_data <- split(chiro, chiro$core_id)

pair_data <- split(pairwise_diss, pairwise_diss$core_id)

one_permutation <- function() {
  
  permuted_pairs <- map2_dfr(
    core_data,
    pair_data,
    function(dat, pairs) {
      
      # randomly reassign the existing sample-class labels
      # within this core
      perm_labels <- sample(as.character(dat$sample_type))
      
      # labels corresponding to the original sample rows
      pairs$perm_sample1 <- perm_labels[pairs$i]
      pairs$perm_sample2 <- perm_labels[pairs$j]
      
      # determine the comparison represented by each distance
      pairs$pair_name <- map2_chr(
        pairs$perm_sample1,
        pairs$perm_sample2,
        get_comparison
      )
      
      pairs$samples_compared <- unname(
        comparison_labels[pairs$pair_name]
      )
      
      pairs
    }
  )
  
  # KW statistic for the permuted data
  kruskal.test(
    dissimilarity ~ samples_compared,
    data = permuted_pairs
  )$statistic
}

set.seed(12)

nperm <- 9999

perm_stats <- replicate(
  nperm,
  one_permutation()
)

p_perm <- (sum(perm_stats >= H_observed) + 1) /
  (nperm + 1)

p_perm

summary(perm_stats)

H_observed
p_perm

diss_plot <- ggplot(pairwise_diss) +
  geom_boxplot(aes(y = dissimilarity, x = samples_compared)) +
  labs(y = "Dissimilarity") +
  geom_text(
    aes(x = Inf, y = Inf, 
        label = paste("p =", round(p_perm, digits = 3))),
    hjust = 1.5, vjust = 2,
    inherit.aes = FALSE) +
  scale_x_discrete(labels = c(
    "post_vs_modern"     = "post-tsunami vs.\n modern",
    "post_vs_pre_bot"    = "post-tsunami vs.\n pre-tsunami bot",
    "post_vs_pre_top"    = "post-tsunami vs.\n pre-tsunami top",
    "pre_bot_vs_modern"  = "pre-tsunami bot vs.\n modern",
    "pre_bot_vs_pre_top" = "pre-tsunami bot vs.\n pre-tsunami top",
    "pre_top_vs_modern"  = "pre-tsunami top vs.\n modern"
  )) +
  labs(x = NULL, title = "(D)")

saveRDS(diss_plot, "figures/dissimilarity_boxplot.rds")
