library(tidyverse)
library(tidypaleo)
theme_bw()

xrf_files <- list.files("data/xrf", full.names = TRUE)

xrf_df_full <- map_dfr(
  xrf_files,
  ~ read_csv(.x, skip = 5) %>%
    mutate(coreID = basename(.x))
) %>%
  rename(depth = "cm...1",
         Ti = "Ti (Value)",
         Si = "Si (Value)",
         Ca = "Ca (Value)",
         Fe = "Fe (Value)") %>% 
  select(depth, Ti, Si, Ca, Fe, coreID) %>% 
  filter(!is.na(Ti)) %>% 
  group_by(coreID) %>% 
  mutate(depth_cor = depth[1],
         depth = depth - depth_cor)


ggplot(xrf_df_full) + 
  geom_line(aes(x = depth, y = Ti)) +
  facet_wrap(.~coreID, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


ggplot(xrf_df_full) + 
  geom_line(aes(x = depth, y = Ca)) +
  facet_wrap(.~coreID, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

ggplot(xrf_df_full) + 
  geom_line(aes(x = depth, y = Fe)) +
  facet_wrap(.~coreID, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

ggplot(xrf_df_full) + 
  geom_line(aes(x = depth, y = Si)) +
  facet_wrap(.~coreID, scales = "free_y", ncol = 1)

