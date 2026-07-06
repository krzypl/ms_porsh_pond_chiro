library(tidyverse)
library(stars)
library(sf)
library(tmap)
library(fields)
library(spData)
library(maps)
library(ggspatial)
library(rnaturalearth)
library(ggmap)
library(patchwork)
library(geodata)
library(grid)
library(rcartocolor)
library(cowplot)
library(metR)
library(OpenStreetMap)
library(ggspatial)
library(terra)
data(world)

#map - panels A and B ------
nf <- gadm("GADM", country="CAN", level=0, resolution = 2)

nf2plot <- st_as_sf(nf)

nf_map <- ggplot(nf2plot) +
  geom_sf() +
  coord_sf(xlim = c(-60, -52),
           ylim = c(46, 52),
           expand = TRUE) +
  theme_bw() +
  #  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering,
                         location = "tl") +
  labs(x = NULL, y = NULL, title = "(A)") +
  annotate("text", x = -56, y = 48.5, label = "Newfoundland") +
  annotate("text", x = -58.4, y = 46.8, label = "Burin Peninsula") +
  annotate("rect", xmin = -55.78, xmax = -55.73, ymin = 46.85, ymax = 46.88, color = "magenta", linewidth = 2) +
  annotate("text", x = -55.84, y = 46.68, label = "(B)") +
  annotate("segment", x = -57, xend = -55.7, y = 46.8, yend = 47, linewidth = 1,
           arrow = arrow(length = unit(0.2, "cm"))) 


inset_map <- ggplot(world) +
  geom_sf() +
  labs(x = NULL, y = NULL) +
  coord_sf(xlim = c(-120, -20),
           ylim = c(20, 60),
           expand = TRUE) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black")) +
  geom_rect(aes(xmin = -60 - 1, xmax = -52 + 1, ymin = 46 - 1, ymax = 52 + 1), color = "red", fill = NA, linewidth = 0.5)




bp2map_prep <- openmap(c(46.8625,-55.75), c(46.8775, -55.775),
                       type="esri-imagery",
                       mergeTiles = TRUE, 
                       #                       zoom = 10
)

bp2map <- openproj(bp2map_prep, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
bp2map_epsg <- openproj(bp2map, projection = "EPSG:3857")

ewbrks <- seq(-55.8, -55.72, by = 0.01)
nsbrks <- seq(46.862, 46.878, by = 0.01)

lake_panel <- tibble(
  xmin = -55.76443,
  xmax = -55.76192,
  ymin = 46.86920,
  ymax = 46.87092
)


wgs84_xbreaks_in_3857 <- function(lon_breaks, lat_for_transform) {
  pts <- st_as_sf(
    data.frame(lon = lon_breaks, lat = lat_for_transform),
    coords = c("lon", "lat"), crs = 4326
  )
  pts_3857 <- st_transform(pts, 3857)
  st_coordinates(pts_3857)[, "X"]
}

wgs84_ybreaks_in_3857 <- function(lat_breaks, lon_for_transform) {
  pts <- st_as_sf(
    data.frame(lon = lon_for_transform, lat = lat_breaks),
    coords = c("lon", "lat"), crs = 4326
  )
  pts_3857 <- st_transform(pts, 3857)
  st_coordinates(pts_3857)[, "Y"]
}

lab_lon <- function(x) paste0(sprintf("%.2f", abs(x)), "°", ifelse(x < 0, "W", "E"))
lab_lat <- function(x) paste0(sprintf("%.3f", abs(x)), "°", ifelse(x < 0, "S", "N"))

ewbrks <- seq(-55.8, -55.72, by = 0.01)
nsbrks <- seq(46.865, 46.878, by = 0.01)

xbrks_3857 <- wgs84_xbreaks_in_3857(ewbrks, lat_for_transform = mean(nsbrks))
ybrks_3857 <- wgs84_ybreaks_in_3857(nsbrks, lon_for_transform = mean(ewbrks))

stopifnot(length(xbrks_3857) == length(ewbrks))
stopifnot(length(ybrks_3857) == length(nsbrks))

lake_panel <- tibble(
  xmin = -55.76443,
  xmax = -55.76192,
  ymin = 46.86920,
  ymax = 46.87092
)

lake_3857 <- lake_panel |>
  rowwise() |>
  mutate(
    ll = st_transform(st_sfc(st_point(c(xmin, ymin)), crs = 4326), 3857),
    ur = st_transform(st_sfc(st_point(c(xmax, ymax)), crs = 4326), 3857),
    xmin_m = st_coordinates(ll)[1,1],
    ymin_m = st_coordinates(ll)[1,2],
    xmax_m = st_coordinates(ur)[1,1],
    ymax_m = st_coordinates(ur)[1,2],
    xmid_m = (xmin_m + xmax_m)/2,
    ymid_m = (ymin_m + ymax_m)/2
  ) |>
  ungroup()

bp_map <-
  autoplot.OpenStreetMap(bp2map_epsg) +
  coord_equal(expand = FALSE) +
  annotation_scale(bar_cols = c("gray", "white"), text_col = "white") +
  annotation_north_arrow(
    location = "tr",
    style = north_arrow_fancy_orienteering(
      text_col = "white",
      line_col = "white",
      fill = c("white")
    )
  ) +
  scale_x_continuous(
    breaks = xbrks_3857,
    labels = function(x) lab_lon(ewbrks)
  ) +
  scale_y_continuous(
    breaks = ybrks_3857,
    labels = function(y) lab_lat(nsbrks),
    position = "right"
  ) +
  annotate("text",
           x = lake_3857$xmid_m, y = lake_3857$ymid_m + 300,
           label = "Porsh Pond", colour = "deepskyblue", size = 4, fontface = "bold") +
  annotate("segment", x = lake_3857$xmid_m, xend = lake_3857$xmid_m,
           y = lake_3857$ymid_m + 250, yend = lake_3857$ymid_m, linewidth = 0.5, color = "deepskyblue",
           arrow = arrow(length = unit(0.2, "cm"))) +
  labs(x = NULL, y = NULL, title = "(B)") +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0))

maps_wrapped <- wrap_plots(
  nf_map,
  bp_map)

final_map <- ggdraw(maps_wrapped) +
  draw_plot(inset_map, x = 0.25, y = 0.76, width = 0.24, height = 0.15)

ggsave(filename="figures/fig_1_panel_a_b.pdf",
       plot = final_map,
       device = pdf,
       width = 11,
       height = 5,
       units = "in")

ggsave(filename="figures/fig_1_panel_a_b.jpg",
       plot = final_map,
       device = jpeg,
       width = 11,
       height = 5,
       units = "in")

ggsave(filename="figures/fig_1_panel_a_b.svg",
       plot = final_map,
       device = svg,
       width = 11,
       height = 5,
       units = "in")

#map, panel C -------
depth_points <- read.csv("data/depth_sounding_data.csv", header = TRUE, sep = ";") |>
  rename(y = ycoord, x = xcoord)

cores_coord <- read_csv("data/sediment_cores_coord.csv")

cores_sf <- st_as_sf(cores_coord, coords = c("x", "y"), crs = 4326)

cores_as_soundings <- cores_coord |>
  transmute(soundingID = core_id, x, y, depth = depth * 100)

soundings <- bind_rows(depth_points, cores_as_soundings)

soundings_sf <- st_as_sf(soundings, coords = c("x", "y"), crs = 4326)

# crs_m <- 32620
# p_m    <- st_transform(p_sf, crs_m)
# lake_m <- st_transform(lake_contour_sf, crs_m)

depth_points <- read.csv("data/depth_sounding_data.csv", head = TRUE, sep = ";") %>% 
  rename(y = ycoord, x = xcoord)

cores_coord_4dp <- cores_coord %>% 
  mutate(depth = depth*100)

names(cores_coord_4dp) <- c("soundingID", "x", "y", "depth")

depth_points <- depth_points %>% 
  add_row(cores_coord_4dp)

depth_points_sf = st_as_sf(depth_points, coords = c("x", "y"))

depth_points_sf <- st_set_crs(depth_points_sf, value = "EPSG:4326")
# 
# bb <- st_bbox(lake_m)
# lake_grid <- st_as_stars(bb, dx = 1, dy = 1)
# st_crs(lake_grid) <- st_crs(lake_m)
# lake_grid <- lake_grid[lake_m]
lake_contour <- st_read("data/tl09_contour.shp")

bb <- st_bbox(lake_contour)
lake_grid <- st_as_stars(bb, dx = 1, dy = 1)
st_crs(lake_grid) <- st_crs(lake_contour)
lake_grid <- lake_grid[lake_contour]

tps <-Tps(st_coordinates(depth_points_sf), depth_points_sf$depth, lambda = 0.0003)

lake_bbox <- st_bbox(lake_contour)
lake_raster <- st_as_stars(lake_bbox, dx = 0.00001, dy = 0.00001)
st_crs(lake_raster) <- st_crs(lake_contour)
lake_raster_clipped <- lake_raster[lake_contour]
lake_raster_clipped$tps_pred <- predict(tps, st_coordinates(lake_raster_clipped))
splain_lake <- lake_raster_clipped[lake_contour]

lake_raster_clipped_filled <- lake_raster_clipped
lake_raster_clipped_filled[[1]] <- lake_raster_clipped$tps_pred

lake_raster_map <- lake_raster_clipped_filled[lake_contour]

lake_raster_map$tps_pred <- (lake_raster_map$tps_pred + 1.675122)/100

breaks_seq <- seq(0, 1.8, by = 0.6)

crs_m <- 32621

cores_utm   <- st_transform(cores_sf, crs_m)
contour_utm <- st_transform(lake_contour, crs_m)
raster_utm <- st_warp(lake_raster_map, crs = st_crs(crs_m), method = "near")

lake_bathy_contours <- 
  tm_shape(raster_utm,
           unit = "m") +
  tm_raster(col = "tps_pred",
            col.scale = tm_scale_intervals(values="brewer.blues", breaks=breaks_seq),
            col.legend = tm_legend(title="Depth (m)")) +
  tm_shape(contour_utm) + tm_borders(lwd = 2, col = "black") +
  tm_shape(cores_utm) +
  tm_symbols(size = 1) +
  tm_shape(cores_utm) +
  tm_text(
    "core_id", 
    col = "black", 
    bgcol = "white",
    bgcol_alpha = 0.6,
    size = 0.8,
    ymod = 1,
    fontface = "bold") +
  tm_scalebar(position = c("left", "bottom"), text.size = 1, breaks = c(0, 50, 100)) +
  tm_layout(legend.position = c("left", "top"))

tmap_save(lake_bathy_contours, filename = "figures/fig_1_panel_c.pdf")
tmap_save(lake_bathy_contours, filename = "figures/fig_1_panel_c.svg")
