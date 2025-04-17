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
#install.packages("package_2_install/maptools_1.1-8.tar.gz", repos = NULL, type = "source")
library(maptools)
#install.packages("package_2_install/ggsn_0.5.0.tar.gz", repos = NULL, type = "source")
library(ggsn)
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

ewbrks <- seq(-55.8, -55.72, by = 0.01)
nsbrks <- seq(46.862, 46.878, by = 0.01)

lake_panel <- tibble(
  xmin = -55.76443,
  xmax = -55.76192,
  ymin = 46.86920,
  ymax = 46.87092
)

bp_map <-
  autoplot.OpenStreetMap(bp2map) +
  scalebar(x.min = -55.774, x.max = -55.763,
           y.min = 46.864, y.max = 46.87,
           dist = 0.4, dist_unit = "km",
           st.dist = 0.05, st.bottom = FALSE,
           st.color = "white", st.size = 3,
           location = "bottomleft", box.color = "white",
           border.size = 0.1,
           transform = TRUE, model = "WGS84") +
  annotation_north_arrow(style = north_arrow_fancy_orienteering,
                         location = "tr") +
  scale_y_latitude(breaks = nsbrks, position = "right") +
  scale_x_longitude(breaks = ewbrks) +
  labs(x = NULL, y = NULL, title = "(B)") +
  theme_bw()

maps_wrapped <- wrap_plots(
  nf_map,
  bp_map)

final_map <- ggdraw(maps_wrapped) +
  draw_plot(inset_map, x = 0.135, y = 0.67, width = 0.25, height = 0.15)

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

depth_points <- read.csv("data/depth_sounding_data.csv", head = TRUE, sep = ";") %>% 
  rename(y = ycoord, x = xcoord)

depth_points_sf = st_as_sf(depth_points, coords = c("x", "y"))
depth_points_sf

depth_points_sf <- st_set_crs(depth_points_sf, value = "EPSG:4326")
plot(depth_points_sf)

depth_zero <- depth_points_sf %>% 
  filter(!depth == 0)

lake_contour <- st_read("data/tl09_contour.shp")

lake_bbox <- st_bbox(lake_contour)

lake_raster <- st_as_stars(lake_bbox, dx = 0.00001, dy = 0.00001)

st_crs(lake_raster) <- st_crs(lake_contour)

lake_raster_clipped <- lake_raster[lake_contour]

tps <- Tps(st_coordinates(depth_points_sf), depth_points_sf$depth)
lake_raster_clipped$tps_pred <- predict(tps, st_coordinates(lake_raster_clipped))
splain_lake <- lake_raster_clipped[lake_contour]

lake_raster_clipped_filled <- lake_raster_clipped
lake_raster_clipped_filled[[1]] <- lake_raster_clipped$tps_pred

lake_raster_map <- lake_raster_clipped_filled[lake_contour]

lake_raster_map$tps_pred <- (lake_raster_map$tps_pred + 1.675122)/100

breaks_seq <- seq(0, 1.8,
                  by = 0.6)

cores_coord <- read_csv("data/sediment_cores_coord.csv")

cores_coord_sf <- st_as_sf(cores_coord, coords = c("x", "y"))

cores_coord_sf <- st_set_crs(cores_coord_sf, value = "EPSG:4326")

lake_raster_map_plot <- 
  tm_shape(lake_raster_map) +
  tm_raster(
    col = "tps_pred",
    col.scale = tm_scale_intervals(
      values = "brewer.blues",
      breaks = breaks_seq
    ),
    col.legend = tm_legend(title = "Depth (m)")
  ) +
  tm_shape(lake_contour) + 
  tm_borders(lwd = 2, col = "black") +
  tm_shape(cores_coord_sf) +
  tm_symbols(size = 1) +
  tm_shape(cores_coord_sf) +
  tm_text(
    "core_id", 
    col = "black", 
    bgcol = "white",
    bgcol_alpha = 0.6,
    size = 0.8,
    ymod = 1,
    fontface = "bold") +
  tm_scalebar(position = c("left", "bottom"), text.size = 1) +
  tm_layout(legend.position = c("left", "top")) #+
#  tm_compass(position = c("right", "top"), text.size = 1, size = 2)

lake_raster_map_plot

tmap_save(lake_raster_map_plot, filename = "figures/fig_1_panel_c.pdf")
tmap_save(lake_raster_map_plot, filename = "figures/fig_1_panel_c.svg")
