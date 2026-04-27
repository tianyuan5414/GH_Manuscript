# This script plots a geomorphological map of the Pyrenees study area, with
# an Europe overview inset and a tight site-zoom inset showing the 6 plants.

# ── Packages ───────────────────────────────────────────────────────────────────
{
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(patchwork)
  library(maps)
  library(metR)
  library(elevatr)
  library(terra)
  library(sf)
  library(osmdata)
  library(dplyr)

  library(showtext)
  font_add(family = "arial", regular = file.path("Packages_R", "Fonts", "arial.ttf"))
  showtext_auto()
}

# ── Site coordinates ───────────────────────────────────────────────────────────
coords    <- read.csv("Data_R/Geo_Coordinates.csv")
coords_sf <- st_as_sf(coords, coords = c("E", "N"), crs = 4326)

# ── Map extents ────────────────────────────────────────────────────────────────
new_zoom_lon <- c(2.21, 2.31)     # ~10 km wide — geomorphological context
new_zoom_lat <- c(42.40, 42.45)
view_lon     <- c(2.208, 2.312)   # slightly padded for label repel room
view_lat     <- c(42.398, 42.452)

# ── Elevation rasters ──────────────────────────────────────────────────────────
{
  # Context DEM: z=12 ≈ 17 m resolution
  wide_bbox_sf     <- st_as_sf(data.frame(x = new_zoom_lon, y = new_zoom_lat),
                               coords = c("x", "y"), crs = 4326)
  elev_wide_raster <- get_elev_raster(locations = wide_bbox_sf, z = 12, clip = "bbox")
  elev_wide_terra  <- rast(elev_wide_raster)
  elev_wide_df     <- as.data.frame(elev_wide_terra, xy = TRUE) |>
    setNames(c("lon", "lat", "elevation"))

  # Hillshade from slope + aspect (sun from NW at 40°)
  slope_r  <- terrain(elev_wide_terra, "slope",  unit = "radians")
  aspect_r <- terrain(elev_wide_terra, "aspect", unit = "radians")
  hill_r   <- shade(slope_r, aspect_r, angle = 40, direction = 315)
  hill_df  <- as.data.frame(hill_r, xy = TRUE) |> setNames(c("lon", "lat", "hillshade"))

  # Site-level DEM: z=13 ≈ 8 m resolution (for site zoom inset)
  zoom_bbox_sf     <- st_as_sf(data.frame(x = c(2.248, 2.265), y = c(42.420, 42.428)),
                               coords = c("x", "y"), crs = 4326)
  elev_zoom_raster <- get_elev_raster(locations = zoom_bbox_sf, z = 13, clip = "bbox")
  elev_zoom_df     <- as.data.frame(rast(elev_zoom_raster), xy = TRUE) |>
    setNames(c("lon", "lat", "elevation"))
}

# ── OSM features ───────────────────────────────────────────────────────────────
{
  bbox_osm <- c(new_zoom_lon[1], new_zoom_lat[1], new_zoom_lon[2], new_zoom_lat[2])

  fetch_osm <- function(key, value) {
    opq(bbox = bbox_osm) |>
      add_osm_feature(key = key, value = value) |>
      osmdata_sf()
  }

  water_osm    <- fetch_osm("natural", "water")
  wood_osm     <- fetch_osm("natural", "wood")
  forest_lu    <- fetch_osm("landuse", "forest")
  glacier_osm  <- fetch_osm("natural", "glacier")
  rock_osm     <- fetch_osm("natural", "bare_rock")
  scree_osm    <- fetch_osm("natural", "scree")
  grassland_osm <- fetch_osm("natural", "grassland")
  heath_osm    <- fetch_osm("natural", "heath")

  forest_polys <- bind_rows(
    wood_osm$osm_polygons  |> select(geometry),
    forest_lu$osm_polygons |> select(geometry)
  )
}

# ── Terrain zone classification ────────────────────────────────────────────────
zone_levels <- c("Snow / Ice", "Alpine Rock", "Subalpine", "Montane Forest", "Valley Forest")

zone_palette <- c(
  "Snow / Ice"     = "#ddeeff",
  "Alpine Rock"    = "#c9bdb5",
  "Subalpine"      = "#d4e8c2",
  "Montane Forest" = "#7ab87a",
  "Valley Forest"  = "#4a8f4a"
)

classify_zones <- function(df) {
  df |> mutate(zone = factor(case_when(
    elevation > 2600 ~ "Snow / Ice",
    elevation > 2300 ~ "Alpine Rock",
    elevation > 2000 ~ "Subalpine",
    elevation > 1750 ~ "Montane Forest",
    TRUE             ~ "Valley Forest"
  ), levels = zone_levels))
}

elev_wide_class <- classify_zones(elev_wide_df)
elev_zoom_class <- classify_zones(elev_zoom_df)

# ── OSM feature palette ────────────────────────────────────────────────────────
osm_palette <- c(
  "Scree"     = "#c4ae97",
  "Bare Rock" = "#8d7b6e",
  "Grassland" = "#c8e6c2",
  "Heath"     = "#b0d49a",
  "Forest"    = "#2e7d32",
  "Water"     = "#4da6db"
)

label_sf <- function(sf_obj, label) sf_obj |> mutate(feature = label)

osm_all <- bind_rows(
  label_sf(scree_osm$osm_polygons,     "Scree"),
  label_sf(rock_osm$osm_polygons,      "Bare Rock"),
  label_sf(grassland_osm$osm_polygons, "Grassland"),
  label_sf(heath_osm$osm_polygons,     "Heath"),
  label_sf(forest_polys,               "Forest"),
  label_sf(water_osm$osm_polygons,     "Water")
)

# ── Panel A: Europe overview inset ────────────────────────────────────────────
europe_data <- map_data("world", region = c("France", "Spain", "Portugal", "Andorra",
                                             "Italy", "Switzerland", "Germany",
                                             "Belgium", "UK", "Morocco", "Algeria",
                                             "Tunisia", "Netherlands", "Luxembourg"))
lon_reg <- c(-2, 4)
lat_reg <- c(41.5, 44.5)

inset_europe <- ggplot() +
  geom_polygon(data = europe_data,
               aes(x = long, y = lat, group = group),
               fill = "grey82", color = "white", linewidth = 0.2) +
  annotate("rect",
           xmin = lon_reg[1], xmax = lon_reg[2],
           ymin = lat_reg[1], ymax = lat_reg[2],
           fill = NA, color = "#c0392b", linewidth = 0.9) +
  coord_fixed(xlim = c(-12, 20), ylim = c(34, 57), ratio = 1.3) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#d6eaf8", color = NA),
        panel.border     = element_rect(fill = NA, color = "grey40", linewidth = 0.6))

# ── Panel B: Site zoom inset (6 labeled plants) ───────────────────────────────
site_inset <- ggplot() +
  geom_raster(data = elev_zoom_class, aes(x = lon, y = lat, fill = zone)) +
  scale_fill_manual(values = zone_palette, guide = "none") +
  geom_contour(data = elev_zoom_df,
               aes(x = lon, y = lat, z = elevation),
               breaks = seq(2200, 2600, by = 25),
               color = "grey30", linewidth = 0.2, alpha = 0.6) +
  geom_text_contour(data = elev_zoom_df,
                    aes(x = lon, y = lat, z = elevation),
                    breaks       = seq(2200, 2600, by = 50),
                    size         = 2.2,
                    color        = "grey20",
                    stroke       = 0.2,
                    stroke.color = "white",
                    skip         = 0) +
  geom_point(data = coords, aes(x = E, y = N),
             shape = 17, color = "#c0392b", size = 2.5) +
  geom_label_repel(data = coords, aes(x = E, y = N, label = Plant),
                   size = 3, fontface = "bold",
                   box.padding = 0.5, point.padding = 0.4,
                   min.segment.length = 0, fill = "white", alpha = 0.9,
                   max.overlaps = Inf) +
  # Tight view: just enough room around the site cluster for labels
  coord_fixed(xlim = c(2.2515, 2.2595), ylim = c(42.4232, 42.4256), ratio = 1.3) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "grey30", linewidth = 0.6))

# ── Crop all layers to exact view extent ──────────────────────────────────────
{
  view_bbox <- st_bbox(c(xmin = view_lon[1], xmax = view_lon[2],
                         ymin = view_lat[1], ymax = view_lat[2]),
                       crs = st_crs(4326))

  osm_all_crop <- st_crop(osm_all, view_bbox)

  elev_wide_class_crop <- elev_wide_class |>
    dplyr::filter(lon >= view_lon[1], lon <= view_lon[2],
                  lat >= view_lat[1], lat <= view_lat[2])

  hill_df_crop <- hill_df |>
    dplyr::filter(lon >= view_lon[1], lon <= view_lon[2],
                  lat >= view_lat[1], lat <= view_lat[2])

  elev_wide_df_crop <- elev_wide_df |>
    dplyr::filter(lon >= view_lon[1], lon <= view_lon[2],
                  lat >= view_lat[1], lat <= view_lat[2])

  # Use the raster's actual data extent as the common clip rectangle so all
  # layers (terrain zone, hillshade, OSM features, contours) share the same border
  raster_xmin <- min(elev_wide_class_crop$lon)
  raster_xmax <- max(elev_wide_class_crop$lon)
  raster_ymin <- min(elev_wide_class_crop$lat)
  raster_ymax <- max(elev_wide_class_crop$lat)

  raster_bbox_sf <- st_bbox(
    c(xmin = raster_xmin, xmax = raster_xmax,
      ymin = raster_ymin, ymax = raster_ymax),
    crs = st_crs(4326)
  )

  osm_all_clip      <- st_crop(osm_all_crop, raster_bbox_sf)
  elev_contour_clip <- elev_wide_df_crop |>
    dplyr::filter(lon >= raster_xmin, lon <= raster_xmax,
                  lat >= raster_ymin, lat <= raster_ymax)
}

# ── Main: Geomorphological map ────────────────────────────────────────────────
map_geo <- ggplot() +
  # Terrain zones
  geom_raster(data = elev_wide_class_crop, aes(x = lon, y = lat, fill = zone)) +
  scale_fill_manual(values = zone_palette, name = "Terrain zone",
                    guide = guide_legend(order = 1, override.aes = list(alpha = 1))) +
  # Hillshade texture
  new_scale_fill() +
  geom_raster(data = hill_df_crop, aes(x = lon, y = lat, fill = hillshade), alpha = 0.35) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  # OSM land cover — clipped to raster rectangle
  new_scale_fill() +
  geom_sf(data = osm_all_clip, aes(fill = feature),
          color = NA, alpha = 0.82, inherit.aes = FALSE) +
  scale_fill_manual(values = osm_palette, name = "Land cover", guide = "none") +
  # Elevation contours at 100 m intervals
  geom_contour(data = elev_contour_clip,
               aes(x = lon, y = lat, z = elevation),
               breaks = seq(1600, 3000, by = 100),
               color = "grey30", linewidth = 0.15, alpha = 0.55) +
  # Elevation labels every 200 m (white stroke keeps them legible over terrain)
  geom_text_contour(data = elev_contour_clip,
                    aes(x = lon, y = lat, z = elevation),
                    breaks       = seq(1600, 3000, by = 200),
                    size         = 2.4,
                    color        = "grey20",
                    stroke       = 0.2,
                    stroke.color = "white",
                    skip         = 0) +
  # Study sites + labels
  geom_sf(data = coords_sf, shape = 17, color = "#c0392b",
          size = 2.5, inherit.aes = FALSE) +
  geom_label_repel(data = coords_sf, aes(label = Plant, geometry = geometry),
                   stat = "sf_coordinates",
                   size = 3.2, fontface = "bold",
                   box.padding = 0.5, point.padding = 0.4,
                   min.segment.length = 0, fill = "white", alpha = 0.9,
                   inherit.aes = FALSE) +
  # Limits set to raster extent — all layers share the same rectangle
  coord_sf(xlim = c(raster_xmin, raster_xmax),
           ylim = c(raster_ymin, raster_ymax),
           expand = FALSE, clip = "on") +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  theme_minimal(base_family = "arial") +
  theme(legend.position  = "right",
        legend.key.size  = unit(0.42, "cm"),
        legend.text      = element_text(size = 8),
        legend.title     = element_text(size = 8.5, face = "bold"),
        legend.spacing.y = unit(0.1, "cm"),
        panel.background = element_rect(fill = NA, color = NA),
        panel.border     = element_rect(fill = NA, color = "grey40", linewidth = 0.5))

# ── Final figure ──────────────────────────────────────────────────────────────
map_pyrenees <- map_geo +
  inset_element(inset_europe,
                left = 0, bottom = 0.72, right = 0.22, top = 1.0,
                align_to = "plot") +
  inset_element(site_inset,
                left = 0, bottom = 0.0, right = 0.38, top = 0.46,
                align_to = "plot")

map_pyrenees

ggsave("Figures_R/Map_Pyrenees.pdf", map_pyrenees, width = 10, height = 6, dpi = 300).
