# PROJECT: Building hare SCR model
# SCRIPT: 02 - Trap locations and state space
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 10 Dec 2025
# COMPLETED: 10 Dec 2025
# LAST MODIFIED: 10 Dec 2025
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# Here I'll read in the trap locations (xj),
# keeping a shapefile and a .csv of x-y coordinates.
# Then, I'll prescribe a rectangular state space
# so that I can use the unif(x, y) prior on ACs si

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

# trap locations
trap.kml <- st_read(dsn = paste0(getwd(), "/Spatial data/", lyr = "final_live.kml"))

# unit boundaries (for visualization)
units <- st_read(dsn = paste0(getwd(), "/Spatial data/units_fixed_utm/", lyr = "units_fixed_utm.shp"))

# ______________________________________________________________________________
# 4. Clean trap location data ----
# ______________________________________________________________________________

traps.sf <- trap.kml %>%
  
  # add site name
  mutate(site = rep(c("1A", "1B", "1C",
                      "2A", "2B", "2C",
                      "3A", "3B", "3C",
                      "4A", "4B", "4C"),
                    each = 36)) %>%
  
  # keep only site and geometry
  dplyr::select(site, geometry) %>%
  
  # transform to UTM
  st_transform(crs = "epsg:32611")

# ______________________________________________________________________________
# 5. Clean unit boundaries ----
# ______________________________________________________________________________

units.sf <- units %>%
  
  # keep only necessary columns
  dplyr::select(name, geometry) %>%
  
  # rename name to site
  rename(site = name) %>%
  
  # arrange by site
  arrange(site)

# ______________________________________________________________________________
# 6. Plot both ----

# just one unit

# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_sf(data = units.sf %>% filter(site == "2A")) +
  
  geom_sf(data = traps.sf %>% filter(site == "2A")) +
  
  coord_sf(datum = st_crs(32611))

# ______________________________________________________________________________
# 7. Rectangular state space ----

# buffering the bounding box of the trap grid seems the most reasonable here
# YES I know that S for 2A and 2B will probably overlap
# deal with it!

# ______________________________________________________________________________

focal.traps.sf <- traps.sf %>% filter(site == "2A")
focal.unit.sf <- units.sf %>% filter(site == "2A")

# bounding box of traps
focal.bbox <- spatialEco::bbox_poly(focal.traps.sf)

plot(st_geometry(focal.bbox))
plot(st_geometry(focal.traps.sf), add = T)

# buffer and extract THAT bbox to make sure it's actually a rectangle
focal.S <- spatialEco::bbox_poly(st_buffer(focal.bbox, dist = 175))

ggplot() +
  
  theme_bw() +
  
  geom_sf(data = focal.S,
          fill = NA) +
  
  geom_sf(data = focal.unit.sf) +
  
  geom_sf(data = focal.traps.sf) +
  
  coord_sf(datum = st_crs(32611))

# 175 m buffer here seems reasonable (yes there will be overlap with 2B)

# ______________________________________________________________________________
# 8. Center S (and traps) on (0, 0) ----

# this is mostly for convenience. Will it affect computation? Doubtful

# ______________________________________________________________________________

# function
center_coords <- function (
    
  S,
  traps.sf
  
) {
  
  # center coordinates (what will end up being (0, 0))
  cent.x = (max(st_coordinates(focal.S)[ , 1]) + min(st_coordinates(focal.S)[ , 1])) / 2
  cent.y = (max(st_coordinates(focal.S)[ , 2]) + min(st_coordinates(focal.S)[ , 2])) / 2
  
  # center S coords
  focal.S.coords.uncenter <- st_coordinates(focal.S)[ , c(1:2)]

  focal.S.coords.center <- matrix(c(focal.S.coords.uncenter[ , 1] - cent.x,
                                    focal.S.coords.uncenter[ , 2] - cent.y),
                                  ncol = 2)
  
  # promote to sf
  S.center <- st_as_sf(as.data.frame(focal.S.coords.center),
                       coords = c("V1", "V2")) %>%
    
    summarize(geometry = st_combine(geometry)) %>%
    
    st_cast("POLYGON")
  
  # center trap coords
  focal.traps.coords.uncenter <- st_coordinates(focal.traps.sf)[ , c(1:2)]
  
  focal.traps.coords.center <- matrix(c(focal.traps.coords.uncenter[ , 1] - cent.x,
                                        focal.traps.coords.uncenter[ , 2] - cent.y),
                                      ncol = 2)
  
  # promote to sf
  traps.center <- st_as_sf(as.data.frame(focal.traps.coords.center),
                           coords = c("V1", "V2"))
  
  # return a list
  return(list(S.center,
              traps.center))
  
}

# use function
centered.sf <- center_coords(focal.S, focal.traps.sf)

# plot
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = centered.sf[[1]],
          fill = NA) +
  
  geom_sf(data = centered.sf[[2]]) +
  
  coord_sf(datum = st_crs(32611))

# ______________________________________________________________________________
# 9. Write to file ----
# ______________________________________________________________________________

# S
st_write(centered.sf[[1]], dsn = paste0(getwd(), "/Data for model/", lyr = "focal_S.shp"))

# traps
# shapefile
st_write(centered.sf[[2]], dsn = paste0(getwd(), "/Data for model/", lyr = "focal_traps.shp"))

# .csv
traps.df <- data.frame(id = 1:nrow(centered.sf[[2]]),
                       x = st_coordinates(centered.sf[[2]])[ , 1],
                       y = st_coordinates(centered.sf[[2]])[ , 2])

write.csv(traps.df, paste0(getwd(), "/Data for model/focal_traps.csv"))
          