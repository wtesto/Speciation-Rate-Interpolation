setwd("~/Desktop")

library(raster)
library(sp)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(scico)
library(smoothr)
library(gstat)

## we'll be working with a behrmann projection -- will use this string later
behrmann<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

##load a world map
world_temp<-ne_countries(returnclass = "sf")
world<-st_transform(world_temp, crs = behrmann)

##read points from CSV
points_temp<-read.csv("Phlegmariurus_occ.csv")

##note that longitude and latitude are in columns 4 and 3 in this CSV - adjust accordingly
points_temp<-st_as_sf(points_temp,coords = c(4,3), crs= '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
points<-st_transform(points_temp, crs = behrmann)


##create a bounding box around occurrence points for use in plotting 
limit <- st_buffer(points, dist = 500000) %>% st_bbox()

##quick check on occurrence points
spGroupColors<-scico(11, palette = "roma") ##make custom color ramp from SciCo palette
occMap <-
  ggplot() +
  geom_sf(data = world, color = "black", size = 0.05, fill = "#f7f7f7")+
  geom_sf(data = points , aes(geometry = geometry, color = SpeciesGroup), size = 0.5) +
  scale_color_manual(values=spGroupColors) +
  geom_sf(data = world, color = "black", size = 0.05, fill = NA)+
  coord_sf(
    xlim = c(limit["xmin"], limit["xmax"]),
    ylim = c(limit["ymin"], limit["ymax"])) +
  labs(x="Longitude", y="Latitude")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(panel.grid.major = element_line(colour = "#c9c9c9", 
                                        linetype = "dashed", 
                                        size = 0.2), 
        panel.background = element_rect(fill = "#f0f8ff"), 
        panel.border = element_rect(fill = NA),
          legend.key = element_blank())
occMap



##let's visualize speciation rates quickly
occSpeciation <-
  ggplot() +
  geom_sf(data = world, color = "black", size = 0.05, fill = "#f7f7f7")+
  geom_sf(data = points , aes(geometry = geometry, color=Rate), size = 0.5) +
  scale_fill_scico_d(direction = 1, palette="lajolla") +
  geom_sf(data = world, color = "black", size = 0.05, fill = NA)+
  coord_sf(
    xlim = c(limit["xmin"], limit["xmax"]),
    ylim = c(limit["ymin"], limit["ymax"])) +
  labs(x="Longitude", y="Latitude", fill = "Speciation Rate")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(panel.grid.major = element_line(colour = "#c9c9c9", 
                                        linetype = "dashed", 
                                        size = 0.2), 
        panel.background = element_rect(fill = "#f0f8ff"), 
        panel.border = element_rect(fill = NA),
        legend.key = element_blank())
occSpeciation



##create a buffer of 100km around each point
##this will be used as a mask for interpolation
buffer<-st_buffer(points$geometry, dist = 50000, nQuadSegs = 50) %>% 
  st_union()%>%
  smooth(method="ksmooth", smoothness=40) %>%  ##smooth the polygon
  st_intersection(world) %>%  ##clip offshore parts
  st_buffer(dist=1) %>%  #additional 1m buffer to overcome imperfect polygons in 'world'
  st_union()


##visualize the buffer
bufferMap <-
  ggplot() +
  geom_sf(data = world, color = "black", size = 0.05, fill = "#f7f7f7")+
  geom_sf(data = buffer , aes(geometry = geometry), size = 0.2, fill = "#e66b55", alpha=0.5) +
  coord_sf(
    xlim = c(limit["xmin"], limit["xmax"]),
    ylim = c(limit["ymin"], limit["ymax"])) +
  labs(x="Longitude", y="Latitude")+
  theme(panel.grid.major = element_line(colour = "#c9c9c9", 
                                        linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "#f0f8ff"), 
        panel.border = element_rect(fill = NA))

bufferMap

##convert buffer to SpatialPolygons format for use in clipping raster
buffer_sp<-as(buffer, "Spatial")


##create irregular grid within limits of buffer
grid <- buffer %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_make_grid(
    cellsize = c(5000, 5000), ##define cell size for interpolation here (in m) -- output size and processing times increase exponentially as cell size decreases
    what = "centers") %>%
  st_as_sf() %>%
  cbind(., st_coordinates(.)) %>% 
  st_drop_geometry() %>% 
  mutate(Z = 0) ##create extra field for interpolation

##rasterize grid
rasterGrid <- grid %>% 
  raster::rasterFromXYZ(
    crs = behrmann)

##create inverse distance weighted interpolation model
IDW_model <- gstat::gstat(
  formula = Rate ~ 1,   ##place colname of field you want to interpolate before "~"
  data = as(points, "Spatial"), ##convert points format
  nmax = 50, nmin = 10, ##nmax is the number of neighboring points used at any given spot; nmin is the minimum
  set = list(idp = 0.5)) #idp is the power of the inverse distance weighting - idp = 0 specify no distance effect - effect strength increases with higher numbers

## carry out interpolation
IDW_output <- interpolate(rasterGrid, IDW_model)

##clip interpolated output to extent of buffer
IDW_output_trimmed<-mask(IDW_output, buffer_sp)

##convert raster to a tibble so it can be plotted with ggplot
IDW_df <- rasterToPoints(IDW_output_trimmed) %>% as_tibble()
colnames(IDW_df) <- c("X", "Y", "Z")


##plot interpolated rates
IDW_rates <-
  ggplot() +
  geom_sf(data = world, color = "black", size = 0.02, fill = "#f7f7f7")+
  geom_sf(data = buffer , aes(geometry = geometry), size = 0.02) +
  geom_tile(data = IDW_df , aes(x = X, y = Y, fill = Z)) +
  scale_fill_scico(direction = 1, palette="lajolla", limits=c(min(IDW_df$Z),max(IDW_df$Z))) +
  geom_sf(data = world, color = "black", size = 0.02, fill = NA)+
  coord_sf(
    xlim = c(limit["xmin"], limit["xmax"]),
    ylim = c(limit["ymin"], limit["ymax"])) +
  labs(x="Longitude", y="Latitude", fill = "Speciation Rate")+
  theme(panel.grid.major = element_line(colour = "#c9c9c9", 
                                        linetype = "dashed", 
                                        size = 0.2), 
        panel.background = element_rect(fill = "#f0f8ff"), 
        panel.border = element_rect(fill = NA))

IDW_rates
