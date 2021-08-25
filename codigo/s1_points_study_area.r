## Get occurrence points, create modelling area #####

FFP
setwd("PATH TO /Tapir_project")

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

library(sf)
options(sf_max.plot=1)
library(raster)

wgs <- 4326
utm17 <- 32717 # https://epsg.io/32717 # mainly here, and area of interest all within this area

## importar puntos
pts.all <- read.csv("data/Table_S1_occ_records.csv")
head(pts.all)
str(pts.all)

# convert to sf
pts.wgs <- st_as_sf(pts.all, coords = c("Longitude", "Latitude"), crs = wgs)
pts.wgs

# convert to UTM 17S
pts.utm <- st_transform(pts.wgs, crs = utm17)

# # write as shape
# st_write(pts.utm, "gis/s_utm/pts_all.shp", delete_layer = T)

## Bring in country shape - guardar shapefile de pais in gis/s_utm
# pais.utm <- st_read("gis/s_utm/co_ec_pe.shp")
# 
# plot(st_geometry(pais.utm))
# plot(st_geometry(pts.utm), add = T, pch = 16, col = "blue")

# cuantos por pais?
table(pts.utm$Country, useNA = "always") #
# ECU PER 
# 73  178 

## Crear area de estudio alrededor de puntos (incluyen dos puntos historicos de Huambo y Cutervo)

# buffer points by 25 k
pts.buff <- st_buffer(st_convex_hull(st_union(pts.utm)), dist = 25000)

# write modelling extent to shp
st_write(pts.buff, "gis/s_utm/pts_buffer.shp", delete_layer = T)

# Ajustar para origin de rasters, 
source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")
buff.ext <- adjExt(raster::extent(as.numeric(st_bbox(pts.buff)[c("xmin", "xmax", "ymin", "ymax")])))
buff.ext

## add 5 km to extent
buff.ext
buff.ext <- raster::extend(buff.ext, 5000)

## create template raster
r <- raster::raster(buff.ext, crs = utm17, res = 1000)
r

## create mask for study area
rmsk <- r
rmsk[] <- 1

rmsk <- raster::mask(rmsk, st_as_sf(pts.buff))
rmsk

plot(rmsk)
plot(pts.buff, add = T)
plot(pts.utm, add = T, pch = 16)

## thin points
source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/SDM/filter_pts.r")

pts.filt <- filter_pts(rmsk, pts.utm)
table(pts.filt$Country) # minus two historic in Peru.
# Ecuador    Peru 
# 73     110 
# total 73+110 = 183

save(pts.buff, pts.utm, pts.all, pts.filt, r, rmsk, file = "data/occ_pts_msk.rdata")

# crear polygon alrededor de puntos en Peru para Global Forest Watch analisis

# load("data/occ_pts_msk.rdata")
# 
# pts.pe <- subset(pts.utm, Country == "Peru")
# 
# plot(pts.buff)
# plot(pts.pe, add = T, pch = 16)
# 
# ## importar biomas (bajar desde pagina web: WWF y guardar)
# biomas <- st_read("K:/GIS_Info/Cover/wwf_ecos/e2017/biomes2017.shp")
# biomas
# unique(biomas$FIRST_BIOM)
# biomas <- st_make_valid(biomas)
# 
# 
# ## intersect biomas and points:
# pts.buff.wgs <- st_transform(pts.buff, crs = wgs)
# biomas.buff <- st_intersection(biomas, pts.buff.wgs)
# 
# # convert back to UTM
# biomas.clip.utm <- st_transform(biomas.buff, crs = utm17)
# 
# plot(biomas.clip.utm[, "FIRST_BIOM"])
# 
# pts.pe.buff <- st_buffer(st_convex_hull(st_union(pts.pe)), dist = 5000)
# 
# biomas.tapir <- st_intersection(biomas.clip.utm, pts.pe.buff)
# 
# plot(st_geometry(biomas.tapir))
# plot(st_geometry(pts.pe), add = T, pch = 16, col = "black")
# 
# st_write(biomas.tapir, "gis/s_utm/biomas_pe_rango.shp", delete_layer = T)
# 
# biomas.tapir.union <- st_union(biomas.tapir)
# 
# plot(biomas.tapir.union)
# plot(st_geometry(pts.pe), add = T, pch = 16, col = "black")
# 
# biomas.tapir.union.wgs <- st_transform(biomas.tapir.union, crs = wgs)
# st_write(biomas.tapir.union.wgs, "gis/biomas_pe_rango_union_wgs.kml", delete_layer = T)
