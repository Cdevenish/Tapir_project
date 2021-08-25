
FFP

## set working directory
setwd("PATH TO /Tapir_project")

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

library(sf)
options(sf_max.plot=1)
library(raster)

wgs <- 4326
utm17 <- 32717 

## set gis folders where climate, etc data is stored
# gis <- "K://GIS_info"


## importar puntos y area de estudio
load("data/occ_pts_msk.rdata") # pts.filt, r, rmsk
rm(pts.buff, pts.all)

# crear un area en wgs para cortar approx los rasters en wgs
# y ## agrandar un poco por si acaso
r.wgs <- extend(projectExtent(r, crs = wgs), 2)
r.wgs



## 1 World clim ####
# PRIMERO BAJAR archivos de worldclim y guardar carpeta... 
dir(file.path(gis, "Clim/wrldclm/v2"))
wc.fn <- list.files(file.path(gis, "Clim/wrldclm/v2"), pattern = "\\.tif$", full.names = T)

wc.stck <- stack(wc.fn)
proj4string(wc.stck)

# crop 
wc.stck <- crop(wc.stck, y = r.wgs)
wc.stck
plot(wc.stck)

beginCluster()
wc.utm <- projectRaster(wc.stck, to = rmsk, method = "bilinear")
endCluster()

names(wc.utm)
names(wc.utm) <- gsub(" ", "0", sprintf("bio%02s", 1:19))

wc.utm$bio01
plot(wc.utm$bio01)

## mask to projected aoi
wc.utm <- mask(wc.utm, rmsk, filename = "gis/r_utm/wc_utm.tif", bylayer =F, overwrite = T) # single file - import as brick later
wc.utm
# change names again
names(wc.utm)
names(wc.utm) <- gsub(" ", "0", sprintf("bio%02s", 1:19))

plot(wc.utm$bio01)

# check in memory
sapply(1:nlayers(wc.utm), function(x) inMemory(wc.utm[[x]]))

# readinto memory
wc.utm <- readAll(wc.utm)

# # to load
# wc.utm <- brick("gis/r_utm/wc_utm.tif")
# names(wc.utm)
# names(wc.utm) <- gsub(" ", "0", sprintf("bio%02s", 1:19))

# save(wc.utm, file = "data/predictors.rdata")

pairs(wc.utm, maxpixels = 500)

# which predictors?
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100) - OUT (similar to 4)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter - OUT
# BIO9 = Mean Temperature of Driest Quarter - OUT
# BIO10 = Mean Temperature of Warmest Quarter - OUT (almost same as 1)
# BIO11 = Mean Temperature of Coldest Quarter - OUT (almost same as 1)
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month - OUT (similar to 12)
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter - OUT (almost same as 13)
# BIO17 = Precipitation of Driest Quarter - OUT  (almost same as 14)
# BIO18 = Precipitation of Warmest Quarter - OUT
# BIO19 = Precipitation of Coldest Quarter - OUT

bio.OUT <- c("bio03", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio13", "bio14", "bio16", "bio17", "bio18", "bio19")

wc <- dropLayer(wc.utm, which(names(wc.utm) %in% bio.OUT))
pairs(wc, maxpixels = 500)


### DHI #####
## combined years. 
## Get rasters in, crop to broad extent in wgs, then project and mask
# Using NDVI, EVI, LAI, GPP
## Bajar estas capas (composite DHI 2003-2014) : http://silvis.forest.wisc.edu/data/DHIs/  y guardar

# names
# f1 Band 1 - cumulative
# f2 Band 2 - minimum
# f3 Band 3 - seasonality


dhi.fn <- list.files(file.path(gis, "cover/DHI"), pattern = "_f\\.tif$", full.names = T)
dhi.fn

dhi.rList <- lapply(dhi.fn, stack) # use stack as these are multilayer rasters
dhi.rList

dhi.crp.rList <- lapply(dhi.rList, crop, y = r.wgs)

beginCluster()
dhi.utm <- lapply(dhi.crp.rList, projectRaster, to = rmsk, method = "bilinear")
endCluster()

dhi.utm.msk <- lapply(dhi.utm, mask, mask = rmsk)

dhi.stck <- stack(dhi.utm.msk)
dhi.stck
dhi.stck[[1]]
plot(dhi.stck)

raster::pairs(dhi.stck, maxpixels = 500)

names(dhi.stck)

dhi.OUT <- c("dhi_eviqa_f.2", "dhi_gppqa_f.2","dhi_ndviqa_f.1","dhi_ndviqa_f.2","dhi_ndviqa_f.3","dhi_lai8qa_f.3", "dhi_lai8qa_f.1")

dhi <- dropLayer(dhi.stck, which(names(dhi.stck) %in% dhi.OUT))
raster::pairs(dhi, maxpixels = 500)

### DEM #####
## BAJAR Y GUARDAR UN DEM PRIMERO....

dem <- raster(file.path(gis, "DEM/SRTM/co_ec_pe.tif"))
dem
3/3600 # 3' 
3/3600 * 100 * 1000 # about 90 m resolution

dem.crop <- crop(dem, r.wgs)

dem.utm <- projectRaster(dem, rmsk, method = "bilinear")
dem.utm
plot(dem.utm)
origin(dem.utm)

terr <- terrain(dem.utm, opt = c("aspect","slope", "TPI"), units = "degrees")

## Make eastness and northness
Nss <- cos(terr$aspect* pi / 180) # "Northness (aspect)" # convert to radians for cos/sin functions
Ess <- sin(terr$aspect* pi / 180)  # eastness

dem <- stack(terr$tpi, terr$slope, Nss) # dem v correlated with bio1
names(dem)
names(dem) <- c("tpi", "slope", "Nss")

pairs(dem, maxpixels = 500)

## put all into stack
allStck <- stack(wc, dhi, dem)
plot(allStck)
pairs(allStck, maxpixels = 500)
allStck

sapply(1:nlayers(allStck), function(x) inMemory(allStck[[x]]))

save(wc, dhi, dem, allStck, file = "data/predictors.rdata")
