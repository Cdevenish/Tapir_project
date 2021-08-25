
### Do model tuning. Final model and predict

FFP
## set working directory
setwd("PATH TO /Tapir_project")


source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

library(sf)
options(sf_max.plot=1)
library(raster)
library(ENMeval)
packageVersion("ENMeval") # [1] '0.3.1'  ## OJO version ahora es 2.0
library(maxnet)

wgs <- 4326
utm17 <- 32717


## Get predictors
load("data/predictors.rdata") # allStck 
rm(wc, dhi, dem)

names(allStck)
plot(allStck[[1]])
plot(allStck)


## Get points 
load("data/occ_pts_msk.rdata") # pts.buff, pts.utm, pts.all, pts.filt, r, rmsk
rm(pts.buff, pts.all, pts.utm)
pts.filt

## remvoe two historic points
pts.filt <- subset(pts.filt, !ID %in% c("PE165", "PE163"))

table(pts.filt$Source) # EC, PE, filtered
table(pts.filt$Country) # EC 73, PE 108. 181 in total 73+108

w.xls(as.data.frame.table(table(data.frame(pts.filt)[,c("Country", "RecordType")]))) # EC, PE, filtered

nCores <- parallel::detectCores()-2

# set.seed(1) # to reproduce results. Results stored

## Do tuning
mxTune <- ENMeval::ENMevaluate(occ = st_coordinates(pts.filt),
                      env = allStck,
                      method = "checkerboard1",
                      aggregation.factor = 3,
                      n.bg = 10000,
                      RMvalues = seq(0.5, 4, 0.5),
                      fc = c("L", "LQ", "LP", "LQP", "LQH", "LQHP", "LQHPT"),
                      #RMvalues = c(1, 2),
                      #fc = c("L", "LQ"),
                      algorithm = "maxnet",
                      parallel = T,
                      numCores = nCores,
                      rasterPreds=FALSE)

str(mxTune, max.level= 2)

head(mxTune@results)

head(mxTune@results[order(mxTune@results[, "avg.test.AUC"], decreasing = TRUE),], 10)
#     settings features rm train.AUC    avg.test.AUC var.test.AUC avg.diff.AUC var.diff.AUC
# 6   LQHP_0.5     LQHP 0.5 0.9417710    0.8929116 4.776094e-05   0.05940637 1.967219e-04
# 5    LQH_0.5      LQH 0.5 0.9264583    0.8923692 4.270917e-06   0.04176913 6.624354e-05
# 13    LQHP_1     LQHP 1.0 0.9311989    0.8899593 7.465609e-05   0.04892041 2.918516e-04

best.res <- mxTune@results[which.max(mxTune@results[, "avg.test.AUC"]),]
best.res

## Do final model

# make full data for absence loop
best.res$features[1] # L Q H P
best.res$rm[1]# 1

dir.exists(file.path(getwd(), "results"))

# same background poitns
head(mxTune@bg.pts)

mxFull <- dismo::maxent(x = allStck, p = st_coordinates(pts.filt), a = mxTune@bg.pts,
                        args = c('beta_lqp=0.5',
                                 "beta_hinge=0.5", # "beta_threshold=1" 
                                 'linear=true',
                                 'quadratic=true',
                                 'product=true',
                                 'threshold=false',
                                 'hinge=true'),
                        path  = file.path(getwd(), "results/mxFull"))

mxFull

mxPred <- dismo::predict(mxFull, x = allStck)

plot(mxPred)


mxFull@results[rownames(mxFull@results) %in% c("X.Training.samples","Training.AUC")]
#[1] 181.0000   0.9116


vars <- paste0(names(allStck), ".contribution")
varsImp <- data.frame(vars = sub(".contribution", "", vars), 
                      imp = mxFull@results[rownames(mxFull@results) %in% vars])

varsImp[order(varsImp[,"imp"], decreasing = T),]

# vars     imp
# 2           bio02 45.0831
# 8   dhi_gppqa_f.1 16.4420
# 5           bio15  9.9770
# 6   dhi_eviqa_f.1  9.2733
# 4           bio12  7.8535
# 9   dhi_gppqa_f.3  4.1706
# 1           bio01  3.2473
# 3           bio04  1.5563
# 7   dhi_eviqa_f.3  0.8562
# 13            Nss  0.4686
# 11            tpi  0.4658
# 10 dhi_lai8qa_f.2  0.4269
# 12          slope  0.1794


vars <- paste0(names(allStck), ".permutation.importance")
varsImp <- data.frame(vars = sub(".permutation.importance", "", vars), 
                      imp = mxFull@results[rownames(mxFull@results) %in% vars])

varsImp[order(varsImp[,"imp"], decreasing = T),]
# vars     imp
# 2           bio02 53.9955
# 5           bio15 13.6845
# 9   dhi_gppqa_f.3  7.0776
# 6   dhi_eviqa_f.1  6.2494
# 4           bio12  6.1809
# 8   dhi_gppqa_f.1  6.0901
# 1           bio01  1.6368
# 7   dhi_eviqa_f.3  1.4571
# 3           bio04  1.3468
# 13            Nss  0.8047
# 12          slope  0.7482
# 10 dhi_lai8qa_f.2  0.6395
# 11            tpi  0.0890


# apply 10% presence threshold.
t10 <- quantile(extract(mxPred, pts.filt), probs = 0.1, na.rm = T)
t10
mxBin <- reclassify(mxPred, rcl = c(0, t10, 0, t10,1, 1))

## do 3x3 filter
mxBin.mf <- focal(mxBin, w = matrix(1, nrow=3, ncol=3), fun = modal)

plot(stack(mxPred, mxBin.mf))

## Convert to Shape
mxBin.poly <- rasterToPolygons(mxBin.mf, fun=function(x){x==1}, n = 8, dissolve = TRUE)
mxBin.sf <- st_as_sf(mxBin.poly)
mxBin.sf

plot(mxBin.sf)

## WRITE AND SAVE
st_write(mxBin.sf, "gis/s_utm/mxBin_final.shp", delete_layer = T)
writeRaster(mxPred, filename = "gis/r_utm/mxPred_fin.tif", datatype = "FLT4S", overwrite = T)
writeRaster(mxBin.mf, filename = "gis/r_utm/mxBin_fin.tif", datatype = "INT1U", overwrite = T)

# Cortar por habitat y peru
hab <- st_read("gis/s_utm/piu-lam-caj_no habitat tapir.shp")
pais.utm <- st_read("gis/s_utm/co_ec_pe.shp")
pais.utm

## cortar para PEru
pe.utm <- subset(pais.utm, FIRST_NAME== "Peru")

plot(st_geometry(pe.utm))
plot(mxBin.sf, add =T)

mxBin.pe <- st_intersection(mxBin.sf, pe.utm)

plot(st_geometry(pe.utm))
plot(mxBin.pe, add = T)

## Cortar por habitat no apto
plot(st_geometry(pe.utm))
plot(hab, add =T) # habitat no apto

# A helper function that erases all of y from x:
st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
mxBin.pe.hab <- st_erase(mxBin.pe, hab)

plot(st_geometry(hab), col = "blue")
plot(mxBin.pe, add = T, col = "green")

plot(mxBin.pe.hab, add = T, col = "red")

st_write(mxBin.pe.hab, "gis/s_utm/mxBin_pe_no_hab.shp", delete_layer = TRUE)

save(mxFull, mxTune, varsImp, mxPred, mxBin, mxBin.mf, mxBin.poly, mxBin.pe.hab, mxBin.pe, 
     file = "results/tapir_mx_final.rdata")

