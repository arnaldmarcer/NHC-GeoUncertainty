# Script to generate a percent tree cover layer for the study area, based on:
# Reference: https://www.science.org/doi/10.1126/science.1244693
# Downloaded from: Global Land Analysis & Discovery group, https://glad.umd.edu

source("scripts/R/init.R")

# FUNCTIONS ====================================================================================== #
# Downloads 2000 tree cover tiles
downloadTreeCoverTiles <- function(region.name){
    filename.prefix <- "Hansen_GFC2014_treecover2000_"
    for(t in tiles[[region.name]]){
        f <- paste0(dir.out, "/", filename.prefix, t, ".tif")
        if(file.exists(f)){
            cat("File", f, "exists. Skipping ...\n")
            next
        }
        file.url <- paste0(url, "tree_cover_2000/", filename.prefix, t, ".tif")
        # dwld.cmd <- paste0("wget --user=", auth$user, " --password=", auth$pwd,
        #                    " -P ", dir.out, " ",
        #                    file.url)
        dwld.cmd <- paste0("wget -P ", dir.out, " ", file.url)
        system(dwld.cmd)
    }
}

# Download 2000 water masks tiles
downloadWaterMaskTiles <- function(region.name){
    filename.prefix <- "Hansen_GFC2014_datamask_"
    for(t in tiles[[region.name]]){
        f <- paste0(dir.out, "/", filename.prefix, t, ".tif")
        if(file.exists(f)){
            cat("File", f, "exists. Skipping ...\n")
            next
        }
        file.url <- paste0(url, "water_mask/", filename.prefix, t, ".tif")
        # dwld.cmd <- paste0("wget --user=", auth$user, " --password=", auth$pwd,
        #                    " -P ", dir.out, " ",
        #                    file.url)
        dwld.cmd <- paste0("wget  -P ", dir.out, " ", file.url)
        system(dwld.cmd)
    }
}

# Apply water mask to tiles and change resolution to study resolution
applyWaterMask <- function(region.name){
    for(i in 1:length(tiles[[region.name]])){
        cat("Processing tree cover tile", tiles[[region.name]][i], "...\n")
        f.out.1km <- paste0(dir.out, "/treecover_2000_", tiles[[region.name]][i], "_1km.tif")
        if(file.exists(f.out.1km)){
            cat("File", f.out.1km, "exists. Skipping ...\n")
            next()
        }
        f.tc <- paste0(dir.out, "/Hansen_GFC2014_treecover2000_", tiles[[region.name]][i], ".tif")
        r.tc <- rast(f.tc)
        f.mk <- paste0(dir.out, "/Hansen_GFC2014_datamask_", tiles[[region.name]][i], ".tif")
        r.mk <- rast(f.mk)
        cat("  creating mask ...\n")
        r.mk[r.mk != 1] <- NA
        cat("  applying mask ...\n")
        r <- r.tc * r.mk
        f.out.30m <- paste0(dir.out, "/treecover_2000_", tiles[[region.name]][i], ".tif")
        cat("  writing 30m masked tile ...\n")
        writeRaster(r, filename = f.out.30m, overwrite = T)
        cat("  changing resolution to 0.008333333 decimal degrees ...\n")
        cmd <- paste("gdalwarp -tr 0.008333333 0.008333333 -r bilinear -overwrite", f.out.30m, f.out.1km)
        system(cmd)
        cat("  ==>", f.out.1km, "\n")
    }
}

# Merge all tiles into a single raster predictor layer for the study area matched to the rest of predictors
generatePredictorLayer <- function(){
    files <- list.files(paste0(dir.out), pattern = "*1km.tif", full.names = T)
    r <- rast(files[1])
    for(i in 2:length(files)){
        cat("Mosaicking. Adding tile", files[i], "...\n")
        r <- merge(r, rast(files[i]))
    }
    f.out <- paste0(data.processed.dir, "predictors/", region.name, "/30s/treecover.tif")
    cat("Projecting to study zone ...\n")
    r <- project(r, r.template, method = "bilinear")
    cat("Writing raster ...\n")
    writeRaster(r, filename = f.out, overwrite = T)
    cat("==>", f.out, "\n")
}

# TREE COVER LAYER CREATION ====================================================================== #
# List of needed tiles
tiles <- c(paste0("00N", paste0("_0", seq(40, 90, 10), "W")),
           c(paste0("10N", paste0("_0", seq(40, 90, 10), "W")), paste0("10N", paste0("_", seq(100, 100, 10), "W"))),
           c(paste0("20N", paste0("_0", seq(50, 90, 10), "W")), paste0("20N", paste0("_", seq(100, 120, 10), "W"))),
           c(paste0("30N", paste0("_0", seq(60, 90, 10), "W")), paste0("30N", paste0("_", seq(100, 130, 10), "W"))),
           c(paste0("40N", paste0("_0", seq(60, 90, 10), "W")), paste0("40N", paste0("_", seq(100, 130, 10), "W"))),
           c(paste0("10S", paste0("_0", seq(40, 90, 10), "W")), paste0("10N", paste0("_", seq(100, 100, 10), "W"))),
           paste0("20S", paste0("_0", seq(40, 90, 10), "W")),
           paste0("30S", paste0("_0", seq(50, 90, 10), "W")),
           paste0("40S", paste0("_0", seq(60, 80, 10), "W")))

url <-  "https://glad.umd.edu/mapdata/"
year <- "2000"
region.name <- "mexico_to_argentina"

dir.out <- paste0(data.raw.dir, "global_forest_cover/", year, "/", region.name)
dir.create(dir.out, showWarnings = F, recursive = T)
f.template <- paste0(data.processed.dir, "predictors/", region.name, "/30s/bio1.tif")
r.template <- rast(f.template)
cat("Downloading tree cover map tiles from", url, "...\n")
downloadTreeCoverTiles(region.name)
cat("Downloading water mask map tiles from", url, "...\n")
downloadWaterMaskTiles(region.name)
cat("Applying water masks to tree cover map tiles ...\n")
applyWaterMask(region.name)
cat("Generating tree cover predictor layer at 30s resolution ...\n")
generatePredictorLayer()

